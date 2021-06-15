#pragma once
#include "defs.h"
#include "kmer_filter.h"
#include "in_reads.h"
#include "reference_reads.h"
#include "parallel_queue.h"
#include "queues_data.h"
#include "params.h"
#include "ref_reads_accepter.h"
#include "hm.h"
#include "hm_compact.h"
#include "hs.h"
#include "murmur64_hash.h"
#include "reference_genome.h"
#include <memory>

#define USE_CINT_HM

#ifdef USE_CINT_HM
using kmers_to_reads_compacted_t = cint_hash_map_lp<uint32_t, MurMur32Hash>;
#else
using kmers_to_reads_compacted_t = hash_map_lp<uint32_t, uint32_t, std::equal_to<uint32_t>, MurMur32Hash>;
#endif

using hash_set_t = hash_set_lp<uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;

class CKmersToReads
{
#ifdef USE_CINT_HM
	uint32_t suffix_len_bits = 31;
#else
	uint32_t suffix_len_bits = 31;
#endif
	uint32_t prefix_len_bits;
	uint64_t suffix_mask = (1ull << suffix_len_bits) - 1;
	std::vector<kmers_to_reads_compacted_t> hash_tables;
public:
	CKmersToReads(uint32_t kmer_len, uint64_t total_kmers, double fill_factor_kmers_to_reads)
	{
		uint32_t kmer_len_bits = kmer_len * 2;
		if (kmer_len_bits > suffix_len_bits)
			prefix_len_bits = kmer_len_bits - suffix_len_bits;
		else
			prefix_len_bits = 0;

		uint32_t n_hash_tables = 1ul << prefix_len_bits;
		uint64_t expected_single_ht_elems = total_kmers / n_hash_tables;
#ifdef USE_CINT_HM
		expected_single_ht_elems = 32; //TODO: in some experiments (joi dataset) it turn out that memory requirements are lower if 16 is set, but it should be inspected 
		hash_tables.assign(n_hash_tables, kmers_to_reads_compacted_t(static_cast<size_t>(expected_single_ht_elems / fill_factor_kmers_to_reads), fill_factor_kmers_to_reads, MurMur32Hash{}));
#else
		expected_single_ht_elems = 16; //TODO: in some experiments (joi dataset) it turn out that memory requirements are lower if 16 is set, but it should be inspected 
		hash_tables.assign(n_hash_tables, kmers_to_reads_compacted_t(std::numeric_limits<uint32_t>::max(), static_cast<size_t>(expected_single_ht_elems / fill_factor_kmers_to_reads), fill_factor_kmers_to_reads, std::equal_to<uint32_t>{}, MurMur32Hash{}));
#endif
	}
	void insert(kmer_type kmer, uint32_t val)
	{
		uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_len_bits);
		uint32_t suffix = static_cast<uint32_t>(kmer & suffix_mask);
		hash_tables[prefix].insert_fast(std::make_pair(suffix, val));
	}

	void prefetch_prefix(kmer_type kmer)
	{
		uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_len_bits);

#ifdef _WIN32
		_mm_prefetch((const char*)(hash_tables.data() + prefix), _MM_HINT_T0);
#else
		__builtin_prefetch(hash_tables.data() + prefix);
#endif
	}

	void prefetch_suffix(kmer_type kmer)
	{
		uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_len_bits);
		uint32_t suffix = static_cast<uint32_t>(kmer & suffix_mask);

		hash_tables[prefix].prefetch(suffix);
	}

	auto find(kmer_type kmer)
	{
		uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_len_bits);
		uint32_t suffix = static_cast<uint32_t>(kmer & suffix_mask);

#ifdef USE_CINT_HM
		return std::make_pair(hash_tables[prefix].find(suffix), hash_tables[prefix].local_end());
#else
		return std::make_pair(hash_tables[prefix].find(suffix), hash_tables[prefix].local_end());
#endif
	}

	void PrintMemoryUsage()
	{
		uint32_t prefix = 0;
		uint32_t n_prefixes = 1ul << prefix_len_bits;
		uint64_t tot_size_bytes{};
		for (; prefix < n_prefixes; ++prefix)
			tot_size_bytes += hash_tables[prefix].allocated_size() * sizeof(kmers_to_reads_compacted_t::value_type);
		std::cerr << "Kmers to reads memory usage: " << (tot_size_bytes / 1024 / 1024) << "MiB\n";
	}
};


#ifdef USE_BETTER_PARALLELIZATION_IN_GRAPH
class CReadsSimilarityGraphInternalThreads;
#endif // USE_BETTER_PARALLELIZATION_IN_GRAPH

class CReadsSimilarityGraph
{
	CParallelQueuePopWaiting<CCompressPack>& compress_queue;
	CReferenceReads& reference_reads;	
	const CKmerFilter& filteredKmers;
	uint32_t kmer_len;
	uint32_t max_candidates;
	uint32_t maxKmerCount;
	uint32_t current_read_id{};

	//kmers_to_reads_t kmers; 
	CKmersToReads kmers;

	uint32_t current_out_elem_id{};	
	CCompressPack current_out_queue_elem;

	std::vector<std::pair<read_t*, std::vector<kmer_type>>> getAcceptedKmers(const read_pack_t& reads_pack);

	void processReadsPack(const read_pack_t& reads_pack);
	void processReadsPackHiFi(const read_pack_t& reads_pack);

	void processReferenceGenome(CReferenceGenome* ref_genome);

	uint32_t id_in_reference{}; //id in reference set, which does not keep reads containing N's. In modes that do not track all reads, not all reads are in reference set

	ReferenceReadsMode referenceReadsMode;

	CRefReadsAccepter& ref_reads_accepter;

	int n_compression_threads;
	DataSource dataSource;

#ifdef USE_BETTER_PARALLELIZATION_IN_GRAPH
	std::unique_ptr<CReadsSimilarityGraphInternalThreads> internalThreads;
#else
#ifdef MEASURE_THREADS_TIMES
	std::vector<double>	twTimes;
#endif
#endif // USE_BETTER_PARALLELIZATION_IN_GRAPH

	
public:
	explicit CReadsSimilarityGraph(CParallelQueue<read_pack_t>& reads_queue,
		CParallelQueuePopWaiting<CCompressPack>& compress_queue,
		CReferenceReads& reference_reads,
		CReferenceGenome* reference_genome,
		const CKmerFilter& filteredKmers, uint32_t kmer_len, uint32_t max_candidates, uint32_t maxKmerCount,
		ReferenceReadsMode referenceReadsMode,
		CRefReadsAccepter& ref_reads_accepter,
		double refReadsFraction,
		int n_compression_threads,
		DataSource dataSource,
		double fill_factor_kmers_to_reads,
		bool verbose);

#ifdef MEASURE_THREADS_TIMES
	std::vector<double>& GetTwTimes();
#endif
	~CReadsSimilarityGraph();
};
