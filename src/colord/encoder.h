/*******************************************************************************
 
 CoLoRd 
 Copyright (C) 2021, M. Kokot, S. Deorowicz, and A. Gudys
 https://github.com/refresh-bio/CoLoRd

 This program is free software: you can redistribute it and/or modify it under 
 the terms of the GNU General Public License as published by the Free Software 
 Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY 
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
 A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this 
 program. If not, see https://www.gnu.org/licenses/.

******************************************************************************/
#pragma once

#include "utils.h"
#include "reads_sim_graph.h"
#include "edit_script.h"
#include "in_reads.h"
#include <vector>
#include <algorithm>
#include "stats_collector.h"
#include "hm.h"
#include "hs.h"
#include "murmur64_hash.h"
#include "reference_reads.h"
#include "parallel_queue.h"
#include <chrono>
#include "../common/libs/libdivide/libdivide.h"

struct Anchor
{
	uint32_t len;
	uint32_t pos_in_enc;
	uint32_t pos_in_ref;

	Anchor(uint32_t len, uint32_t pos_in_enc, uint32_t pos_in_ref) : len(len), pos_in_enc(pos_in_enc), pos_in_ref(pos_in_ref) {};
};


struct Candidate
{
	bool shouldReverse = false;
	uint32_t ref_read_id;
	std::vector<Anchor> anchors;
	uint32_t tot_anchor_len;
};

struct EncodeCandidates
{	
	std::vector<Candidate> candidates;
};

enum class AnalyseRefReadRes { EmptyIntersection, TooManyMatches, TooLowAnchors, Accept };

enum class AnalyseRefReadWithKmersRes { NoAnchors, CorespondingKmersIncompatibile, Accept };

//using CMmers = CMmersSTL;
class CMmersHashMapLP;
class CKmersHashMapLP;
class CKmersHashSetLP;
class CMmersHashMapDuplicateOptimizedLP;
//using CMmers = CMmersHashMapLP;
using CMmers = CMmersHashMapDuplicateOptimizedLP;
using CKmers = CKmersHashMapLP;

class CBloomFilter
{
#define BLOOM_1_FUNC
	uint64_t* data;
	uint64_t* raw_data;
	uint64_t size_in_bits;
	uint64_t size_mask;
	uint64_t word_mask;
	const uint64_t sub_bits = 6;
	const uint64_t sub_mask = (1ull << sub_bits) - 1ull;

	void allocate(uint64_t no_items)
	{
		if (no_items < 64)
			no_items = 64;

		size_in_bits = no_items * 32;

		// Round up to the power of 2
		if ((size_in_bits & (size_in_bits - 1)))
		{
			while ((size_in_bits & (size_in_bits - 1)))
				size_in_bits &= size_in_bits - 1;
			size_in_bits *= 2;
		}
		
		uint64_t size_in_words = size_in_bits / 64;

		size_mask = size_in_bits - 1ull;
		word_mask = size_mask >> 6;

		raw_data = new uint64_t[size_in_words + 8];
		data = raw_data + 8 - (((uint64_t)raw_data) / 8) % 8;

//		std::fill_n(data, size_in_words, 0ull);
		std::fill_n(data, size_in_words, ~0ull);
	}

	void deallocate()
	{
		if (raw_data)
		{
			delete[] raw_data;
			raw_data = nullptr;
		}
	}

	uint64_t basic_hash(uint64_t h)
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
/*		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;*/
	
		return h;
	}

	void set_bit(uint64_t hash_val)
	{
		data[hash_val & word_mask] &= ~(1ull << (hash_val >> 58));
	}

	void set_bits(uint64_t hash_val)
	{
		uint64_t m = (1ull << (hash_val & 63)) + (1ull << (hash_val >> 58));
		data[(hash_val >> 6) & word_mask] &= ~m;
	}

	bool test_bit(uint64_t hash_val)
	{
		return (data[hash_val & word_mask] & (1ull << (hash_val >> 58))) == 0;
	}

	bool test_bits(uint64_t hash_val)
	{
		uint64_t m = (1ull << (hash_val & 63)) + (1ull << (hash_val >> 58));
		return (data[(hash_val >> 6) & word_mask] & m) == 0;
	}

	void insert(uint64_t x)
	{
		uint64_t h = basic_hash(x);

#ifdef BLOOM_1_FUNC
		set_bit(h);
#else
		set_bits(h);
#endif
	}

public:
	explicit CBloomFilter(const read_t& read, uint32_t m, bool canonical = false) : raw_data(nullptr)
	{
		if (read_len(read) < m)
			return;

		allocate(read_len(read) - m + 1);

		anchor_type mask = (1ull << (2 * m)) - 1;

		if (canonical)
		{
			anchor_type mmer{}, rev{};
			uint32_t pos = 0;
			for (; pos < m - 1; ++pos)
			{
				assert(read[pos] < 4); //only ACGT
				mmer <<= 2;
				mmer += read[pos];

				rev >>= 2;
				rev += ((uint64_t)(3 - read[pos])) << (2 * (m - 1));
			}

			for (; pos < read_len(read); ++pos)
			{
				assert(read[pos] < 4); //only ACGT
				mmer <<= 2;
				mmer += read[pos];
				mmer &= mask;

				rev >>= 2;
				rev += ((uint64_t)(3 - read[pos])) << (2 * (m - 1));
				auto can = mmer < rev ? mmer : rev;

				insert(can);
			}
		}
		else
		{
			anchor_type mmer{};
			uint32_t pos = 0;
			for (; pos < m - 1; ++pos)
			{
				assert(read[pos] < 4); //only ACGT
				mmer <<= 2;
				mmer += read[pos];
			}

			for (; pos < read_len(read); ++pos)
			{
				assert(read[pos] < 4); //only ACGT
				mmer <<= 2;
				mmer += read[pos];
				mmer &= mask;
				insert(mmer);
			}
		}
	}

/*	uint64_t prefetch(uint64_t x)
	{
		uint64_t h = basic_hash(x);
		uint64_t h_main = h & size_mask_cl;

#ifdef _WIN32
		_mm_prefetch((const char*)(data + h_main), _MM_HINT_T0);
#else
			__builtin_prefetch(data + h_main);
#endif

		return h;
	}*/

//	bool __declspec(noinline) test(uint64_t x)
	bool test(uint64_t x)
	{
		uint64_t h = basic_hash(x);

#ifdef BLOOM_1_FUNC
		return test_bit(h);
#else
		return test_bits(h);
#endif
	}

/*	bool test_hint(uint64_t h)
	{
		uint64_t h_main = h & size_mask_cl;

#ifdef BLOOM_2_FUNC
		h >>= 32;
		uint64_t b0 = h & 0x1ff;
		h >>= 9;
		uint64_t b1 = h & 0x1ff;

		return test_bit(h_main, b0) && test_bit(h_main, b1);
#else
		h >>= 32;
		uint64_t b0 = h & 0x1ff;
		h >>= 9;
		uint64_t b1 = h & 0x1ff;
		h >>= 9;
		uint64_t b2 = h & 0x1ff;

		return test_bit(h_main, b0) && test_bit(h_main, b1) && test_bit(h_main, b2);
#endif
	}*/

	~CBloomFilter()
	{
		deallocate();
	}
};



class CEncoder
{
	CStatsCollector stats;

	CParallelQueuePopWaiting<CCompressPack>& compress_queue;
	const CReferenceReads& reference_reads;	
	CParallelPriorityQueue<std::vector<es_t>>& compressed_queue;
	CParallelPriorityQueue<std::vector<es_t>>& edit_script_for_qual_queue;
	std::vector<es_t> current_encoded_reads;
	uint32_t anchor_len;
	double minFractionOfMmersInEncodeToAlwaysEncode;
	double minFractionOfMmersInEncode;
	double maxMatchesMultiplier;
	double editScriptCostMultiplier;
	uint32_t minPartLenToConsiderAltRead;
	uint32_t maxRecurence;	
	uint32_t minAnchors;
	bool is_fastq;
	uint32_t modulo;
	libdivide::divider<uint64_t> modulo_div;
	uint32_t kmerLen;
	DataSource dataSource;

	CEntropyEstimator entropyEstimator;
	AnalyseRefReadRes AnalyseRefRead(CMmers& encode_mmers, CBloomFilter& bloom_mmers, const read_t& enc_read, const read_t& ref_read, Candidate& candidate, int decision);
	AnalyseRefReadWithKmersRes AnalyseRefReadWithKmers(CKmers& enc_kmers, CBloomFilter& bloom_kmers, const read_t& enc_read, const read_t& ref_read, Candidate& candidate, CKmersHashSetLP& common_kmers);

	bool KmerBasedAnchors(CKmersHashMapLP& enc_kmers, CBloomFilter& bloom_mmers, 
		const read_t& enc_read, const read_t& ref_read,
		const read_t& rev_compl_ref_read, Candidate& candidate, Candidate& rev_compl_candidate,
		const std::vector<kmer_type>& common_kmers, EncodeCandidates& res);

	void MmerBasedAnchors(CMmers& encode_mmers, CBloomFilter& bloom_mmers,
		const read_t& enc_read, const read_t& ref_read, const read_t& rev_compl_ref_read,
		Candidate& candidate, Candidate& rev_compl_candidate, int decision, EncodeCandidates& res,
		bool& refused_too_many_matches, bool& refules_too_low_anchors);
	
	std::vector<std::tuple<anchor_type, uint32_t, uint32_t>> get_aligned_mmers_LIS(
		const std::vector<std::pair<uint32_t, anchor_type>>& sorted_mmers_enc_vec,
		const std::vector<std::pair<uint32_t, anchor_type>>& sorted_mmers_ref_vec,
		const std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& mmers_ref_read_vec
	);

	void AddPlainRead(const read_t& read);
	void AddPlainReadWithN(const read_t& read);

	//second must be sorted
	uint64_t get_number_of_matches(std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& first, std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& second);

	std::vector<std::pair<uint32_t, anchor_type>> Convert(const std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& in, size_t len);

	uint32_t MergeAnchors(const std::vector<std::tuple<anchor_type, uint32_t, uint32_t>>& inputAnchors, std::vector<Anchor>& res);

	uint32_t AdjustAnchors(std::vector<Anchor>& anchors, uint32_t new_enc_start, uint32_t new_enc_end);

	EncodeCandidates prepareEncodeCandidates(const read_t& enc_read, const std::vector<uint32_t>& neighbours);
	
	EncodeCandidates prepareEncodeCandidatesHiFi(const read_t& enc_read, const std::vector<uint32_t>& neighbours, const std::vector<std::vector<kmer_type>>& common_kmers);
	
	EditDistRes GetEditDist(read_view refPart, read_view encPart, uint32_t frag_no, uint32_t n_fragments);

	uint32_t CountDeletions(std::string_view es);

	std::string_view GetEditScriptEntropyInput(std::string_view es);

	bool EncodeWithEditScript(const EditDistRes& ed, read_view refPart, read_view encPart, uint32_t frag_no, uint32_t n_fragments);

	bool EncodeWithAlternativeRead(const std::vector<Candidate>& candidates, uint32_t level, read_view encPart, uint32_t cur_pos_in_encode_read, uint32_t end_enc,
		std::vector<Candidate>& alt_candidates);

	void StoreFrag(std::string& currentBigEditScirpt, uint32_t level, es_t& encoded_read, uint32_t reference_read_id, uint32_t main_ref_id,
		uint32_t& last_pos_in_ref, uint32_t cur_pos_in_ref_read, bool& first, bool shouldReverse);

	void EncodePart(uint32_t level, read_view encode_read, uint32_t frag_no, uint32_t n_fragments,
		uint32_t end_enc, uint32_t end_ref, std::vector<Candidate>& candidates,
		uint32_t reference_read_id, 
		read_t &ref_read,
		uint32_t main_ref_id, uint32_t cur_pos_in_ref_read, uint32_t cur_pos_in_encode_read,
		std::string& currentBigEditScirpt, es_t& encoded_read,
		uint32_t& last_pos_in_ref, bool& first);

	void AddEncodedReadWithCandidates(const read_view encode_read, std::vector<Candidate>& candidates, uint32_t level, es_t& encoded_read, uint32_t main_ref_id, bool& first);

	void fixOverlapingInEncodeRead(std::vector<Anchor>& anchors);
	void fixOverlapingInReferenceRead(std::vector<Anchor>& anchors);

	void processComprElem(const CCompressElem& compr_elem);

	bool pop(CCompressPack& compress_pack);
	std::chrono::duration<double> waitOnQueue;
public:
	explicit CEncoder(bool verbose,
		CParallelQueuePopWaiting<CCompressPack>& compress_queue,
		const CReferenceReads& reference_reads, 
		CParallelPriorityQueue<std::vector<es_t>>& compressed_queue,
		CParallelPriorityQueue<std::vector<es_t>>& edit_script_for_qual_queue,
		uint32_t anchor_len,
		double minFractionOfMmersInEncodeToAlwaysEncode,
		double minFractionOfMmersInEncode,
		double maxMatchesMultiplier,
		double editScriptCostMultiplier,
		uint32_t minPartLenToConsiderAltRead,
		uint32_t maxRecurence,	
		uint32_t minAnchors,
		bool is_fastq,
		uint32_t modulo,
		uint32_t kmerLen,
		DataSource dataSource) :

		stats(verbose),
		compress_queue(compress_queue),
		reference_reads(reference_reads),		
		compressed_queue(compressed_queue),
		edit_script_for_qual_queue(edit_script_for_qual_queue),
		anchor_len(anchor_len),
		minFractionOfMmersInEncodeToAlwaysEncode(minFractionOfMmersInEncodeToAlwaysEncode),
		minFractionOfMmersInEncode(minFractionOfMmersInEncode),
		maxMatchesMultiplier(maxMatchesMultiplier),
		editScriptCostMultiplier(editScriptCostMultiplier),
		minPartLenToConsiderAltRead(minPartLenToConsiderAltRead),
		maxRecurence(maxRecurence),		
		minAnchors(minAnchors),
		is_fastq(is_fastq),
		modulo(modulo),
		kmerLen(kmerLen),
		dataSource(dataSource)
	{
		modulo_div = libdivide::divider<uint64_t>(modulo);
	}

	void Encode();

	double GetWaitOnQueueTime()
	{
		return waitOnQueue.count();
	}

};
