#pragma once
#include "defs.h"
#include <string>
#include <iostream>

//#include <unordered_map>
#include "hs.h"
#include "murmur64_hash.h"
#include "libs/libdivide/libdivide.h"

class CCompactedKmers
{
	using hash_set_t = hash_set_lp<uint32_t, std::equal_to<uint32_t>, MurMur32Hash>;
	uint32_t suffix_len_bits = 31;
	uint32_t prefix_len_bits;
	uint64_t suffix_mask = (1ull << suffix_len_bits) - 1;
	mutable std::vector<hash_set_t> hash_tables;

#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
	std::vector<uint64_t> counts_per_prefix;
#endif

public:
	CCompactedKmers(uint32_t kmer_len, uint64_t total_kmers, double fill_factor_filtered_kmers)
	{
		uint32_t kmer_len_bits = kmer_len * 2;
		if (kmer_len_bits > suffix_len_bits)
			prefix_len_bits = kmer_len_bits - suffix_len_bits;
		else
			prefix_len_bits = 0;

		uint32_t n_hash_tables = 1ul << prefix_len_bits;
		uint64_t expected_sinle_ht_elems = total_kmers / n_hash_tables;
		if (!expected_sinle_ht_elems)
			expected_sinle_ht_elems = 16;
		hash_tables.assign(n_hash_tables, hash_set_t(std::numeric_limits<uint32_t>::max(), static_cast<size_t>(expected_sinle_ht_elems / fill_factor_filtered_kmers), fill_factor_filtered_kmers, std::equal_to<uint32_t>{}, MurMur32Hash{}));
#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
		counts_per_prefix.resize(n_hash_tables);
#endif
	}
	void insert(kmer_type kmer
#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
		, uint32_t count
#endif
	)
	{
		uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_len_bits);
		uint32_t suffix = static_cast<uint32_t>(kmer & suffix_mask);
//		hash_tables[prefix].insert(suffix);
		hash_tables[prefix].insert_fast(suffix);
#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
		counts_per_prefix[prefix] += count;
#endif
	}

	auto check(kmer_type kmer) const
	{
		uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_len_bits);
		uint32_t suffix = static_cast<uint32_t>(kmer & suffix_mask);
		return hash_tables[prefix].check(suffix);
	}


#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
	uint64_t GetCountPerPrefix(uint32_t prefix) const
	{
		return counts_per_prefix[prefix];
	}
	void ReleaseCountsPerPrefix()
	{
		counts_per_prefix.clear();
		counts_per_prefix.shrink_to_fit();
	}
#endif

	uint64_t GetMemoryUsage() const
	{
		uint32_t prefix = 0;
		uint32_t n_prefixes = 1ul << prefix_len_bits;
		uint64_t tot_size_bytes{};
		for (; prefix < n_prefixes; ++prefix)
			tot_size_bytes += hash_tables[prefix].allocated_size() * sizeof(hash_set_t::key_type);
		return tot_size_bytes + sizeof(std::vector<hash_set_t>) * hash_tables.size();
	}

	void PrintMemoryUsage() const
	{
		std::cerr << "Filtered k-mers memory usage: " << (GetMemoryUsage() / 1024 / 1024) << "MiB\n";
	}

	uint32_t GetNHashTables() const
	{
		return hash_tables.size();
	}

};

class CKmerFilter
{
	kmer_type hash_mm(kmer_type x) const;
	CCompactedKmers kmers;
	uint64_t total_count_filtered{};
	uint32_t modulo;
	libdivide::divider<uint64_t> div;


public:
	explicit CKmerFilter(const std::string& kmcDbPath, uint32_t modulo, uint32_t kmer_len, uint64_t tot_uniqe_kmers, double fill_factor_filtered_kmers, bool verbose);

	bool Possible(kmer_type kmer) const
	{
		auto h = hash_mm(kmer);
		return h - modulo * (h / div) == 0;
	}

	bool Check(kmer_type kmer) const
	{
/*		switch (modulo)
		{
		case 1: break;
		case 2: if (hash_mm(kmer) % 2) return false; break;
		case 3: if (hash_mm(kmer) % 3) return false; break;
		case 4: if (hash_mm(kmer) % 4) return false; break;
		case 5: if (hash_mm(kmer) % 5) return false; break;
		case 6: if (hash_mm(kmer) % 6) return false; break;
		case 7: if (hash_mm(kmer) % 7) return false; break;
		case 8: if (hash_mm(kmer) % 8) return false; break;
		case 9: if (hash_mm(kmer) % 9) return false; break;
		case 10: if (hash_mm(kmer) % 10) return false; break;
		case 11: if (hash_mm(kmer) % 11) return false; break;
		case 12: if (hash_mm(kmer) % 12) return false; break;
		case 13: if (hash_mm(kmer) % 13) return false; break;
		case 14: if (hash_mm(kmer) % 14) return false; break;
		case 15: if (hash_mm(kmer) % 15) return false; break;
		case 16: if (hash_mm(kmer) % 16) return false; break;
		case 17: if (hash_mm(kmer) % 17) return false; break;
		case 18: if (hash_mm(kmer) % 18) return false; break;
		case 19: if (hash_mm(kmer) % 19) return false; break;
		case 20: if (hash_mm(kmer) % 20) return false; break;
		default:
			if (hash_mm(kmer) % modulo)	return false;
		}*/

/*		if (hash_mm(kmer) % modulo)
			return false;*/

		return kmers.check(kmer);
	}	
	
	uint64_t GetTotalKmers() const
	{
		return total_count_filtered;
	}

#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
	uint64_t GetCountsPerPrefix(uint32_t prefix) const
	{
		return kmers.GetCountPerPrefix(prefix);
	}
	void ReleaseCountsPerPrefix()
	{
		kmers.ReleaseCountsPerPrefix();
	}
#endif

	uint64_t GetMemoryUsage() const
	{
		return kmers.GetMemoryUsage();
	}

	void PrintMemoryUsage() const
	{
		kmers.PrintMemoryUsage();
	}

	uint32_t GetNHashTables() const
	{
		return kmers.GetNHashTables();
	}
};
