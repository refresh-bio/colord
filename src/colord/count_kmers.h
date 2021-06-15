#pragma once
#include "defs.h"
#include "params.h"
#include <string>

class CKmerCounter
{
	uint32_t n_reads;
	uint64_t tot_kmers;
	uint64_t n_unique_counted_kmers;
public:
	explicit CKmerCounter(uint32_t k, uint32_t ci, uint32_t cs, uint32_t n_threads, uint32_t modulo, const std::string& inputPath, const std::string& outPath, const std::string& tmpPath, bool is_fastq, bool verbose);
	uint32_t GetNReads() const { return n_reads; }
	uint64_t GetTotKmers() const { return tot_kmers; }
	uint64_t GetNUniqueCounted() const { return n_unique_counted_kmers; }
};
