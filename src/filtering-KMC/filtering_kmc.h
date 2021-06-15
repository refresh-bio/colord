#pragma once
#include <cinttypes>
#include <string>

struct CFilteringParams
{
	uint32_t kmerLen;
	uint32_t cutoffMin;
	uint32_t maxCount;
	uint32_t nThreads;
	uint32_t modulo;
	std::string inputPath;
	std::string outputPath;
	std::string tmpPath;	
	std::string statsFile;
	bool is_fasta;
};
int run_filtering_kmc(const CFilteringParams& params);
