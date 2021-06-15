#pragma once
#include "utils.h"
#include <vector>

struct CCompressElem
{
	uint32_t read_id; //probably not used, kept for debug purposes
	bool hasN;
	read_t read;
	std::vector<uint32_t> ref_reads;
	std::vector<std::vector<kmer_type>> common_kmers;
};

struct CCompressPack
{
	uint32_t id; //for priorities
	std::vector<CCompressElem> data;
};