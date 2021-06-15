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
#include "kmer_filter.h"
#include "kmc_file.h"
#include "utils.h"


kmer_type CKmerFilter::hash_mm(kmer_type x) const
{
	x ^= x >> 33;
	x *= 0xff51afd7ed558ccduLL;
	x ^= x >> 33;
	x *= 0xc4ceb9fe1a85ec53uLL;
	x ^= x >> 33;
	return x;
}

std::string kmer_to_str(kmer_type kmer, uint32_t len)
{
	std::string res{};
	for (int32_t i = len - 1; i >= 0; --i)
	{
		auto s = (kmer >> (i * 2)) & 3;
		res.push_back("ACGT"[s]);
	}
	return res;
}

CKmerFilter::CKmerFilter(const std::string& kmcDbPath, uint32_t modulo, uint32_t kmer_len, uint64_t tot_uniqe_kmers, double fill_factor_filtered_kmers, bool verbose):
	kmers(kmer_len, tot_uniqe_kmers, fill_factor_filtered_kmers),
	modulo(modulo), 
	div(modulo)
{
	CKMCFile kmc_file;

	if (!kmc_file.OpenForListing(kmcDbPath))
	{
		std::cerr << "Cannot open kmc database\n";
		exit(1);
	}
	
	CKmerAPI kmer(kmc_file.KmerLength());
	uint32_t count;
	std::vector<uint64> kmer_v;
	
	CKMCFileInfo info;
	kmc_file.Info(info);	
	CPercentProgress progress(info.total_kmers);	
	
//	div = libdivide::divider(modulo);

	while (kmc_file.ReadNextKmer(kmer, count))
	{	
		progress.LogIter();
		kmer.to_long(kmer_v);
		uint64_t lowest = kmer_v.back();
//		if (hash_mm(lowest) % modulo == 0) // KMC was responsible for filtering :-)
		{
			total_count_filtered += count;
			//this->kmers.emplace(kmer_type(lowest), count);
			this->kmers.insert(lowest
#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
				, count
#endif
			);
		}
	}
	
	if(verbose)		
		std::cerr << "Total count filtered: " << total_count_filtered << "\n";	
}