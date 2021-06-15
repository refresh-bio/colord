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
#include "count_kmers.h"
#include "libs/count_kmers/filtering_kmc.h"
#include "utils.h"
#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>

using namespace std;
CKmerCounter::CKmerCounter(uint32_t k, uint32_t ci, uint32_t cs, uint32_t n_threads, uint32_t modulo, const std::string& inputPath, const std::string& outPath, const std::string& tmpPath, bool is_fastq, bool verbose)
{
	std::cerr << "Counting k-mers.\n";
	//if (fileExists(outPath + ".kmc_pre") && fileExists(outPath + ".kmc_suf"))
	//	return;
#ifdef _WIN32
	const std::string kmc_path = "kmer_counter.exe"; 
#else
	const std::string kmc_path = "./kmc";
#endif
	//bool isFasta = !isFastq(inputPath);
	//std::string inType = "";
	//if (isFasta)
	//	inType = " -fa";
	//std::ostringstream stream;
	std::string kmc_stats_file = (std::filesystem::path(tmpPath) / "kmc-stats.json").string();

	//stream << kmc_path << " -k" << k << " -ci" << ci << " -cs" << cs << " -t" << n_threads << " -j" << kmc_stats_file <<" " << " -f" << modulo << inType << " \"" << inputPath << "\" \"" << outPath << "\" \"" << tmpPath << "\"";
	//auto command = stream.str();
	//if(verbose)
	//	std::cerr << "kmc command:\n" << command << "\n";
	//int stat = system(command.c_str());
	//if (stat != 0)
	//{
	//	std::cerr << "Error while counting k-mers. The command was: " << command << "\n";
	//	exit(1);
	//}


	CFilteringParams params;
	params.kmerLen = k;
	params.cutoffMin = ci;
	params.maxCount = cs;
	params.nThreads = n_threads;
	params.modulo = modulo;
	params.inputPath = inputPath;
	params.outputPath = outPath;
	params.tmpPath = tmpPath;
	params.statsFile = kmc_stats_file;
	params.is_fasta = !is_fastq;
	int stat = run_filtering_kmc(params);	
	if (stat != 0)
	{
		std::cerr << "Error while counting k-mers\n";
		exit(1);
	}

	auto in = inOpenOrDie(kmc_stats_file);
	if (!in)
	{
		std::cerr << "Error reading k-mer counter stats\n";
		exit(1);
	}
	std::string line;
	while (std::getline(in, line))
	{
		if (auto pos = line.find("#Unique_counted_k-mers"); pos != std::string::npos)
			std::istringstream(line) >> line >> n_unique_counted_kmers;
		else if (auto pos = line.find("#Total no. of k-mers"); pos != std::string::npos)
			std::istringstream(line) >> line >> line >> line >> line >> tot_kmers;
		else if (auto pos = line.find("#Total_reads"); pos != std::string::npos)
			std::istringstream(line) >> line >> n_reads;
	}
	in.close();
}
