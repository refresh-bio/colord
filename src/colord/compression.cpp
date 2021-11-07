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
#include "defs.h"
#include "compression.h"
#include "count_kmers.h"
#include "in_reads.h"
#include "reads_sim_graph.h"
#include "encoder.h"
#include "parallel_queue.h"
#include "queues_data.h"
#include "reference_reads.h"
#include "entr_read.h"
#include "entr_qual.h"
#include "entr_header.h"
#include "archive.h"
#include "ref_reads_accepter.h"
#include "reference_genome.h"
#include <thread>
#include <memory>
#include <iostream>
#include <filesystem>
#include "timer.h"

static std::mutex mtx_cerr;

void adjustKmerAndAnchorLen(uint32_t& kmerLen, uint32_t& anchorLen, bool is_gzip_input, bool is_fastq, const std::string& filePath)
{
	if (kmerLen && anchorLen) 
		return;

	ifstream tmp(filePath, ios::binary);
	const auto begin = tmp.tellg();
	tmp.seekg(0, ios::end);
	uint64_t base_count = static_cast<uint64_t>(tmp.tellg() - begin);

	if (is_gzip_input)
		if (is_fastq)
			base_count = static_cast<uint64_t>(2.08 * base_count);
		else
			base_count = static_cast<uint64_t>(3.98 * base_count);
	else
		if (is_fastq)
			base_count = static_cast<uint64_t>(0.49 * base_count);
		else
			base_count = static_cast<uint64_t>(0.98 * base_count);


	if (base_count < 1'000'000'000ull)
	{
		kmerLen = 20;
		anchorLen = 16;
	}
	else if (base_count < 4'000'000'000ull)
	{
		kmerLen = 21;
		anchorLen = 18;
	}
	else if (base_count < 16'000'000'000ull)
	{
		kmerLen = 23;
		anchorLen = 21;
	}
	else if (base_count < 48'000'000'000ull)
	{
		kmerLen = 24;
		anchorLen = 22;
	}
	else if (base_count < 128'000'000'000ull)
	{
		kmerLen = 25;
		anchorLen = 22;
	}
	else
	{
		kmerLen = 26;
		anchorLen = 23;
	}
}

struct CTimeCollector
{
	bool is_fastq;

#ifdef MEASURE_THREADS_TIMES
	double inputReads{};
	double similarity_finder;
	std::vector<double> encoders;
	std::vector<double> encodersWaitOnQueueTime;
	double entropyCoder{};
	double entropyCoderQual{};
	double entropyCoderHeader{};
	double storeResult{};
	std::mutex mtxEncoders;
	std::vector<double> similarity_finder_internal;

	void RegisterNewEncoderTime(double time)
	{
		std::lock_guard<std::mutex> lck(mtxEncoders);
		encoders.push_back(time);
	}
#endif

private:
	void printVec(std::vector<double>& v, const std::string& name)
	{
		std::sort(v.begin(), v.end());
		double sum{};
		std::cerr << "\t" << name << " (sorted):";
		for (auto tw : v)
		{
			std::cerr << " " << tw;
			sum += tw;
		}
		std::cerr << "\n";
		std::cerr << "\tmean "<<name << ": " << sum / v.size() << "\n";
	}
public:
	void PrintReport()
	{
#ifdef MEASURE_THREADS_TIMES
	
#ifndef __APPLE__
		std::cerr << "Threads times summary:\n";
		std::cerr << "\treader: " << inputReads << "\n";
		std::cerr << "\tsimilarity finder: " << similarity_finder << "\n";

		printVec(similarity_finder_internal, "similarity finder internal");
		printVec(encoders, "encoders");
		printVec(encodersWaitOnQueueTime, "encoders wait on queue");

		std::cerr << "\tentr. compr.: " << entropyCoder << "\n";
		if(is_fastq)
			std::cerr << "\tentr. compr. qual.: " << entropyCoderQual << "\n";

		std::cerr << "\tentr. compr. header: " << entropyCoderHeader << "\n";
		std::cerr << "\tstore result: " << storeResult << "\n";

#endif // !__APPLE__

#endif
	}

	CTimeCollector(bool is_fastq) :is_fastq(is_fastq)
	{

	}
};

void PrintParams(const CCompressorParams& params, uint32_t kmerLen, uint32_t anchorLen)
{
	std::cerr << " * * * * * * * * * * * * Parameters * * * * * * * * * * * * \n";
	//std::cerr << "compression level: " << params.compressionLevel <<"\n";
	std::cerr << "\t" << "input file path: " << params.inputFilePath << "\n";
	std::cerr << "\t" << "output file path: " << params.outputFilePath << "\n";

	std::cerr << "\t" << "number of threads: " << params.nThreads << "\n";

	std::cerr << "\t" << "k-mer length: " << kmerLen << "\n";
	std::cerr << "\t" << "anchor length: " << anchorLen << "\n";
	
	std::cerr << "\t" << "data source type: " << dataSourceToString(params.dataSource) <<"\n";
	std::cerr << "\t" << "multipier for predicted cost of storing read part as edit script: " << params.editScriptCostMultiplier << "\n";
	std::cerr << "\t" << "filter modulo: " << params.filterHashModulo << "\n";
	std::cerr << "\t" << "header compression mode: " << headerComprModeToString(params.headerComprMode) << "\n";
				 
	std::cerr << "\t" << "max candidates: " << params.maxCandidates << "\n";
	std::cerr << "\t" << "min k-mer count: " << params.minKmerCount << "\n";
	std::cerr << "\t" << "max k-mer count: " << params.maxKmerCount << "\n";
	std::cerr << "\t" << "max matches multiplier: " << params.maxMatchesMultiplier << "\n";
	std::cerr << "\t" << "max recurence: " << params.maxRecurence << "\n";
	std::cerr << "\t" << "min anchors: " << params.minAnchors << "\n";
	std::cerr << "\t" << "min fraction of m-mers in encode: " << params.minFractionOfMmersInEncode << "\n";
	std::cerr << "\t" << "min fraction of m-mers in encode to always encode: " << params.minFractionOfMmersInEncodeToAlwaysEncode << "\n";	
	std::cerr << "\t" << "min part length to consider alternative reference read: " << params.minPartLenToConsiderAltRead << "\n";	
	std::cerr << "\t" << "compression priority: " << compressionPriorityToString(params.priority) << "\n";
	std::cerr << "\t" << "quality compression mode: " << qualityComprModeToString(params.qualityComprMode) << "\n";
	std::cerr << "\t" << "quality thresholds: " << vec_to_string(params.qualityFwdThresholds) << "\n";
	std::cerr << "\t" << "quality values: " << vec_to_string(params.qualityRevThresholds) << "\n";
	std::cerr << "\t" << "reference reads mode: " << referenceReadsModeToString(params.referenceReadsMode) << "\n";
	std::cerr << "\t" << "sparse mode exponent: " << params.sparseMode_exponent << "\n";
	std::cerr << "\t" << "sparse mode range: " << params.sparseMode_range_symbols << "\n";

	std::cerr << "\t" << "fill factor filtered k-mers: " << params.fillFactorFilteredKmers << "\n";
	std::cerr << "\t" << "fill factor k-mers to reads: " << params.fillFactorKmersToReads << "\n";

	if (params.refGenomePath != "")
	{
		std::cerr << "\t" << "reference genome path: " << params.refGenomePath << "\n";
		std::cerr << "\t" << "store reference genome in compressed file: " << (params.storeRefGenome ? "true" : "false") << "\n";
	}
}

uint32_t CalcApproxNRefReads(uint32_t tot_n_reads, uint32_t n_ref_genome_pseudo_reads, uint32_t sparseMode_range, double sparseMode_exponent)
{
	uint32_t n_ranges = (tot_n_reads + sparseMode_range - 1) / sparseMode_range;
	double res{};
	for (uint32_t rno = 0; rno < n_ranges; ++rno)
		res += pow(1.0 / (rno + 1), sparseMode_exponent);
	return sparseMode_range * res + n_ref_genome_pseudo_reads;
}

struct ApproxSizes
{
	uint64_t ref_reads_bytes;
	uint64_t kmers_to_reads_bytes;
};
ApproxSizes getApproxSizes(const CCompressorParams& params,
	const CKmerFilter& filtered_kmers,
	uint64_t tot_n_reads,
	uint32_t n_ref_genome_pseudo_reads,
	uint64_t mean_read_len,
	uint32_t sparseMode_range,
	double sparseMode_exponent,
	double fill_factor_kmers_to_reads)
{
	ApproxSizes res;
	auto approx_ref_reads = CalcApproxNRefReads(tot_n_reads, n_ref_genome_pseudo_reads, sparseMode_range, sparseMode_exponent);

	auto ref_reads_expected_size_bytes = approx_ref_reads * ((mean_read_len + 3) / 4 + 1 + sizeof(std::vector<read_t>)); //+ 1 for guard

	uint32_t n_hash_tables = filtered_kmers.GetNHashTables();

	uint64_t kmers_to_reads_expected_size_bytes = 0;
	
#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
	auto n_reads_factor = double(approx_ref_reads) / (tot_n_reads + n_ref_genome_pseudo_reads);
	for (uint32_t i = 0; i < n_hash_tables; ++i)
	{
		auto kmers_to_reads_expected_elems = filtered_kmers.GetCountsPerPrefix(i) * n_reads_factor;
		kmers_to_reads_expected_elems = round_to_pow_of_2(kmers_to_reads_expected_elems / fill_factor_kmers_to_reads);
		kmers_to_reads_expected_size_bytes += sizeof(kmers_to_reads_compacted_t::value_type) * kmers_to_reads_expected_elems;
	}
#else
	uint64_t tot_kmers_counts = filtered_kmers.GetTotalKmers();
	auto kmers_to_reads_expected_elems = tot_kmers_counts * double(approx_ref_reads) / (tot_n_reads + n_ref_genome_pseudo_reads);
	kmers_to_reads_expected_elems = round_to_pow_of_2(kmers_to_reads_expected_elems * (1.0 / fill_factor_kmers_to_reads));
	kmers_to_reads_expected_size_bytes = sizeof(kmers_to_reads_compacted_t::value_type) * kmers_to_reads_expected_elems;
#endif

	kmers_to_reads_expected_size_bytes += sizeof(std::vector<kmers_to_reads_compacted_t>) * n_hash_tables;

	res.ref_reads_bytes = ref_reads_expected_size_bytes;
	res.kmers_to_reads_bytes = kmers_to_reads_expected_size_bytes;
	return res;
}

uint64_t calcQueuesSize(bool is_fastq, uint32_t edit_script_for_qual_queue_size, uint32_t compressed_queue_size, uint32_t n_compression_threads, uint64_t mean_read_len, uint32_t maxCandidates, bool verbose)
{
	uint64_t reads_queue_bytes = (reads_queue_size // queue size
		+ 1 // one pack in reader
		+ 1 // one pack in consumer -> similarity finder
		) * reads_pack_size;

	uint64_t quals_queue_bytes = (quals_queue_size // queue size
		+ 1 // one pack in readers
		+ 1 // one pack in consumer -> qual entr
		) * reads_pack_size; //pack size for reads and quals is the same


	uint64_t headers_queue_bytes = (headers_queue_size //queue size
		+ 1 // one pack in reader
		+ 1 // one pack in consumer -> header entr
		) * headers_pack_size;


	uint64_t compress_queue_bytes = (compress_queue_size // queue size
		+ 1 // one pack in similarity finder
		+ n_compression_threads // n_compression_threads in CEncoder
		) * (sizeof(CCompressPack) + sizeof(CCompressElem) + sizeof(uint32_t) * maxCandidates + reads_pack_size);

	//I assume a tuple size is 1.25 bytes, and for single read we need read_len tuples. It is not very accurate because there are anchors, skips etc.
	double avg_compressed_pack_size = 1.25;

	uint64_t es_queue_bytes = (compressed_queue_size + // queue size
		+ n_compression_threads // n_compression_threads in CEncoder
		+ 1 // one in consumer -> read entr
		) * avg_compressed_pack_size * reads_pack_size;

	uint64_t es_queue_for_qual = 0;
	if(is_fastq)
		es_queue_for_qual = (edit_script_for_qual_queue_size + // queue size
			+n_compression_threads // n_compression_threads in CEncoder
			+ 1 // one in consumer -> read entr
			) * avg_compressed_pack_size * reads_pack_size;

	if (verbose)
	{
		std::cerr << "reads queue approx. size: " << reads_queue_bytes / 1024 << "KiB\n";
		std::cerr << "quals queue approx. size: " << quals_queue_bytes / 1024 << "KiB\n";
		std::cerr << "headers queue approx. size: " << headers_queue_bytes / 1024 << "KiB\n";
		std::cerr << "compress queue approx. size: " << compress_queue_bytes / 1024 << "KiB\n";
		std::cerr << "es queue approx. size: " << es_queue_bytes / 1024 << "KiB\n";
		std::cerr << "es queue for qual approx. size: " << es_queue_for_qual / 1024 << "KiB\n";
	}

	return reads_queue_bytes + quals_queue_bytes + headers_queue_bytes + compress_queue_bytes + es_queue_bytes + es_queue_for_qual;
}

void adjustMemorySparseMode(
	const CCompressorParams& params,
	const CKmerFilter& filtered_kmers,
	uint64_t tot_kmers_counts,
	uint64_t tot_n_reads,
	uint32_t n_ref_genome_pseudo_reads,
	uint64_t mean_read_len,
	uint64_t queuesApproxSize,
	
	uint32_t& sparseMode_range, double& sparseMode_exponent,
	double& fill_factor_kmers_to_reads)
{
	uint64_t filtered_kmers_bytes = filtered_kmers.GetMemoryUsage();
	
	ApproxSizes ref_reads_and_graph = getApproxSizes(params, filtered_kmers, tot_n_reads, n_ref_genome_pseudo_reads, mean_read_len, sparseMode_range, sparseMode_exponent, fill_factor_kmers_to_reads);

	uint64_t approx_total_memory = queuesApproxSize + filtered_kmers_bytes + ref_reads_and_graph.ref_reads_bytes + ref_reads_and_graph.kmers_to_reads_bytes;

	if (params.verbose)
	{
		std::cerr << "queues approx size: " << queuesApproxSize / 1024 / 1024 << "MiB\n";
		std::cerr << "filtered kmers size: " << filtered_kmers_bytes / 1024 / 1024 << "MiB\n";
		std::cerr << "ref reads expected size: " << ref_reads_and_graph.ref_reads_bytes / 1024 / 1024 << "MiB\n";
		std::cerr << "kmers to reads expected size: " << ref_reads_and_graph.kmers_to_reads_bytes / 1024 / 1024 << "MiB\n";
		std::cerr << "approx total memory: " << approx_total_memory / 1024 / 1024 << "MiB\n";
	}

}

void runCompression(const CCompressorParams& params, CInfo& info)
{
	info.version_major = version_major;
	info.version_minor = version_minor;
	info.version_patch = version_patch;
	Timer total_timer;
	total_timer.Start();

	CArchive archive(false);

	if (!archive.Open(params.outputFilePath))
	{
		std::cerr << "Error: cannot open archive: " << params.outputFilePath << "\n";
		exit(1);
	}
	
	int s_meta = archive.RegisterStream("meta");

	bool is_gzip_input = izGzipFile(params.inputFilePath);
	bool is_fastq = isFastq(params.inputFilePath);

	if (params.verbose)
	{
		if (is_gzip_input)
			std::cerr << "input is gzipped\n";
		else
			std::cerr << "input is not gzipped\n";
	}
	uint32_t input_reader_threads = is_gzip_input;
	uint32_t entropy_compr_threads = 1 + (is_fastq && params.qualityComprMode != QualityComprMode::None);

	uint32_t similarity_threads = 1; 

	int n_compression_threads = params.nThreads - input_reader_threads - entropy_compr_threads - similarity_threads;
	if (n_compression_threads < 1)
		n_compression_threads = 1;

	if (params.nThreads < 20)
	{
		if (is_fastq)
			n_compression_threads += 2;
		else
			n_compression_threads += 1;
	}


	uint32_t kmerLen = params.kmerLen;
	uint32_t anchorLen = params.anchorLen;

	adjustKmerAndAnchorLen(kmerLen, anchorLen, is_gzip_input, is_fastq, params.inputFilePath);
	if (params.verbose)	
		PrintParams(params, kmerLen, anchorLen);	

	auto tmp_dir_path = create_tmp_dir(std::filesystem::path(params.outputFilePath).remove_filename().string());

	std::string kmersDbPath = (std::filesystem::path(tmp_dir_path) / std::filesystem::path(params.inputFilePath).filename()).string() + "." + std::to_string(kmerLen) + "mers";
	
	bool ref_genome_available = params.refGenomePath != "";
	
	std::string kmcInputPath = params.inputFilePath;
	std::unique_ptr<CReferenceGenome> ref_genome;
	
	uint32_t ref_genome_overlap_size = (kmerLen - 1) * 10;
	uint32_t ref_genome_read_len{};
	if (ref_genome_available)
	{
		bool calc_checksum = !params.storeRefGenome;
		ref_genome = std::make_unique<CReferenceGenome>(params.refGenomePath, ref_genome_overlap_size, calc_checksum, params.verbose);
		std::string refGenomeKmcPath = "refGen.fa";
		if(is_fastq)
			refGenomeKmcPath = "refGen.fq";

		refGenomeKmcPath = (std::filesystem::path(tmp_dir_path) / std::filesystem::path(refGenomeKmcPath)).string();
		ref_genome->Store(refGenomeKmcPath, is_fastq);
		if(params.storeRefGenome)
			ref_genome->Store(archive);
		std::string newKmcInputPath = (std::filesystem::path(tmp_dir_path) / std::filesystem::path("kmc_file_list.txt")).string();
		std::ofstream newKmcInput(newKmcInputPath);
		if (!newKmcInput)
		{
			std::cerr << "Error: cannot open file: " << newKmcInputPath << "\n";
			exit(1);
		}
		newKmcInput << kmcInputPath << "\n";
		newKmcInput << refGenomeKmcPath << "\n";
		kmcInputPath = "@" + newKmcInputPath;
	}

	CKmerCounter kmer_counter(kmerLen, params.minKmerCount, params.maxKmerCount, params.nThreads, params.filterHashModulo, kmcInputPath, kmersDbPath, tmp_dir_path, is_fastq, params.verbose);
	auto tot_n_reads = kmer_counter.GetNReads();
	
	auto tot_kmers = kmer_counter.GetTotKmers();
	auto n_uniq_counted_kmers = kmer_counter.GetNUniqueCounted();
	std::cerr << "\n";
	if(params.verbose)
	{
		std::cerr << "tot k-mers: " << tot_kmers << "\n";
		std::cerr << "n uniq counted: " << n_uniq_counted_kmers << "\n";
	}
	uint64_t mean_read_len = static_cast<uint64_t>((double(tot_kmers * params.filterHashModulo) / tot_n_reads + kmerLen - 1));

	if (ref_genome_available) //correct stats if ref. genome is used
	{		
		mean_read_len = double(mean_read_len * tot_n_reads - ref_genome->GetTotSeqsLen()) / (tot_n_reads - ref_genome->GetTotNSeqs());
		tot_n_reads -= ref_genome->GetTotNSeqs();

		ref_genome_read_len = 20 * mean_read_len;

		ref_genome->SetReadLen(ref_genome_read_len);
	}

	info.total_reads = tot_n_reads;

	if (params.verbose)
		std::cerr << "approx. avg. read len: " << mean_read_len << "\n";

	Timer timer;

	std::cerr << "Filtering k-mers.\n";
	timer.Start();
	CKmerFilter filtered_kmers(kmersDbPath, params.filterHashModulo, kmerLen, n_uniq_counted_kmers, params.fillFactorFilteredKmers, params.verbose);
	std::error_code ec;
	std::filesystem::remove_all(tmp_dir_path, ec);
	if (ec)
	{
		std::cerr << "Warning: cannot remove tmp dir: " << tmp_dir_path << "\n";
	}
	//std::cerr << "Done.\n";
	if(params.verbose)
		timer.Log(std::cerr);

	uint32_t edit_script_for_qual_queue_size = 2 * n_compression_threads;
	uint32_t compressed_queue_size = 2 * n_compression_threads;

	CQueueMonitor queue_monitor(std::cerr, false, true);		// output_stream, single_line, report
//	CQueueMonitor queue_monitor(std::cerr, false, false);		// output_stream, single_line, report

	queue_monitor.register_queue(0, "reads_queue", reads_queue_size);
	queue_monitor.register_queue(1, "quals_queue", quals_queue_size);
	queue_monitor.register_queue(2, "headers_queue", headers_queue_size);
	queue_monitor.register_queue(3, "edit_script_for_qual_queue", edit_script_for_qual_queue_size);
	queue_monitor.register_queue(4, "compress_queue", compress_queue_size);
	queue_monitor.register_queue(5, "compressed_queue", compressed_queue_size);

	CParallelQueue<read_pack_t> reads_queue(reads_queue_size, 1, &queue_monitor, 0);
	CParallelQueue<qual_pack_t> quals_queue(quals_queue_size, 1, &queue_monitor, 1);
	CParallelQueue<header_pack_t> headers_queue(headers_queue_size, 1, &queue_monitor, 2);

//	CParallelQueue<std::vector<es_t>> edit_script_for_qual_queue(edit_script_for_qual_queue_size, 1, &queue_monitor, 3);
	CParallelPriorityQueue<std::vector<es_t>> edit_script_for_qual_queue(edit_script_for_qual_queue_size, n_compression_threads, &queue_monitor, 3);

	CParallelQueuePopWaiting<CCompressPack> compress_queue(compress_queue_size, &queue_monitor, 4);
	
	CParallelPriorityQueue<std::vector<es_t>> compressed_queue(compressed_queue_size, n_compression_threads, &queue_monitor, 5);

	uint32_t tot_ref_reads = tot_n_reads;
	
	uint32_t sparseMode_range = static_cast<uint32_t>((params.sparseMode_range_symbols * n_uniq_counted_kmers * params.filterHashModulo) / mean_read_len);
	if (sparseMode_range == 0)
		sparseMode_range = 1;
	double sparseMode_exponent = params.sparseMode_exponent;

	double fill_factor_kmers_to_reads = params.fillFactorKmersToReads;

	uint32_t n_ref_genome_pseudo_reads = ref_genome_available ? ref_genome->GetNPseudoReads() : 0;

	if (params.verbose && ref_genome_available)
		std::cerr << "# ref genome pseudo reads: " << n_ref_genome_pseudo_reads << "\n";

	tot_ref_reads += n_ref_genome_pseudo_reads;

	if (params.referenceReadsMode == ReferenceReadsMode::Sparse)
	{
		auto queuesApproxSize = calcQueuesSize(is_fastq, edit_script_for_qual_queue_size, compressed_queue_size, n_compression_threads,
			mean_read_len, params.maxCandidates, params.verbose);
		adjustMemorySparseMode(params, filtered_kmers, filtered_kmers.GetTotalKmers(), tot_n_reads, n_ref_genome_pseudo_reads, mean_read_len, queuesApproxSize,
			sparseMode_range, sparseMode_exponent, fill_factor_kmers_to_reads);
#ifdef ESTIMATE_MEMORY_WITH_COUNTS_PER_PREFIX
		filtered_kmers.ReleaseCountsPerPrefix();
#endif
		if (params.verbose)
		{
			std::cerr << "adjusted memory related settings:\n";
			std::cerr << "\tfill_factor_kmers_to_reads: " << fill_factor_kmers_to_reads << "\n";
			std::cerr << "\tsparseMode_range: " << sparseMode_range << "\n";
			std::cerr << "\tsparseMode_exponent: " << sparseMode_exponent << "\n";

		}
	}

	CRefReadsAccepter ref_reads_accepter(sparseMode_range, sparseMode_exponent, n_ref_genome_pseudo_reads);
	if (params.referenceReadsMode == ReferenceReadsMode::Sparse)
	{
		if(params.verbose)
			std::cerr << "sparse mode range in reads: " << sparseMode_range << "\n";
		tot_ref_reads = ref_reads_accepter.GetNAccepted(tot_n_reads);
	}

	CReferenceReads reference_reads(tot_ref_reads);

	CTimeCollector tc(is_fastq);

	uint64_t total_symb_header;
	std::thread reader([&params, &reads_queue, &quals_queue, &headers_queue, &tc, &info, &total_symb_header]{
#ifdef MEASURE_THREADS_TIMES
		CThreadWatch tw;
		tw.startTimer();
#endif

		CInputReads reads(params.verbose, params.inputFilePath, reads_queue, quals_queue, headers_queue);
		reads.GetStats(info.total_bytes, info.total_bases, total_symb_header);

#ifdef MEASURE_THREADS_TIMES
		tw.stopTimer();
		tc.inputReads = tw.getElapsedTime();
#endif

	});
	

	std::thread similarity_finder([&reads_queue, &compress_queue, &reference_reads, &ref_genome, &filtered_kmers, &params, kmerLen,
		&ref_reads_accepter, tot_n_reads, tot_ref_reads, &tc, n_compression_threads, &fill_factor_kmers_to_reads]{
#ifdef MEASURE_THREADS_TIMES
		CThreadWatch tw;
		tw.startTimer();
#endif

		//replace with nothing for debug purposes :)
		//read_pack_t p;
		//while (reads_queue.Pop(p))
		//	;
		//compress_queue.MarkCompleted();

		CReadsSimilarityGraph sim_graph(reads_queue, compress_queue, reference_reads, ref_genome.get(), filtered_kmers, 
			kmerLen, params.maxCandidates, params.maxKmerCount, params.referenceReadsMode, ref_reads_accepter,
			(double)tot_ref_reads/tot_n_reads, n_compression_threads, params.dataSource, fill_factor_kmers_to_reads, params.verbose);

#ifdef MEASURE_THREADS_TIMES
		tw.stopTimer();
		tc.similarity_finder = tw.getElapsedTime();
		tc.similarity_finder_internal = std::move(sim_graph.GetTwTimes());
#endif
	});
	
	std::vector<thread> encoders;
#ifdef MEASURE_THREADS_TIMES
	tc.encodersWaitOnQueueTime.resize(n_compression_threads);
#endif
	for (uint32_t i = 0; i < static_cast<uint32_t>(n_compression_threads); ++i)
		encoders.emplace_back([&compress_queue, &reference_reads, &compressed_queue, &edit_script_for_qual_queue, &params, anchorLen, &tc, i, is_fastq, kmerLen] {
#ifdef MEASURE_THREADS_TIMES
			CThreadWatch tw;
			tw.startTimer();
#endif

			//replace with nothing for debug purposes :)
			//CCompressPack p;
			//while (compress_queue.Pop(p))
			//	;
			//compressed_queue.MarkCompleted();

			CEncoder encoder(params.verbose, compress_queue, reference_reads, compressed_queue, edit_script_for_qual_queue, anchorLen,
			params.minFractionOfMmersInEncodeToAlwaysEncode, 
			params.minFractionOfMmersInEncode,
			params.maxMatchesMultiplier,
			params.editScriptCostMultiplier,
			params.minPartLenToConsiderAltRead,
			params.maxRecurence,			
			params.minAnchors,
			is_fastq,
			params.filterHashModulo,
			kmerLen,
			params.dataSource);
			
			encoder.Encode();

#ifdef MEASURE_THREADS_TIMES
			tc.encodersWaitOnQueueTime[i] = encoder.GetWaitOnQueueTime();
			tw.stopTimer();
			tc.RegisterNewEncoderTime(tw.getElapsedTime());
#endif
	});



	std::thread entropy_compressor([&reference_reads, &params, &compressed_queue, &archive, tot_n_reads, mean_read_len, n_ref_genome_pseudo_reads, &tc] {
#ifdef MEASURE_THREADS_TIMES
		CThreadWatch tw;
		tw.startTimer();
#endif

		//replace with nothing for debug purposes :)
		//std::vector<es_t> p;
		//while (compressed_queue.Pop(p))
		//	;

		CEntrComprReads compr{ compressed_queue, reference_reads, params.verbose, params.maxCandidates,
			params.compressionLevel, tot_n_reads * mean_read_len, archive, tot_n_reads, n_ref_genome_pseudo_reads };
		compr.Compress();

#ifdef MEASURE_THREADS_TIMES
		tw.stopTimer();
		tc.entropyCoder = tw.getElapsedTime();
#endif
	});

	std::thread entr_compr_qual;
	
	if (is_fastq)
	{		
		entr_compr_qual = std::thread([&quals_queue, &archive, &params, &edit_script_for_qual_queue, tot_n_reads, mean_read_len, &tc] {
#ifdef MEASURE_THREADS_TIMES
			CThreadWatch tw;
			tw.startTimer();
#endif

			CEntrComprQuals compr{ quals_queue, archive, params.qualityComprMode, params.qualityFwdThresholds, params.qualityRevThresholds, params.verbose, params.compressionLevel, tot_n_reads * mean_read_len, edit_script_for_qual_queue, params.dataSource };
			compr.Compress();

#ifdef MEASURE_THREADS_TIMES
			tw.stopTimer();
			tc.entropyCoderQual = tw.getElapsedTime();
#endif
		});
	}

	
	std::thread entr_compr_header([&headers_queue, &archive, &params, &tc] {
#ifdef MEASURE_THREADS_TIMES
		CThreadWatch tw;
		tw.startTimer();
#endif


		//replace with nothing for debug purposes :)
		//header_pack_t p;
		//while (headers_queue.Pop(p))
		//	;
		CEntrComprHeaders compr{ headers_queue, archive, params.headerComprMode, params.compressionLevel, params.verbose };
		compr.Compress();		

#ifdef MEASURE_THREADS_TIMES
		tw.stopTimer();
		tc.entropyCoderHeader = tw.getElapsedTime();
#endif
	});

#ifdef MEASURE_THREADS_TIMES
	CThreadWatch tw;
	tw.startTimer();
#endif

	if (is_fastq)
		entr_compr_qual.join();
	entr_compr_header.join();
	similarity_finder.join();
	reader.join();
	for (auto& th : encoders)
		th.join();
	entropy_compressor.join();

	std::vector<uint8_t> _config;
	StoreLittleEndian(_config, tot_ref_reads);
	StoreLittleEndian(_config, params.maxCandidates);
	StoreLittleEndian(_config, params.compressionLevel);
	StoreLittleEndian(_config, static_cast<uint8_t>(params.dataSource));

	StoreLittleEndian(_config, tot_n_reads * mean_read_len);
	if(is_fastq)
	{
		_config.push_back(static_cast<uint8_t>(params.qualityComprMode));

		switch (params.qualityComprMode)
		{
		case QualityComprMode::None:
			assert(params.qualityRevThresholds.size() == 1);
			StoreLittleEndian(_config, params.qualityRevThresholds[0]);
			break;
		case QualityComprMode::Original:
			break;
		case QualityComprMode::BinaryThreshold:
			assert(params.qualityRevThresholds.size() == 2);
			StoreLittleEndian(_config, params.qualityRevThresholds[0]);
			StoreLittleEndian(_config, params.qualityRevThresholds[1]);
			break;
		case QualityComprMode::QuadThreshold:
			assert(params.qualityRevThresholds.size() == 4);
			StoreLittleEndian(_config, params.qualityRevThresholds[0]);
			StoreLittleEndian(_config, params.qualityRevThresholds[1]);
			StoreLittleEndian(_config, params.qualityRevThresholds[2]);
			StoreLittleEndian(_config, params.qualityRevThresholds[3]);
			break;
		case QualityComprMode::QuinaryThreshold:
			assert(params.qualityRevThresholds.size() == 5);
			StoreLittleEndian(_config, params.qualityRevThresholds[0]);
			StoreLittleEndian(_config, params.qualityRevThresholds[1]);
			StoreLittleEndian(_config, params.qualityRevThresholds[2]);
			StoreLittleEndian(_config, params.qualityRevThresholds[3]);
			StoreLittleEndian(_config, params.qualityRevThresholds[4]);
			break;
		case QualityComprMode::Average:
			break;
		case QualityComprMode::BinaryAverage:
			break;
		case QualityComprMode::QuadAverage:
			break;
		default:
			break;
		}
	}
	_config.push_back(static_cast<uint8_t>(params.headerComprMode));
	_config.push_back(static_cast<uint8_t>(params.referenceReadsMode));

	if (params.referenceReadsMode == ReferenceReadsMode::Sparse)
	{
		StoreLittleEndian(_config, sparseMode_range);
		StoreLittleEndian(_config, sparseMode_exponent);
	}

	_config.push_back(static_cast<uint8_t>(ref_genome_available));
	if (ref_genome_available)
	{
		_config.push_back(static_cast<uint8_t>(params.storeRefGenome));
		StoreLittleEndian(_config, ref_genome_read_len);
		StoreLittleEndian(_config, ref_genome_overlap_size);
		StoreLittleEndian(_config, n_ref_genome_pseudo_reads);
		
		if (!params.storeRefGenome)
		{
			auto checksum = ref_genome->GetChecksum();
			for (auto c : checksum)
				StoreLittleEndian(_config, c);
		}
	}

	archive.AddPart(s_meta, _config, 0);

	int s_info = archive.RegisterStream("info");
	std::vector<uint8_t> _info = info.Serialize();
	archive.AddPart(s_info, _info, 0ull);

	archive.Close();

#ifdef MEASURE_THREADS_TIMES
	tw.stopTimer();
	tc.storeResult = tw.getElapsedTime();

	if (params.verbose) 
	{
		tc.PrintReport();
	}
#endif
	if (params.verbose)
	{
		cerr << "Input stream sizes (before compression):\n";
		cerr << "\tDNA size     : " << info.total_bases << "\n";
		cerr << "\tHeader size  : " << total_symb_header << "\n";
	}
	cerr << "DNA size        : " << archive.GetStreamPackedSize(archive.GetStreamId("dna")) << endl;
	cerr << "Quality size    : " << archive.GetStreamPackedSize(archive.GetStreamId("qual")) << endl;
	cerr << "Header size     : " << archive.GetStreamPackedSize(archive.GetStreamId("header")) << endl;
	cerr << "Meta size       : " << archive.GetStreamPackedSize(archive.GetStreamId("meta")) << endl;
	cerr << "Info size       : " << archive.GetStreamPackedSize(archive.GetStreamId("info")) << endl;
	if(params.storeRefGenome)
		cerr << "Ref genome size : " << archive.GetStreamPackedSize(archive.GetStreamId("ref-genome")) << endl;

	if(params.verbose)
	{
		reference_reads.PrintMemoryUsage();
	}

	cerr << "Total time      : " << total_timer.GetTimeInSec() << "s\n";
}
