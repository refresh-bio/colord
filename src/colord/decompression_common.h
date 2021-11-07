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
#include <string>
#include <sstream>
#include "parallel_queue.h"
#include "utils.h"
#include "archive.h"
#include "reference_reads.h"
#include "reference_genome.h"
#include "ref_reads_accepter.h"
#include "params.h"

class ILogger
{
public:
	virtual void Log(const std::string& msg) = 0;
	virtual ~ILogger() = default;
};

class IErrorHandler
{
public:
	virtual void LogError(const std::string& msg) = 0;
	virtual ~IErrorHandler() = default;
};

class IDecompressedStreamConsumer
{
public:
	virtual void ConsumeFasta(CParallelQueue<decomp_read_pack_t>& read_decompr_queue, CParallelQueue<header_pack_t>& header_decompr_queue) = 0;
	virtual void ConsumeFastq(CParallelQueue<decomp_read_pack_t>& read_decompr_queue, CParallelQueue<header_pack_t>& header_decompr_queue, CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue) = 0;
	virtual ~IDecompressedStreamConsumer() = default;
};

struct CMetaData
{
	int32_t compressionLevel{};
	DataSource dataSource{};
	QualityComprMode qualityComprMode{};
	HeaderComprMode headerComprMode;
};
class CDecmpressionModule
{
	const std::string& inputFilePath;
	const std::string& refGenomePath;
	ILogger& logger;
	bool verbose;
	IErrorHandler& error_handler;
	IDecompressedStreamConsumer& decompressed_stream_consumer;
	CArchive archive;
	bool hide_progress;

	CInfo info;
	CMetaData meta_data;
	bool is_fastq;
	
	vector<uint32_t> qualityFwdThresholds;
	vector<uint32_t> qualityRevThresholds;

	std::unique_ptr<CReferenceReads> ref_reads;
	std::unique_ptr<CReferenceGenome> ref_genome;
	std::unique_ptr<CRefReadsAccepter> ref_reads_accepter;

	std::unique_ptr<CParallelQueue<decomp_read_pack_t>> read_decompr_queue;

	std::vector<CParallelQueue<decomp_read_pack_t>*> read_decompr_queues;

	std::unique_ptr<CParallelQueue<decomp_read_pack_t>> read_decompr_queue_for_qual;
	std::unique_ptr<CParallelQueue<decomp_qual_pack_t>> qual_decompr_queue;
	std::unique_ptr<CParallelQueue<header_pack_t>> header_decompr_queue;

	std::vector<std::thread> running_threads;

	bool checkVersion(CInfo& info)
	{
		return info.version_major == version_major;
	}
public:
	CDecmpressionModule(const std::string& inputFilePath,
		const std::string& refGenomePath,
		ILogger& logger,
		bool verbose,
		IErrorHandler& error_handler,
		IDecompressedStreamConsumer& decompressed_stream_consumer,
		bool hide_progress = false)
		:
		inputFilePath(inputFilePath),
		refGenomePath(refGenomePath),
		logger(logger),
		verbose(verbose),
		error_handler(error_handler),
		decompressed_stream_consumer(decompressed_stream_consumer),
		archive(true),
		hide_progress(hide_progress)
	{
		//std::cerr << "Running decompression.\n";
		logger.Log("Running decompression.\n");
		if (!archive.Open(inputFilePath))
		{
			std::ostringstream oss;
			oss << "cannot open archive : " << inputFilePath;
			error_handler.LogError(oss.str());
			//std::cerr << "Error: cannot open archive: " << params.inputFilePath << "\n";
			//exit(1);
		}
	}

	const CInfo& GetInfo() const { return info; }
	const CMetaData& GetMetaData() const { return meta_data; }
	const vector<uint32_t>& GetQualityRevThresholds() const { return qualityRevThresholds; }
	bool GetIsFastq() const { return is_fastq; }
	void Run();

	void Cancel()
	{
		for (auto& q : read_decompr_queues)
			q->Cancel();

		if (qual_decompr_queue)
			qual_decompr_queue->Cancel();

		header_decompr_queue->Cancel();
		WaitForThreads();
	}

	void WaitForThreads()
	{
		for (auto& th : running_threads)
			th.join();
	}
	~CDecmpressionModule()
	{		
		archive.Close();
	}
};
