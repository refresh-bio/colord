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
#include "decompression_common.h"
#include "params.h"
#include "reference_genome.h"
#include "reference_reads.h"
#include "entr_read.h"
#include "entr_header.h"
#include "entr_qual.h"

void CDecmpressionModule::Run()
{
	int s_info = archive.GetStreamId("info");
	size_t meta;
	std::vector<uint8_t> _info;
	archive.GetPart(s_info, _info, meta);
	info.Deserialize(_info);

	if (!checkVersion(info))
	{
		std::ostringstream oss;
		oss << "incompatibile archive version";
		error_handler.LogError(oss.str());
		//std::cerr << "Error: incompatibile archive version\n";
		//exit(1);
	}

	int s_qual = archive.GetStreamId("qual");
	int s_meta = archive.GetStreamId("meta");

	is_fastq = s_qual != -1;

	std::vector<uint8_t> _params;

	archive.GetPart(s_meta, _params, meta);
	const uint8_t* ptr_params = _params.data();

	uint32_t tot_ref_reads{};
	LoadLittleEndian(ptr_params, tot_ref_reads);
	ptr_params += sizeof(tot_ref_reads);
	
	uint32_t maxCandidates{};
	LoadLittleEndian(ptr_params, maxCandidates);
	ptr_params += sizeof(maxCandidates);

	LoadLittleEndian(ptr_params, meta_data.compressionLevel);
	ptr_params += sizeof(meta_data.compressionLevel);

	uint8_t tmp;
	LoadLittleEndian(ptr_params, tmp);
	ptr_params += sizeof(tmp);
	meta_data.dataSource = static_cast<DataSource>(tmp);

	uint64_t approx_stream_size; //approximate size of DNA (or qual, cause the are the same) BEFORE compression (or after decompression)
	LoadLittleEndian(ptr_params, approx_stream_size);
	ptr_params += sizeof(approx_stream_size);

	if (is_fastq)
	{
		LoadLittleEndian(ptr_params, tmp);
		ptr_params += sizeof(tmp);
		meta_data.qualityComprMode = static_cast<QualityComprMode>(tmp);

		switch (meta_data.qualityComprMode)
		{
		case QualityComprMode::None:
			qualityRevThresholds.emplace_back();
			LoadLittleEndian(ptr_params, qualityRevThresholds.back());
			ptr_params += sizeof(qualityRevThresholds.back());
			break;
		case QualityComprMode::Original:
			break;
		case QualityComprMode::BinaryThreshold:
			for (uint32_t i = 0; i < 2; ++i)
			{
				qualityRevThresholds.emplace_back();
				LoadLittleEndian(ptr_params, qualityRevThresholds.back());
				ptr_params += sizeof(qualityRevThresholds.back());
			}
			break;
		case QualityComprMode::QuadThreshold:
			for (uint32_t i = 0; i < 4; ++i)
			{
				qualityRevThresholds.emplace_back();
				LoadLittleEndian(ptr_params, qualityRevThresholds.back());
				ptr_params += sizeof(qualityRevThresholds.back());
			}
			break;
		case QualityComprMode::QuinaryThreshold:
			for (uint32_t i = 0; i < 5; ++i)
			{
				qualityRevThresholds.emplace_back();
				LoadLittleEndian(ptr_params, qualityRevThresholds.back());
				ptr_params += sizeof(qualityRevThresholds.back());
			}
			break;
		case QualityComprMode::Average:
			break;
		case QualityComprMode::BinaryAverage:
			break;
		case QualityComprMode::QuadAverage:
			break;
		case QualityComprMode::QuinaryAverage:
			break;
		default:
			break;
		}

		if (verbose)
		{
			std::cerr << "qualityRevThresholds: ";
			for (auto v : qualityRevThresholds)
				std::cerr << v << " ";
			std::cerr << "\n";
			switch (meta_data.qualityComprMode)
			{
			case QualityComprMode::None:
				std::cerr << "qualityComprMode: None\n";
				break;

			case QualityComprMode::Original:
				std::cerr << "qualityComprMode: Original\n";
				break;
			case QualityComprMode::BinaryThreshold:
				std::cerr << "qualityComprMode: BinaryThreshold\n";
				break;
			case QualityComprMode::QuadThreshold:
				std::cerr << "qualityComprMode: QuadThreshold\n";
				break;
			case QualityComprMode::QuinaryThreshold:
				std::cerr << "qualityComprMode: QuinaryThreshold\n";
				break;
			case QualityComprMode::Average:
				std::cerr << "qualityComprMode: Average\n";
				break;
			case QualityComprMode::BinaryAverage:
				std::cerr << "qualityComprMode: BinaryAverage\n";
				break;
			case QualityComprMode::QuadAverage:
				std::cerr << "qualityComprMode: QuadAverage\n";
				break;
			case QualityComprMode::QuinaryAverage:
				std::cerr << "qualityComprMode: QuinaryAverage\n";
				break;
			default:
				break;
			}
		}
	}
	LoadLittleEndian(ptr_params, tmp);
	ptr_params += sizeof(tmp);
	meta_data.headerComprMode = static_cast<HeaderComprMode>(tmp);

	if (verbose)
	{
		switch (meta_data.headerComprMode)
		{
		case HeaderComprMode::Original:
			std::cerr << "headerComprMode: Original\n";
			break;
		case HeaderComprMode::Main:
			std::cerr << "headerComprMode: Main\n";
			break;
		case HeaderComprMode::None:
			std::cerr << "headerComprMode: None\n";
			break;
		default:
			break;
		}
	}

	LoadLittleEndian(ptr_params, tmp);
	ptr_params += sizeof(tmp);
	ReferenceReadsMode referenceReadsMode = static_cast<ReferenceReadsMode>(tmp);

	if (verbose)
	{
		switch (referenceReadsMode)
		{
		case ReferenceReadsMode::All:
			std::cerr << "reference reads mode: All\n";
			break;
		case ReferenceReadsMode::Sparse:
			std::cerr << "reference reads mode: Sparse\n";
			break;
		default:
			break;
		}
	}

	uint32_t sparseRange = 0;
	double sparseExponent = 0;
	if (referenceReadsMode == ReferenceReadsMode::Sparse)
	{
		LoadLittleEndian(ptr_params, sparseRange);
		ptr_params += sizeof(sparseRange);

		LoadLittleEndian(ptr_params, sparseExponent);
		ptr_params += sizeof(sparseExponent);
		if (verbose)
		{
			std::cerr << "sparse range: " << sparseRange << "\n";
			std::cerr << "sparse exponent: " << sparseExponent << "\n";
		}
	}

	bool ref_genome_available;
	LoadLittleEndian(ptr_params, ref_genome_available);
	ptr_params += sizeof(ref_genome_available);

	bool ref_genome_in_arch;

	uint32_t n_ref_genome_pseudo_reads{};

	if (ref_genome_available)
	{
		std::vector<uint8_t> ref_genome_checksum;
		//uint64_t mean_read_len{};
		//uint32_t kmerLen{};

		uint32_t ref_genome_overlap_size;
		uint32_t ref_genome_read_len;

		LoadLittleEndian(ptr_params, ref_genome_in_arch);
		ptr_params += sizeof(ref_genome_in_arch);

		LoadLittleEndian(ptr_params, ref_genome_read_len);
		ptr_params += sizeof(ref_genome_read_len);

		LoadLittleEndian(ptr_params, ref_genome_overlap_size);
		ptr_params += sizeof(ref_genome_overlap_size);

		LoadLittleEndian(ptr_params, n_ref_genome_pseudo_reads);
		ptr_params += sizeof(n_ref_genome_pseudo_reads);

		if (!ref_genome_in_arch)
		{
			ref_genome_checksum.resize(16);
			for (auto& c : ref_genome_checksum)
			{
				LoadLittleEndian(ptr_params, c);
				ptr_params += sizeof(c);
			}
		}

		if (ref_genome_in_arch)
			ref_genome = std::make_unique<CReferenceGenome>(archive, ref_genome_overlap_size, verbose);
		else
		{
			if (refGenomePath == "")
			{
				std::ostringstream oss;
				oss << "compressed file was created without -s switch, reference genome is required for decompression";
				error_handler.LogError(oss.str());
				//std::cerr << "Error: compressed file was created without -s switch, reference genome is required for decompression";
				//exit(1);
			}
			ref_genome = std::make_unique<CReferenceGenome>(refGenomePath, ref_genome_overlap_size, true, verbose);
			if (ref_genome->GetChecksum() != ref_genome_checksum)
			{
				std::cerr << "Error: different reference genome was used during compression. Decompression impossible.\n";
				exit(1);
			}
		}

		ref_genome->SetReadLen(ref_genome_read_len);
	}

	ref_reads = std::make_unique<CReferenceReads>(tot_ref_reads);

	if (ref_genome_available)
	{
		read_pack_t reads = ref_genome->GetPseudoReads();
		for (const auto& [tmp, read] : reads)
			ref_reads->Add(read);
	}

	ref_reads_accepter = std::make_unique<CRefReadsAccepter>(sparseRange, sparseExponent, n_ref_genome_pseudo_reads);

	if (verbose)
	{
		std::cerr << "maxCandidates: " << maxCandidates << "\n";
		std::cerr << "compressionLevel: " << meta_data.compressionLevel << "\n";
		std::cerr << "approx stream size: " << approx_stream_size << "\n";
	}

	read_decompr_queue = std::make_unique<CParallelQueue<decomp_read_pack_t>>(read_decompress_queue_size);
	read_decompr_queues.push_back(read_decompr_queue.get());

	
	header_decompr_queue = std::make_unique< CParallelQueue<header_pack_t>>(header_decompress_queue_size);

	if (is_fastq)//second queue for quality decompression
	{
		read_decompr_queue_for_qual = std::make_unique<CParallelQueue<decomp_read_pack_t>>(read_decompress_queue_size);
		read_decompr_queues.push_back(read_decompr_queue_for_qual.get());
		qual_decompr_queue = std::make_unique<CParallelQueue<decomp_qual_pack_t>>(qual_decompress_queue_size);
	}

	running_threads.emplace_back([&archive = archive, &read_decompr_queues = read_decompr_queues, verbose = verbose, maxCandidates, referenceReadsMode, &ref_reads = *ref_reads.get(), &ref_reads_accepter = *ref_reads_accepter.get(), compressionLevel = meta_data.compressionLevel, approx_stream_size, n_ref_genome_pseudo_reads, &info = info, hide_progress = hide_progress]{
		CEntropyDecomprReads decompr{ archive, read_decompr_queues, verbose, maxCandidates, compressionLevel, approx_stream_size, info.total_reads, referenceReadsMode, ref_reads, ref_reads_accepter, n_ref_genome_pseudo_reads, hide_progress };
		decompr.Decompress();
		});	

	running_threads.emplace_back([&archive = archive, &header_decompr_queue = header_decompr_queue, headerComprMode = meta_data.headerComprMode, compressionLevel = meta_data.compressionLevel, verbose = verbose]{
		CEntrDecomprHeaders decompr{ archive, *header_decompr_queue.get(), headerComprMode, compressionLevel, verbose };
		decompr.Decompress();
		});
	
	if (is_fastq)	
		running_threads.emplace_back([&archive = archive, &read_decompr_queues = read_decompr_queues, &qual_decompr_queue = *qual_decompr_queue.get(), verbose = verbose, qualityComprMode = meta_data.qualityComprMode, dataSource = meta_data.dataSource, &qualityFwdThresholds = qualityFwdThresholds, &qualityRevThresholds = qualityRevThresholds, compressionLevel = meta_data.compressionLevel, approx_stream_size]{
			CEntrDecomprQuals decompr{ archive, *read_decompr_queues[1], qual_decompr_queue, verbose, qualityComprMode, dataSource, qualityFwdThresholds, qualityRevThresholds, compressionLevel, approx_stream_size };
			decompr.Decompress();
			});

	if (is_fastq)
	{
		decompressed_stream_consumer.ConsumeFastq(*read_decompr_queue.get(), *header_decompr_queue.get(), *qual_decompr_queue.get());
		//CFastqWriter{ params.outputFilePath, read_decompr_queue, header_decompr_queue, *qual_decompr_queue.get() }.Write();			
	}
	else
	{
		decompressed_stream_consumer.ConsumeFasta(*read_decompr_queue.get(), *header_decompr_queue.get());
		//CFastaWriter{ params.outputFilePath, read_decompr_queue, header_decompr_queue }.Write();
	}
}