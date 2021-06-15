#pragma once

#include "utils.h"
#include "parallel_queue.h"
#include "dna_coder.h"
#include "archive.h"
#include "ref_reads_accepter.h"

class CEntrComprReads
{
	CParallelPriorityQueue<std::vector<es_t>>& compressed_queue;
	entropy_coder::CDNACoder dna_coder;
	CArchive& archive;
	int s_dna;	
	uint32_t tot_reads;
	uint32_t n_ref_genome_pseudo_reads;

public:
	CEntrComprReads(CParallelPriorityQueue<std::vector<es_t>>& compressed_queue,
		CReferenceReads& ref_reads, bool verbose,
		uint32_t maxCandidates,
		int32_t compression_level,
		uint64_t approx_input_stream_size,
		CArchive& archive,
		uint32_t tot_reads,
		uint32_t n_ref_genome_pseudo_reads
		) :
		compressed_queue(compressed_queue),
		dna_coder(ref_reads, verbose),	
		archive(archive),
		s_dna(archive.RegisterStream("dna")),
		tot_reads(tot_reads),
		n_ref_genome_pseudo_reads(n_ref_genome_pseudo_reads)
	{
		dna_coder.Init(true, maxCandidates, compression_level, approx_input_stream_size, n_ref_genome_pseudo_reads);
	}

	void Compress()
	{
		std::vector<es_t> encoded_redas_part;
		vector<uint8_t> v_output;
		std::cerr << "Running compression.\n";
		CPercentProgress progress(tot_reads);
		while (compressed_queue.Pop(encoded_redas_part))
		{	
			for (const auto& encoded_read : encoded_redas_part)
			{
				dna_coder.Encode(encoded_read);
			}

			dna_coder.Finish();

			dna_coder.GetOutput(v_output);
			dna_coder.Restart();
			
			size_t n_reads = encoded_redas_part.size();

			archive.AddPart(s_dna, v_output, n_reads);
			v_output.clear();			

			progress.LongNIters(encoded_redas_part.size());
		}
	}
};

class CEntropyDecomprReads
{
	CArchive& archive;
	int s_dna;
	std::vector<CParallelQueue<decomp_read_pack_t>*>& read_decompr_queues; //first for output, second is optional for quality
	entropy_coder::CDNACoder coder;
	
	uint64_t n_reads{};
	uint64_t maxCandidates;
	int32_t compression_level;
	uint64_t output_stream_size; // approx
	uint32_t total_reads;

	ReferenceReadsMode referenceReadsMode;

	CRefReadsAccepter& ref_reads_accepter;

	uint32_t n_ref_genome_pseudo_reads;

	void packToQueues(decomp_read_pack_t& pack)
	{		
		for (uint32_t i = 1; i < read_decompr_queues.size(); ++i)
		{
			auto copy = pack;
			read_decompr_queues[i]->Push(std::move(copy));
		}

		// Mask flags
		if(compression_level > 1)
			for (auto& r : pack)
				for (auto& c : r)
					c &= entropy_coder::base_mask_flags;

		read_decompr_queues[0]->Push(std::move(pack));
	}
public:
	CEntropyDecomprReads(CArchive& archive, std::vector<CParallelQueue<decomp_read_pack_t>*>& read_decompr_queues, bool verbose, uint32_t maxCandidates,
		int32_t compression_level,
		uint64_t output_stream_size,
		uint32_t total_reads,
		ReferenceReadsMode referenceReadsMode, CReferenceReads& ref_reads, CRefReadsAccepter& ref_reads_accepter,
		uint32_t n_ref_genome_pseudo_reads) :
		archive(archive),
		s_dna(archive.GetStreamId("dna")),
		read_decompr_queues(read_decompr_queues),
		coder(ref_reads, verbose),
		maxCandidates(maxCandidates),
		compression_level(compression_level),
		output_stream_size(output_stream_size),
		total_reads(total_reads),
		referenceReadsMode(referenceReadsMode),
		ref_reads_accepter(ref_reads_accepter),
		n_ref_genome_pseudo_reads(n_ref_genome_pseudo_reads)
	{
	}

	void Decompress()
	{
		vector<uint8_t> v_input;
		size_t meta;
		archive.GetPart(s_dna, v_input, meta);
		uint64_t local_n_reads = meta;

		n_reads += local_n_reads;

		coder.SetInput(v_input);
// !!! TODO		coder.Init(false, maxCandidates, compression_level, input_stream_size);
		coder.Init(false, static_cast<int>(maxCandidates), compression_level, 0, n_ref_genome_pseudo_reads);
		decomp_read_pack_t read_pack;

		//CPercentProgress progress(output_stream_size);
		CPercentProgress progress(total_reads);
		while (true)
		{
			read_t read;
			for (uint32_t i = 0; i < local_n_reads; ++i)
			{
				coder.Decode(read, ref_reads_accepter, referenceReadsMode == ReferenceReadsMode::All);
				//progress.LongNIters(read_len(read));
				progress.LogIter();
				read_pack.emplace_back(std::move(read));
			}

			packToQueues(read_pack);

			if (!archive.GetPart(s_dna, v_input, meta))
				break;
			local_n_reads = meta;

			n_reads += local_n_reads;

			coder.SetInput(v_input);
			coder.Restart();
		}

		//progress.ForceFinish();

		for (auto q : read_decompr_queues)
			q->MarkCompleted();
	}
};