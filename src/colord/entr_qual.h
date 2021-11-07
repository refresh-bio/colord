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

#include "utils.h"
#include "parallel_queue.h"
#include "archive.h"
#include "quality_coder.h"

class CEntrComprQuals
{
	CParallelQueue<qual_pack_t>& quals_queue;	
	CArchive& archive;
	int s_qual;
	entropy_coder::CQualityCoder quality_coder;

	size_t q_size = 0;

	CParallelPriorityQueue<std::vector<es_t>>& edit_script_for_qual_queue;
	
	qual_pack_t quals_pack;
	uint32_t quals_pack_pos{};
	std::vector<es_t> es_pack;
	uint32_t es_pack_pos{};

	DataSource dataSource;

	bool nextReadAndQual(qual_elem_t& read_and_qual)
	{
		if (quals_pack_pos == quals_pack.size())
		{
			quals_pack_pos = 0;
			if (!quals_queue.Pop(quals_pack))
				return false;
		}
		read_and_qual = std::move(quals_pack[quals_pack_pos++]);
		return true;
	}

	bool nextEditScript(es_t& es)
	{
		if (es_pack_pos == es_pack.size())
		{
			es_pack_pos = 0;
			if (!edit_script_for_qual_queue.Pop(es_pack))
				return false;
		}
		es = std::move(es_pack[es_pack_pos++]);
		return true;
	}

	void storeCurEncoded(std::vector<uint8_t>& v_output)
	{
		quality_coder.Finish();

		quality_coder.GetOutput(v_output);
		archive.AddPart(s_qual, v_output, 0);

		q_size += v_output.size();

		quality_coder.Restart();
	}
public:
	CEntrComprQuals(CParallelQueue<qual_pack_t>& quals_queue,
		CArchive& archive,
		QualityComprMode qualityComprMode,
		const std::vector<uint32_t>& qualityFwdThresholds,
		const std::vector<uint32_t>& qualityRevThresholds,
		bool verbose,
		int32_t compression_level,
		uint64_t approx_input_stream_size,
		CParallelPriorityQueue<std::vector<es_t>>& edit_script_for_qual_queue,
		DataSource dataSource) :
		quals_queue(quals_queue),		
		archive(archive),
		s_qual(archive.RegisterStream("qual")),
		quality_coder(verbose),
		edit_script_for_qual_queue(edit_script_for_qual_queue),
		dataSource(dataSource)
	{
		quality_coder.Init(true, qualityComprMode, dataSource, qualityFwdThresholds, qualityRevThresholds, compression_level, approx_input_stream_size);
	}
	
	void Compress()
	{		
		std::vector<uint8_t> v_output;

		qual_elem_t read_and_qual;
		es_t es;
		bool hasNextReadAndQual = nextReadAndQual(read_and_qual);
		bool hasNextEditScirpt = nextEditScript(es);

		assert(quals_pack_pos == es_pack_pos);

		while (hasNextReadAndQual && hasNextEditScirpt)
		{
			const auto& [read, qual] = read_and_qual;
			quality_coder.Encode(read, qual, es); 

			if(quals_pack_pos == quals_pack.size())
				storeCurEncoded(v_output);

			hasNextReadAndQual = nextReadAndQual(read_and_qual);
			hasNextEditScirpt = nextEditScript(es);

			assert(quals_pack_pos == es_pack_pos);
		}
		if (quals_pack_pos)
			storeCurEncoded(v_output);

		if (hasNextReadAndQual != hasNextEditScirpt)
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}

		//cout << "Quality size: " << q_size << endl;
	}
};

class CEntrDecomprQuals
{	
	CArchive& archive;
	int s_qual;
	CParallelQueue<decomp_read_pack_t>& read_decompr_queue;
	CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue;

	entropy_coder::CQualityCoder quality_coder;

	std::vector<uint8_t> v_archive;

	decomp_read_pack_t read_pack;
	uint32_t read_pack_pos{};
	decomp_qual_pack_t qual_pack;

	bool nextCompressedQual(std::vector<uint8_t>& read, std::vector<uint8_t>& qual, bool& cancel)
	{
		if (read_pack_pos == 0)
			return false;

		if (read_pack_pos == 1)
		{
			size_t meta;
			if (!archive.GetPart(s_qual, v_archive, meta))
				return false;

			if (!qual_pack.empty())
			{
				if (!qual_decompr_queue.PushOrCancel(std::move(qual_pack)))
					cancel = true;
				qual_pack.clear();
			}

			quality_coder.SetInput(v_archive);
			quality_coder.Restart();
		}

		qual.clear();

/*		string q;
		
		quality_coder.Decode(read, q);

		qual.assign(q.begin(), q.end());*/

		quality_coder.Decode(read, qual);

		return true;
	}

	bool nextRead(read_t& read)
	{
		if (read_pack_pos == read_pack.size())
		{
			read_pack_pos = 0;
			if (!read_decompr_queue.Pop(read_pack))
				return false;
		}
		read = std::move(read_pack[read_pack_pos++]);
		return true;
	}
public:
	CEntrDecomprQuals(CArchive& archive,
		CParallelQueue<decomp_read_pack_t>& read_decompr_queue,
		CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue,
		bool verbose,
		QualityComprMode qualityComprMode,
		DataSource dataSource,
		const std::vector<uint32_t>& qualityFwdThresholds,
		const std::vector<uint32_t>& qualityRevThresholds,
		int32_t compression_level,
		uint64_t approx_output_stream_size) :
		archive(archive),
		s_qual(archive.GetStreamId("qual")),
		read_decompr_queue(read_decompr_queue),
		qual_decompr_queue(qual_decompr_queue),
		quality_coder(verbose)
	{
		quality_coder.Init(false, qualityComprMode, dataSource, qualityFwdThresholds, qualityRevThresholds, compression_level, approx_output_stream_size);
	}

	void Decompress()
	{
		std::vector<uint8_t> qual;
		read_t read;
		bool cancel = false;
		bool hasNextRead = nextRead(read);
		bool hasNextQual = nextCompressedQual(read, qual, cancel);

		while (!cancel && hasNextRead && hasNextQual)
		{
			qual_pack.emplace_back(qual.begin(), qual.end());
//			for (auto c : qual)
//				qual_pack.back().push_back(c);

/*			if (qual_pack.size() == 1024)
			{
				qual_decompr_queue.Push(std::move(qual_pack));
				qual_pack.clear();
			}*/

			hasNextRead = nextRead(read);
			hasNextQual = nextCompressedQual(read, qual, cancel);
		}

		if (qual_pack.size())
			if (!qual_decompr_queue.PushOrCancel(std::move(qual_pack)))
				cancel = true;


		if (!cancel && hasNextQual != hasNextRead)
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}

		qual_decompr_queue.MarkCompleted();
	}
};