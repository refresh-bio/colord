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
#include "decompression.h"
#include "decompression_common.h"
#include "utils.h"
#include "reference_genome.h"
#include "dna_coder.h"
#include "archive.h"
#include "parallel_queue.h"
#include "entr_read.h"
#include "entr_header.h"
#include "entr_qual.h"
#include <cstring>
#include <sstream>

// *****************************************************************
void store_read(ofstream& out, read_t& read)
{
	for (uint32_t i = 0; i < read.size() - 1; ++i)
		read[i] = "ACGTN"[read[i]];

	read.back() = '\n';

	out.write((char*)read.data(), read.size());
}

class CBufferedWrite
{
	std::ofstream out;
	const size_t write_buff_size = 1ull << 26;
	std::vector<uint8_t> write_buf;
	size_t pos{};
public:
	CBufferedWrite(const std::string& outputFilePath) :
		out(outOpenOrDie(outputFilePath, std::ios::binary)),
		write_buf(write_buff_size)
	{

	}
	void write(const char* data, size_t size)
	{		
		while (size)
		{
			if (pos + size <= write_buff_size)
			{
				memcpy(write_buf.data() + pos, data, size);
				pos += size;
				return;
			}
			else
			{
				size_t part_size = write_buff_size - pos;
				memcpy(write_buf.data() + pos, data, part_size);
				pos += part_size;
				out.write(reinterpret_cast<char*>(write_buf.data()), pos);
				data += part_size;				
				size -= part_size;				
				pos = 0;
			}
		}
	}
	~CBufferedWrite()
	{
		out.write(reinterpret_cast<char*>(write_buf.data()), pos);
	}
};

class CFastaWriter
{
	CBufferedWrite out;
	CParallelQueue<decomp_read_pack_t>& read_decompr_queue;
	CParallelQueue<header_pack_t>& header_decompr_queue;	

	decomp_read_pack_t read_pack;
	uint32_t read_pack_pos{};

	header_pack_t header_pack;
	uint32_t header_pack_pos{};

	
	template<typename T_ELEM, typename PACK_T>
	bool nextImpl(T_ELEM& elem, CParallelQueue<PACK_T>& queue, PACK_T& pack, uint32_t& pos)
	{
		if (pos == pack.size())
		{
			pos = 0;
			if (!queue.Pop(pack))
				return false;
		}
		elem = std::move(pack[pos++]);
		return true;
	}

	bool nextRead(read_t& read)
	{
		return nextImpl(read, read_decompr_queue, read_pack, read_pack_pos);
	}
	
	bool nextHeader(header_elem_t& header)
	{
		return nextImpl(header, header_decompr_queue, header_pack, header_pack_pos);
	}
public:
	CFastaWriter(const std::string& outputFilePath,
		CParallelQueue<decomp_read_pack_t>& read_decompr_queue,
		CParallelQueue<header_pack_t>& header_decompr_queue
		) :
		out(outputFilePath),
		read_decompr_queue(read_decompr_queue),
		header_decompr_queue(header_decompr_queue)		
	{

	}

	void Write()
	{
		read_t read;
		qual_t qual;
		header_elem_t header;

		bool hasNextRead = nextRead(read);		
		bool hasNextHeader = nextHeader(header);
		while (hasNextRead && hasNextHeader)
		{
			out.write(">", 1);
			out.write(header.first.c_str(), header.first.size());
			out.write("\n", 1);

			for (uint32_t i = 0; i < read.size() - 1; ++i)
				read[i] = "ACGTN"[read[i]];
			read.back() = '\n';
			out.write((char*)read.data(), read.size());
			
			hasNextRead = nextRead(read);			
			hasNextHeader = nextHeader(header);
		}
		if (hasNextRead + hasNextHeader != 0)
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}
	}
};

class CFastqWriter
{
	CBufferedWrite out;
	CParallelQueue<decomp_read_pack_t>& read_decompr_queue;
	CParallelQueue<header_pack_t>& header_decompr_queue;
	CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue;

	decomp_read_pack_t read_pack;
	uint32_t read_pack_pos{};

	header_pack_t header_pack;
	uint32_t header_pack_pos{};

	decomp_qual_pack_t qual_pack;
	uint32_t qual_pack_pos{};
	
	template<typename T_ELEM, typename PACK_T>
	bool nextImpl(T_ELEM& elem, CParallelQueue<PACK_T>& queue, PACK_T& pack, uint32_t& pos)
	{
		if (pos == pack.size())
		{
			pos = 0;
			if (!queue.Pop(pack))
				return false;
		}
		elem = std::move(pack[pos++]);
		return true;
	}

	bool nextRead(read_t& read)
	{
		return nextImpl(read, read_decompr_queue, read_pack, read_pack_pos);		
	}
	bool nextQual(qual_t& qual)
	{
		return nextImpl(qual, qual_decompr_queue, qual_pack, qual_pack_pos);
	}
	bool nextHeader(header_elem_t& header)
	{
		return nextImpl(header, header_decompr_queue, header_pack, header_pack_pos);
	}
	void write_via_buff(uint8_t* data, size_t size)
	{

	}
public:
	CFastqWriter(const std::string& outputFilePath,
		CParallelQueue<decomp_read_pack_t>& read_decompr_queue,
		CParallelQueue<header_pack_t>& header_decompr_queue,
		CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue) :
		out(outputFilePath),
		read_decompr_queue(read_decompr_queue),
		header_decompr_queue(header_decompr_queue),
		qual_decompr_queue(qual_decompr_queue)
		
	{

	}

	void Write()
	{
		read_t read;
		qual_t qual;
		header_elem_t header;

		bool hasNextRead = nextRead(read);
		bool hasNextQual = nextQual(qual);
		bool hasNextHeader = nextHeader(header);
		while (hasNextRead && hasNextQual && hasNextHeader)
		{
			out.write("@", 1);
			out.write(header.first.c_str(), header.first.size());
			out.write("\n", 1);

			for (uint32_t i = 0; i < read.size()-1; ++i)
				read[i] = read[i] = "ACGTN"[read[i]];
			read.back() = '\n';
			out.write((char*)read.data(), read.size());

			out.write("+", 1);
			if(header.second == qual_header_type::eq_read_header)
				out.write(header.first.c_str(), header.first.size());
			out.write("\n", 1);
			
			out.write((char*)qual.data(), qual.size());
			out.write("\n", 1);

			hasNextRead = nextRead(read);
			hasNextQual = nextQual(qual);
			hasNextHeader = nextHeader(header);
		}
		if (hasNextRead + hasNextQual + hasNextHeader != 0)
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}
	}
};


class CToFileDecompressedStreamConsumer : public IDecompressedStreamConsumer
{
	const std::string& outputFilePath;
public:
	CToFileDecompressedStreamConsumer(const std::string& outputFilePath)
		:
		outputFilePath(outputFilePath)
	{

	}

	virtual void ConsumeFasta(CParallelQueue<decomp_read_pack_t>& read_decompr_queue, CParallelQueue<header_pack_t>& header_decompr_queue) override
	{
		CFastaWriter{ outputFilePath, read_decompr_queue, header_decompr_queue }.Write();
	}
	virtual void ConsumeFastq(CParallelQueue<decomp_read_pack_t>& read_decompr_queue, CParallelQueue<header_pack_t>& header_decompr_queue, CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue) override
	{
		CFastqWriter{ outputFilePath, read_decompr_queue, header_decompr_queue, qual_decompr_queue }.Write();
	}
};

class CCerrLogger : public ILogger
{
public:
	virtual void Log(const std::string& msg) override
	{
		std::cerr << msg;
	}
};

class CExitErrorHandler : public IErrorHandler
{
public:
	virtual void LogError(const std::string& msg) override
	{
		std::cerr <<"Error: " << msg << "\n";
		exit(1);
	}
};

void runDecompression(const CDecompressorParams& params)
{
	CCerrLogger logger;
	CExitErrorHandler error_handler;
	CToFileDecompressedStreamConsumer consumer(params.outputFilePath);

	CDecmpressionModule decmpression_module(params.inputFilePath,
		params.refGenomePath,
		logger,
		params.verbose,
		error_handler,
		consumer
	);
	decmpression_module.Run();
	decmpression_module.WaitForThreads();
}
