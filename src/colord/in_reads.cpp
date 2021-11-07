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
#include "in_reads.h"
#include <iostream>

using namespace std;

read_t CInputReads::to_read_t(const std::string& str, bool& hasN)
{
	read_t res(str.size() + 1);
	hasN = false;
	for (size_t i = 0; i < str.length(); ++i)
	{
		int8_t code = SymbToBinMap[(uint8_t)str[i]];
		if (code == -1)
		{
			std::cerr << "Only ACGTN symbols supported inside a read\n";
			exit(1);
		}
		hasN |= code == 4;
		
		res[i] = code;
	}
	res[str.size()] = 255; //guard
	return res;
}

void CInputReads::addReadHeader()
{
	total_symb_header += currentLine.size();
	currentLine.erase(currentLine.begin()); //remove '>'
	if (!is_fastq)
	{		
		current_header_bytes += currentLine.size();
		headers.emplace_back(std::move(currentLine), qual_header_type::empty);
		if (current_header_bytes >= headers_pack_size)
		{
			current_header_bytes = 0;
			headers_queue.Push(std::move(headers));
		}
	}
	else
		last_read_header = std::move(currentLine);
}

void CInputReads::addRead()
{
	bool hasN;
	last_read = to_read_t(currentLine, hasN);
	total_bases += currentLine.length();
	current_reads_bytes += last_read.size();

	reads.emplace_back(hasN, last_read);

	stats.LogRead(read_len(reads.back().second));
	if (current_reads_bytes >= reads_pack_size)
	{
		current_reads_bytes = 0;
		reads_queue.Push(std::move(reads));
	}
}

void CInputReads::addQualHeader()
{
	qual_header_type type = qual_header_type::empty;
	total_symb_header += currentLine.size();
	currentLine.erase(currentLine.begin()); //remove '+'
	if (currentLine.size() > 0)
	{
		if (currentLine != last_read_header)
		{
			std::cerr << "Error: quality header not empty but different than read header\n";
			exit(1);
		}
		type = qual_header_type::eq_read_header;
	}
	current_header_bytes += last_read_header.size();

	headers.emplace_back(std::move(last_read_header), type);

	if (current_header_bytes >= headers_pack_size)
	{
		current_header_bytes = 0;
		headers_queue.Push(std::move(headers));
	}
}

void CInputReads::addQual()
{
	qual_t qual(currentLine.begin(), currentLine.end());

//	quals.emplace_back(std::move(last_read), std::move(currentLine));
	quals.emplace_back(std::move(last_read), std::move(qual));
	if(current_reads_bytes == 0) //if reads was just added then qual should be also, because the number of records should be the same for quals and for reads	
		quals_queue.Push(std::move(quals));
}

void CInputReads::porcessFastaOrMultiFasta(std::vector<uint8_t>& buff, uint64_t readed, uint32_t buf_size, gzFile gzfile)
{
	enum class FastaReadState { Header, EOLsAfterHeader, Read, EOLsAfterOrInsideRead };
	FastaReadState fastaReadState = FastaReadState::Header;
	while (readed)
	{
		for(uint32_t pos = 0; pos < readed; ++pos)
		{
			auto symb = buff[pos];
			bool is_eol = symb == '\n' || symb == '\r';
			switch (fastaReadState)
			{
				case FastaReadState::Header:
					if (is_eol)
					{
						currentLine.shrink_to_fit();
						addReadHeader();
						currentLine.clear();
						fastaReadState = FastaReadState::EOLsAfterHeader;
					}
					else
						currentLine.push_back(symb);
					break;
				case FastaReadState::EOLsAfterHeader:
					if (!is_eol)
					{
						currentLine.push_back(symb);
						fastaReadState = FastaReadState::Read;
					}
					break;
				case FastaReadState::Read:
					if (is_eol)
					{
						fastaReadState = FastaReadState::EOLsAfterOrInsideRead;
					}
					else
						currentLine.push_back(symb);
					break;
				case FastaReadState::EOLsAfterOrInsideRead:
					if (!is_eol)
					{
						if (symb == '>')
						{
							currentLine.shrink_to_fit();
							addRead();
							currentLine.clear();
							fastaReadState = FastaReadState::Header;
						}
						else
							fastaReadState = FastaReadState::Read;
						currentLine.push_back(symb);
					}
					break;
				default:
					break;
			}
		}
		readed = gzfread(buff.data(), 1, buf_size, gzfile);
		total_bytes += readed;
	}
	currentLine.shrink_to_fit();
	addRead();
	currentLine.clear();
}

void CInputReads::processFastq(std::vector<uint8_t>& buff, uint64_t readed, uint32_t buf_size, gzFile gzfile)
{
	enum class WhereInRead { read_header, read, qual_header, qual };
	WhereInRead whereInRead = WhereInRead::read_header;

	uint32_t record_lines = 4; //4 lines form fastq record

	while (readed)
	{
		uint32_t pos = 0;
		while (pos < readed)
		{
			if (buff[pos] == '\n' || buff[pos] == '\r') // EOL reached
			{
				if (currentLine.empty()) //we are skipping windows EOL
					++pos;
				else
				{
					currentLine.shrink_to_fit();
					switch (whereInRead)
					{
					case WhereInRead::read_header:
						addReadHeader();
						break;
					case WhereInRead::read:
						addRead();
						break;
					case WhereInRead::qual_header:
						addQualHeader();
						break;
					case WhereInRead::qual:
						addQual();
						break;
					default:
						break;
					}
					whereInRead = (WhereInRead)(((int)whereInRead + 1) % record_lines);
					currentLine.clear();
					++pos;
				}
			}
			else
				currentLine.push_back(buff[pos++]);
		}
		readed = gzfread(buff.data(), 1, buf_size, gzfile);
		total_bytes += readed;
	}
}


CInputReads::CInputReads(bool verbose, const std::string& path, CParallelQueue<read_pack_t>& reads_queue, CParallelQueue<qual_pack_t>& quals_queue, CParallelQueue<header_pack_t>& headers_queue) :
	stats(verbose),
	reads_queue(reads_queue),
	quals_queue(quals_queue),
	headers_queue(headers_queue)
{
	auto gzfile = gzopen(path.c_str(), "rb");
	if (!gzfile)
	{
		cerr << "Error: cannot open file: " << path << "\n";
		exit(1);
	}

	const uint32_t buf_size = 1ul << 25;
	std::vector<uint8_t> buff(buf_size);

	uint64_t readed = gzfread(buff.data(), 1, buf_size, gzfile);
	total_bytes += readed;
	if (!readed)
	{
		int code;
		auto errmsg = gzerror(gzfile, &code);
		if (code < 0)
		{
			std::cerr << "zblib error: " << errmsg << "\n";
			exit(1);
		}
		std::cerr << "Error: file " << path << " is empty\n";
		exit(1);
	}
	if (buff[0] != '@' && buff[0] != '>')
	{
		std::cerr << "Error: unknown file format\n";
		exit(1);
	}
	is_fastq = buff[0] == '@';

	//FASTA or multi-fasta
	if(!is_fastq)
		porcessFastaOrMultiFasta(buff, readed, buf_size, gzfile);
	else
		processFastq(buff, readed, buf_size, gzfile);

	int code;
	auto errmsg = gzerror(gzfile, &code);
	if (code < 0)
	{
		std::cerr << "zblib error: " << errmsg << "\n";
		exit(1);
	}

	if (!currentLine.empty())
	{
		std::cerr << "Error: something went wrong during input reading\n";
		exit(1);
	}
	gzclose(gzfile);

	
	if (reads.size())
		reads_queue.Push(std::move(reads));
	if (quals.size())
		quals_queue.Push(std::move(quals));
	if (headers.size())
		headers_queue.Push(std::move(headers));
	
	reads_queue.MarkCompleted();
	quals_queue.MarkCompleted();
	headers_queue.MarkCompleted();
}