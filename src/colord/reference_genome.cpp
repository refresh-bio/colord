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
#include "reference_genome.h"
#include "zlib.h"
#include "archive.h"
#include "dna_coder.h"
#include "utils.h"
#include "md5_wrapper.h"
#include <cassert>
#include <iostream>
#include <fstream>

void CReferenceGenome::packSeq(std::vector<uint8_t>& seq)
{
	//seq.shrink_to_fit();
		//return;
	size_t out_len = (seq.size() + 3) / 4 + 1; //one for number of symbols in last byte
	std::vector<uint8_t> packed(out_len);

	auto full_bytes = seq.size() / 4; //number of fully filled bytes in output
	uint8_t symbols_in_last_byte = seq.size() % 4;

	uint32_t in_pos = 0;
	for (uint32_t i = 0; i < full_bytes; ++i)
	{
		packed[i] = seq[in_pos++] << 6;
		packed[i] += seq[in_pos++] << 4;
		packed[i] += seq[in_pos++] << 2;
		packed[i] += seq[in_pos++];
	}

	switch (symbols_in_last_byte)
	{
	case 3:
		packed[out_len - 2] += seq[in_pos++] << 6;
		packed[out_len - 2] += seq[in_pos++] << 4;
		packed[out_len - 2] += seq[in_pos++] << 2;
		break;
	case 2:
		packed[out_len - 2] += seq[in_pos++] << 6;
		packed[out_len - 2] += seq[in_pos++] << 4;
		break;
	case 1:
		packed[out_len - 2] += seq[in_pos++] << 6;
		break;
	}
	assert(in_pos == seq.size());
	packed[out_len - 1] = symbols_in_last_byte;

	seq = std::move(packed);
}

std::vector<uint8_t> CReferenceGenome::unpackSeq(const std::vector<uint8_t>& seq)
{
	//return seq;
	auto symbols_in_last_byte = seq.back();
	auto full_bytes = seq.size() - 2 + !symbols_in_last_byte;
	auto out_len = full_bytes * 4 + symbols_in_last_byte;

	std::vector<uint8_t> res(out_len);
	uint32_t out_pos{};
	for (uint32_t i = 0; i < full_bytes; ++i)
	{
		res[out_pos++] = seq[i] >> 6;
		res[out_pos++] = (seq[i] >> 4) & 3;
		res[out_pos++] = (seq[i] >> 2) & 3;
		res[out_pos++] = seq[i] & 3;
	}
	switch (symbols_in_last_byte)
	{
	case 3:
		res[out_pos++] = seq[full_bytes] >> 6;
		res[out_pos++] = (seq[full_bytes] >> 4) & 3;
		res[out_pos++] = (seq[full_bytes] >> 2) & 3;
		break;
	case 2:
		res[out_pos++] = seq[full_bytes] >> 6;
		res[out_pos++] = (seq[full_bytes] >> 4) & 3;
		break;
	case 1:
		res[out_pos++] = seq[full_bytes] >> 6;
		break;
	default:
		break;
	}

	return res;
}

CReferenceGenome::CReferenceGenome(const std::string& path, uint32_t overlap_size, bool calc_checksum, bool verbose):
	overlap_size(overlap_size),
	verbose(verbose)
{
	if (verbose)
		std::cerr << "Reading reference genome...";
	Timer timer;
	timer.Start();

	auto gzfile = gzopen(path.c_str(), "rb");
	if (!gzfile)
	{
		std::cerr << "Error: cannot open file: " << path << "\n";
		exit(1);
	}

	const uint32_t buf_size = 1ul << 25;
	std::vector<uint8_t> buff(buf_size);

	uint64_t readed = gzfread(buff.data(), 1, buf_size, gzfile);

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

	if (buff[0] != '>')
	{
		std::cerr << "Error: wrong reference genome file format, multi fasta expected\n";
		exit(1);
	}

	enum class states { header, seq, skipping_eol_header, skipping_eol_seq };
	states state = states::header;
	newSeq();
	while (readed)
	{
		uint32_t pos = 0;
		while (pos < readed)
		{
			char c = buff[pos++];
			switch (state)
			{
			case states::header:
				if (c == '\n' || c == '\r')
					state = states::skipping_eol_header;
				break;
			case states::seq:
				if (c == '\n' || c == '\r')				
					state = states::skipping_eol_seq;				
				else
					addSymb(c);
				break;
			case states::skipping_eol_seq:
				if (c == '\n' || c == '\r')
					;
				else if (c == '>')
				{
					state = states::header;
					seqEnd();
					newSeq();
				}
				else
				{
					state = states::seq;
					addSymb(c);
				}
				break;
			case states::skipping_eol_header:
				if (c == '\n' || c == '\r')
					;
				else
				{
					state = states::seq;
					addSymb(c);
				}
			}
		
		}

		readed = gzfread(buff.data(), 1, buf_size, gzfile);
	}
	seqEnd();

	if (verbose)
	{
		std::cerr << "Done.\n";
		std::cerr << "Time: " << timer.GetTimeInSec() << "\n";
	}
	timer.Start();
	
	if (calc_checksum)
	{
		if (verbose)
			std::cerr << "Calculating checksum...";

		CMD5 md5;
		for (const auto& seq : sequences)
			md5.Update(seq.data(), seq.size());
		checksum = md5.Get();

		if (verbose)
		{
			std::cerr << "Done.\n";
			std::cerr << "Time: " << timer.GetTimeInSec() << "\n";
		}
	}
	if (verbose)
	{		
		std::cerr << "total sequences in reference genome file: " << sequences.size() << "\n";
		if (calc_checksum)
		{
			std::cerr << "md5: ";
			for (auto c : checksum)
				std::cerr << std::hex << (int)c;
			std::cerr << std::dec << "\n";
		}

	}
}

CReferenceGenome::CReferenceGenome(CArchive& archive, uint32_t overlap_size, bool verbose)
	:
	overlap_size(overlap_size),
	verbose(verbose)
{
	Timer timer;
	timer.Start();
	if (verbose)
		std::cerr << "Reading reference genome from archive...";
	int s_ref_genome = archive.GetStreamId("ref-genome");
	
	//commented part is without entropy compression, just 2bits/symbol
	//std::vector<uint8_t> v_input;
	//size_t meta;
	//while (archive.GetPart(s_ref_genome, v_input, meta))
	//	sequences.emplace_back(std::move(v_input));

	vector<uint8_t> v_input;
	size_t n_seqs;
	archive.GetPart(s_ref_genome, v_input, n_seqs);

	CRefReadsAccepter fake_ref_reads_accepter(1, 100000.0, 0);
	CReferenceReads fake_ref_reads(fake_ref_reads_accepter.GetNAccepted(n_seqs));
	entropy_coder::CDNACoder coder(fake_ref_reads, verbose);

	coder.SetInput(v_input);

	coder.Init(false, 0, 9, 0, 0);

	read_t read;
	for (uint32_t i = 0; i < n_seqs; ++i)
	{			
		coder.Decode(read, fake_ref_reads_accepter, false);		
			
		read.pop_back();
		sequences.emplace_back(std::move(read));
		packSeq(sequences.back());
	}
	if (verbose)
	{
		std::cerr << "Done.\n";
		std::cerr << "Time: " << timer.GetTimeInSec() << "\n";
	}

}

void CReferenceGenome::Store(const std::string& path, bool as_fastq)
{
	Timer timer;
	timer.Start();
	if (verbose)
		std::cerr << "Store reference sequences as reads in file " << path << "...";
	std::ofstream out(path);
	if (!out)
	{
		std::cerr << "Error: cannot open file: " << path << "\n";
		exit(1);
	}

	std::string fake_qual;
	for (const auto& _seq : sequences)
	{
		auto seq = unpackSeq(_seq);
		for (auto& c : seq)
			c = "ACGT"[c];
		auto to_store = seq.size();
		if (as_fastq)
			out.write("@\n", 2);
		else
			out.write(">\n", 2);

		out.write(reinterpret_cast<const char*>(seq.data()), to_store);
		out.write("\n", 1);
		if (as_fastq)
		{
			out.write("+\n", 2);
			if (fake_qual.size() < to_store)
				fake_qual.assign(to_store, 'I');
			out.write(reinterpret_cast<const char*>(fake_qual.data()), to_store);
			out.write("\n", 1);
		}
	}

	if (verbose)
	{
		std::cerr << "Done.\n";
		std::cerr << "Time: " << timer.GetTimeInSec() << "\n";
	}
}

void CReferenceGenome::Store(CArchive& archive)
{	
	Timer timer;
	timer.Start();
	if (verbose)
		std::cerr << "Store reference genome in archive (entropy compression)...";
	int s_ref_genome = archive.RegisterStream("ref-genome");

	//commented part is without entropy compression, just 2bits/symbol
	//for(std::vector<uint8_t>& seq : sequences)
	//	archive.AddPart(s_ref_genome, seq, 0);

	CReferenceReads fake_ref_reads(0);
	entropy_coder::CDNACoder coder(fake_ref_reads, verbose);

	coder.Init(true, 0, 9, tot_seqs_len, 0);

	es_t es_seq;
	for (std::vector<uint8_t>& _seq : sequences)
	{
		read_t seq = unpackSeq(_seq);
		//seq.emplace_back(255);
	
		es_seq.append(tuple_types::start_plain, 0, 0);
		for (auto c : seq)
			es_seq.append(tuple_types::plain, c, 1);

		coder.Encode(es_seq);
		es_seq.clear();
	}

	coder.Finish();

	vector<uint8_t> v_output;

	coder.GetOutput(v_output);
	//coder.Restart();

	archive.AddPart(s_ref_genome, v_output, sequences.size());

	if (verbose)
	{
		std::cerr << "Done.\n";
		std::cerr << "Time: " << timer.GetTimeInSec() << "\n";
	}
}

uint32_t CReferenceGenome::GetNPseudoReads()
{
	
	assert(read_len != 0);	
	uint32_t res{};
	for (const auto& seq : sequences)
	{
		auto symbols_in_last_byte = seq.back();
		auto full_bytes = seq.size() - 2 + !symbols_in_last_byte;
		auto seq_len = full_bytes * 4 + symbols_in_last_byte;

		uint32_t n_reads = (seq_len + (read_len - overlap_size) - 1) / (read_len - overlap_size);

		res += n_reads;
	}

	return res;
}

read_pack_t CReferenceGenome::GetPseudoReads()
{
	Timer timer;
	timer.Start();
	if (verbose)
		std::cerr << "Preparing pseudoreads...";
	assert(read_len != 0);	
	read_pack_t res;
	for (const auto& _seq : sequences)
	{
		auto seq = unpackSeq(_seq);
		uint32_t start = 0;
		
		while (start < seq.size())
		{
			uint32_t end = start + read_len;
			if (end > seq.size())
				end = seq.size();

			read_t read;
			read.reserve(end - start + 1);
			read.assign(seq.begin() + start, seq.begin() + end);
			read.push_back(255); //guard
			res.emplace_back(false, read);

			start += read_len - overlap_size;
		}

	}
	if (verbose)
	{
		std::cerr << "Done.\n";
		std::cerr << "Time: " << timer.GetTimeInSec() << "\n";
		std::cerr << "overlap size: " << overlap_size << "\n";
		std::cerr << "pseudoread len: " << read_len << "\n";
	}

	return res;
}

void CReferenceGenome::Release()
{
	sequences.clear();
	sequences.shrink_to_fit();
	read_len = 0;
}