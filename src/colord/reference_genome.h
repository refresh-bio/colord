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
#include "archive.h"
#include <string>
#include <vector>

class CReferenceGenome
{
	uint64_t tot_seqs_len{};
	uint32_t overlap_size;
	uint32_t read_len{};
	std::vector<uint8_t> checksum;

	bool verbose;
	std::vector<std::vector<uint8_t>> sequences;

	void packSeq(std::vector<uint8_t>& seq);

	std::vector<uint8_t> unpackSeq(const std::vector<uint8_t>& seq);

	void newSeq()
	{
		sequences.emplace_back();
	}

	void seqEnd()
	{
		tot_seqs_len += sequences.back().size();
		packSeq(sequences.back());
	}

	void addSymb(char c)
	{
		c = toupper(c);
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			sequences.back().push_back(SymbToBinMap[(uint8_t)c]);
	}

public:
	explicit CReferenceGenome(const std::string& path, uint32_t overlap_size, bool calc_checksum, bool verbose);
	explicit CReferenceGenome(CArchive& archive, uint32_t overlap_size, bool verbose);

	void Store(const std::string& path, bool as_fastq);

	void Store(CArchive& archive);

	uint64_t GetTotSeqsLen() const
	{
		return tot_seqs_len;
	}

	uint32_t GetTotNSeqs() const
	{
		return sequences.size();
	}

	void SetReadLen(uint32_t val)
	{
		this->read_len = val;
	}

	uint32_t GetNPseudoReads();

	read_pack_t GetPseudoReads();

	void Release();

	std::vector<uint8_t> GetChecksum() const
	{
		return checksum;
	}
};
