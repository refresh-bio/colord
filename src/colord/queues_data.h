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
#include <vector>

struct CCompressElem
{
	uint32_t read_id; //probably not used, kept for debug purposes
	bool hasN;
	read_t read;
	std::vector<uint32_t> ref_reads;
	std::vector<std::vector<kmer_type>> common_kmers;
};

struct CCompressPack
{
	uint32_t id; //for priorities
	std::vector<CCompressElem> data;
};