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
#include "params.h"
#include "id_coder.h"

class CEntrComprHeaders
{
	CParallelQueue<header_pack_t>& headers_queue;
	CArchive& archive;
	int s_header;
	HeaderComprMode headerComprMode;	
	int32_t compression_level;

	entropy_coder::CIDCoder id_coder;

public:
	CEntrComprHeaders(CParallelQueue<header_pack_t>& headers_queue, 		
		CArchive& archive,
		HeaderComprMode headerComprMode,
		int32_t compression_level,
		bool verbose) :

		headers_queue(headers_queue),		
		archive(archive),
		s_header(archive.RegisterStream("header")),
		headerComprMode(headerComprMode),
		compression_level(compression_level),
		id_coder(verbose)
	{
	}

	~CEntrComprHeaders()
	{
	}

	void Compress();
};


class CEntrDecomprHeaders
{
	CArchive& archive;
	CParallelQueue<header_pack_t>& header_decompr_queue;
	int s_header;
	HeaderComprMode headerComprMode;
	int32_t compression_level;
	entropy_coder::CIDCoder id_coder;

public:
	CEntrDecomprHeaders(CArchive& archive,
		CParallelQueue<header_pack_t>& header_decompr_queue,
		HeaderComprMode headerComprMode,
		int32_t compression_level,
		bool verbose) :

		archive(archive),
		header_decompr_queue(header_decompr_queue),
		s_header(archive.GetStreamId("header")),
		headerComprMode(headerComprMode),
		compression_level(compression_level),
		id_coder(verbose)
	{
	}

	~CEntrDecomprHeaders()
	{
	}

	void Decompress();
};