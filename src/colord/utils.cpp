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
#define _CRT_NONSTDC_NO_WARNINGS
#include "utils.h"
#include "zlib.h"
#include <string>
#include <iostream>
#include <vector>
#include <string_view>
#include <cmath>
#include <random>
#include <filesystem>
#ifdef _WIN32
#include <process.h>
#else // linux and mac os
#include <unistd.h>
#endif
#include "archive.h"
#include "params.h"


std::ifstream inOpenOrDie(const std::string& path, std::ios_base::openmode mode)
{
	std::ifstream in(path, mode);
	if (!in)
	{
		std::cerr << "Error: cannot open input file: " << path << "\n";
		exit(1);
	}
	return in;
}

std::ofstream outOpenOrDie(const std::string& path, std::ios_base::openmode mode)
{
	std::ofstream out(path, mode);
	if (!out)
	{
		std::cerr << "Error: cannot open output file: " << path << "\n";
		exit(1);
	}
	return out;
}

bool isFastq(const std::string& path)
{
	auto f = gzopen(path.c_str(), "rb");
	if (!f)
	{
		std::cerr << "Error: cannot open file " << path << "\n";
		exit(1);
	}
	char c;
	if (gzread(f, &c, 1) != 1)
	{
		std::cerr << "Error: file " << path << " is empty\n";
		exit(1);
	}
	gzclose(f);
	return c == '@';
}

bool izGzipFile(const std::string& path)
{
	auto in = inOpenOrDie(path, std::ios::binary);
	uint8_t magic[2];
	in.read((char*)magic, 2);
	return magic[0] == 0x1f && magic[1] == 0x8b;
}


bool fileExists(const std::string& path)
{
	std::ifstream in(path);
	return !!in;
}


kmer_type rev_compl(const kmer_type& kmer, uint32_t len)
{
	uint64_t res{};
	for (uint32_t i = 0; i < len; ++i)
	{
		auto symb = (kmer >> (2 * i)) & 3;
		res <<= 2;
		res += 3 - symb;
	}
	return res;
}

kmer_type canonical(kmer_type kmer, uint32_t len)
{
	auto rc = rev_compl(kmer, len);
	return kmer < rc ? kmer : rc;	
}


// Code fast_upper_bound z https://assets.ctfassets.net/s72atsk5w5jo/Y42yE3L0u4ksUwYc2sCaS/160c2f1da122861ea0c4dfa5bf107566/blog.cpp
inline int fast_upper_bound(std::vector<std::pair<int, int>>& vec, int size, int value)
{
	int low = 0;

	while (size >= 8) {
		int half = size / 2;
		int other_half = size - half;
		int probe = low + half;
		int other_low = low + other_half;
		int v = vec[probe].first;
		size = half;
		low = value > v ? other_low : low;

		half = size / 2;
		other_half = size - half;
		probe = low + half;
		other_low = low + other_half;
		v = vec[probe].first;
		size = half;
		low = value > v ? other_low : low;

		half = size / 2;
		other_half = size - half;
		probe = low + half;
		other_low = low + other_half;
		v = vec[probe].first;
		size = half;
		low = value > v ? other_low : low;
	}

	while (size > 0) {
		int half = size / 2;
		int other_half = size - half;
		int probe = low + half;
		int other_low = low + other_half;
		int v = vec[probe].first;
		size = half;
		low = value > v ? other_low : low;
	};

	return low;
}

void LIS(std::vector<int>& v_in, std::vector<int>& v_out)
{
	std::vector<int> v_pred(v_in.size(), -1);
	std::vector<std::pair<int, int>> v_tmp;

	v_tmp.clear();

	if (v_in.empty())
		return;

	v_tmp.reserve(v_in.size());
	v_tmp.emplace_back(v_in[0], 0);

	int in_size = static_cast<int>(v_in.size());
	int out_len = 1;

	for (int i = 1; i < in_size; ++i)
	{
		int pos;
		int x = v_in[i];

		if (v_tmp[out_len - 1].first < x)
			pos = out_len;
		else
			pos = fast_upper_bound(v_tmp, out_len, x);

		if (pos == out_len)
		{
			v_tmp.emplace_back(x, i);
			++out_len;
		}
		else
		{
			v_tmp[pos].first = x;
			v_tmp[pos].second = i;
		}

		if (pos > 0)
			v_pred[i] = v_tmp[pos - 1].second;
		else
			v_pred[i] = -1;
	}

	v_out.resize(v_tmp.size());

	int cur = v_tmp.back().second;

	for (int i = (int)v_tmp.size() - 1; i >= 0; --i)
	{
		v_out[i] = v_in[cur];
		cur = v_pred[cur];
	}
}

std::string encode_RLE(std::string_view str)
{
	std::string res;
	uint32_t start = 0;
	const uint32_t minimum_symbols_to_store_compressed = 3;
	for (uint32_t i = 1; i < str.size(); ++i)
	{
		if (str[i] != str[start]) //new symbol
		{
			auto len = i - start;

			if (len >= minimum_symbols_to_store_compressed)
				res += str[start] + std::to_string(len);
			else
				for (uint32_t i = 0; i < len; ++i)
					res.push_back(str[start]);

			start = i;
		}
	}
	//last
	auto len = str.size() - start;

	if (len >= minimum_symbols_to_store_compressed)
		res += str[start] + std::to_string(len);
	else
		for (uint32_t i = 0; i < len; ++i)
			res.push_back(str[start]);

	return res;
}


std::string decode_RLE(std::string_view str)
{
	if (str.length() == 0)
		return "";
	std::string res;
	char symb = str[0];
	for (uint32_t i = 1; i < str.length(); ++i)
	{
		if (std::isdigit(str[i]))
		{
			uint32_t len = str[i] - '0';
			++i;
			while (i < str.length() && std::isdigit(str[i]))
			{
				len *= 10;
				len += str[i] - '0';
				++i;
			}

			for (uint32_t j = 0; j < len; ++j)
				res.push_back(symb);

			if (i < str.length())
				symb = str[i];
		}
		else
		{
			res.push_back(symb);
			symb = str[i];
		}
	}
	//last 
	if (!std::isdigit(str[str.length() - 1]))
		res.push_back(symb);

	return res;
}

std::string rand_string(const int length) {

	std::string res;
	static const char symbols[] =
		"0123456789"
		"abcdefghijklmnopqrstuvwxyz"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ";

	static std::uniform_int_distribution<unsigned> dist(0, sizeof(symbols) - 2); //-2 is because '\0' and distribution gets max allowed value
	static std::random_device eng;

	for (int i = 0; i < length; ++i)
		res.push_back(symbols[dist(eng)]);

	return res;
}

std::string get_random_name()
{
	std::string res = ".CoLoRd-" + std::to_string(getpid()) + "-" + rand_string(7);
	return res;
}


std::string create_tmp_dir(const std::string& where)
{
	while (true)
	{
		auto name = get_random_name();
		auto path = std::filesystem::path(where) / name;
		if (!std::filesystem::exists(path))
		{
			std::error_code ec;
			std::filesystem::create_directory(path, ec);
			if (ec)
			{
				std::cerr << "Error creating tmp directory (random path is " << path << ")";
				exit(1);
			}
			return path.string();
		}
	}
}

std::vector<uint8_t> CInfo::Serialize() const
{
	std::vector<uint8_t> res;
	StoreLittleEndian(res, version_major);
	StoreLittleEndian(res, version_minor);
	StoreLittleEndian(res, version_patch);

	StoreLittleEndian(res, total_bytes);
	StoreLittleEndian(res, total_bases);
	StoreLittleEndian(res, total_reads);
	StoreLittleEndian(res, time);
	StoreLittleEndian(res, static_cast<uint32_t>(full_command_line.size()));
	for (char c : full_command_line)
		res.push_back(static_cast<uint8_t>(c));

	return res;
}



void CInfo::Deserialize(const std::vector<uint8_t>& data)
{
	const uint8_t* ptr = data.data();	
	
	loadSingle(ptr, version_major);
	loadSingle(ptr, version_minor);
	loadSingle(ptr, version_patch);

	loadSingle(ptr, total_bytes);
	loadSingle(ptr, total_bases);
	loadSingle(ptr, total_reads);
	loadSingle(ptr, time);
	uint32_t cmd_size;
	loadSingle(ptr, cmd_size);
	full_command_line.reserve(cmd_size);
	for (uint32_t i = 0; i < cmd_size; ++i)
		full_command_line.push_back((char)*ptr++);
}

