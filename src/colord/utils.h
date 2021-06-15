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

#include "defs.h"
#include <fstream>
#include <vector>
#include <iostream>

#include <cassert>
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <type_traits>
#include <array>
#include <numeric>

using read_t = std::vector<uint8_t>;

//using qual_t = std::string;
using qual_t = std::vector<uint8_t>;
enum class qual_header_type { empty, eq_read_header };

enum class quality_mode_t { none, binary, lossless };

using qual_elem_t = std::pair<read_t, qual_t>;
using header_elem_t = std::pair<std::string, qual_header_type>;

using read_pack_t = std::vector<std::pair<bool, read_t>>;
using qual_pack_t = std::vector<qual_elem_t>;
using header_pack_t = std::vector<header_elem_t>;

using decomp_qual_pack_t = std::vector<qual_t>;
using decomp_read_pack_t = std::vector<read_t>;

const int min_anchor_len = 15;
const int min_skip_len = 10;

enum class tuple_types {
	insertion = 0, deletion = 1, match = 2, substitution = 3, anchor = 4, skip = 5, alt_id = 6, main_ref = 7,
	plain = 8, start_plain = 9, start_es = 10, start_plain_with_Ns = 11, none = 12
};

struct tuple_t {
	tuple_types type;
	int value;
	int no_rep;

	tuple_t(tuple_types _type = tuple_types::none, int _value = 0, int _no_rep = 0) : type(_type), value(_value), no_rep(_no_rep) {};
};

class CompactES {
private:
	std::vector<uint8_t> data;
	uint32_t no_tuples;
/*	uint32_t pos;
	uint32_t max_pos;*/

	std::vector<uint8_t>::iterator ptr;
	std::vector<uint8_t>::iterator end_ptr;

	void pack4(tuple_types type) {
		data.emplace_back(((uint8_t)type) << 4);
	}

	void pack4_4(tuple_types type, uint8_t val) {
		data.emplace_back((((uint8_t)type) << 4) + val);
	}

	void pack4_28(tuple_types type, uint32_t val) {
		data.emplace_back((((uint8_t)type) << 4) + (uint8_t) (val >> 24));
		data.emplace_back((val >> 16) & 0xff);
		data.emplace_back((val >> 8) & 0xff);
		data.emplace_back(val & 0xff);
	}

	void pack32(uint32_t val) {
		data.emplace_back(val >> 24);
		data.emplace_back((val >> 16) & 0xff);
		data.emplace_back((val >> 8) & 0xff);
		data.emplace_back(val & 0xff);
	}

	void unpack4()
	{
		++ptr;
	}

	uint32_t unpack4_4()
	{
		return *ptr++ & 0xf;
	}

	uint32_t unpack4_28()
	{
		uint32_t x = ((uint32_t) (*ptr++ & 0xf)) << 24;
		x += ((uint32_t)*ptr++) << 16;
		x += ((uint32_t)*ptr++) << 8;
		x += (uint32_t)*ptr++;

		return x;
	}

	uint32_t unpack32()
	{
		uint32_t x = ((uint32_t) *ptr++) << 24;
		x += ((uint32_t)*ptr++) << 16;
		x += ((uint32_t)*ptr++) << 8;
		x += (uint32_t)*ptr++;

		return x;
	}

	tuple_types decode_type() {
		return (tuple_types) (*ptr >> 4);
	}

public:
	CompactES() : no_tuples(0)
	{}

	void clear() {
		data.clear();
		no_tuples = 0;
	}

	uint32_t size() {
		return no_tuples;
	}

	uint32_t raw_size() {
		return static_cast<uint32_t>(data.size());
	}

	void restart_reading() {
		ptr = data.begin();
		end_ptr = data.end();
	}

	bool eof() {
		return ptr == end_ptr;
	}

	FORCE_INLINE void append(tuple_types type, int val1, int val2)
	{
		switch (type)
		{
		case tuple_types::insertion:
			pack4_4(tuple_types::insertion, val1);
			break;
		case tuple_types::deletion:
			pack4(tuple_types::deletion);
			break;
		case tuple_types::match:
			pack4(tuple_types::match);
			break;
		case tuple_types::substitution:
			pack4_4(tuple_types::substitution, val1);
			break;
		case tuple_types::anchor:
			pack4_28(tuple_types::anchor, val2);
			break;
		case tuple_types::skip:
			pack4_28(tuple_types::skip, val2);
			break;
		case tuple_types::alt_id:
			pack4_4(tuple_types::alt_id, val2);
			pack32(val1);
			break;
		case tuple_types::main_ref:
			pack4(tuple_types::main_ref);
			break;
		case tuple_types::plain:
			pack4_4(tuple_types::plain, val1);
			break;
		case tuple_types::start_plain:
			pack4(tuple_types::start_plain);
			break;
		case tuple_types::start_es:
			pack4_4(tuple_types::start_es, val2);
			pack32(val1);
			break;
		case tuple_types::start_plain_with_Ns:
			pack4(tuple_types::start_plain_with_Ns);
			break;
		case tuple_types::none:
			break;  //silence g++ warnings
		}

		++no_tuples;
	}

	FORCE_INLINE bool load(tuple_types& type, uint32_t& val1, uint32_t& val2)
	{
		if (ptr == end_ptr)
			return false;

		type = decode_type();

		switch (type)
		{
		case tuple_types::insertion:
			val1 = unpack4_4();
			break;
		case tuple_types::deletion:
			unpack4();
			break;
		case tuple_types::match:
			unpack4();
			break;
		case tuple_types::substitution:
			val1 = unpack4_4();
			break;
		case tuple_types::anchor:
			val2 = unpack4_28();
			break;
		case tuple_types::skip:
			val2 = unpack4_28();
			break;
		case tuple_types::alt_id:
			val2 = unpack4_4();
			val1 = unpack32();
			break;
		case tuple_types::main_ref:
			unpack4();
			break;
		case tuple_types::plain:
			val1 = unpack4_4();
			break;
		case tuple_types::start_plain:
			unpack4();
			break;
		case tuple_types::start_es:
			val2 = unpack4_4();
			val1 = unpack32();
			break;
		case tuple_types::start_plain_with_Ns:
			unpack4();
			break;
		case tuple_types::none:
			break;  //silence g++ warnings
		}

		return true;
	}

	void reserve(uint32_t res_size)
	{
		data.reserve(res_size);
	}

	void compact()
	{
		data.shrink_to_fit();
	}
};

//using es_t = std::vector<tuple_t>;
using es_t = CompactES;

using context_t = uint64_t;


struct CMissmatchCoder
{
#ifdef PLAIN_SUBSTITUTION
	static const char* get_missmatch_codes()
	{
		static char codes[256];
		std::fill_n(codes, 256, -1);
		codes['a'] = 0;
		codes['c'] = 1;
		codes['g'] = 2;
		codes['t'] = 3;
		return codes;
	}

#else
	static const char* get_symb_codes()
	{
		static char codes[256];
		std::fill_n(codes, 256, -1);
		codes[static_cast<unsigned char>('A')] = 0;
		codes[static_cast<unsigned char>('C')] = 1;
		codes[static_cast<unsigned char>('G')] = 2;
		codes[static_cast<unsigned char>('T')] = 3;
		return codes;
	}

	static const char* get_missmatch_codes()
	{
		static char codes[256];
		std::fill_n(codes, 256, -1);
		codes[static_cast<unsigned char>('X')] = 0;
		codes[static_cast<unsigned char>('Y')] = 1;
		codes[static_cast<unsigned char>('Z')] = 2;
		return codes;
	}
#endif

public:
	static int GetMissmatchCode(char missmatch_symb)
	{
		static const char* enc_codes = get_missmatch_codes();
		return enc_codes[(uint8_t)missmatch_symb];
	}

#ifdef PLAIN_SUBSTITUTION
	static char encode_missmatch_symb(char ref_symb, char new_symb)
	{
		assert(ref_symb != new_symb);
		return std::tolower(new_symb);
	}
	static char decode_mismatch_symb(char ref_symb, char encoded)
	{
		assert(ref_symb != std::toupper(encoded));
		return std::toupper(encoded);
	}
	static bool is_missmatch_es_symb(char symb)
	{
		return symb == 'a' || symb == 'c' || symb == 'g' || symb == 't';
	}
#else
	static char encode_missmatch_symb(char ref_symb, char new_symb)
	{
		static const char* codes = get_symb_codes();
		//Matches on diagnal important for plain read part when store as big edit script
		static const char mm[4][4] = {
			{'M', 'X', 'Y', 'Z'}, //input symbol is A
			{'X', 'M', 'Y', 'Z'}, //input symbol is C
			{'X', 'Y', 'M', 'Z'}, //input symbol is G
			{'X', 'Y', 'Z', 'M'}  //input symbol is T
		};
		return mm[(uint8_t)codes[(uint8_t)ref_symb]][(uint8_t)codes[(uint8_t)new_symb]];
	}
	static char decode_mismatch_symb(char ref_symb, char encoded)
	{
		static const char* codes = get_symb_codes();
		static const char* enc_codes = get_missmatch_codes();
		static const char mm[4][3] = {
			{'C', 'G', 'T'}, //input symbol is A
			{'A', 'G', 'T'}, //input symbol is C
			{'A', 'C', 'T'}, //input symbol is G
			{'A', 'C', 'G'}  //input symbol is T
		};
		return mm[(uint8_t)codes[(uint8_t)ref_symb]][(uint8_t)enc_codes[(uint8_t)encoded]];
	}
	static bool is_missmatch_es_symb(char symb)
	{
		return symb == 'X' || symb == 'Y' || symb == 'Z';
	}
#endif
};

inline size_t read_len(const read_t& read)
{
	assert(read.size());
	return read.size() - 1;
}


inline read_t get_rev_compl(const read_t& read)
{
	auto len = read_len(read);
	read_t res(len + 1);
	res[len] = 255; //guard

	for (uint32_t i = 0; i < len; ++i)
		res[i] = 3 - read[len - i - 1];
	return res;
}


class read_view
{
	const uint8_t* _data;
	size_t len;
public:
	read_view(const read_t& read) :
		_data(read.data()),
		len(read_len(read))
	{

	}
	read_view(const uint8_t* _data, size_t len) :
		_data(_data), len(len) {}

	size_t length() const
	{
		return len;
	}
	size_t size() const
	{
		return len;
	}
	bool empty() const
	{
		return len == 0;
	}
	read_view substr(size_t start, size_t len = std::string::npos)
	{
		len = std::min(len, this->len - start);
		return read_view(_data + start, len);
	}

	bool operator==(const read_view& rhs) const
	{
		if (len != rhs.len)
			return false;
		for (size_t i = 0; i < len; ++i)
			if (_data[i] != rhs._data[i])
				return false;
		return true;
	}
	bool operator!=(const read_view& rhs) const
	{
		return !operator==(rhs);
	}
	const uint8_t& operator[](int32_t idx) const
	{
		return _data[idx];
	}
	const uint8_t* data() const
	{
		return _data;
	}

	const uint8_t* begin() const
	{
		return _data;
	}
	const uint8_t* end() const
	{
		return _data + len;
	}

	read_t get_reversed() const
	{
		read_t res(len + 1);

		std::reverse_copy(_data, _data + len, res.begin());
		res[len] = 255;		// guard

		return res;
	}
};

template<typename T> constexpr T min3(T a, T b, T c) {
	return (a < b) ? (a < c ? a : c) : (b < c ? b : c);
}

static const char* SymbToBinMap = []() {
	static char _m[256];
	std::fill(std::begin(_m), std::end(_m), -1);
	//_m['A'] = _m['a'] = 0;
	//_m['C'] = _m['c'] = 1;
	//_m['G'] = _m['g'] = 2;
	//_m['T'] = _m['t'] = 3;
	_m[static_cast<unsigned char>('A')] = 0;
	_m[static_cast<unsigned char>('C')] = 1;
	_m[static_cast<unsigned char>('G')] = 2;
	_m[static_cast<unsigned char>('T')] = 3;
	_m[static_cast<unsigned char>('N')] = 4;
	return _m;
}();

template<typename T>
void StoreLittleEndian(uint8_t* buff, const T& data)
{	
	for (uint32_t b = 0; b < sizeof(data); ++b)	
		*buff++ = data >> (8 * b);
}

template<typename T>
void StoreLittleEndian(std::vector<uint8_t>& buff, const T& data)
{
	for (uint32_t b = 0; b < sizeof(data); ++b)
		buff.push_back(static_cast<uint8_t>(data >> (8 * b)));
}

template<>
inline void StoreLittleEndian<double>(std::vector<uint8_t>& buff, const double& data)
{
	static bool isLittleEndian = [] {
		int x = 1;
		return *reinterpret_cast<uint8_t*>(&x) == 1;
	}();
	if (isLittleEndian)
	{
		auto ptr = reinterpret_cast<const uint8_t*>(&data);
		for (uint32_t b = 0; b < sizeof(data); ++b)
			buff.push_back(*ptr++);
	}
	else
	{
		auto ptr = reinterpret_cast<const uint8_t*>(&data) + sizeof(data) - 1;
		for (uint32_t b = 0; b < sizeof(data); ++b)
			buff.push_back(*ptr--);
	}
}

template<typename T>
void StoreLittleEndian(std::ostream& buff, const T& data)
{
	for (uint32_t b = 0; b < sizeof(data); ++b)
	{
		uint8_t byte = data >> (8 * b);
		buff.write(reinterpret_cast<char*>(&byte), 1);
	}		
}

template<typename T>
void LoadLittleEndian(const uint8_t* buff, T& data)
{
	data = T{};
	for (uint32_t b = 0; b < sizeof(data); ++b)
	{
		uint8_t byte = *buff++;
		data += (T)byte << (8 * b);
	}
}

template<>
inline void LoadLittleEndian<double>(const uint8_t* buff, double& data)
{
	static bool isLittleEndian = [] {
		int x = 1;
		return *reinterpret_cast<uint8_t*>(&x) == 1;
	}();
	if (isLittleEndian)
	{
		auto ptr = reinterpret_cast<uint8_t*>(&data);
		for (uint32_t b = 0; b < sizeof(data); ++b)
			*ptr++ = *buff++;
	}
	else
	{
		auto ptr = reinterpret_cast<int8_t*>(&data) + sizeof(data) - 1;
		for (uint32_t b = 0; b < sizeof(data); ++b)
			*ptr-- = *buff++;
	}
}

template<typename T>
void LoadLittleEndian(const std::vector<uint8_t>& buff, T& data)
{	
	LoadLittleEndian(buff.data(), data);	
}

template<typename T>
void LoadLittleEndian(std::istream& buff, T& data)
{
	data = T{};
	for (uint32_t b = 0; b < sizeof(data); ++b)
	{
		uint8_t byte;
		buff.read((char*)&byte, 1);
		data += (uint64_t)byte << (8 * b);	
	}	
}

std::ifstream inOpenOrDie(const std::string& path, std::ios_base::openmode mode = std::ios_base::in);

std::ofstream outOpenOrDie(const std::string& path, std::ios_base::openmode mode = std::ios_base::out);

bool isFastq(const std::string& path);

bool izGzipFile(const std::string& path);

bool fileExists(const std::string& path);

kmer_type rev_compl(const kmer_type& kmer, uint32_t len);

kmer_type canonical(kmer_type kmer, uint32_t len);

void LIS(std::vector<int>& v_in, std::vector<int>& v_out);

std::string encode_RLE(std::string_view str);

std::string decode_RLE(std::string_view str);

class CPercentProgress
{
	uint64_t max_iter;
	uint64_t cur_iter = 0;
	uint32_t cur_percent = std::numeric_limits<uint32_t>::max();
	bool finished = false;
	void update()
	{
		if (finished)
			return;
		uint32_t new_percent = static_cast<uint32_t>(cur_iter * 100 / max_iter);
		if (new_percent != cur_percent)
		{
			std::cerr << new_percent << "%\r";
			cur_percent = new_percent;
		}
		if (new_percent == 100)
		{
			finished = true;
			std::cerr << "\n";
		}
	}
public:
	/*
	Max iter may be approximate, but one needs to call ForceFinish in such a case
	*/
	CPercentProgress(uint64_t max_iter):max_iter(max_iter)
	{

	}
	void LongNIters(uint64_t n)
	{
		cur_iter += n;
		update();
	}
	void LogIter()
	{
		LongNIters(1u);
	}

	void ForceFinish()
	{
		cur_iter = max_iter;
		update();
	}
};

class Timer
{
	//std::chrono::steady_clock::time_point start;
	std::chrono::time_point<std::chrono::high_resolution_clock> start;
public:
	void Start()
	{
		start = std::chrono::high_resolution_clock::now();
	}
	void Log(std::ostream& log)
	{
		auto end = std::chrono::high_resolution_clock::now();
		auto ts = std::chrono::duration<double>(end - start).count();
		//auto tm = std::chrono::duration<double, std::ratio<60>>(end - start).count();
		//log << "Time: " << ts << "sec = " << tm << " min\n";
		log << "Time: " << ts << "s\n";
	}

	auto GetTimeInSec()
	{
		auto end = std::chrono::high_resolution_clock::now();
		return std::chrono::duration<double>(end - start).count();		
	}
};


std::string create_tmp_dir(const std::string& where=".");


struct CInfo
{
	uint32_t version_major;
	uint32_t version_minor;
	uint64_t total_bytes;
	uint64_t total_bases;
	uint32_t total_reads;
	uint64_t time;
	std::string full_command_line;
private:
	template<typename T>
	void loadSingle(const uint8_t*& ptr, T& out)
	{
		LoadLittleEndian(ptr, out);
		ptr += sizeof(T);
	}
public:
	std::vector<uint8_t> Serialize() const;	
	void Deserialize(const std::vector<uint8_t>& data);
};

class CEntropy
{
	static const uint32_t plain_num = 4;
	static const uint32_t es_num = 11;

	static inline const uint8_t es_sym[] = { 'A', 'C', 'D', 'G', 'M', 'T', 'X', 'Y', 'Z', 'S', 'R' };

public:
	template<typename T>
	static double entropy_dna(T str)
	{
		uint32_t histo[plain_num]{};

		for (auto c : str)
			++histo[c];

		double sum{};
		for (uint32_t i = 0; i < plain_num; ++i)
			sum += histo[i];

		double sum_rec = 1.0 / sum;

		double entr{};
		for (uint32_t c = 0; c < plain_num; ++c)
			if (histo[c])
			{
				double p = (double)histo[c] * sum_rec;
				entr += log2(p) * p;
			}

		return -entr;
	}

	template<typename T>
	static double entropy_es(T str)
	{
		uint32_t histo[128]{};

		for (auto c : str)
			++histo[static_cast<typename std::make_unsigned<decltype(c)>::type>(c)]; //signed subscribt warnings in g++...

		double sum{};
		for (uint32_t i = 0; i < es_num; ++i)
			sum += histo[es_sym[i]];

		double sum_rec = 1.0 / sum;

		double entr{};
		for (uint32_t c = 0; c < es_num; ++c)
			if (histo[es_sym[c]])
			{
				double p = (double)histo[es_sym[c]] * sum_rec;
				entr += log2(p) * p;
			}

		return -entr;
	}
};


class CEntropyEstimator
{
	using dna_stats_t = std::array<uint32_t, 4>;
	using es_stats_t = std::array<uint32_t, 12>;
	using decision_stats_t = std::array<uint32_t, 2>;

	dna_stats_t dna_stats;
	es_stats_t es_stats;
	decision_stats_t decision_stats;

	std::array<double, 4> dna_logs;
	std::array<double, 12> es_logs;
	std::array<double, 2> decision_logs;

	uint32_t dna_sum;
	uint32_t es_sum;
	uint32_t decision_sum;

	const uint32_t dna_sum_max = 1 << 20;
	const uint32_t es_sum_max = 1 << 20;
	const uint32_t decision_sum_max = 1 << 20;

	std::array<uint32_t, 128> es_codes;
	std::vector<uint32_t> es_lens;

	template<typename STATS>
	void rescale(STATS& arr, uint32_t &sum, uint32_t max_sum)
	{
		while (sum > max_sum)
		{
			sum = 0;

			for (auto& x : arr)
			{
				x = (x + 1) / 2;
				sum += x;
			}
		}
	}

	template<typename STATS, typename LOGS>
	void calc_logs(STATS& arr_stat, LOGS& arr_log, uint32_t sum)
	{
		double sum_rec = 1.0 / sum;

		for (size_t i = 0; i < arr_stat.size(); ++i)
			if (arr_stat[i] != 0)
				arr_log[i] = -log2((double)arr_stat[i] * sum_rec);
			else
				arr_log[i] = 0.0;
	}

	template<typename STATS>
	void clear(STATS& arr, uint32_t &sum, uint32_t init_value = 0)
	{
		std::fill(arr.begin(), arr.end(), init_value);
		sum = static_cast<uint32_t>(arr.size()) * init_value;
	}

	template<typename STATS>
	void add_stats(STATS &arr_global, STATS &arr_local, uint32_t &sum_global, uint32_t sum_local)
	{
		for (size_t i = 0; i < arr_global.size(); ++i)
			arr_global[i] += arr_local[i];

		sum_global += sum_local;
	}

	void analyze_es(std::string_view edit_script, es_stats_t& loc_es_stats, uint32_t& loc_es_sum, es_stats_t& read_es_stats, std::vector<uint32_t> &es_lens)
	{
		char c = ' ';
		uint32_t len = 0;

		es_lens.clear();

		for(size_t i = 0; i <= edit_script.size(); ++i)
		{
			char e = (i < edit_script.size()) ? edit_script[i] : ' ';

			if (e == c)
				++len;
			else
			{
				if (c == 'D')
				{
					if (len >= min_skip_len)
					{
						++loc_es_stats[es_codes['S']];
						++loc_es_sum;
						++read_es_stats[es_codes['S']];
						es_lens.emplace_back(len);
					}
					else
					{
						loc_es_stats[es_codes['D']] += len;
						loc_es_sum += len;
						read_es_stats[es_codes['D']] += len;
					}
				}
				else if (c == 'M')
				{
					if (len >= min_anchor_len)
					{
						++loc_es_stats[es_codes['R']];
						++loc_es_sum;
						++read_es_stats[es_codes['R']];
						es_lens.emplace_back(len);
					}
					else
					{
						loc_es_stats[es_codes['M']] += len;
						loc_es_sum += len;
						read_es_stats[es_codes['M']] += len;
					}
				}
				else if(c != ' ')
				{
					++loc_es_stats[es_codes[c]];
					++loc_es_sum;
					++read_es_stats[es_codes[c]];
				}

				c = e;
				len = 1;
			}
		}
	}

	void reset()
	{
/*		clear(dna_stats, dna_sum, 1u);
		clear(es_stats, es_sum, 1u);
		clear(decision_stats, decision_sum, 1u);*/
		clear(dna_stats, dna_sum, 1u);
		clear(es_stats, es_sum, 1u);
		clear(decision_stats, decision_sum, 1u);

		calc_logs(dna_stats, dna_logs, dna_sum);
		calc_logs(es_stats, es_logs, es_sum);
		calc_logs(decision_stats, decision_logs, decision_sum);
	}
	
	// *****************************************************************
	constexpr uint64_t ilog2(uint64_t x)
	{
		uint64_t r = 0;

		for (; x; ++r)
			x >>= 1;

		return r;
	}

public:
	CEntropyEstimator() 
	{
		std::fill(es_codes.begin(), es_codes.end(), 11);
		es_codes['A'] = 0;
		es_codes['C'] = 1;
		es_codes['G'] = 2;
		es_codes['T'] = 3;
		es_codes['D'] = 4;
		es_codes['M'] = 5;
		es_codes['X'] = 6;
		es_codes['Y'] = 7;
		es_codes['Z'] = 8;
		es_codes['S'] = 9;		// skip
		es_codes['R'] = 10;		// anchor
//		es_codes[''] = 9;		// reserved for anchors

		reset();
	}

	//to bym wywolywal jak przychodzi nowa paczka z kolejki
	void Reset()
	{
		reset();
	}

	//czy dla readów ca³kiem plain tez wywoa³ywac
	//czy dla readów z Nkami tez wywo³ywaæ
	//mozliwe symbole read: 
	//{0, 1, 2, 3} 
	//lub
	//{-1, 0, 1, 2, 3} jezeli wrzucac tutaj tez ready z Nkami   
	//to bym wywo³ywa³ przed rozpoczeciem kodowania danego reada
	void LogRead(read_view read)
	{
		for (auto c : read)
			++dna_stats[c];

		dna_sum += static_cast<uint32_t>(read.size());

		rescale(dna_stats, dna_sum, dna_sum_max);
		calc_logs(dna_stats, dna_logs, dna_sum);
	}

#if 0
	//to by moglo zwracac true jezeli nalezy kodowac z edit scriptem, false jesli plain
	//i to by tez od razu sobie zbieralo statystyki dla edit scripta
	// mozliwe symbole edit_script: 
	//static inline const uint8_t es_sym[] = { 'A', 'C', 'D', 'G', 'M', 'T', 'X', 'Y', 'Z' }
	//mozliwe symbole plain: 
	//{0, 1, 2, 3}
	//to bym wywolywal w trakcie kodowania reada, dla kawalkow miedzykotwicznych do sprawdzenia czego uzywac
	bool EncodeWithEditScript0(std::string_view edit_script, read_view plain, size_t ref_len)
	{
//		return CEntropy::entropy_es(edit_script) < CEntropy::entropy_dna(plain);

		es_stats_t loc_es_stats;
		es_stats_t loc_plain_stats;
		es_stats_t read_es_stats;
		es_stats_t read_plain_stats;
		uint32_t loc_es_sum;
		uint32_t loc_plain_sum;
		uint32_t read_es_sum = 0;
		uint32_t read_plain_sum = 0;

		clear(read_es_stats, read_es_sum);
		clear(read_plain_stats, read_plain_sum);

		loc_es_stats = es_stats;
		loc_plain_stats = es_stats;
		loc_es_sum = es_sum;
		loc_plain_sum = es_sum;

		double es_cost = 0;
		double plain_cost = 0;

		for (auto e : edit_script)
		{
			++loc_es_stats[es_codes[e]];
			++loc_es_sum;
			++read_es_stats[es_codes[e]];
			++read_es_sum;
		}

		for (auto c : plain)
		{
			++loc_plain_stats[c];
			++loc_plain_sum;
			++read_plain_stats[c];
			++read_plain_sum;
		}

		uint32_t d_cost = ilog2(ref_len + 1);
//		uint32_t i_cost = ilog2(plain.size() + 1);
//		d_cost = std::min(plain.size(), 4ull);
		d_cost = 3;

		loc_plain_stats[es_codes['D']] += d_cost;
		loc_plain_sum += d_cost;
/*		read_plain_stats[es_codes['D']] += ref_len;
		read_plain_sum += ref_len;*/
		read_plain_stats[es_codes['D']] += d_cost;
		read_plain_sum += d_cost;
/*		read_plain_stats[9] += i_cost;
		read_plain_sum += i_cost;*/

		std::array<double, 10> loc_es_logs;

		calc_logs(loc_es_stats, loc_es_logs, loc_es_sum);
//		calc_logs(read_es_stats, loc_es_logs, read_es_sum);

		std::array<double, 10> loc_plain_logs;
//		calc_logs(loc_plain_stats, loc_plain_logs, loc_plain_sum);

		auto tmp_plain_stats = loc_plain_stats;
		fill(tmp_plain_stats.begin() + 5, tmp_plain_stats.end(), 0u);
		auto tmp_plain_sum = std::accumulate(tmp_plain_stats.begin(), tmp_plain_stats.begin() + 5, 0u);

		calc_logs(tmp_plain_stats, loc_plain_logs, tmp_plain_sum);

		for (size_t i = 0; i < es_stats.size(); ++i)
			es_cost += read_es_stats[i] * loc_es_logs[i];

		for (size_t i = 0; i < es_stats.size(); ++i)
			plain_cost += read_plain_stats[i] * loc_plain_logs[i];

		if (plain_cost < es_cost)
		{
			es_stats = loc_plain_stats;
			es_sum = loc_plain_sum;

			rescale(es_stats, es_sum, es_sum_max);

			return false;
		}
		else
		{
			es_stats = loc_es_stats;
			es_sum = loc_es_sum;

			rescale(es_stats, es_sum, es_sum_max);

			return true;
		}
	}
#endif

	bool EncodeWithEditScript(std::string_view edit_script, read_view plain, size_t ref_len)
	{
		es_stats_t loc_es_stats;
		es_stats_t read_es_stats;
		dna_stats_t read_plain_stats;
		uint32_t loc_es_sum;

		uint32_t read_es_sum = 0;
		uint32_t read_plain_sum = 0;

		clear(read_es_stats, read_es_sum);
		clear(read_plain_stats, read_plain_sum);

		loc_es_stats = es_stats;
		loc_es_sum = es_sum;

		double es_cost = decision_logs[0];
		double plain_cost = decision_logs[1];

		analyze_es(edit_script, loc_es_stats, loc_es_sum, read_es_stats, es_lens);

/*		for(auto c : edit_script)
		{
			++loc_es_stats[es_codes[c]];
			++loc_es_sum;
			++read_es_stats[es_codes[c]];
		}*/

		for (auto c : plain)
		{
			++read_plain_stats[c];
			++read_plain_sum;
		}

//		std::array<double, 10> loc_es_logs;

		calc_logs(loc_es_stats, es_logs, loc_es_sum);

		for (size_t i = 0; i < es_stats.size(); ++i)
			es_cost += read_es_stats[i] * es_logs[i];
		for (auto x : es_lens)
			es_cost += ilog2(x) + 1;

		for (size_t i = 0; i < dna_stats.size(); ++i)
			plain_cost += read_plain_stats[i] * dna_logs[i];
//		plain_cost += 30;		// for deletions after insert
		plain_cost += ilog2(ref_len) + 1;

		bool choose_plain = plain_cost < es_cost;

		if (choose_plain)
		{
			++decision_stats[1];

			rescale(es_stats, es_sum, es_sum_max);
		}
		else
		{
			++decision_stats[0];
			es_stats = loc_es_stats;
			es_sum = loc_es_sum;

			rescale(es_stats, es_sum, es_sum_max);
		}

		++decision_sum;
		rescale(decision_stats, decision_sum, decision_sum_max);
		calc_logs(decision_stats, decision_logs, decision_sum);

		return !choose_plain;
	}
};

inline uint64_t round_to_pow_of_2(uint64_t to_round)
{
	if ((to_round & (to_round - 1)))
	{
		while ((to_round & (to_round - 1)))
			to_round &= to_round - 1;
		to_round *= 2;
	}
	return to_round;
}
