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

#include "basic_coder.h"
#include "params.h"

#include <array>
#include <string>

namespace entropy_coder
{

// *****************************************************************
class CIDCoder : public CBasicCoder
{
	bool verbose = false;
	enum class token_type_t { literal = 0, numeric = 1, hexadecimal = 2 };
	array<bool, 256> a_numeric, a_literal, a_hexadecimal;

	typedef vector<tuple<token_type_t, uint8_t, uint32_t, uint32_t>> v_tokens_t;
	v_tokens_t v_tokens_prev, v_tokens_cur;
	vector<int64_t> v_deltas;
	string id_prev;
	uint32_t size_prev;

	HeaderComprMode header_mode;
	int compression_level;

	context_t ctx_flags{};

	//	using rcmfs_quality_original_t = CRangeCoderModelFixedSize<CVectorIOStream, 96, 1 << 20, 32>;

	using rcmfs_plus_id_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 15, 1>;
	using rcmfs_flags_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 15, 1>;
	using rcmfs_numeric_size_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 19, 32>;
	using rcmfs_numeric_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 19, 32>;
	using rcmfs_numeric_small_t = CRangeCoderModelFixedSize<CVectorIOStream, 4, 1 << 15, 1>;
	using rcmfs_hexa_t = CRangeCoderModelFixedSize<CVectorIOStream, 16, 1 << 15, 1>;
	using rcmfs_literal_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 20, 64>;
	using rcmfs_literal_same_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 15, 1>;
	using rcmfs_literal_same_length_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 15, 1>;
	using rcmfs_plain_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 128, 1 << 19, 32>;

	using ctx_map_plus_id_t = CContextHM<rcmfs_plus_id_t>;
	using ctx_map_flags_t = CContextHM<rcmfs_flags_t>;
	using ctx_map_numeric_size_t = CContextHM<rcmfs_numeric_size_t>;
	using ctx_map_numeric_t = CContextHM<rcmfs_numeric_t>;
	using ctx_map_numeric_small_t = CContextHM<rcmfs_numeric_small_t>;
	using ctx_map_hexa_t = CContextHM<rcmfs_hexa_t>;
	using ctx_map_literal_t = CContextHM<rcmfs_literal_t>;
	using ctx_map_literal_same_t = CContextHM<rcmfs_literal_same_t>;
	using ctx_map_literal_same_length_t = CContextHM<rcmfs_literal_same_length_t>;
	using ctx_map_plain_t = CContextHM<rcmfs_plain_t>;


	ctx_map_plus_id_t m_ctx_rc_plus_id;
	ctx_map_flags_t m_ctx_rc_flags;
	ctx_map_numeric_size_t m_ctx_rc_numeric_size;
	ctx_map_numeric_t m_ctx_rc_numeric;
	ctx_map_numeric_small_t m_ctx_rc_numeric_small;
	ctx_map_hexa_t m_ctx_rc_hexa;
	ctx_map_literal_t m_ctx_rc_literal;
	ctx_map_literal_same_t m_ctx_rc_literal_same;
	ctx_map_literal_same_length_t m_ctx_rc_literal_same_length;
	ctx_map_plain_t m_ctx_rc_plain;

	rcmfs_plus_id_t* tpl_ctx_rc_plus_id;
	rcmfs_flags_t* tpl_ctx_rc_flags;
	rcmfs_numeric_size_t* tpl_ctx_rc_numeric_size;
	rcmfs_numeric_t* tpl_ctx_rc_numeric;
	rcmfs_numeric_small_t* tpl_ctx_rc_numeric_small;
	rcmfs_hexa_t* tpl_ctx_rc_hexa;
	rcmfs_literal_t* tpl_ctx_rc_literal;
	rcmfs_literal_same_t* tpl_ctx_rc_literal_same;
	rcmfs_literal_same_length_t* tpl_ctx_rc_literal_same_length;
	rcmfs_plain_t* tpl_ctx_rc_plain;

	string extract_instrument(const string& id);

	uint32_t tokenize(const string &id, v_tokens_t& v_tokens);

	void compress_lossless(const bool plus_id, const string &id);
	void compress_instrument(const bool plus_id, const string& id);

	void decompress_none(bool &plus_id, string &id);
	void decompress_lossless(bool &plus_id, string& id);
	void decompress_instrument(bool &plus_id, string& id);

	bool is_numeric(uint8_t c)
	{
		return a_numeric[c];
	}

	bool is_literal(uint8_t c)
	{
		return a_literal[c];
	}

	bool is_hexadecimal(uint8_t c)
	{
		return a_hexadecimal[c];
	}

	void init_symbol_classes();

	bool token_types_same(const v_tokens_t& v1, const v_tokens_t& v2)
	{
		if (v1.size() != v2.size())
			return false;

		for (size_t i = 0; i < v1.size(); ++i)
			if (get<0>(v1[i]) != get<0>(v2[i]) || get<1>(v1[i]) != get<1>(v2[i]))
				return false;

		return true;
	}

	int64_t get_int(const string &p, const uint32_t begin, const uint32_t end)
	{
		int64_t r = 0;

		for (uint32_t i = begin; i < end; ++i)
			r = r * 10 + (int64_t)(p[i] - '0');

		return r;
	}

	int32_t store_int(string &id, int64_t val)
	{
		int n_dig = 0;
		char p[16];

		if (val < 10ll)						n_dig = 1;
		else if (val < 100ll)				n_dig = 2;
		else if (val < 1000ll)				n_dig = 3;
		else if (val < 10000ll)				n_dig = 4;
		else if (val < 100000ll)			n_dig = 5;
		else if (val < 1000000ll)			n_dig = 6;
		else if (val < 10000000ll)			n_dig = 7;
		else if (val < 100000000ll)			n_dig = 8;
		else if (val < 1000000000ll)		n_dig = 9;
		else if (val < 10000000000ll)		n_dig = 10;
		else if (val < 100000000000ll)		n_dig = 11;
		else if (val < 1000000000000ll)		n_dig = 12;
		else if (val < 10000000000000ll)	n_dig = 13;
		else if (val < 100000000000000ll)	n_dig = 14;
		else if (val < 1000000000000000ll)	n_dig = 15;
		else
		{
			cerr << "Error in store_int(): " << val << endl;
		}


		for (int i = n_dig - 1; i >= 0; --i)
		{
			p[i] = '0' + (val % 10);
			val /= 10;
		}

		id.append(p, n_dig);

		return n_dig;
	}

public:
	CIDCoder(bool verbose) : CBasicCoder(), verbose(verbose)
	{
		tpl_ctx_rc_plus_id = nullptr;
		tpl_ctx_rc_flags = nullptr;
		tpl_ctx_rc_numeric_size = nullptr;
		tpl_ctx_rc_numeric = nullptr;
		tpl_ctx_rc_numeric_small = nullptr;
		tpl_ctx_rc_hexa = nullptr;
		tpl_ctx_rc_literal = nullptr;
		tpl_ctx_rc_literal_same = nullptr;
		tpl_ctx_rc_literal_same_length = nullptr;
		tpl_ctx_rc_plain = nullptr;
	}

	~CIDCoder()
	{
		delete tpl_ctx_rc_plus_id;
		delete tpl_ctx_rc_flags;
		delete tpl_ctx_rc_numeric_size;
		delete tpl_ctx_rc_numeric;
		delete tpl_ctx_rc_numeric_small;
		delete tpl_ctx_rc_hexa;
		delete tpl_ctx_rc_literal;
		delete tpl_ctx_rc_literal_same;
		delete tpl_ctx_rc_literal_same_length;
		delete tpl_ctx_rc_plain;

		if(verbose)
		{
			cout << "Plus id size: " << m_ctx_rc_plus_id.get_size() << endl;
			cout << "Flags size: " << m_ctx_rc_flags.get_size() << endl;
			cout << "Numeric size size: " << m_ctx_rc_numeric_size.get_size() << endl;
			cout << "Numeric small size: " << m_ctx_rc_numeric_small.get_size() << endl;
			cout << "Numeric size: " << m_ctx_rc_numeric.get_size() << endl;
			cout << "Hexadecimal size: " << m_ctx_rc_hexa.get_size() << endl;
			cout << "Literal size: " << m_ctx_rc_literal.get_size() << endl;
			cout << "Literal same size: " << m_ctx_rc_literal_same.get_size() << endl;
			cout << "Literal same length size: " << m_ctx_rc_literal_same_length.get_size() << endl;
			cout << "Plain size: " << m_ctx_rc_plain.get_size() << endl;
		}
	}

	void Init(bool is_compressing, HeaderComprMode _header_mode, int _compression_level);

	void Restart();

	void Finish();

	void Encode(bool plus_id, string id);

	void Decode(bool &plus_id, string& id);

	void ResetReadPrev();
};
} //namespace entropy_coder

  // EOF