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
#include "id_coder.h"

namespace entropy_coder
{
// *****************************************************************
void CIDCoder::Init(bool _is_compressing, HeaderComprMode _header_mode, int _compression_level)
{
	is_compressing = _is_compressing;
	header_mode = _header_mode;
	compression_level = _compression_level;

	v_vios_io = new CVectorIOStream(v_io);

	init_symbol_classes();

	CBasicRangeCoder<CVectorIOStream>* rc;

	if (is_compressing)
	{
		rc = rce = new CRangeEncoder<CVectorIOStream>(*v_vios_io);
		rce->Start();
	}
	else
	{
		rc = rcd = new CRangeDecoder<CVectorIOStream>(*v_vios_io);
		rcd->Start();
	}

	delete tpl_ctx_rc_plus_id;
	tpl_ctx_rc_plus_id = new rcmfs_plus_id_t(rc, nullptr, is_compressing);

	// *** Flags for read ids
	delete tpl_ctx_rc_flags;
	tpl_ctx_rc_flags = new rcmfs_flags_t(rc, nullptr, is_compressing);

	delete tpl_ctx_rc_numeric_size;
	tpl_ctx_rc_numeric_size = new rcmfs_numeric_size_t(rc, nullptr, is_compressing);

	delete tpl_ctx_rc_numeric;
	tpl_ctx_rc_numeric = new rcmfs_numeric_t(rc, nullptr, is_compressing);

	delete tpl_ctx_rc_numeric_small;
	tpl_ctx_rc_numeric_small = new rcmfs_numeric_small_t(rc, nullptr, is_compressing);
	
	delete tpl_ctx_rc_hexa;
	tpl_ctx_rc_hexa = new rcmfs_hexa_t(rc, nullptr, is_compressing);
	
	delete tpl_ctx_rc_literal;
	tpl_ctx_rc_literal = new rcmfs_literal_t(rc, nullptr, is_compressing);

	delete tpl_ctx_rc_literal_same;
	tpl_ctx_rc_literal_same = new rcmfs_literal_same_t(rc, nullptr, is_compressing);

	delete tpl_ctx_rc_literal_same_length;
	tpl_ctx_rc_literal_same_length = new rcmfs_literal_same_length_t(rc, nullptr, is_compressing);
	
	delete tpl_ctx_rc_plain;
	tpl_ctx_rc_plain = new rcmfs_plain_t(rc, nullptr, is_compressing);
}

// *****************************************************************
void CIDCoder::Restart()
{
	if (is_compressing)
		rce->Start();
	else
	{
		v_vios_io->RestartRead();
		rcd->Start();
	}
	ctx_flags = 0;
}

// *****************************************************************
void CIDCoder::Finish()
{
	if (is_compressing)
		rce->End();
	else
		rcd->End();
}

// *****************************************************************
void CIDCoder::Encode(bool plus_id, string id)
{
	if (header_mode == HeaderComprMode::None)
		return;
	else if (header_mode == HeaderComprMode::Original)
		compress_lossless(plus_id, id);
	else
		compress_instrument(plus_id, id);
}

// *****************************************************************
void CIDCoder::Decode(bool &plus_id, string& id)
{
	if (header_mode == HeaderComprMode::None)
		decompress_none(plus_id, id);
	else if (header_mode == HeaderComprMode::Original)
		decompress_lossless(plus_id, id);
	else
		decompress_instrument(plus_id, id);
}

// *****************************************************************
void CIDCoder::init_symbol_classes()
{
	fill_n(a_literal.begin(), 256, false);
	fill_n(a_numeric.begin(), 256, false);
	fill_n(a_hexadecimal.begin(), 256, false);

	for (int i = '0'; i <= '9'; ++i)
//		a_numeric[i] = a_literal[i] = a_hexadecimal[i] = true;
		a_literal[i] = a_hexadecimal[i] = true;

	for (int i = 'a'; i <= 'f'; ++i)
		a_hexadecimal[i] = true;

	for (int i = 'A'; i <= 'Z'; ++i)
		a_literal[i] = true;
	for (int i = 'a'; i <= 'z'; ++i)
		a_literal[i] = true;
	a_literal['@'] = true;
}

//*****************************************************************************************************
void CIDCoder::ResetReadPrev()
{
	id_prev.clear();
	size_prev = 0;
	v_tokens_cur.clear();
	v_tokens_prev.clear();
	v_deltas.clear();

	ctx_flags = 0;
}

// *****************************************************************
string CIDCoder::extract_instrument(const string &id)
{
	auto size = id.size();

	for (uint32_t i = 0; i < size; ++i)
		if (id[i] == '.' || id[i] == ' ' || id[i] == ':')
			return string(id.begin(), id.begin()+i);

	return id;
}

// *****************************************************************
uint32_t CIDCoder::tokenize(const string &id, v_tokens_t& v_tokens)
{
	uint32_t token_start_pos = 0;
//	token_type_t token_type = token_type_t::numeric;
	token_type_t token_type = token_type_t::literal;
	auto size = id.size();

	v_tokens.clear();

	for (uint32_t i = 0; i < size; ++i)
	{
		if (!is_literal(id[i]))
		{
			if (token_type == token_type_t::numeric && ((i - token_start_pos >= 11) || (i == token_start_pos)))
				token_type = token_type_t::literal;
			v_tokens.emplace_back(make_tuple(token_type, id[i], token_start_pos, i));

//			token_type = token_type_t::numeric;
			token_type = token_type_t::literal;
			token_start_pos = i + 1;
		}
/*		else if (!is_hexadecimal(id[i]))
			token_type = token_type_t::literal;
		else if (!is_numeric(id[i]))
			token_type = token_type_t::hexadecimal;*/
		else if (!is_numeric(id[i]))
			token_type = token_type_t::literal;
	}

	v_tokens.emplace_back(make_tuple(token_type, 0, token_start_pos, size));

/*	get<0>(v_tokens[16]) = token_type_t::numeric;
	get<0>(v_tokens[18]) = token_type_t::numeric;
	get<0>(v_tokens[21]) = token_type_t::numeric;
	get<0>(v_tokens[22]) = token_type_t::numeric;
	get<0>(v_tokens[24]) = token_type_t::numeric;*/

	return (uint32_t)v_tokens.size();
}

// *****************************************************************
void CIDCoder::compress_lossless(const bool plus_id, const string &id)
{
	auto n_tokens = tokenize(id, v_tokens_cur);

	auto rc_plus_id = find_rc_context(m_ctx_rc_plus_id, 0, tpl_ctx_rc_plus_id);
	rc_plus_id->Encode((int) plus_id);

	auto rc_flags = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
	if (token_types_same(v_tokens_cur, v_tokens_prev))
	{
		rc_flags->Encode(1);
		ctx_flags = ((ctx_flags << 1) + 1) & 0xff;

		for (uint32_t i = 0; i < n_tokens; ++i)
		{
			if (get<0>(v_tokens_cur[i]) == token_type_t::literal)
			{
				bool same_length = (get<3>(v_tokens_cur[i]) - get<2>(v_tokens_cur[i])) == (get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]));
				bool same = false;
				if (same_length)
					same = equal(id.begin() + get<2>(v_tokens_cur[i]), id.begin() + get<3>(v_tokens_cur[i]), id_prev.begin() + get<2>(v_tokens_prev[i]));

				auto rc_same = find_rc_context(m_ctx_rc_literal_same, i, tpl_ctx_rc_literal_same);

				if (same)
					rc_same->Encode(1);
				else
				{
					rc_same->Encode(0);
					auto rc_same_length = find_rc_context(m_ctx_rc_literal_same_length, i, tpl_ctx_rc_literal_same_length);
					if (same_length)
					{
						rc_same_length->Encode(1);
						bool prev_eq = true;
						for (uint32_t j = 0; j < get<3>(v_tokens_cur[i]) - get<2>(v_tokens_cur[i]); ++j)
						{
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + (1ll << 32) + j + (((context_t) i) << 40) + (prev_eq ? (1ull << 60) : 0ull), tpl_ctx_rc_literal);
							if (id[j + get<2>(v_tokens_cur[i])] == id_prev[j + get<2>(v_tokens_prev[i])])
							{
								rc->Encode(0);
//								prev_eq = true;
							}
							else
							{
								rc->Encode(id[j + get<2>(v_tokens_cur[i])]);
//								prev_eq = false;
							}
						}
					}
					else
					{
						rc_same_length->Encode(0);
						for (uint32_t j = get<2>(v_tokens_cur[i]); j < get<3>(v_tokens_cur[i]); ++j)
						{
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + j - get<2>(v_tokens_cur[i]) + (1ll << 32) + (((context_t)i) << 40), tpl_ctx_rc_literal);
							rc->Encode(id[j]);
						}

						auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + (get<3>(v_tokens_cur[i]) - get<2>(v_tokens_cur[i])) + (1ll << 32) + (((context_t)i) << 40), tpl_ctx_rc_literal);
						rc->Encode(0);
					}
				}
			}
			else
			{
				int64_t v_prev = get_int(id_prev, get<2>(v_tokens_prev[i]), get<3>(v_tokens_prev[i]));
				int64_t v_cur = get_int(id, get<2>(v_tokens_cur[i]), get<3>(v_tokens_cur[i]));
				int64_t delta = v_cur - v_prev;

				context_t ctx_type = (uint64_t)i << 40;
				ctx_type += (uint64_t)ilog2(abs(v_deltas[i])) << 31;
				ctx_type += (uint64_t)(v_deltas[i] < 0) << 30;
				auto rc_numeric_size = find_rc_context(m_ctx_rc_numeric_size, ctx_type, tpl_ctx_rc_numeric_size);
				auto rc_numeric_small = find_rc_context(m_ctx_rc_numeric_small, ctx_type, tpl_ctx_rc_numeric_small);

				v_deltas[i] = delta;

				if (delta >= -1 && delta <= 1)
				{
					rc_numeric_small->Encode((int)(delta + 1));
				}
				else
				{
					rc_numeric_small->Encode(3);
					int n_bytes = 0;
					if (delta >= -123 && delta <= 123)
						rc_numeric_size->Encode(((uint32_t)(delta + 123)) & 0xff);
					else if (delta > 0 && delta < 0x10000ll)
					{
						rc_numeric_size->Encode(247);
						n_bytes = 2;
						ctx_type += 0x10;
					}
					else if (delta > 0 && delta < 0x1000000ll)
					{
						rc_numeric_size->Encode(248);
						n_bytes = 3;
						ctx_type += 0x20;
					}
					else if (delta > 0 && delta < 0x100000000)
					{
						rc_numeric_size->Encode(249);
						n_bytes = 4;
						ctx_type += 0x30;
					}
					else if (delta > 0)
					{
						rc_numeric_size->Encode(250);
						n_bytes = 8;
						ctx_type += 0x40;
					}
					else if (delta < 0 && delta > -0x10000ll)
					{
						rc_numeric_size->Encode(251);
						delta = -delta;
						n_bytes = 2;
						ctx_type += 0x50;
					}
					else if (delta < 0 && delta > -0x1000000ll)
					{
						rc_numeric_size->Encode(252);
						delta = -delta;
						n_bytes = 3;
						ctx_type += 0x60;
					}
					else if (delta < 0 && delta > -0x100000000ll)
					{
						rc_numeric_size->Encode(253);
						delta = -delta;
						n_bytes = 4;
						ctx_type += 0x70;
					}
					else if (delta < 0)
					{
						rc_numeric_size->Encode(254);
						delta = -delta;
						n_bytes = 8;
						ctx_type += 0x80;
					}

					for (int j = 0; j < n_bytes; ++j)
					{
						auto rc = find_rc_context(m_ctx_rc_numeric_size, ctx_type + j, tpl_ctx_rc_numeric_size);
						rc->Encode(((uint64_t)delta >> (8 * j)) & 0xff);
					}
				}
			}
		}
	}
	else
	{
		rc_flags->Encode(0);
		ctx_flags = ((ctx_flags << 1) + 0) & 0xff;

		auto size = id.size();

		// Encode id plain
		for (uint32_t i = 0; i < size; ++i)
		{
			auto rc = find_rc_context(m_ctx_rc_plain, i, tpl_ctx_rc_plain);
			rc->Encode(id[i]);
		}

		auto rc = find_rc_context(m_ctx_rc_plain, size, tpl_ctx_rc_plain);
		rc->Encode(0);

		v_deltas.clear();
		v_deltas.resize(v_tokens_cur.size(), 0u);
	}

	v_tokens_prev.swap(v_tokens_cur);

	id_prev = id;
	size_prev = static_cast<uint32_t>(id.size());
}

// *****************************************************************
void CIDCoder::compress_instrument(const bool plus_id, const string& id)
{

}

// *****************************************************************
void CIDCoder::decompress_none(bool &plus_id, string& id)
{
	id = "@";
}

// *****************************************************************
void CIDCoder::decompress_lossless(bool& plus_id, string& id)
{
	id.clear();

	auto rc_plus_id = find_rc_context(m_ctx_rc_plus_id, 0, tpl_ctx_rc_plus_id);
	plus_id = (bool) rc_plus_id->Decode();

	auto rc_flags = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
	
	if (rc_flags->Decode() == 1)	// tokens of the same type
	{
		ctx_flags = ((ctx_flags << 1) + 1) & 0xff;
		auto n_tokens = v_tokens_prev.size();

		for (uint32_t i = 0; i < n_tokens; ++i)
		{
			if (get<0>(v_tokens_prev[i]) == token_type_t::literal)
			{
				auto rc_same = find_rc_context(m_ctx_rc_literal_same, i, tpl_ctx_rc_literal_same);

				if (rc_same->Decode() == 1)		// same token
					id.append(id_prev, get<2>(v_tokens_prev[i]), get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]));
				else
				{
					auto rc_same_length = find_rc_context(m_ctx_rc_literal_same_length, i, tpl_ctx_rc_literal_same_length);
					
					if (rc_same_length->Decode() == 1)	// same lengths
					{
						bool prev_eq = true;
						for (uint32_t j = 0; j < get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]); ++j)
						{
							// !!! Fix mixing literal contexts
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + (1ll << 32) + j + (((context_t)i) << 40) + (prev_eq ? (1ull << 60) : 0ull), tpl_ctx_rc_literal);
							auto d = rc->Decode();

							if (d == 0)
							{
								id.push_back(id_prev[j + get<2>(v_tokens_prev[i])]);
//								prev_eq = true;
							}
							else
							{
								id.push_back(d);
//								prev_eq = false;
							}
						}
					}
					else
					{
						for(uint32_t j = 0; ; ++j)
						{
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + j + (1ll << 32) + (((context_t)i) << 40), tpl_ctx_rc_literal);
							auto d = rc->Decode();
							if (d)
								id.push_back(d);
							else
								break;
						}
					}
				}
			}
			else
			{
				int64_t v_prev = get_int(id_prev, get<2>(v_tokens_prev[i]), get<3>(v_tokens_prev[i]));
				int64_t v_cur;
				int64_t delta;

				context_t ctx_type = (uint64_t)i << 40;
				ctx_type += (uint64_t)ilog2(abs(v_deltas[i])) << 31;
				ctx_type += (uint64_t)(v_deltas[i] < 0) << 30;
				auto rc_numeric_size = find_rc_context(m_ctx_rc_numeric_size, ctx_type, tpl_ctx_rc_numeric_size);
				auto rc_numeric_small = find_rc_context(m_ctx_rc_numeric_small, ctx_type, tpl_ctx_rc_numeric_small);
				
				int64_t d = rc_numeric_small->Decode();

				if (d < 3)
				{
					delta = d - 1ll;
				}
				else
				{
					int n_bytes = 0;
					delta = 0;

					d = rc_numeric_size->Decode();
					if (d <= 246)
					{
						delta = d - 123ll;
					}
					else if (d == 247)
					{
						n_bytes = 2;
						ctx_type += 0x10;
					}
					else if (d == 248)
					{
						n_bytes = 3;
						ctx_type += 0x20;
					}
					else if (d == 249)
					{
						n_bytes = 4;
						ctx_type += 0x30;
					}
					else if (d == 250)
					{
						n_bytes = 8;
						ctx_type += 0x40;
					}
					else if (d == 251)
					{
						n_bytes = 2;
						ctx_type += 0x50;
					}
					else if (d == 252)
					{
						n_bytes = 3;
						ctx_type += 0x60;
					}
					else if (d == 253)
					{
						n_bytes = 4;
						ctx_type += 0x70;
					}
					else if (d == 254)
					{
						n_bytes = 8;
						ctx_type += 0x80;
					}

					for (int j = 0; j < n_bytes; ++j)
					{
						auto rc = find_rc_context(m_ctx_rc_numeric_size, ctx_type + j, tpl_ctx_rc_numeric_size);
						delta += (int64_t) rc->Decode() << (8 * j);
					}
					
					if (d >= 251)
						delta = -delta;
				}

				v_deltas[i] = delta;

				v_cur = delta + v_prev;

				store_int(id, v_cur);
			}

			id.push_back(get<1>(v_tokens_prev[i]));
		}

		if (id.back() == 0)
			id.pop_back();

		tokenize(id, v_tokens_cur);
	}
	else
	{
		ctx_flags = ((ctx_flags << 1) + 0) & 0xff;

		// Decode id plain
		for (uint32_t i = 0; ; ++i)
		{
			auto rc = find_rc_context(m_ctx_rc_plain, i, tpl_ctx_rc_plain);
			id.push_back(rc->Decode());

			if (id.back() == 0)
				break;
			if (id.back() == 0xA)
			{
//				id.push_back(0);
				break;
			}
		}

		if (id.back() == 0)
			id.pop_back();

		tokenize(id, v_tokens_cur);

		v_deltas.clear();
		v_deltas.resize(v_tokens_cur.size(), 0u);
	}

	v_tokens_prev.swap(v_tokens_cur);
	id_prev = id;
//	size_prev = size;
}

// *****************************************************************
void CIDCoder::decompress_instrument(bool& plus_id, string& id)
{

}
}
// EOF
