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

#include <unordered_map>
#include "basic_coder.h"
#include "ref_reads_accepter.h"
#include "reference_reads.h"

namespace entropy_coder
{

// *****************************************************************
class CDNACoder : public CBasicCoder
{
	bool verbose;
	const uint32_t max_enc_anchor_len = 24;
	const uint32_t max_enc_skip_len = 256;
	//	const int max_enc_skip_len = 64;

	const int subst_to_code[4][4] = { {1, 0, 0, 0}, {2, 2, 1, 1}, {3, 3, 3, 2}, {3, 3, 3, 3} };

	constexpr uint64_t bm_type(tuple_types t) { return 1ull << (int)t; };

	const uint64_t exc_mask_anchor = bm_type(tuple_types::anchor) + bm_type(tuple_types::match);
	const uint64_t exc_mask_match = bm_type(tuple_types::anchor);
	const uint64_t exc_mask_deletion = bm_type(tuple_types::skip);
	const uint64_t exc_mask_skip = bm_type(tuple_types::deletion) + bm_type(tuple_types::skip);
	const uint64_t exc_mask_alt_id = bm_type(tuple_types::alt_id) + bm_type(tuple_types::main_ref);
	const uint64_t exc_mask_main_ref = bm_type(tuple_types::alt_id) + bm_type(tuple_types::main_ref);

	using rcmfs_rev_comp_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 15, 1>;
	using rcmfs_read_type_t = CRangeCoderModelFixedSize<CVectorIOStream, 3, 1 << 15, 1>;
	using rcmfs_seen_read_id_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 15, 1>;
	using rcmfs_read_len_no_bits_t = CRangeCoderModelFixedSize<CVectorIOStream, 32, 1 << 18, 8>;
	using rcmfs_read_len_data_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 18, 8>;
	using rcmfs_symbols_t = CRangeCoderModelFixedSize<CVectorIOStream, 4, 1 << 10, 1>;
	using rcmfs_symbols_with_Ns_t = CRangeCoderModelFixedSize<CVectorIOStream, 5, 1 << 10, 1>;
	using rcmfs_read_id_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 13, 1>;
	using rcmfs_skip_len_distant_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 15, 1>;
	using rcmfs_tuple_type_t = CRangeCoderModelFixedSize<CVectorIOStream, 8, 1 << 15, 1>;
	using rcm_read_id_short_t = CRangeCoderModel<CVectorIOStream, 1 << 13, 1>;
	using rcm_anchor_len_t = CRangeCoderModel<CVectorIOStream, 1 << 15, 1>;
	using rcm_skip_len_local_t = CRangeCoderModel<CVectorIOStream, 1 << 15, 1>;

	using ctx_vece_rev_comp_t = CContextVecEmb<rcmfs_rev_comp_t>;
	using ctx_vece_read_type_t = CContextVecEmb<rcmfs_read_type_t>;
	using ctx_vece_seen_read_id_t = CContextVecEmb<rcmfs_seen_read_id_t>;
	using ctx_vece_read_len_no_bits_t = CContextVecEmb<rcmfs_read_len_no_bits_t>;
	using ctx_vec_read_len_data_t = CContextVec<rcmfs_read_len_data_t>;
	using ctx_vece_symbols_t = CContextVecEmb<rcmfs_symbols_t>;
	using ctx_vece_symbols_with_Ns_t = CContextVecEmb<rcmfs_symbols_with_Ns_t>;
	using ctx_vec_read_id_t = CContextVec<rcmfs_read_id_t>;
	using ctx_vec_skip_len_distant_t = CContextVec<rcmfs_skip_len_distant_t>;
	using ctx_vece_tuple_type_t = CContextVecEmb<rcmfs_tuple_type_t>;
	using ctx_vece_read_id_short_t = CContextVecEmb<rcm_read_id_short_t>;
	using ctx_vece_anchor_len_t = CContextVecEmb<rcm_anchor_len_t>;
	using ctx_vece_skip_len_local_t = CContextVecEmb<rcm_skip_len_local_t>;

	CReferenceReads& ref_reads;

	const context_t ctx_mask_read_type = 0xff;
	const context_t ctx_mask_rev_comp = 0xf;

/*	const context_t ctx_mask_tuple_type = 0xffffull;
	const int no_symbols_in_mask = 8;
	const context_t ctx_mask_symbol = (1ull << (2 * no_symbols_in_mask)) - 1ull;*/
	int no_tuples_in_mask;
	context_t ctx_mask_tuple_type;
	int no_symbols_in_mask;
	context_t ctx_mask_symbol;

	unordered_map<int, bool> uo_rev_comp;

	int compression_level;
	uint64_t input_stream_size;
	int max_no_alt_refs;
	int cur_read_id;
	int cur_ref_delta;

	context_t ctx_read_type;
	context_t ctx_rev_comp;
	context_t ctx_tuple_type;
	context_t ctx_symbol;

	ctx_vece_rev_comp_t m_ctx_rc_rev_comp;
	ctx_vece_read_type_t m_ctx_rc_read_type;
	ctx_vece_seen_read_id_t m_ctx_rc_seen_read_id;
	ctx_vece_symbols_t m_ctx_rc_symbols;
	ctx_vece_symbols_with_Ns_t m_ctx_rc_symbols_with_Ns;
	ctx_vece_tuple_type_t m_ctx_rc_tuple_type;
	ctx_vece_read_len_no_bits_t m_ctx_rc_read_len_no_bits;
	ctx_vec_read_len_data_t m_ctx_rc_read_len_data;

	ctx_vec_read_id_t m_ctx_rc_read_id;
	ctx_vece_read_id_short_t m_ctx_rc_read_id_short;
	ctx_vece_anchor_len_t m_ctx_rc_anchor_len;
	ctx_vece_skip_len_local_t m_ctx_rc_skip_len_local;
	ctx_vec_skip_len_distant_t m_ctx_rc_skip_len_distant;

	rcmfs_rev_comp_t* tpl_ctx_rc_rev_comp;
	rcmfs_read_type_t* tpl_ctx_rc_read_type;
	rcmfs_seen_read_id_t* tpl_ctx_rc_seen_read_id;
	rcmfs_read_len_no_bits_t* tpl_ctx_rc_read_len_no_bits;
	rcmfs_read_len_data_t* tpl_ctx_rc_read_len_data;

	rcmfs_symbols_t* tpl_ctx_rc_symbols;
	rcmfs_symbols_with_Ns_t* tpl_ctx_rc_symbols_with_Ns;
	rcmfs_tuple_type_t* tpl_ctx_rc_tuple_type;

	rcmfs_read_id_t* tpl_ctx_rc_read_id;
	rcm_read_id_short_t* tpl_ctx_rc_read_id_short;
	rcm_anchor_len_t* tpl_ctx_rc_anchor_len;
	rcm_skip_len_local_t* tpl_ctx_rc_skip_len_local;
	rcmfs_skip_len_distant_t* tpl_ctx_rc_skip_len_distant;

	void encode_rev_comp_flag(int read_id, bool is_rev_comp);
	void encode_read_flag(tuple_types type);
	void encode_read_len(uint32_t len);
	void encode_read_id(int id, int cur_read_id);
	bool encode_alt_read_id(int id, int cur_read_id, map<int, int>& m_alt_ids);
	void encode_symbol_plain(uint32_t symbol);
	void encode_symbol_plain_with_Ns(uint32_t symbol);
//	void __declspec(noinline) encode_tuple_type(tuple_types type, int ref_symbol, tuple_types last_flag, bool is_first_tuple = false);
	void encode_tuple_type(tuple_types type, uint32_t ref_symbol, tuple_types last_flag, bool is_first_tuple = false);
	void encode_anchor_len(uint32_t len);
	void encode_insertion(uint32_t symbol, uint32_t base);
	void encode_substitution(uint32_t symbol, uint32_t base);
	void encode_skip_len(uint32_t len, bool local);

	void prefetch_insertion(uint32_t base);
	void prefetch_tuple_type(tuple_types type, uint32_t ref_symbol, tuple_types last_flag, bool is_first_tuple = false);


	bool decode_rev_comp_flag(int read_id);
	tuple_types decode_read_flag();
	uint32_t decode_read_len();
	uint32_t decode_symbol_plain();
	uint32_t decode_symbol_plain_with_Ns();
	int decode_read_id(int cur_read_id);
	int decode_alt_read_id(int cur_read_id, vector<int>& v_alt_ids);
	tuple_types decode_tuple_type(uint32_t ref_symbol, tuple_types last_flag, bool is_first_tuple = false);
	uint32_t decode_anchor_len();
	uint32_t decode_insertion(uint32_t base);
	uint32_t decode_substitution(uint32_t base);
	uint32_t decode_skip_len(bool local);

/*	const int subst_to_code(int c, int base)
	{
		if (c == 0 && base == 0)	return 1;
		if (c == 1 && base == 0)	return 2;
		if (c == 2 && base == 0)	return 3;
		if (c == 0 && base == 1)	return 0;
		if (c == 1 && base == 1)	return 2;
		if (c == 2 && base == 1)	return 3;
		if (c == 0 && base == 2)	return 0;
		if (c == 1 && base == 2)	return 1;
		if (c == 2 && base == 2)	return 3;
		return c;
	}*/

public:
	// Remarks:
	//		* symols stored as 0, 1, 2, 3 (no ACGT)
	//		* at the end of a read there should ba a guard: 255 (only for more efficient encoding)	
	CDNACoder(CReferenceReads& _ref_reads, bool verbose) : CBasicCoder(), verbose(verbose), ref_reads(_ref_reads)
	{
		max_no_alt_refs = 1;

		tpl_ctx_rc_rev_comp = nullptr;
		tpl_ctx_rc_read_type = nullptr;
		tpl_ctx_rc_seen_read_id = nullptr;
		tpl_ctx_rc_read_len_no_bits = nullptr;
		tpl_ctx_rc_read_len_data = nullptr;
		tpl_ctx_rc_symbols = nullptr;
		tpl_ctx_rc_symbols_with_Ns = nullptr;
		tpl_ctx_rc_read_id = nullptr;
		tpl_ctx_rc_read_id_short = nullptr;
		tpl_ctx_rc_tuple_type = nullptr;
		tpl_ctx_rc_anchor_len = nullptr;
		tpl_ctx_rc_skip_len_local = nullptr;
		tpl_ctx_rc_skip_len_distant = nullptr;
	}

	~CDNACoder()
	{
		delete tpl_ctx_rc_rev_comp;
		delete tpl_ctx_rc_read_type;
		delete tpl_ctx_rc_seen_read_id;
		delete tpl_ctx_rc_read_len_no_bits;
		delete tpl_ctx_rc_read_len_data;
		delete tpl_ctx_rc_symbols;
		delete tpl_ctx_rc_symbols_with_Ns;
		delete tpl_ctx_rc_read_id;
		delete tpl_ctx_rc_read_id_short;
		delete tpl_ctx_rc_tuple_type;
		delete tpl_ctx_rc_anchor_len;
		delete tpl_ctx_rc_skip_len_local;
		delete tpl_ctx_rc_skip_len_distant;
		if(verbose)
		{
			cout << "Read type size: " << m_ctx_rc_read_type.get_size() << endl;
			cout << "Rev. comp. size: " << m_ctx_rc_rev_comp.get_size() << endl;
			cout << "Seen read id size: " << m_ctx_rc_seen_read_id.get_size() << endl;
			cout << "Symbols size: " << m_ctx_rc_symbols.get_size() << endl;
			cout << "Symbols with Ns size: " << m_ctx_rc_symbols_with_Ns.get_size() << endl;
			cout << "Read len no. bits size: " << m_ctx_rc_read_len_no_bits.get_size() << endl;
			cout << "Read len data size: " << m_ctx_rc_read_len_data.get_size() << endl;
			cout << "Read id size: " << m_ctx_rc_read_id.get_size() << endl;
			cout << "Read id short size: " << m_ctx_rc_read_id_short.get_size() << endl;
			cout << "Tuple type size: " << m_ctx_rc_tuple_type.get_size() << endl;
			cout << "Anchor len size: " << m_ctx_rc_anchor_len.get_size() << endl;
			cout << "Skip len local size: " << m_ctx_rc_skip_len_local.get_size() << endl;
			cout << "Skip len distant size: " << m_ctx_rc_skip_len_distant.get_size() << endl;
		}
	}

	// Initialization - required before (de)compression
	void Init(bool is_compressing, int _max_no_alt_refs, int _compression_level, uint64_t _input_stream_size, uint32_t start_read_id);

	// Compression restart without model cleanup
	void Restart();

	// Should be called after (de)compression
	void Finish();

	// Encoding read in an edit script form
	// Plain reads:
	//		(start_plain, 0, 0), (plain, symbol, 1), [(plain, symbol, 1)...]
	// Reads as edit script
	//		(start_es, ref_read_id, 0), data_tuple, [data_uple, ...]
	// Types of data tuples:
	//		(anchor, 0, length) - anchor and its length
	//		(match, 0, 1) - matches, that are not a part of anchor, given always as a single match (no RLE); match cannot be directly before/after anchor
	//		(insertion, symbol, 1) - insertion of symbol (no RLE)
	//		(deletion, 0, 1) - deletion of symbol (no RLE)
	//		(substitution, symbol, 1) - substitution of symbol (no RLE);
	//		(skip, 0, length) - long deletion (min. len. 10, but it may be parametrized); 
	//								before/after skip deletion is not allowed; 
	//		(alt_id, id_reada, 0) - switch to alternative read; 
	//		(main_ref, 0, 0) - return to main reference read
	// Remarks:
	//		* symbols represented as 0, 1, 2, 3 (no ACGT)
	void Encode(es_t es);

	// Read decoding
	void Decode(read_t& read, CRefReadsAccepter& ref_read_accepter, bool acceptAll);
};
} //namespace entropy_coder
// EOF