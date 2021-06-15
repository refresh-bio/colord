#include "dna_coder.h"

#include <cassert>

namespace entropy_coder
{
// *****************************************************************
void CDNACoder::Encode(es_t es)
{
	ctx_tuple_type = ctx_mask_tuple_type;
	ctx_symbol = ctx_mask_symbol;

	ctx_rev_comp = ctx_mask_rev_comp;
	uo_rev_comp.clear();

	tuple_types tuple_type = tuple_types::none;
	uint32_t tuple_val1{};
	uint32_t tuple_val2{};

	es.restart_reading();

	es.load(tuple_type, tuple_val1, tuple_val2);

	encode_read_flag(tuple_type);
	encode_read_len(static_cast<uint32_t>(es.size() - 1u));

	if (tuple_type == tuple_types::start_plain)
	{
		while (es.load(tuple_type, tuple_val1, tuple_val2))
			encode_symbol_plain(tuple_val1);

		++cur_read_id;

		return;
	}

	if (tuple_type == tuple_types::start_plain_with_Ns)
	{
		while (es.load(tuple_type, tuple_val1, tuple_val2))
			encode_symbol_plain_with_Ns(tuple_val1);

		++cur_read_id;

		return;
	}

	int ref_read_id = -1;
	int alt_read_id = -1;

	map<int, int> m_alt_ids;
	map<int, int> m_alt_pos;
	map<int, read_t> m_alt_read;
	
	int ref_pos = 0;
	int alt_pos = 0;
	int cur_pos = 0;
	tuple_types last_tuple_type = tuple_types::none;

	cur_ref_delta = 0;

	bool is_main_ref = true;

	encode_read_id(tuple_val1, cur_read_id);
	ref_read_id = tuple_val1;

	encode_rev_comp_flag(tuple_val1, (bool) tuple_val2);
	read_t ref_read = ref_reads.GetRefRead(ref_read_id, tuple_val2);

	read_t alt_ref_read;
	bool is_first_tuple = true;

	tuple_types last_flag = tuple_types::none;

	while(es.load(tuple_type, tuple_val1, tuple_val2))
	{
		int ref_symbol = is_main_ref ? ref_read[ref_pos] : alt_ref_read[alt_pos];

//		prefetch_insertion(ref_symbol);
		encode_tuple_type(tuple_type, ref_symbol, last_flag, is_first_tuple);
		is_first_tuple = false;
		last_flag = tuple_type;

		if (tuple_type == tuple_types::alt_id)
		{
			// Remember last position in reference read
			if (!is_main_ref)
				m_alt_pos[alt_read_id] = alt_pos;

			bool is_new_alt = encode_alt_read_id(tuple_val1, cur_read_id, m_alt_ids);
			alt_read_id = tuple_val1;

			if (is_new_alt)
				m_alt_read[alt_read_id] = ref_reads.GetRefRead(alt_read_id, tuple_val2);
			encode_rev_comp_flag(tuple_val1, (bool) tuple_val2);

			alt_ref_read = m_alt_read[alt_read_id];
			alt_pos = 0;
			is_main_ref = false;
			//			no_alt_encoded_symbols = 0;
			cur_ref_delta = 0;
		}
		else if (tuple_type == tuple_types::anchor)
		{
			encode_anchor_len(tuple_val2);
			is_main_ref ? (ref_pos += tuple_val2) : (alt_pos += tuple_val2);
			cur_pos += tuple_val2;

			for (int i = no_symbols_in_mask; i > 0; --i)
				ctx_symbol = (ctx_symbol << 2) + (is_main_ref ? ref_read : alt_ref_read)[is_main_ref ? ref_pos - i : alt_pos - i]; 
			ctx_symbol &= ctx_mask_symbol;

			/*			if(is_main_ref)
							last_anchor_in_ref_end_pos = cur_pos;*/
			cur_ref_delta = 0;
		}
		else if (tuple_type == tuple_types::match)
		{
			//assert(v_reads[cur_read_id][cur_pos] == ref_symbol);

			ctx_symbol = ((ctx_symbol << 2) + (ref_symbol)) & ctx_mask_symbol;
			is_main_ref ? ++ref_pos : ++alt_pos;
			++cur_pos;
		}
		else if (tuple_type == tuple_types::insertion)
		{
			encode_insertion(tuple_val1, ref_symbol);
			ctx_symbol = ((ctx_symbol << 2) + (tuple_val1)) & ctx_mask_symbol;
			++cur_pos;

			++cur_ref_delta;
		}
		else if (tuple_type == tuple_types::deletion)
		{
			is_main_ref ? ++ref_pos : ++alt_pos;

			--cur_ref_delta;
		}
		else if (tuple_type == tuple_types::substitution)
		{
			int symbol = subst_to_code[tuple_val1][ref_symbol];
	
			encode_substitution(symbol, ref_symbol);

			ctx_symbol = ((ctx_symbol << 2) + (symbol)) & ctx_mask_symbol;
			is_main_ref ? ++ref_pos : ++alt_pos;
			++cur_pos;
		}
		else if (tuple_type == tuple_types::skip)
		{
			int skip_len = tuple_val2;

			/*			if (!is_main_ref && last_tuple_type == tuple_types::alt_id)
						{
							skip_len -= alt_pos;
						}*/

						/*			if (no_alt_encoded_symbols)
									{
										if (skip_len < no_alt_encoded_symbols)
											skip_len = 2 * (no_alt_encoded_symbols - skip_len) - 1;
										else if (skip_len - no_alt_encoded_symbols < no_alt_encoded_symbols)
											skip_len = 2 * (skip_len - no_alt_encoded_symbols);
										else
											skip_len += 2 * no_alt_encoded_symbols;
									}
						*/
			//no_alt_encoded_symbols = 0;

			cur_ref_delta -= skip_len;
			//			encode_skip_len(p->no_rep, last_tuple_type != tuple_types::alt_id && last_tuple_type != tuple_types::none);
			if (!is_main_ref && last_tuple_type == tuple_types::alt_id)
			{
				int mod_skip_len = skip_len - m_alt_pos[alt_read_id];

				//				if(skip_len != mod_skip_len)
				//					cout << skip_len << " " << mod_skip_len << endl;

				if (mod_skip_len > 0)
					encode_skip_len(static_cast<uint32_t>(mod_skip_len), false);
				else
				{
					//					cout << skip_len << " " << mod_skip_len << endl;
					encode_skip_len(0, false);
					encode_skip_len(static_cast<uint32_t>(-mod_skip_len), false);
				}
			}
			else
				encode_skip_len(static_cast<uint32_t>(skip_len), last_tuple_type != tuple_types::alt_id && last_tuple_type != tuple_types::none);

			is_main_ref ? (ref_pos += tuple_val2) : (alt_pos += tuple_val2);
		}
		else if (tuple_type == tuple_types::main_ref)
		{
			is_main_ref = true;

			// Remember last position in reference read
			m_alt_pos[alt_read_id] = alt_pos;

			/*			if (last_anchor_in_ref_end_pos)
						{
							no_alt_encoded_symbols = cur_pos - last_anchor_in_ref_end_pos;
							last_anchor_in_ref_end_pos = 0;
						}*/
			cur_ref_delta = 0;
		}
		else
			cerr << "Unknown flag: " << (int)tuple_type << endl;

		last_tuple_type = tuple_type;
	}

	++cur_read_id;
}

// *****************************************************************
void CDNACoder::Decode(read_t& read, CRefReadsAccepter& ref_read_accepter, bool acceptAll)
{
	read.clear();

	ctx_tuple_type = ctx_mask_tuple_type;
	ctx_symbol = ctx_mask_symbol;

	ctx_rev_comp = ctx_mask_rev_comp;
	uo_rev_comp.clear();

	auto read_flag = decode_read_flag();
	auto read_len = decode_read_len();

	read.reserve(read_len + 1u);

	bool acceptReadAsRef = read_flag != tuple_types::start_plain_with_Ns;
	if (!acceptAll)
	{
		acceptReadAsRef &= ref_read_accepter.ShouldAddToReference(cur_read_id);
	}

	if (read_flag == tuple_types::start_plain)
	{
		for (uint32_t i = 0; i < read_len; ++i)
			read.emplace_back(decode_symbol_plain());

		read.emplace_back(255);
		++cur_read_id;

		if (acceptReadAsRef)
			ref_reads.Add(read);

		return;
	}

	if (read_flag == tuple_types::start_plain_with_Ns)
	{
		for (uint32_t i = 0; i < read_len; ++i)
			read.emplace_back(decode_symbol_plain_with_Ns());

		read.emplace_back(255);
		++cur_read_id;

		if (acceptReadAsRef)
			ref_reads.Add(read);

		return;
	}

	bool is_main_ref = true;
	int ref_read_id = decode_read_id(cur_read_id);

	bool is_rev_comp = decode_rev_comp_flag(ref_read_id);
	read_t ref_read = ref_reads.GetRefRead(ref_read_id, is_rev_comp);

	read_t alt_ref_read;

	int alt_read_id = -1;
	vector<int> v_alt_ids;
	map<int, int> m_alt_pos;

	tuple_types last_tuple_type = tuple_types::none;
	read_t read_without_flags;

	int ref_pos = 0;
	int alt_pos = 0;
	int cur_pos = 0;

	cur_ref_delta = 0;
	uint8_t effective_base_flag_match = (compression_level > 1) ? base_flag_match : 0u;
	uint8_t effective_base_flag_anchor = (compression_level > 1) ? base_flag_anchor : 0u;

	for (uint32_t cur_es_len = 0; cur_es_len < read_len; ++cur_es_len)
	{
		int ref_symbol = 0;

		if (is_main_ref && ref_pos < static_cast<int>(ref_read.size()))
			ref_symbol = ref_read[ref_pos];
		else if (!is_main_ref && alt_pos < static_cast<int>(alt_ref_read.size()))
			ref_symbol = alt_ref_read[alt_pos];
		//		= is_main_ref ? v_reads[ref_read_id][ref_pos] : v_reads[alt_read_id][alt_pos];

		auto t = decode_tuple_type(ref_symbol, last_tuple_type, cur_es_len == 0);

		if (t == tuple_types::alt_id)
		{
			// Remember last position in reference read
			if (!is_main_ref)
				m_alt_pos[alt_read_id] = alt_pos;

			alt_read_id = decode_alt_read_id(cur_read_id, v_alt_ids);

			bool is_rev_comp = decode_rev_comp_flag(alt_read_id);
			alt_ref_read = ref_reads.GetRefRead(alt_read_id, is_rev_comp);

			alt_pos = 0;
			is_main_ref = false;
			cur_ref_delta = 0;
		}
		else if (t == tuple_types::anchor)
		{
			uint32_t anchor_len = decode_anchor_len();

			if (is_main_ref)
				for (uint32_t i = 0; i < anchor_len; ++i)
				{
					read.emplace_back(ref_read[ref_pos + i] | effective_base_flag_anchor);
					read_without_flags.emplace_back(ref_read[ref_pos + i]);
				}
			else
				for (uint32_t i = 0; i < anchor_len; ++i)
				{
					read.emplace_back(alt_ref_read[alt_pos + i] | effective_base_flag_anchor);
					read_without_flags.emplace_back(alt_ref_read[alt_pos + i]);
				}

			is_main_ref ? (ref_pos += anchor_len) : (alt_pos += anchor_len);
			cur_pos += anchor_len;
			cur_ref_delta = 0;

			for (int i = no_symbols_in_mask; i > 0; --i)
				ctx_symbol = (ctx_symbol << 2) + read_without_flags[cur_pos - i];
			ctx_symbol &= ctx_mask_symbol;
		}
		else if (t == tuple_types::match)
		{
			read.emplace_back(((uint8_t) ref_symbol) | effective_base_flag_match);
			read_without_flags.emplace_back(((uint8_t) ref_symbol));

			ctx_symbol = ((ctx_symbol << 2) + (ref_symbol)) & ctx_mask_symbol;
			is_main_ref ? ++ref_pos : ++alt_pos;
			++cur_pos;
		}
		else if (t == tuple_types::insertion)
		{
			uint32_t x = decode_insertion(ref_symbol);

			read.emplace_back(x);
			read_without_flags.emplace_back(x);

			ctx_symbol = ((ctx_symbol << 2) + (x)) & ctx_mask_symbol;
			++cur_pos;
			++cur_ref_delta;
		}
		else if (t == tuple_types::deletion)
		{
			is_main_ref ? ++ref_pos : ++alt_pos;
			--cur_ref_delta;
		}
		else if (t == tuple_types::substitution)
		{
			uint32_t symbol = decode_substitution(ref_symbol);

			read.emplace_back(symbol);
			read_without_flags.emplace_back(symbol);

			ctx_symbol = ((ctx_symbol << 2) + (symbol)) & ctx_mask_symbol;
			is_main_ref ? ++ref_pos : ++alt_pos;
			++cur_pos;
		}
		else if (t == tuple_types::skip)
		{
			//			int skip_len = decode_skip_len(last_tuple_type != tuple_types::alt_id && last_tuple_type != tuple_types::none);
			int skip_len;

			if (!is_main_ref && last_tuple_type == tuple_types::alt_id)
			{
				int mod_skip_len = static_cast<int>(decode_skip_len(false));

				if (mod_skip_len == 0)
					mod_skip_len = -static_cast<int>(decode_skip_len(false));

				skip_len = mod_skip_len + m_alt_pos[alt_read_id];

				//cout << skip_len << " " << mod_skip_len << endl;
			}
			else
				skip_len = decode_skip_len(last_tuple_type != tuple_types::alt_id && last_tuple_type != tuple_types::none);

			cur_ref_delta -= skip_len;
			is_main_ref ? (ref_pos += skip_len) : (alt_pos += skip_len);
		}
		else if (t == tuple_types::main_ref)
		{
		// Remember last position in reference read
			m_alt_pos[alt_read_id] = alt_pos;

			is_main_ref = true;
			cur_ref_delta = 0;
		}
		else
			cerr << "Unknown flag: " << (int)t << endl;

		last_tuple_type = t;
	}

	read.emplace_back(255);
	read_without_flags.emplace_back(255);

	if (acceptReadAsRef)
		ref_reads.Add(read_without_flags);

	++cur_read_id;
}

// *****************************************************************
void CDNACoder::encode_read_flag(tuple_types type)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	auto p = find_rce_context(m_ctx_rc_read_type, ctx_read_type, tpl_ctx_rc_read_type);

	uint32_t flag;

	if (type == tuple_types::start_plain)
		flag = 0;
	else if (type == tuple_types::start_plain_with_Ns)
		flag = 1;
	else
		flag = 2;

	p->Encode(flag);

	ctx_read_type <<= 2;
	ctx_read_type += flag;
	ctx_read_type &= ctx_mask_read_type;
}

// *****************************************************************
tuple_types CDNACoder::decode_read_flag()
{
	auto p = find_rce_context(m_ctx_rc_read_type, ctx_read_type, tpl_ctx_rc_read_type);

	auto flag = p->Decode();

	tuple_types type;
	
	if (flag == 0)
		type = tuple_types::start_plain;
	else if(flag == 1)
		type = tuple_types::start_plain_with_Ns;
	else
		type = tuple_types::start_es;

	ctx_read_type <<= 2;
	ctx_read_type += flag;
	ctx_read_type &= ctx_mask_read_type;

	return type;
}

// *****************************************************************
void CDNACoder::encode_rev_comp_flag(int read_id, bool is_rev_comp)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	// Read orientation already stored for the read_id
	if (uo_rev_comp.count(read_id))
		return;

	auto p = find_rce_context(m_ctx_rc_rev_comp, ctx_rev_comp, tpl_ctx_rc_rev_comp);

	p->Encode((int) is_rev_comp);

	uo_rev_comp[read_id] = is_rev_comp;

	ctx_rev_comp <<= 2;
	ctx_rev_comp += (int) is_rev_comp;
	ctx_rev_comp &= ctx_mask_rev_comp;
}

// *****************************************************************
bool CDNACoder::decode_rev_comp_flag(int read_id)
{
	// Read orientation already stored for the read_id
	
	auto q = uo_rev_comp.find(read_id);
		
	if (q != uo_rev_comp.end())
		return q->second;

	auto p = find_rce_context(m_ctx_rc_rev_comp, ctx_rev_comp, tpl_ctx_rc_rev_comp);

	auto flag = p->Decode();

	uo_rev_comp[read_id] = (bool) flag;

	ctx_rev_comp <<= 2;
	ctx_rev_comp += flag;
	ctx_rev_comp &= ctx_mask_rev_comp;

	return (bool) flag;
}

// *****************************************************************
void CDNACoder::encode_read_id(int id, int cur_read_id)
{
	int n = static_cast<int>(no_bytes(cur_read_id));

#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	for (int i = n - 1; i >= 0; --i)
	{
		context_t add = (i == n - 2) ? ((id >> (8 * (n - 1))) & 0xff) : 0;
		auto p = find_rc_context(m_ctx_rc_read_id, i + (add << 3), tpl_ctx_rc_read_id);

		p->Encode((id >> (8 * i)) & 0xff);
	}
}

// *****************************************************************
int CDNACoder::decode_read_id(int cur_read_id)
{
	int n = static_cast<int>(no_bytes(cur_read_id));
	uint32_t id = 0;

	for (int i = n - 1; i >= 0; --i)
	{
		context_t add = (i == n - 2) ? id : 0;
		auto p = find_rc_context(m_ctx_rc_read_id, i + (add << 3), tpl_ctx_rc_read_id);

		id <<= 8;
		id += p->Decode();
	}

	return static_cast<int>(id);
}

// *****************************************************************
bool CDNACoder::encode_alt_read_id(int id, int cur_read_id, map<int, int>& m_alt_ids)
{
	if (m_alt_ids.empty())
	{
		encode_read_id(id, cur_read_id);
		m_alt_ids[id] = 0;

		return true;
	}

#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return true;
#endif

	bool r = false;
	int alt_read_id_short = -1;
	int no_seen_ids = static_cast<int>(m_alt_ids.size());

	auto q = m_alt_ids.find(id);
	if (q == m_alt_ids.end())
	{
		int new_id_short = static_cast<int>(m_alt_ids.size());
		m_alt_ids[id] = new_id_short;
		r = true;
	}
	else
		alt_read_id_short = q->second;

	context_t ctx_seen_read_id = no_seen_ids;
	auto p_seen_read_id = find_rce_context(m_ctx_rc_seen_read_id, ctx_seen_read_id, tpl_ctx_rc_seen_read_id);
	p_seen_read_id->Encode(alt_read_id_short >= 0);

	if (alt_read_id_short < 0)
		encode_read_id(id, cur_read_id);
	else
	{
		context_t ctx_read_id_short = no_seen_ids;
		auto p_read_id_short = find_rce_context(m_ctx_rc_read_id_short, ctx_read_id_short, tpl_ctx_rc_read_id_short);
		p_read_id_short->Encode(alt_read_id_short);
	}

	return r;
}

// *****************************************************************
int CDNACoder::decode_alt_read_id(int cur_read_id, vector<int>& v_alt_ids)
{
	if (v_alt_ids.empty())
	{
		auto id = decode_read_id(cur_read_id);
		v_alt_ids.push_back(id);

		return id;
	}

	int no_seen_ids = static_cast<int>(v_alt_ids.size());

	context_t ctx_seen_read_id = no_seen_ids;
	auto p_seen_read_id = find_rce_context(m_ctx_rc_seen_read_id, ctx_seen_read_id, tpl_ctx_rc_seen_read_id);

	if (!p_seen_read_id->Decode())
	{
		int id = decode_read_id(cur_read_id);
		v_alt_ids.push_back(id);

		return id;
	}
	else
	{
		context_t ctx_read_id_short = no_seen_ids;
		auto p_read_id_short = find_rce_context(m_ctx_rc_read_id_short, ctx_read_id_short, tpl_ctx_rc_read_id_short);

		return v_alt_ids[p_read_id_short->Decode()];
	}
}

// *****************************************************************
//void __declspec(noinline) CDNACoder::encode_tuple_type(tuple_types type, int ref_symbol, tuple_types last_flag, bool is_first_tuple)
void CDNACoder::encode_tuple_type(tuple_types type, uint32_t ref_symbol, tuple_types last_flag, bool is_first_tuple)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	context_t ctx = 0;
	uint32_t shift = 0;

	ctx = ctx_tuple_type;
	shift += 3 * no_tuples_in_mask;

	ctx += ((ctx_symbol & 0xfull) << shift);
	shift += 4;

	ctx += ((uint64_t)ref_symbol) << shift;
	shift += 2;

/*	if (cur_ref_delta < -100)
		ctx += 5ull << shift;
	else if (cur_ref_delta < -10)
		ctx += 1ull << shift;
	else if (cur_ref_delta < -1)
		ctx += 2ull << shift;
	else if (cur_ref_delta > 100)
		ctx += 6ull << shift;
	else if (cur_ref_delta > 10)
		ctx += 3ull << shift;
	else if (cur_ref_delta > 1)
		ctx += 4ull << shift;
		*/
	if (cur_ref_delta < -10)
		ctx += 1ull << shift;
	else if (cur_ref_delta < -1)
		ctx += 2ull << shift;
	else if (cur_ref_delta > 10)
		ctx += 3ull << shift;
	else if (cur_ref_delta > 1)
		ctx += 4ull << shift;

	auto p = find_rce_context(m_ctx_rc_tuple_type, ctx, tpl_ctx_rc_tuple_type);
	int flag = (int)type;

	if(is_first_tuple)
		p->Encode(flag);
	else if (last_flag == tuple_types::insertion || last_flag == tuple_types::substitution)
		p->Encode(flag);
	else if (last_flag == tuple_types::match)
		p->EncodeExcluding(flag, (uint32_t)tuple_types::anchor);
	else if (last_flag == tuple_types::deletion)
		p->EncodeExcluding(flag, (uint32_t)tuple_types::skip);
	else if (last_flag == tuple_types::anchor)
		p->EncodeExcluding(flag, (uint32_t) tuple_types::anchor, (uint32_t) tuple_types::match);
	else if (last_flag == tuple_types::skip)
		p->EncodeExcluding(flag, (uint32_t) tuple_types::deletion, (uint32_t) tuple_types::skip);
	else if (last_flag == tuple_types::main_ref)
		p->EncodeExcluding(flag, (uint32_t) tuple_types::alt_id, (uint32_t) tuple_types::main_ref);
	else if (last_flag == tuple_types::alt_id)
		p->EncodeExcluding(flag, (uint32_t) tuple_types::alt_id, (uint32_t) tuple_types::main_ref);
/*	else
		p->Encode(flag);*/

	ctx_tuple_type <<= 3;
	ctx_tuple_type += flag;
	ctx_tuple_type &= ctx_mask_tuple_type;
}

// *****************************************************************
tuple_types CDNACoder::decode_tuple_type(uint32_t ref_symbol, tuple_types last_flag, bool is_first_tuple)
{
	context_t ctx = 0;
	uint32_t shift = 0;

	ctx = ctx_tuple_type;
	shift = 3 * no_tuples_in_mask;

	ctx += ((ctx_symbol & 0xfull) << shift);
	shift += 4;

	ctx += ((uint64_t)ref_symbol) << shift;
	shift += 2;

	if (cur_ref_delta < -10)
		ctx += 1ull << shift;
	else if (cur_ref_delta < -1)
		ctx += 2ull << shift;
	else if (cur_ref_delta > 10)
		ctx += 3ull << shift;
	else if (cur_ref_delta > 1)
		ctx += 4ull << shift;

	auto p = find_rce_context(m_ctx_rc_tuple_type, ctx, tpl_ctx_rc_tuple_type);

	uint32_t flag{};

	if(is_first_tuple)
		flag = p->Decode();
	else if (last_flag == tuple_types::insertion || last_flag == tuple_types::substitution)
		flag = p->Decode();
	else if (last_flag == tuple_types::match)
		flag = p->DecodeExcluding((uint32_t)tuple_types::anchor);
	else if (last_flag == tuple_types::deletion)
		flag = p->DecodeExcluding((uint32_t)tuple_types::skip);
	else if (last_flag == tuple_types::anchor)
		flag = p->DecodeExcluding((uint32_t)tuple_types::anchor, (uint32_t)tuple_types::match);
	else if (last_flag == tuple_types::skip)
		flag = p->DecodeExcluding((uint32_t)tuple_types::deletion, (uint32_t)tuple_types::skip);
	else if (last_flag == tuple_types::main_ref)
		flag = p->DecodeExcluding((uint32_t)tuple_types::alt_id, (uint32_t)tuple_types::main_ref);
	else if (last_flag == tuple_types::alt_id)
		flag = p->DecodeExcluding((uint32_t)tuple_types::alt_id, (uint32_t)tuple_types::main_ref);

	ctx_tuple_type <<= 3;
	ctx_tuple_type += flag;
	ctx_tuple_type &= ctx_mask_tuple_type;

	return (tuple_types) flag;
}

// *****************************************************************
void CDNACoder::encode_insertion(uint32_t symbol, uint32_t base)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	uint32_t shift = 0;
	context_t ctx = 2;		// insertion marker
	shift += 2;
		
	if (compression_level == 1)
	{
		ctx += (ctx_symbol & 0xffull) << shift;
		shift += 8;
	}
	else if (compression_level == 2)
	{
		ctx += (ctx_symbol & 0x3ffull) << shift;
		shift += 10;
	}
	else // (compression_level == 3)
	{
		ctx += (ctx_symbol & 0x3ffull) << shift;
		shift += 10;
		ctx += ((uint64_t)(((ctx_symbol >> 10) & 0x3u) == ((ctx_symbol >> 8) & 0x3u))) << shift;
		++shift;
	}

	ctx += ((context_t)base) << shift;
	shift += 2;
	ctx += (ctx_tuple_type & 0777ull) << shift;

	auto p = find_rce_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);

	assert(symbol < 4);
	assert(base < 4 || base == 255);

	p->Encode(symbol);
}

// *****************************************************************
void CDNACoder::prefetch_insertion(uint32_t base)
{
	uint32_t shift = 0;
	context_t ctx = 2;		// insertion marker
	shift += 2;

	if (compression_level == 1)
	{
		ctx += (ctx_symbol & 0xffull) << shift;
		shift += 8;
	}
	else if (compression_level == 2)
	{
		ctx += (ctx_symbol & 0x3ffull) << shift;
		shift += 10;
	}
	else // (compression_level == 3)
	{
		ctx += (ctx_symbol & 0x3ffull) << shift;
		shift += 10;
		ctx += ((uint64_t)(((ctx_symbol >> 10) & 0x3u) == ((ctx_symbol >> 8) & 0x3u))) << shift;
		++shift;
	}

	ctx += ((context_t)base) << shift;
	shift += 2;

	ctx += ((uint64_t)tuple_types::insertion) << shift;
	shift += 3;

	ctx += (ctx_tuple_type & 077ull) << shift;

	//prefetch_rc_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);
}

// *****************************************************************
uint32_t CDNACoder::decode_insertion(uint32_t base)
{
	uint32_t shift = 0;
	context_t ctx = 2;		// insertion marker
	shift += 2;

	if (compression_level == 1)
	{
		ctx += (ctx_symbol & 0xffull) << shift;
		shift += 8;
	}
	else if (compression_level == 2)
	{
		ctx += (ctx_symbol & 0x3ffull) << shift;
		shift += 10;
	}
	else // (compression_level == 3)
	{
		ctx += (ctx_symbol & 0x3ffull) << shift;
		shift += 10;
		ctx += ((uint64_t)(((ctx_symbol >> 10) & 0x3u) == ((ctx_symbol >> 8) & 0x3u))) << shift;
		++shift;
	}

	ctx += ((context_t)base) << shift;
	shift += 2;
	ctx += (ctx_tuple_type & 0777ull) << shift;

	auto p = find_rce_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);

	uint32_t symbol = p->Decode();

	assert(symbol < 4);
	assert(base < 4 || base == 255);

	return symbol;
}

// *****************************************************************
void CDNACoder::encode_substitution(uint32_t symbol, uint32_t base)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif
	
	uint32_t shift = 0;

	context_t ctx = 1;		// substitution marker
	shift += 2;

	ctx += (ctx_symbol & 0x3full) << shift;
	shift += 6;

	if (compression_level == 3)
	{
		ctx += ((uint64_t)(((ctx_symbol >> 6) & 0x3u) == ((ctx_symbol >> 4) & 0x3u))) << shift;
		shift++;
	}

	ctx += ((uint64_t) base) << shift;
	shift += 2;

	ctx += (ctx_tuple_type & 07777ull) << shift;

	auto p = find_rce_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);

	assert(symbol < 4);
	assert(base < 4);
	assert(symbol != base);

	p->EncodeExcluding(symbol, base);
}

// *****************************************************************
uint32_t CDNACoder::decode_substitution(uint32_t base)
{
	uint32_t shift = 0;

	context_t ctx = 1;		// substitution marker
	shift += 2;

	ctx += (ctx_symbol & 0x3full) << shift;
	shift += 6;

	if (compression_level == 3)
	{
		ctx += ((uint64_t)(((ctx_symbol >> 6) & 0x3u) == ((ctx_symbol >> 4) & 0x3u))) << shift;
		shift++;
	}

	ctx += ((uint64_t)base) << shift;
	shift += 2;

	ctx += (ctx_tuple_type & 07777ull) << shift;

	auto p = find_rce_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);

	uint32_t symbol = p->DecodeExcluding(base);

	assert(symbol < 4);
	assert(base < 4);
	assert(symbol != base);

	return symbol;
}

// *****************************************************************
void CDNACoder::encode_anchor_len(uint32_t len)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	for (int i_part = 0; len; ++i_part)
	{
		auto p = find_rce_context(m_ctx_rc_anchor_len, i_part, tpl_ctx_rc_anchor_len);

		if (len < max_enc_anchor_len - 1)
		{
			p->Encode(len);
			break;
		}

		p->Encode(max_enc_anchor_len - 1);
		len -= max_enc_anchor_len - 2;
	}
}

// *****************************************************************
uint32_t CDNACoder::decode_anchor_len()
{
	uint32_t len = 0;

	for (int i_part = 0; ; ++i_part)
	{
		auto p = find_rce_context(m_ctx_rc_anchor_len, i_part, tpl_ctx_rc_anchor_len);

		uint32_t x = p->Decode();

		if (x < max_enc_anchor_len - 1)
		{
			len += x;
			break;
		}

		len += max_enc_anchor_len - 2;
	}

	return len;
}

// *****************************************************************
void CDNACoder::encode_read_len(uint32_t len)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	context_t ctx = 0u;
	int no_bits = static_cast<int>(ilog2(len));

	auto p_no_bits = find_rce_context(m_ctx_rc_read_len_no_bits, 0, tpl_ctx_rc_read_len_no_bits);
	p_no_bits->Encode(no_bits);

	if (no_bits < 2)
		return;

	ctx = ((uint64_t)no_bits) << 3;

	len -= 1u << (no_bits - 1);

	uint32_t prefix;
	uint32_t suffix;

	if (no_bits <= 9)
	{
		prefix = len;
		suffix = 0;
	}
	else
	{
		prefix = len >> (no_bits - 9);
		suffix = len - (prefix << (no_bits - 9));
	}

	auto p_prefix = find_rc_context(m_ctx_rc_read_len_data, ctx, tpl_ctx_rc_read_len_data);
	p_prefix->Encode(prefix);

	if (no_bits <= 9)
		return;

	no_bits -= 9;

	ctx += 1ull << 2;

	for (; no_bits > 0; no_bits -= 8)
	{
		auto p = find_rc_context(m_ctx_rc_read_len_data, ctx, tpl_ctx_rc_read_len_data);

		p->Encode(suffix & 0xff);
		suffix >>= 8;
		++ctx;
	}
}

// *****************************************************************
uint32_t CDNACoder::decode_read_len()
{
	uint32_t len = 0;

	auto p_no_bits = find_rce_context(m_ctx_rc_read_len_no_bits, 0, tpl_ctx_rc_read_len_no_bits);
	int no_bits = p_no_bits->Decode();

	if (no_bits < 2)
		return no_bits;

	context_t ctx = ((uint64_t)no_bits) << 3;

	len = 1u << (no_bits - 1);

	auto p_prefix = find_rc_context(m_ctx_rc_read_len_data, ctx, tpl_ctx_rc_read_len_data);
	uint32_t prefix = p_prefix->Decode();

	if (no_bits <= 9)
		return len + prefix;

	len += prefix << (no_bits - 9);

	no_bits -= 9;

	ctx += 1ull << 2;

	uint32_t shift = 0;

	for (; no_bits > 0; no_bits -= 8)
	{
		auto p = find_rc_context(m_ctx_rc_read_len_data, ctx, tpl_ctx_rc_read_len_data);

		uint32_t r = p->Decode();
		len += r << shift;

		shift += 8;
		++ctx;
	}

	return len;
}

// *****************************************************************
void CDNACoder::encode_skip_len(uint32_t len, bool local)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	if (local)
		for (int i_part = 0; len; ++i_part)
		{
			auto p = find_rce_context(m_ctx_rc_skip_len_local, i_part, tpl_ctx_rc_skip_len_local);

			if (len < max_enc_skip_len - 1)
			{
				p->Encode(len);
				break;
			}

			p->Encode(max_enc_skip_len - 1);
			len -= max_enc_skip_len - 2;
		}
	else
	{
		uint32_t encoded_part = 0;
		for (int i = 3; i >= 0; --i)
		{
			int x = (len >> (8 * i)) & 0xff;

			context_t ctx = (uint32_t)i * 64ull + ilog2(encoded_part);

			auto p = find_rc_context(m_ctx_rc_skip_len_distant, ctx, tpl_ctx_rc_skip_len_distant);
			p->Encode(x);

			encoded_part = (encoded_part << 8) + x;
		}
	}
}

// *****************************************************************
uint32_t CDNACoder::decode_skip_len(bool local)
{
	uint32_t len = 0;

	if (local)
		for (int i_part = 0; ; ++i_part)
		{
			auto p = find_rce_context(m_ctx_rc_skip_len_local, i_part, tpl_ctx_rc_skip_len_local);

			uint32_t x = p->Decode();

			if (x < max_enc_skip_len - 1)
			{
				len += x;
				break;
			}

			len += max_enc_skip_len - 2;
		}
	else
	{
		//uint32_t decoded_part = 0;
		for (int i = 3; i >= 0; --i)
		{
			context_t ctx = (uint32_t)i * 64ull + ilog2(len);

			auto p = find_rc_context(m_ctx_rc_skip_len_distant, ctx, tpl_ctx_rc_skip_len_distant);
			uint32_t x = p->Decode();

			len = (len << 8) + x;
		}
	}

	return len;
}

// *****************************************************************
void CDNACoder::encode_symbol_plain(uint32_t symbol)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	context_t ctx = ctx_symbol << 2;		// plain marker

//	cout << ctx << endl;

	auto p = find_rce_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);

	assert(symbol < 4);

	p->Encode(symbol);

	ctx_symbol = ((ctx_symbol << 2) + symbol) & ctx_mask_symbol;
}

// *****************************************************************
uint32_t CDNACoder::decode_symbol_plain()
{
	context_t ctx = ctx_symbol << 2;		// plain marker

	auto p = find_rce_context(m_ctx_rc_symbols, ctx, tpl_ctx_rc_symbols);

	auto symbol = p->Decode();

	ctx_symbol = ((ctx_symbol << 2) + symbol) & ctx_mask_symbol;

	return symbol;
}

// *****************************************************************
void CDNACoder::encode_symbol_plain_with_Ns(uint32_t symbol)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	auto p = find_rce_context(m_ctx_rc_symbols_with_Ns, ctx_symbol, tpl_ctx_rc_symbols_with_Ns);

	assert(symbol < 5);

	p->Encode(symbol);

	ctx_symbol = ((ctx_symbol << 4) + symbol) & ctx_mask_symbol;
}

// *****************************************************************
uint32_t CDNACoder::decode_symbol_plain_with_Ns()
{
	auto p = find_rce_context(m_ctx_rc_symbols_with_Ns, ctx_symbol, tpl_ctx_rc_symbols_with_Ns);

	auto symbol = p->Decode();

	ctx_symbol = ((ctx_symbol << 4) + symbol) & ctx_mask_symbol;

	return symbol;
}

// *****************************************************************
void CDNACoder::Init(bool _is_compressing, int _max_no_alt_refs, int _compression_level, uint64_t _input_stream_size, uint32_t start_read_id)
{
	is_compressing = _is_compressing;
	max_no_alt_refs = _max_no_alt_refs;
	compression_level = _compression_level;
	input_stream_size = _input_stream_size;

	v_vios_io = new CVectorIOStream(v_io);

	CBasicRangeCoder<CVectorIOStream>* rc;

	if (compression_level == 3)
	{
		no_tuples_in_mask = 4;
		ctx_mask_tuple_type = 07777ull;
		no_symbols_in_mask = 8;
		ctx_mask_symbol = (1ull << (2 * no_symbols_in_mask)) - 1ull;
	}
	else if (compression_level == 2)
	{
		no_tuples_in_mask = 3;
		ctx_mask_tuple_type = 0777ull;
		no_symbols_in_mask = 7;
		ctx_mask_symbol = (1ull << (2 * no_symbols_in_mask)) - 1ull;
	}
	else if (compression_level == 1)
	{
		no_tuples_in_mask = 2;
		ctx_mask_tuple_type = 077ull;
		no_symbols_in_mask = 5;
		ctx_mask_symbol = (1ull << (2 * no_symbols_in_mask)) - 1ull;
	}
	else
	{
		no_tuples_in_mask = 1;
		ctx_mask_tuple_type = 07ull;
		no_symbols_in_mask = 1;
		ctx_mask_symbol = (1ull << (2 * no_symbols_in_mask)) - 1ull;
	}

//	m_ctx_rc_symbols.resize(ctx_mask_symbol + 1ull);

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

	// *** Flags for read types
	delete tpl_ctx_rc_read_type;
	tpl_ctx_rc_read_type = new rcmfs_read_type_t(rc, nullptr, is_compressing);
	ctx_read_type = 0;

	// *** Flags for marking rev. comp. reference reads
	delete tpl_ctx_rc_rev_comp;
	tpl_ctx_rc_rev_comp = new rcmfs_rev_comp_t(rc, nullptr, is_compressing);
	ctx_rev_comp = 0;

	// *** Flags telling whether alt_read_id was recently seen in a read
	delete tpl_ctx_rc_seen_read_id;
	tpl_ctx_rc_seen_read_id = new rcmfs_seen_read_id_t(rc, nullptr, is_compressing);

	// *** Read length
	delete tpl_ctx_rc_read_len_no_bits;
	tpl_ctx_rc_read_len_no_bits = new rcmfs_read_len_no_bits_t(rc, nullptr, is_compressing);

	delete tpl_ctx_rc_read_len_data;
	tpl_ctx_rc_read_len_data = new rcmfs_read_len_data_t(rc, nullptr, is_compressing);

	// *** Plain symbols
	delete tpl_ctx_rc_symbols;
	tpl_ctx_rc_symbols = new rcmfs_symbols_t(rc, nullptr, is_compressing);

	// *** Plain symbols with Ns
	delete tpl_ctx_rc_symbols_with_Ns;
	tpl_ctx_rc_symbols_with_Ns = new rcmfs_symbols_with_Ns_t(rc, nullptr, is_compressing);

	// *** Read ids
	delete tpl_ctx_rc_read_id;
	tpl_ctx_rc_read_id = new rcmfs_read_id_t(rc, nullptr, is_compressing);

	// *** Read ids short (for read ids allready seen in a read)
	delete tpl_ctx_rc_read_id_short;
	tpl_ctx_rc_read_id_short = new rcm_read_id_short_t(rc, max_no_alt_refs, nullptr, is_compressing);

	// *** Anchor length
	delete tpl_ctx_rc_anchor_len;
	tpl_ctx_rc_anchor_len = new rcm_anchor_len_t(rc, max_enc_anchor_len, nullptr, is_compressing);

	// *** Skip length - local
	delete tpl_ctx_rc_skip_len_local;
	tpl_ctx_rc_skip_len_local = new rcm_skip_len_local_t(rc, max_enc_skip_len, nullptr, is_compressing);

	// *** Skip length - distant
	delete tpl_ctx_rc_skip_len_distant;
	tpl_ctx_rc_skip_len_distant = new rcmfs_skip_len_distant_t(rc, nullptr, is_compressing);

	// *** Flags for read types
	delete tpl_ctx_rc_tuple_type;
	tpl_ctx_rc_tuple_type = new rcmfs_tuple_type_t(rc, nullptr, is_compressing);

	cur_read_id = start_read_id;
}

// *****************************************************************
void CDNACoder::Restart()
{
	if (is_compressing)
		rce->Start();
	else
	{
		v_vios_io->RestartRead();
		rcd->Start();
	}
}

// *****************************************************************
void CDNACoder::Finish()
{
	if (is_compressing)
		rce->End();
	else
		rcd->End();
}

} // namespace entropy_coder

// EOF
