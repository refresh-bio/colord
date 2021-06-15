#include "quality_coder.h"
#include <algorithm>

namespace entropy_coder
{
// *****************************************************************
void CQualityCoder::analyze_es(es_t& es, uint32_t read_size)
{
	v_base_flags.clear();

	tuple_types tuple_type = tuple_types::none;
	uint32_t tuple_val1;
	uint32_t tuple_val2;

	es.restart_reading();
	es.load(tuple_type, tuple_val1, tuple_val2);

	if (tuple_type == tuple_types::start_plain || tuple_type == tuple_types::start_plain_with_Ns)
	{
		v_base_flags.resize(read_size, 'P');

		return;
	}

	v_base_flags.resize(read_size, ' ');
	auto ptr = v_base_flags.begin();

	while (es.load(tuple_type, tuple_val1, tuple_val2))
	{
		if (tuple_type == tuple_types::anchor)
		{
			fill_n(ptr, tuple_val2, 'A');
			ptr += tuple_val2;
		}
		else if (tuple_type == tuple_types::match)
		{
			*ptr++ = 'M';
		}
		else if (tuple_type == tuple_types::insertion)
		{
			++ptr;
		}
		else if (tuple_type == tuple_types::deletion)
		{
		}
		else if (tuple_type == tuple_types::substitution)
		{
			++ptr;
		}
		else if (tuple_type == tuple_types::skip)
		{
		}
		else if (tuple_type == tuple_types::main_ref)
		{
		}
	}
}
	
// *****************************************************************
void CQualityCoder::encode_original(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (compression_level == 3)
		{
			if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
			shift += 2;
		}
		else
		{
			if (i > 1)					ctx += ((context_t)(read[i - 2] == read[i - 1])) << shift;
			shift += 1;
		}

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;

		/*			if (i > 1)					ctx += ((context_t)(read[i - 2] == read[i-1])) << shift;
					shift += 1;
					if (i + 2u < read_size)		ctx += ((context_t)(read[i + 1u] == read[i + 2u])) << shift;
					shift += 1;*/

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;
			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rc_context(m_quality_original_rc, ctx, tpl_quality_original_rc);

		assert(quality[i] >= 33);

		auto q = quality_code_map_fwd[quality[i] - 33u];

		p_rc->Encode(q);

		update_context(context, quality_code_quantize[q]);
	}
}

// *****************************************************************
void CQualityCoder::encode_quinary_average(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	array<pair<double, uint32_t>, 5> q_stats{ make_pair(0.0, 0), make_pair(0.0, 0), make_pair(0.0, 0), make_pair(0.0, 0), make_pair(0.0, 0) };
	array<uint32_t, 128> q_stats128{};

	for (uint32_t i = 0; i < read_size; ++i)
		++q_stats128[quality[i]];

	for (uint32_t i = 33; i < 128; ++i)
	{
		q_stats[quality_code_map_fwd[i - 33u]].first += (double)(i - 33u) * q_stats128[i];
		q_stats[quality_code_map_fwd[i - 33u]].second += q_stats128[i];
	}

	context_t ctx_p = 0;

	for (uint32_t i = 0; i < 5; ++i)
	{
		double avg = q_stats[i].second != 0 ? q_stats[i].first / q_stats[i].second : 0.0;

		encode_avg((1ull << 30) + ((context_t)i << 24) + (ctx_p << 16), avg);
		ctx_p = (context_t)avg;
	}

	context_t dna_ctx = ((context_t)read[0]);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		dna_ctx <<= 2;
		if (i + 1 < read_size)
			dna_ctx += (context_t)read[i + 1];
		dna_ctx &= 0xff;

		ctx += dna_ctx << shift;
		shift += 8;

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;
			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rce_context(m_quality_quinary_rc, ctx, tpl_quality_quinary_rc);
		auto q = quality_code_map_fwd[quality[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}
}

// *****************************************************************
void CQualityCoder::encode_quad_average(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	array<pair<double, uint32_t>, 4> q_stats{ make_pair(0.0, 0), make_pair(0.0, 0), make_pair(0.0, 0), make_pair(0.0, 0) };
	array<uint32_t, 128> q_stats128{};

	for (uint32_t i = 0; i < read_size; ++i)
		++q_stats128[quality[i]];

	for (uint32_t i = 33; i < 128; ++i)
	{
		q_stats[quality_code_map_fwd[i - 33u]].first += (double)(i - 33u) * q_stats128[i];
		q_stats[quality_code_map_fwd[i - 33u]].second += q_stats128[i];
	}

	context_t ctx_p = 0;

	for (uint32_t i = 0; i < 4; ++i)
	{
		double avg = q_stats[i].second != 0 ? q_stats[i].first / q_stats[i].second : 0.0;

		encode_avg((1ull << 30) + ((context_t)i << 24) + (ctx_p << 16), avg);
		ctx_p = (context_t)avg;
	}

	context_t dna_ctx = ((context_t)read[0]);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		dna_ctx <<= 2;
		if (i + 1 < read_size)
			dna_ctx += (context_t)read[i + 1];
		dna_ctx &= 0xff;

		ctx += dna_ctx << shift;
		shift += 8;

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;
			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rce_context(m_quality_quad_rc, ctx, tpl_quality_quad_rc);
		auto q = quality_code_map_fwd[quality[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}
}

// *****************************************************************
void CQualityCoder::encode_binary_average(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	array<pair<double, uint32_t>, 2> q_stats{ make_pair(0.0, 0), make_pair(0.0, 0) };
	array<uint32_t, 128> q_stats128{};

	for (uint32_t i = 0; i < read_size; ++i)
		++q_stats128[quality[i]];

	for (uint32_t i = 33; i < 128; ++i)
	{
		q_stats[quality_code_map_fwd[i - 33u]].first += (double)(i - 33u) * q_stats128[i];
		q_stats[quality_code_map_fwd[i - 33u]].second += q_stats128[i];
	}

	context_t ctx_p = 0;

	for (uint32_t i = 0; i < 2; ++i)
	{
		double avg = q_stats[i].second != 0 ? q_stats[i].first / q_stats[i].second : 0.0;

		encode_avg((1ull << 30) + ((context_t)i << 24) + (ctx_p << 16), avg);
		ctx_p = (context_t)avg;
	}

	context_t dna_ctx = ((context_t)read[0]);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		dna_ctx <<= 2;
		if (i + 1 < read_size)
			dna_ctx += (context_t)read[i + 1];
		dna_ctx &= 0xff;

		ctx += dna_ctx << shift;
		shift += 8;

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;
			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rce_context(m_quality_binary_rc, ctx, tpl_quality_binary_rc);
		auto q = quality_code_map_fwd[quality[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}
}

// *****************************************************************
void CQualityCoder::encode_quinary_threshold(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
		shift += 2;

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;
		//			if (i > 2)					ctx += ((context_t)read[i - 3]) << 44;

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;

			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rce_context(m_quality_quinary_rc, ctx, tpl_quality_quinary_rc);
		auto q = quality_code_map_fwd[quality[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}
}

// *****************************************************************
void CQualityCoder::encode_quad_threshold(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
		shift += 2;

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;
		//			if (i > 2)					ctx += ((context_t)read[i - 3]) << 44;

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;

			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rce_context(m_quality_quad_rc, ctx, tpl_quality_quad_rc);
		auto q = quality_code_map_fwd[quality[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}
}

// *****************************************************************
void CQualityCoder::encode_binary_threshold(const read_t& read, const qual_t& quality, es_t& es)
{
	context_t context;
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
		shift += 2;

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;

		if (compression_level > 1)
		{
			ctx += ((context_t)(v_base_flags[i] == 'M')) << shift;
			++shift;
			ctx += ((context_t)(v_base_flags[i] == 'A')) << shift;
		}

		auto p_rc = find_rce_context(m_quality_binary_rc, ctx, tpl_quality_binary_rc);
		auto q = quality_code_map_fwd[quality[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}
}

// *****************************************************************
void CQualityCoder::encode_average(const read_t& read, const qual_t& quality, es_t& es)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	double avg_qual = 0.0;

	for (uint32_t i = 0; i < read_size; ++i)
		avg_qual += quality[i] - 33u;

	avg_qual /= read_size;

	encode_avg(0ull, avg_qual);
}

// *****************************************************************
void CQualityCoder::decode_original(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	quality.resize(read_size);

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (compression_level == 3)
		{
			if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
			shift += 2;
		}
		else
		{
			if (i > 1)					ctx += ((context_t)(read[i - 2] == read[i - 1])) << shift;
			shift += 1;
		}

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
			/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
						ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/
		}

		auto p_rc = find_rc_context(m_quality_original_rc, ctx, tpl_quality_original_rc);
		auto d = p_rc->Decode();
		auto q = quality_code_map_rev[d];
//		quality.push_back(q + 33);
		quality[i] = q + 33;

		update_context(context, quality_code_quantize[q]);
	}
}

// *****************************************************************
void CQualityCoder::decode_quinary_average(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	reset_context(context);

	array<double, 5> avg;
	context_t ctx_p = 0;

	for (uint32_t i = 0; i < 5; ++i)
	{
		avg[i] = decode_avg((1ull << 30) + ((context_t)i << 24) + (ctx_p << 16));

		ctx_p = (context_t)avg[i];
	}

	array<double, 5> avg_sum = { 0, 0, 0, 0, 0 };
	array<double, 5> qual_sum = { 0, 0, 0, 0, 0 };

	context_t dna_ctx = ((context_t)read[0]);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		dna_ctx <<= 2;
		if (i + 1 < read_size)
			dna_ctx += (context_t)read[i + 1];
		dna_ctx &= 0xff;

		ctx += dna_ctx << shift;
		shift += 8;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
			/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
						ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/
		}

		auto p_rc = find_rce_context(m_quality_quinary_rc, ctx, tpl_quality_quinary_rc);
		auto d = p_rc->Decode();

		avg_sum[d] += avg[d];
		auto q = (uint32_t)(avg_sum[d] - qual_sum[d]);
		qual_sum[d] += q;
		quality.push_back(q + 33);

		update_context(context, d);
	}
}

// *****************************************************************
void CQualityCoder::decode_quad_average(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	reset_context(context);

	array<double, 4> avg;
	context_t ctx_p = 0;

	for (uint32_t i = 0; i < 4; ++i)
	{
		avg[i] = decode_avg((1ull << 30) + ((context_t)i << 24) + (ctx_p << 16));

		ctx_p = (context_t)avg[i];
	}

	array<double, 4> avg_sum = { 0, 0, 0, 0 };
	array<double, 4> qual_sum = { 0, 0, 0, 0 };

	context_t dna_ctx = ((context_t)read[0]);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		dna_ctx <<= 2;
		if (i + 1 < read_size)
			dna_ctx += (context_t)read[i + 1];
		dna_ctx &= 0xff;

		ctx += dna_ctx << shift;
		shift += 8;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
			/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
						ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/
		}

		auto p_rc = find_rce_context(m_quality_quad_rc, ctx, tpl_quality_quad_rc);
		auto d = p_rc->Decode();

		avg_sum[d] += avg[d];
		auto q = (uint32_t)(avg_sum[d] - qual_sum[d]);
		qual_sum[d] += q;
		quality.push_back(q + 33);

		update_context(context, d);
	}
}

// *****************************************************************
void CQualityCoder::decode_binary_average(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	reset_context(context);

	array<double, 2> avg;
	context_t ctx_p = 0;

	for (uint32_t i = 0; i < 2; ++i)
	{
		avg[i] = decode_avg((1ull << 30) + ((context_t)i << 24) + (ctx_p << 16));

		ctx_p = (context_t)avg[i];
	}

	array<double, 2> avg_sum = { 0, 0 };
	array<double, 2> qual_sum = { 0, 0 };

	context_t dna_ctx = ((context_t)read[0]);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		dna_ctx <<= 2;
		if (i + 1 < read_size)
			dna_ctx += (context_t)read[i + 1];
		dna_ctx &= 0xff;

		ctx += dna_ctx << shift;
		shift += 8;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
			/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
						ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/
		}

		auto p_rc = find_rce_context(m_quality_binary_rc, ctx, tpl_quality_binary_rc);
		auto d = p_rc->Decode();

		avg_sum[d] += avg[d];
		auto q = (uint32_t)(avg_sum[d] - qual_sum[d]);
		qual_sum[d] += q;
		quality.push_back(q + 33);

		update_context(context, d);
	}
}

// *****************************************************************
void CQualityCoder::decode_quinary_threshold(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
		shift += 2;

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
			/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
						ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/
		}

		auto p_rc = find_rce_context(m_quality_quinary_rc, ctx, tpl_quality_quinary_rc);
		auto d = p_rc->Decode();
		auto q = quality_code_map_rev[d];
		quality.push_back(q + 33);

		update_context(context, d);
	}
}

// *****************************************************************
void CQualityCoder::decode_quad_threshold(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
		shift += 2;

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
			/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
						ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/
		}

		auto p_rc = find_rce_context(m_quality_quad_rc, ctx, tpl_quality_quad_rc);
		auto d = p_rc->Decode();
		auto q = quality_code_map_rev[d];
		quality.push_back(q + 33);

		update_context(context, d);
	}
}

// *****************************************************************
void CQualityCoder::decode_binary_threshold(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;
	context_t context;

	reset_context(context);

	for (uint32_t i = 0; i < read_size; ++i)
	{
		context_t ctx = context;
		uint32_t shift = ctx_size_in_bits;

		ctx += ((context_t)read[i]) << shift;
		shift += 2;

		if (i > 0)					ctx += ((context_t)read[i - 1]) << shift;
		shift += 2;

		if (i > 1)					ctx += ((context_t)read[i - 2]) << shift;
		shift += 2;

		if (i + 1u < read_size)		ctx += ((context_t)read[i + 1u]) << shift;
		shift += 2;

		if (compression_level > 1)
		{
			if (v_base_flags[i])
				ctx += ((context_t)v_base_flags[i]) << shift;
		}
		/*			ctx += ((context_t)(v_flags[i] == 'M')) << 45;
					ctx += ((context_t)(v_flags[i] == 'A')) << 46;*/

		auto p_rc = find_rce_context(m_quality_binary_rc, ctx, tpl_quality_binary_rc);
		auto d = p_rc->Decode();
		auto q = quality_code_map_rev[d];
		quality.push_back(q + 33);

		update_context(context, d);
	}
}

// *****************************************************************
void CQualityCoder::decode_average(read_t& read, qual_t& quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	double avg = decode_avg(0ull);

	double avg_sum = 0.0;
	double qual_sum = 0.0;
	uint32_t q;

	for (uint32_t i = 0; i < read_size; ++i)
	{
		avg_sum += avg;
		q = (uint32_t)(avg_sum - qual_sum);
		qual_sum += q;
		quality.push_back(q + 33);
	}
}


// *****************************************************************
void CQualityCoder::encode_avg(context_t ctx_base, double x)
{
	uint32_t a = (uint32_t)(x * 256);
	uint32_t a1 = a >> 8;
	uint32_t a2 = a & 0xffu;

	context_t ctx = ctx_base;
	auto p_rc1 = find_rc_context(m_quality_byte_rc, ctx, tpl_quality_byte_rc);
	p_rc1->Encode(a1);

	ctx = a1 + 0x100ull;
	auto p_rc2 = find_rc_context(m_quality_byte_rc, ctx, tpl_quality_byte_rc);
	p_rc2->Encode(a2);
}

// *****************************************************************
double CQualityCoder::decode_avg(context_t ctx_base)
{
	context_t ctx = ctx_base;
	auto p_rc1 = find_rc_context(m_quality_byte_rc, ctx, tpl_quality_byte_rc);
	uint32_t a1 = p_rc1->Decode();

	ctx = a1 + 0x100ull;
	auto p_rc2 = find_rc_context(m_quality_byte_rc, ctx, tpl_quality_byte_rc);
	uint32_t a2 = p_rc2->Decode();

	uint32_t a = (a1 << 8) + a2;
	return (double)a / 256.0;
}

}