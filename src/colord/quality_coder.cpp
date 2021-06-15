#include "quality_coder.h"
#include <algorithm>

namespace entropy_coder
{

// *****************************************************************
void CQualityCoder::Init(bool _is_compressing, QualityComprMode _quality_mode, DataSource _data_source, const std::vector<uint32_t>& _quality_fwd_thrs, const std::vector<uint32_t>& _quality_rev_thrs, 
	int _compression_level, uint64_t _input_stream_size)
{
	is_compressing = _is_compressing;
	quality_mode = _quality_mode;
	data_source = _data_source;
	quality_fwd_thrs = _quality_fwd_thrs;
	quality_rev_thrs = _quality_rev_thrs;
	compression_level = _compression_level;
	input_stream_size = _input_stream_size;

	v_vios_io = new CVectorIOStream(v_io);

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

	delete tpl_quality_binary_rc;
	delete tpl_quality_quad_rc;
	delete tpl_quality_quinary_rc;
	delete tpl_quality_original_rc;
	delete tpl_quality_byte_rc;

	if (data_source == DataSource::ONT)
	{
		switch (quality_mode)
		{
		case QualityComprMode::Original:
			tpl_quality_original_rc = new rcmfs_quality_original_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_lossless();
			no_bits_per_symbol = 4;		// quantized values
			no_ctx_symbols = 2;
			break;

		case QualityComprMode::QuinaryAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_quinary_rc = new rcmfs_quality_quinary_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_quinary();
			no_bits_per_symbol = 3;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::QuadAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_quad_rc = new rcmfs_quality_quad_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_quad();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::BinaryAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_binary_rc = new rcmfs_quality_binary_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_binary();
			no_bits_per_symbol = 1 + 1;
			no_ctx_symbols = 6;
			break;

		case QualityComprMode::QuinaryThreshold:
			tpl_quality_quinary_rc = new rcmfs_quality_quinary_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_quinary();
			no_bits_per_symbol = 3;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::QuadThreshold:
			tpl_quality_quad_rc = new rcmfs_quality_quad_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_quad();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::BinaryThreshold:
			tpl_quality_binary_rc = new rcmfs_quality_binary_t(rc, nullptr, is_compressing);
			adjust_quality_map_ONT_binary();
			no_bits_per_symbol = 1 + 1;
			no_ctx_symbols = 6;
			break;

		case QualityComprMode::Average:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			no_bits_per_symbol = 8;
			no_ctx_symbols = 2;
			break;
		default:
			break;
		}
	}
	else if (data_source == DataSource::PBRaw)
	{
		switch (quality_mode)
		{
		case QualityComprMode::Original:
			tpl_quality_original_rc = new rcmfs_quality_original_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_lossless();
			no_bits_per_symbol = 4;		// quantized values
			no_ctx_symbols = 2;
			break;
		
		case QualityComprMode::QuinaryAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_quinary_rc = new rcmfs_quality_quinary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_quinary();
			no_bits_per_symbol = 3;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::QuadAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_quad_rc = new rcmfs_quality_quad_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_quad();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::BinaryAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_binary_rc = new rcmfs_quality_binary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_binary();
			no_bits_per_symbol = 1 + 1;
			no_ctx_symbols = 6;
			break;

		case QualityComprMode::QuinaryThreshold:
			tpl_quality_quinary_rc = new rcmfs_quality_quinary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_quinary();
			no_bits_per_symbol = 3;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::QuadThreshold:
			tpl_quality_quad_rc = new rcmfs_quality_quad_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_quad();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::BinaryThreshold:
			tpl_quality_binary_rc = new rcmfs_quality_binary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBRaw_binary();
			no_bits_per_symbol = 1 + 1;
			no_ctx_symbols = 6;
			break;

		case QualityComprMode::Average:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			no_bits_per_symbol = 8;
			no_ctx_symbols = 2;
			break;
		default:
			break;
		}
	}
	else if (data_source == DataSource::PBHiFi)
	{
		switch (quality_mode)
		{
		case QualityComprMode::Original:
			tpl_quality_original_rc = new rcmfs_quality_original_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_lossless();
			no_bits_per_symbol = 4;		// quantized values
			no_ctx_symbols = 2;
			break;

		case QualityComprMode::QuinaryAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_quinary_rc = new rcmfs_quality_quinary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_quinary();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::QuadAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_quad_rc = new rcmfs_quality_quad_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_quad();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::BinaryAverage:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			tpl_quality_binary_rc = new rcmfs_quality_binary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_binary();
			no_bits_per_symbol = 1 + 1;
			no_ctx_symbols = 6;
			break;

		case QualityComprMode::QuinaryThreshold:
			tpl_quality_quinary_rc = new rcmfs_quality_quinary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_quinary();
			no_bits_per_symbol = 3;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::QuadThreshold:
			tpl_quality_quad_rc = new rcmfs_quality_quad_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_quad();
			no_bits_per_symbol = 2 + 1;
			no_ctx_symbols = 3;
			break;
		case QualityComprMode::BinaryThreshold:
			tpl_quality_binary_rc = new rcmfs_quality_binary_t(rc, nullptr, is_compressing);
			adjust_quality_map_PBHiFi_binary();
			no_bits_per_symbol = 1 + 1;
			no_ctx_symbols = 6;
			break;

		case QualityComprMode::Average:
			tpl_quality_byte_rc = new rcmfs_quality_byte_t(rc, nullptr, is_compressing);
			no_bits_per_symbol = 8;
			no_ctx_symbols = 2;
			break;
		default:
			break;
		}
	}
	else
		assert(0);		// Cannot be here

	ctx_size_in_bits = no_bits_per_symbol * no_ctx_symbols;
	ctx_mask = (1ull << ctx_size_in_bits) - 1ull;
	ctx_inc = 1ull << 48;
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_symbols(uint32_t n)
{
	assert(n >= 2);

	if (!quality_fwd_thrs.empty())
	{
		for (uint32_t i = 0; i < quality_fwd_thrs[0]; ++i)
			quality_code_map_fwd[i] = 0;

		for(uint32_t bin = 1; bin < n - 1u; ++bin)
			for (uint32_t i = quality_fwd_thrs[bin-1u]; i < quality_fwd_thrs[bin]; ++i)
				quality_code_map_fwd[i] = bin;

		for (uint32_t i = quality_fwd_thrs[n - 2u]; i < 96u; ++i)
			quality_code_map_fwd[i] = n - 1u;
	}

	if (!quality_rev_thrs.empty())
		for(uint32_t i = 0; i < n; ++i)
			quality_code_map_rev[i] = quality_rev_thrs[i];
}

#ifndef __GNUC__
#pragma region(adjust_ONT)
#endif
//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_ONT_lossless()
{
	for (int i = 0; i < 96; ++i)
	{
		quality_code_map_fwd[i] = i;
		quality_code_map_rev[i] = i;
	}

	auto p = quality_code_quantize.begin();

	if (compression_level == 3)
	{
		quality_code_quantize[0] = 0;
		quality_code_quantize[1] = 1;
		fill(p + 2, p + 4, 2);
		fill(p + 4, p + 7, 3);
		fill(p + 7, p + 11, 4);
		fill(p + 11, p + 16, 5);
		fill(p + 16, p + 22, 6);
		fill(p + 22, p + 29, 7);
		fill(p + 29, p + 37, 8);
		fill(p + 37, p + 46, 9);
		fill(p + 46, p + 56, 10);
		fill(p + 56, p + 67, 11);
		fill(p + 67, p + 79, 12);
		fill(p + 79, p + 90, 13);
		fill(p + 90, p + 96, 14);
	}
	else if (compression_level == 2)
	{
		quality_code_quantize[0] = 0;
		quality_code_quantize[1] = 1;
		fill(p + 2, p + 5, 2);
		fill(p + 5, p + 10, 3);
		fill(p + 10, p + 15, 4);
		fill(p + 15, p + 20, 5);
		fill(p + 20, p + 25, 6);
		fill(p + 25, p + 35, 7);
		fill(p + 35, p + 50, 8);
		fill(p + 50, p + 70, 9);
		fill(p + 70, p + 96, 10);
	}
	else if (compression_level == 1)
	{
		quality_code_quantize[0] = 0;
		quality_code_quantize[1] = 1;
		fill(p + 2, p + 5, 2);
		fill(p + 5, p + 10, 3);
		fill(p + 10, p + 15, 4);
		fill(p + 15, p + 20, 5);
		fill(p + 20, p + 25, 6);
		fill(p + 25, p + 35, 7);
		fill(p + 35, p + 50, 8);
		fill(p + 50, p + 70, 9);
		fill(p + 70, p + 96, 10);
	}
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_ONT_binary()
{
	adjust_quality_map_symbols(2);
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_ONT_quad()
{
	adjust_quality_map_symbols(4);
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_ONT_quinary()
{
	adjust_quality_map_symbols(5);
}
#ifndef __GNUC__
#pragma endregion
#endif

#ifndef __GNUC__
#pragma region(adjust_PBRaw)
#endif
//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBRaw_lossless()
{
	for (int i = 0; i < 96; ++i)
	{
		quality_code_map_fwd[i] = i;
		quality_code_map_rev[i] = i;
	}

	auto p = quality_code_quantize.begin();

	if (compression_level == 3)
	{
		quality_code_quantize[0] = 0;
		fill(p + 1, p + 10, 1);
		fill(p + 10, p + 20, 2);
		fill(p + 20, p + 30, 3);
		fill(p + 30, p + 39, 4);
		fill(p + 39, p + 45, 5);
		fill(p + 45, p + 51, 6);
		fill(p + 51, p + 57, 7);
		fill(p + 57, p + 63, 8);
		fill(p + 63, p + 69, 9);
		fill(p + 69, p + 75, 10);
		fill(p + 75, p + 81, 11);
		fill(p + 81, p + 87, 12);
		fill(p + 87, p + 93, 13);

		quality_code_quantize[93] = 14;
	}
	else if (compression_level == 2)
	{
		quality_code_quantize[0] = 0;
		fill(p + 1, p + 15, 1);
		fill(p + 15, p + 29, 2);
		fill(p + 29, p + 41, 3);
		fill(p + 41, p + 53, 4);
		fill(p + 53, p + 63, 5);
		fill(p + 63, p + 72, 6);
		fill(p + 72, p + 80, 7);
		fill(p + 80, p + 87, 8);
		fill(p + 87, p + 93, 9);

		quality_code_quantize[93] = 10;
	}
	else if (compression_level == 1)
	{
		quality_code_quantize[0] = 0;
		fill(p + 1, p + 15, 1);
		fill(p + 15, p + 29, 2);
		fill(p + 29, p + 41, 3);
		fill(p + 41, p + 53, 4);
		fill(p + 53, p + 63, 5);
		fill(p + 63, p + 72, 6);
		fill(p + 72, p + 80, 7);
		fill(p + 80, p + 87, 8);
		fill(p + 87, p + 93, 9);

		quality_code_quantize[93] = 10;
	}
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBRaw_binary()
{
	adjust_quality_map_symbols(2);
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBRaw_quad()
{
	adjust_quality_map_symbols(4);
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBRaw_quinary()
{
	adjust_quality_map_symbols(5);
}
#ifndef __GNUC__
#pragma endregion
#endif

#ifndef __GNUC__
#pragma region(adjust_PBHiFi)
#endif
//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBHiFi_lossless()
{
	for (int i = 0; i < 96; ++i)
	{
		quality_code_map_fwd[i] = i;
		quality_code_map_rev[i] = i;
	}

	auto p = quality_code_quantize.begin();

	if (compression_level == 3)
	{
		quality_code_quantize[0] = 1;
		fill(p + 1, p + 10, 2);
		fill(p + 10, p + 20, 3);
		fill(p + 20, p + 30, 4);
		fill(p + 30, p + 39, 5);
		fill(p + 39, p + 45, 6);
		fill(p + 45, p + 51, 7);
		fill(p + 51, p + 57, 8);
		fill(p + 57, p + 63, 9);
		fill(p + 63, p + 69, 10);
		fill(p + 69, p + 75, 11);
		fill(p + 75, p + 81, 12);
		fill(p + 81, p + 87, 13);
		fill(p + 87, p + 93, 14);

		quality_code_quantize[93] = 0;
	}
	else if (compression_level == 2)
	{
		quality_code_quantize[0] = 1;
		fill(p + 1, p + 15, 2);
		fill(p + 15, p + 29, 3);
		fill(p + 29, p + 41, 4);
		fill(p + 41, p + 53, 5);
		fill(p + 53, p + 63, 6);
		fill(p + 63, p + 72, 7);
		fill(p + 72, p + 80, 8);
		fill(p + 80, p + 87, 9);
		fill(p + 87, p + 93, 10);

		quality_code_quantize[93] = 0;
	}
	else if (compression_level == 1)
	{
		quality_code_quantize[0] = 1;
		fill(p + 1, p + 15, 2);
		fill(p + 15, p + 29, 3);
		fill(p + 29, p + 41, 4);
		fill(p + 41, p + 53, 5);
		fill(p + 53, p + 63, 6);
		fill(p + 63, p + 72, 7);
		fill(p + 72, p + 80, 8);
		fill(p + 80, p + 87, 9);
		fill(p + 87, p + 93, 10);

		quality_code_quantize[93] = 0;
	}
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBHiFi_binary()
{
	adjust_quality_map_symbols(2);
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBHiFi_quad()
{
	adjust_quality_map_symbols(4);
}

//*****************************************************************************************************
void CQualityCoder::adjust_quality_map_PBHiFi_quinary()
{
	adjust_quality_map_symbols(5);
}
#ifndef __GNUC__
#pragma endregion
#endif

//*****************************************************************************************************
void CQualityCoder::reset_context(context_t& ctx)
{
	ctx = ctx_mask;
}

//*****************************************************************************************************
void CQualityCoder::update_context(context_t& ctx, uint8_t x)
{
	ctx = ((ctx << no_bits_per_symbol) + (context_t)x) & ctx_mask;
}
// *****************************************************************
void CQualityCoder::Restart()
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
void CQualityCoder::Finish()
{
	if (is_compressing)
		rce->End();
	else
		rcd->End();
}

// *****************************************************************
void CQualityCoder::Encode(const read_t& read, const qual_t& quality, es_t &es)
{
#ifdef ALLOW_ZERO_COMPRESSION_MODE
	if (compression_level == 0)
		return;
#endif

	if (quality_mode == QualityComprMode::None)
		return;

	if(compression_level > 1)
		analyze_es(es, static_cast<uint32_t>(read.size()) - 1u);

	switch (quality_mode)
	{
	case QualityComprMode::Original:
		encode_original(read, quality, es); 
		break;
	case QualityComprMode::QuinaryAverage:
		encode_quinary_average(read, quality, es); 
		break;
	case QualityComprMode::QuadAverage:
		encode_quad_average(read, quality, es); 
		break;
	case QualityComprMode::BinaryAverage:
		encode_binary_average(read, quality, es); 
		break;
	case QualityComprMode::QuinaryThreshold:
		encode_quinary_threshold(read, quality, es); 
		break;
	case QualityComprMode::QuadThreshold:
		encode_quad_threshold(read, quality, es); 
		break;
	case QualityComprMode::BinaryThreshold:
		encode_binary_threshold(read, quality, es); 
		break;
	case QualityComprMode::Average:
		encode_average(read, quality, es); 
		break;
	case QualityComprMode::None:
		;
	}
}

// *****************************************************************
void CQualityCoder::Decode(read_t& read, qual_t &quality)
{
	uint32_t read_size = static_cast<uint32_t>(read.size()) - 1u;

	quality.reserve(read.size());

	if (quality_mode == QualityComprMode::None)
	{
		for (uint32_t i = 0; i < read_size; ++i)
			quality.push_back(33 + quality_rev_thrs[0]);

		return;
	}

	v_base_flags.clear();

	if(compression_level > 1)
		for (uint32_t i = 0; i < read_size; ++i)
		{
			v_base_flags.push_back(read[i] & ~base_mask_flags);
			read[i] &= base_mask_flags;
		}

	switch (quality_mode)
	{
	case QualityComprMode::Original:
		decode_original(read, quality);
		break;
	case QualityComprMode::QuinaryAverage:
		decode_quinary_average(read, quality);
		break;
	case QualityComprMode::QuadAverage:
		decode_quad_average(read, quality);
		break;
	case QualityComprMode::BinaryAverage:
		decode_binary_average(read, quality);
		break;
	case QualityComprMode::QuinaryThreshold:
		decode_quinary_threshold(read, quality);
		break;
	case QualityComprMode::QuadThreshold:
		decode_quad_threshold(read, quality);
		break;
	case QualityComprMode::BinaryThreshold:
		decode_binary_threshold(read, quality);
		break;
	case QualityComprMode::Average:
		decode_average(read, quality);
		break;
	case QualityComprMode::None:
		;
	}
}

}

// EOF
