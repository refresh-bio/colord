#pragma once
#include "basic_coder.h"
#include "params.h"

#include <array>

namespace entropy_coder
{
class CQualityCoder : public CBasicCoder
{
	bool verbose;
	QualityComprMode quality_mode;
	DataSource data_source;
	vector<uint32_t> quality_fwd_thrs;
	vector<uint32_t> quality_rev_thrs;

//	using rcmfs_quality_original_t = CRangeCoderModelFixedSize<CVectorIOStream, 96, 1 << 20, 32>;
	using rcmfs_quality_original_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 96, 1 << 20, 32>;
	using rcmfs_quality_quad_t = CRangeCoderModelFixedSize<CVectorIOStream, 4, 1 << 18, 8>;
	using rcmfs_quality_quinary_t = CRangeCoderModelFixedSize<CVectorIOStream, 5, 1 << 18, 8>;
	using rcmfs_quality_binary_t = CRangeCoderModelFixedSize<CVectorIOStream, 2, 1 << 18, 8>;
//	using rcmfs_quality_byte_t = CRangeCoderModelFixedSize<CVectorIOStream, 256, 1 << 18, 8>;
	using rcmfs_quality_byte_t = CRangeCoderFenwickTreeFixedSize<CVectorIOStream, 256, 1 << 18, 8>;

	using ctx_map_quality_original_t = CContextVec<rcmfs_quality_original_t>;
	using ctx_map_quality_byte_t = CContextHM<rcmfs_quality_byte_t>;
	using ctx_vece_quality_binary_t = CContextVecEmb<rcmfs_quality_binary_t>;
	using ctx_vece_quality_quad_t = CContextVecEmb<rcmfs_quality_quad_t>;
	using ctx_vece_quality_quinary_t = CContextVecEmb<rcmfs_quality_quinary_t>;

	array<uint32_t, 96> quality_code_map_fwd;
	array<uint32_t, 96> quality_code_map_rev;
	array<uint32_t, 96> quality_code_quantize;

	uint32_t no_ctx_symbols;
	uint32_t no_bits_per_symbol;
	uint32_t ctx_size_in_bits;
	int compression_level;
	uint64_t input_stream_size;

	context_t ctx_mask;
	context_t ctx_inc;

	ctx_vece_quality_binary_t m_quality_binary_rc;
	ctx_vece_quality_quad_t m_quality_quad_rc;
	ctx_vece_quality_quinary_t m_quality_quinary_rc;
	ctx_map_quality_original_t m_quality_original_rc;
	ctx_map_quality_byte_t m_quality_byte_rc;

	vector<uint8_t> v_base_flags;

	rcmfs_quality_binary_t* tpl_quality_binary_rc;
	rcmfs_quality_quad_t* tpl_quality_quad_rc;
	rcmfs_quality_quinary_t* tpl_quality_quinary_rc;
	rcmfs_quality_original_t* tpl_quality_original_rc;
	rcmfs_quality_byte_t* tpl_quality_byte_rc;

	void reset_context(context_t& ctx);
	void update_context(context_t& ctx, uint8_t x);

	void adjust_quality_map_symbols(uint32_t n);

	void adjust_quality_map_ONT_lossless();
	void adjust_quality_map_ONT_quad();
	void adjust_quality_map_ONT_quinary();
	void adjust_quality_map_ONT_binary();
	void adjust_quality_map_PBRaw_lossless();
	void adjust_quality_map_PBRaw_quad();
	void adjust_quality_map_PBRaw_quinary();
	void adjust_quality_map_PBRaw_binary();
	void adjust_quality_map_PBHiFi_lossless();
	void adjust_quality_map_PBHiFi_quad();
	void adjust_quality_map_PBHiFi_quinary();
	void adjust_quality_map_PBHiFi_binary();

	void encode_original(const read_t& read, const qual_t& quality, es_t& es);
	void encode_quinary_average(const read_t& read, const qual_t& quality, es_t& es);
	void encode_quad_average(const read_t& read, const qual_t& quality, es_t& es);
	void encode_binary_average(const read_t& read, const qual_t& quality, es_t& es);
	void encode_quinary_threshold(const read_t& read, const qual_t& quality, es_t& es);
	void encode_quad_threshold(const read_t& read, const qual_t& quality, es_t& es);
	void encode_binary_threshold(const read_t& read, const qual_t& quality, es_t& es);
	void encode_average(const read_t& read, const qual_t& quality, es_t& es);

	void decode_original(read_t& read, qual_t& quality);
	void decode_quinary_average(read_t& read, qual_t& quality);
	void decode_quad_average(read_t& read, qual_t& quality);
	void decode_binary_average(read_t& read, qual_t& quality);
	void decode_quinary_threshold(read_t& read, qual_t& quality);
	void decode_quad_threshold(read_t& read, qual_t& quality);
	void decode_binary_threshold(read_t& read, qual_t& quality);
	void decode_average(read_t& read, qual_t& quality);

	void analyze_es(es_t& es, uint32_t read_size);

	constexpr context_t quantize_diff(uint8_t q1, uint8_t q2) {
		int qd = (int)q1 - (int)q2;

		if (qd == 0)
			return 0;
		if (qd == 1)
			return 1;
		if (qd == -1)
			return 2;
		if (qd >= 2 && qd <= 4)
			return 3;
		if (-qd >= 2 && -qd <= 4)
			return 4;
		if (qd >= 5 && qd <= 10)
			return 5;
		if (-qd >= 5 && -qd <= 10)
			return 6;
		if (qd >= 5)
			return 7;
		return 8;
	}

	void encode_avg(context_t ctx_base, double x);
	double decode_avg(context_t ctx_base);

public:
	CQualityCoder(bool verbose) : CBasicCoder(), verbose(verbose), quality_mode(QualityComprMode::Original), data_source(DataSource::ONT), quality_fwd_thrs({ 7 }), quality_rev_thrs({ 1, 13 })
	{
		tpl_quality_binary_rc = nullptr;
		tpl_quality_quad_rc = nullptr;
		tpl_quality_original_rc = nullptr;
		tpl_quality_byte_rc = nullptr;
		tpl_quality_quinary_rc = nullptr;

		for (int i = 0; i < 96; ++i)
		{
			quality_code_map_fwd[i] = 0;
			quality_code_map_rev[i] = 0;
		}
	}

	~CQualityCoder()
	{
		if(verbose)
		{
			cout << "Quality quality_binary size  : " << m_quality_binary_rc.get_size() << endl;
			cout << "Quality quality_quad size    : " << m_quality_quad_rc.get_size() << endl;
			cout << "Quality quality_quinary size : " << m_quality_quinary_rc.get_size() << endl;
			cout << "Quality quality_original size: " << m_quality_original_rc.get_size() << endl;
			cout << "Quality quality_byte size    : " << m_quality_byte_rc.get_size() << endl;
		}

		delete tpl_quality_binary_rc;
		delete tpl_quality_quinary_rc;
		delete tpl_quality_quad_rc;
		delete tpl_quality_original_rc;
		delete tpl_quality_byte_rc;
	}

	// Initialization - required before (de)compression
	void Init(bool _is_compressing, QualityComprMode _quality_mode, DataSource _data_source, const std::vector<uint32_t>& _quality_fwd_thrs, const std::vector<uint32_t>& _quality_rev_thrs, int _compression_level, uint64_t _input_stream_size);

	// compression restart without model clenup
	void Restart();

	// Should be called after (de)compression
	void Finish();

	void Encode(const read_t &read, const qual_t &quality, es_t& es);

	// Read quality decoding
	void Decode(read_t& read, qual_t &quality);
};

}

// EOF