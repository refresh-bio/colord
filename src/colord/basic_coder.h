#pragma once

#include<vector>
#include <map>

#include "utils.h"
#include "rc.h"
#include "sub_rc.h"
#include "context_hm.h"

using namespace std;


namespace entropy_coder
{
const uint8_t base_flag_anchor = 0b10000000;
const uint8_t base_flag_match = 0b01000000;
const uint8_t base_mask_flags = 0b00111111;

// *****************************************************************
constexpr uint64_t ilog2(uint64_t x)
{
	uint64_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

// *****************************************************************
constexpr uint64_t no_bytes(uint64_t x)
{
	uint64_t r = 1;

	x >>= 8;

	for (; x; ++r)
		x >>= 8;

	return r;
}

// *******************************************************************************************
// Class for storage of range coder compressed data
class CVectorIOStream
{
	vector<uint8_t>& v;
	size_t read_pos;

public:
	CVectorIOStream(vector<uint8_t>& _v) : v(_v), read_pos(0)
	{}

	void RestartRead()
	{
		read_pos = 0;
	}

	bool Eof() const
	{
		return read_pos >= v.size();
	}

	uint8_t GetByte()
	{
		return v[read_pos++];
	}

	void PutByte(uint8_t x)
	{
		v.push_back(x);
	}

	template<typename T>
	void PutByte(T x)
	{
		v.push_back(static_cast<uint8_t>(x));
	}

	size_t Size()
	{
		return v.size();
	}
};


// *******************************************************************************************
class CBasicCoder
{
protected:
	CRangeEncoder<CVectorIOStream>* rce;
	CRangeDecoder<CVectorIOStream>* rcd;

	CVectorIOStream* v_vios_io;
	vector<uint8_t> v_io;

	template<typename T, typename V>
	V* find_rc_context(T& m_ctx_rc, context_t ctx, V* tpl)
	{
		auto p = m_ctx_rc.find(ctx);

		if (p == nullptr)
			m_ctx_rc.insert(ctx, p = new V(*tpl));

		return p;
	}

	template<typename T, typename V>
	V* find_rce_context(T& m_ctx_rc, context_t ctx, V* tpl)
	{
		auto p = m_ctx_rc.find(ctx);

		if (p == nullptr)
			p = m_ctx_rc.insert(ctx, *tpl);

		return p;
	}

	bool is_compressing;

public:
	CBasicCoder()
	{
		rce = nullptr;
		rcd = nullptr;

		v_vios_io = nullptr;
	}

	virtual ~CBasicCoder()
	{
		delete rce;
		delete rcd;

		delete v_vios_io;
	}

	// Reading compression output (after Finish)
	void GetOutput(vector<uint8_t>& _v_io)
	{
		_v_io = move(v_io);
	}

	// Setting compressed input (before Init)
	void SetInput(vector<uint8_t>& _v_io)
	{
		v_io = move(_v_io);
	}
};

}

// EOF