#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// This file is based on the code from the MSAC project: https://sun.aei.polsl.pl/REFRESH/msac
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "defs.h"
#include "sub_rc.h"
#include "meta_switch.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <numeric>

//#define STATS_MODE

using namespace std;

// *******************************************************************************************
//
// *******************************************************************************************
#define CASE(n)		case n: left_freq += stats[n-1]

namespace entropy_coder
{
template<unsigned MAX_TOTAL, unsigned ADDER>
class CSimpleModel
{
	uint32_t n_symbols;
	uint32_t* stats;
	uint32_t total;

	void rescale()
	{
		while (total >= MAX_TOTAL)
		{
			total = 0;
			for (uint32_t i = 0; i < n_symbols; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public:
	CSimpleModel() : n_symbols(0), stats(nullptr), total(0)
	{
	};

	~CSimpleModel()
	{
		if (stats)
			delete[] stats;
	};

	CSimpleModel(const CSimpleModel& c)
	{
		n_symbols = c.n_symbols;
		total = c.total;
		stats = new uint32_t[n_symbols];

		for (uint32_t i = 0; i < n_symbols; ++i)
			stats[i] = c.stats[i];
	}

	CSimpleModel& operator=(const CSimpleModel&c)
	{
		if (this != &c)
		{
			if (stats)
				delete[] stats;

			n_symbols = c.n_symbols;
			total = c.total;

			stats = new uint32_t[n_symbols];

			for (uint32_t i = 0; i < n_symbols; ++i)
				stats[i] = c.stats[i];
		}

		return *this;
	}

	void Init(uint32_t _n_symbols, uint32_t* _init_stats)
	{
		if (stats)
		{
			if (n_symbols != _n_symbols)
			{
				delete[] stats;
				n_symbols = _n_symbols;
				stats = new uint32_t[n_symbols];
			}
		}
		else
		{
			n_symbols = _n_symbols;
			stats = new uint32_t[n_symbols];
		}

		if (_init_stats)
			for (uint32_t i = 0; i < n_symbols; ++i)
				stats[i] = _init_stats[i];
		else
			std::fill_n(stats, n_symbols, 1);

		total = std::accumulate(stats, stats + n_symbols, 0u);
		rescale();
	}

	void Init(const CSimpleModel& c)
	{
		n_symbols = c.n_symbols;

		if (stats)
			delete[] stats;

		stats = new uint32_t[n_symbols];
		std::copy_n(c.stats, n_symbols, stats);
		total = std::accumulate(stats, stats + n_symbols, 0u);
	}

	void GetFreq(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf)
	{
//		left_freq = 0;

//		if(symbol >= 100)
			left_freq = accumulate(stats, stats + symbol, 0u);
/*		else
		{
			switch (symbol)
			{
				CASE(99);	CASE(98);	CASE(97);	CASE(96);	CASE(95);	CASE(94);	CASE(93);	CASE(92);	CASE(91);	CASE(90);
				CASE(89);	CASE(88);	CASE(87);	CASE(86);	CASE(85);	CASE(84);	CASE(83);	CASE(82);	CASE(81);	CASE(80);
				CASE(79);	CASE(78);	CASE(77);	CASE(76);	CASE(75);	CASE(74);	CASE(73);	CASE(72);	CASE(71);	CASE(70);
				CASE(69);	CASE(68);	CASE(67);	CASE(66);	CASE(65);	CASE(64);	CASE(63);	CASE(62);	CASE(61);	CASE(60);
				CASE(59);	CASE(58);	CASE(57);	CASE(56);	CASE(55);	CASE(54);	CASE(53);	CASE(52);	CASE(51);	CASE(50);
				CASE(49);	CASE(48);	CASE(47);	CASE(46);	CASE(45);	CASE(44);	CASE(43);	CASE(42);	CASE(41);	CASE(40);
				CASE(39);	CASE(38);	CASE(37);	CASE(36);	CASE(35);	CASE(34);	CASE(33);	CASE(32);	CASE(31);	CASE(30);
				CASE(29);	CASE(28);	CASE(27);	CASE(26);	CASE(25);	CASE(24);	CASE(23);	CASE(22);	CASE(21);	CASE(20);
				CASE(19);	CASE(18);	CASE(17);	CASE(16);	CASE(15);	CASE(14);	CASE(13);	CASE(12);	CASE(11);	CASE(10);
				CASE(9);	CASE(8);	CASE(7);	CASE(6);	CASE(5);	CASE(4);	CASE(3);	CASE(2);	CASE(1);

			case 0: break;
			default:
				for (int i = 0; i < symbol; ++i)
					left_freq += stats[i];
			}
		}*/

		sym_freq = stats[symbol];
		totf = total;
	}

	void GetFreqExc(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf, const uint32_t exc)
	{
		left_freq = 0;

		for (uint32_t i = 0; i < symbol; ++i)
			if (i != exc)
				left_freq += stats[i];

		sym_freq = stats[symbol];
		totf = total - stats[exc];
	}

	void Update(uint32_t symbol)
	{
		stats[symbol] += ADDER;
		total += ADDER;

		if (total >= MAX_TOTAL)
			rescale();
	}

	uint32_t GetSym(uint32_t left_freq)
	{
		uint32_t t = 0;

		for (uint32_t i = 0; i < n_symbols; ++i)
		{
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return ~0u;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	uint32_t GetTotalExc(const uint32_t exc)
	{
		return total - stats[exc];
	}

	uint32_t* GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t* stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < n_symbols; ++i)
			total += stats[i] = stats_to_set[i];
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template <unsigned N_SYMBOLS, unsigned MAX_TOTAL, unsigned ADDER> class CSimpleModelFixedSize
{
	uint32_t stats[N_SYMBOLS];
	uint32_t total;
#ifdef STATS_MODE
	size_t no_updates;
#endif

	void rescale()
	{
		while (total >= MAX_TOTAL)
		{
			total = 0;
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public:
#ifdef STATS_MODE
	CSimpleModelFixedSize() : total(0), no_updates(0)
#else
	CSimpleModelFixedSize() : total(0)
#endif
	{
	};

	~CSimpleModelFixedSize()
	{
	};

	CSimpleModelFixedSize(const CSimpleModelFixedSize& c)
	{
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			stats[i] = c.stats[i];
		total = c.total;
	}

	CSimpleModelFixedSize& operator=(const CSimpleModelFixedSize& c)
	{
		if (this != &c)
		{
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				stats[i] = c.stats[i];
			total = c.total;
		}

		return *this;
	}

	void Init(const uint32_t* _init_stats)
	{
		if (_init_stats)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				stats[i] = _init_stats[i];
		else
			std::fill_n(stats, N_SYMBOLS, 1);

		total = std::accumulate(stats, stats + N_SYMBOLS, 0u);
		rescale();

#ifdef STATS_MODE
		no_updates = 0;
#endif
	}

	void Init(const CSimpleModelFixedSize& c)
	{
		std::copy_n(c.stats, N_SYMBOLS, stats);
		total = std::accumulate(stats, stats + N_SYMBOLS, 0u);

#ifdef STATS_MODE
		no_updates = 0;
#endif
	}

	void get_freq(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf)
	{
		left_freq = 0;

		if constexpr (N_SYMBOLS > 10)
			/*			for (uint32_t i = 0; i < symbol; ++i)
									left_freq += stats[i];*/
			left_freq = accumulate(stats, stats + symbol, 0u);
		else
			SwitchImpl<N_SYMBOLS>::get_freq(symbol, left_freq, stats);

		sym_freq = stats[symbol];
		totf = total;
	}

	void GetFreq(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf)
	{
		get_freq(symbol, sym_freq, left_freq, totf);
	}

	void GetFreqExc(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf, const uint32_t exc)
	{
		get_freq(symbol, sym_freq, left_freq, totf);

		totf -= stats[exc];

		if (exc < symbol)
			left_freq -= stats[exc];
	}

	void GetFreqExc(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf, const uint32_t exc1, const uint32_t exc2)
	{
		get_freq(symbol, sym_freq, left_freq, totf);

		totf -= stats[exc1];
		totf -= stats[exc2];

		if (exc1 < symbol)
			left_freq -= stats[exc1];
		if (exc2 < symbol)
			left_freq -= stats[exc2];
	}

	void Update(uint32_t symbol)
	{
		stats[symbol] += ADDER;
		total += ADDER;

		if (total >= MAX_TOTAL)
			rescale();

#ifdef STATS_MODE
		++no_updates;
#endif
	}

	uint32_t GetSym(uint32_t left_freq)
	{
		uint32_t t = 0;

		if constexpr (N_SYMBOLS > 10)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			{
				t += stats[i];
				if (t > left_freq)
					return i;
			}
		else
		{
			return IfImpl<N_SYMBOLS>::get_sym(left_freq, stats);
		}

		return ~0u;
	}

	void GetSymFreqAndUpdate(uint32_t& left_freq, uint32_t& sym, uint32_t& sym_freq)
	{
		if constexpr (N_SYMBOLS > 10)
		{
			uint32_t t = 0;
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			{
				t += stats[i];
				if (t > left_freq)
				{
					sym = i;
					left_freq = t - stats[i];
					break;
				}
			}
		}
		else
		{
			sym = IfImpl<N_SYMBOLS>::get_sym_lf(left_freq, stats);
		}

		sym_freq = stats[sym];
		Update(sym);
	}

	uint32_t GetSymExc(uint32_t left_freq, uint32_t exc)
	{
		uint32_t t = 0;

		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
		{
			if (i != exc)
				t += stats[i];
			if (t > left_freq)
				return i;
		}

		return ~0u;
	}

	// !!! TODO: can be without a loop
	uint32_t GetSymExc(uint32_t left_freq, uint32_t exc1, uint32_t exc2)
	{
		uint32_t t = 0;

		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
		{
			if (i != exc1 && i != exc2)
				t += stats[i];
			if (t > left_freq)
				return i;
		}

		return ~0u;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	uint32_t GetTotalExc(const uint32_t exc)
	{
		return total - stats[exc];
	}

	uint32_t GetTotalExc(const uint32_t exc1, const uint32_t exc2)
	{
		return total - stats[exc1] - stats[exc2];
	}

	uint32_t* GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t* stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			total += stats[i] = stats_to_set[i];
	}

#ifdef STATS_MODE
	void GetLogStats(size_t& _no_updates, std::vector<float>& _v_freq)
	{
		_no_updates = no_updates;

		_v_freq.resize(N_SYMBOLS, 0u);
		auto sum = accumulate(stats, stats + N_SYMBOLS, 0u);

		if (sum)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				_v_freq[i] = (float)stats[i] / sum;
	}

	size_t GetNoUpdates()
	{
		return no_updates;
	}
#endif
};


// *******************************************************************************************
//
// *******************************************************************************************
// Max value of N_SYMBOLS = 256
template <unsigned N_SYMBOLS, unsigned MAX_TOTAL, unsigned ADDER> class CFenwickTreeModelFixedSize
{
	uint32_t arr[N_SYMBOLS];
	uint32_t total;

	static constexpr uint32_t LSB[] = { 0, 1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8,
		16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
		64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
		64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
		128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
		128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
		128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 256 };

	uint32_t prefix_sum(uint32_t symbol)
	{
		uint32_t sum = 0;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];
		symbol &= symbol - 1;

		if (!symbol)		return sum;
		sum += arr[symbol - 1];

		return sum;
/*		for (; symbol > 0; symbol &= symbol - 1)
			sum += arr[symbol - 1];*/

		return sum;
	}

	void increment(uint32_t symbol)
	{
		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		symbol |= symbol + 1;

		if (symbol >= N_SYMBOLS)	return;
		arr[symbol] += ADDER;
		
/*		for (; symbol < N_SYMBOLS; symbol |= symbol + 1)
			arr[symbol] += delta;*/
	}

	// Convert counts in arr into Fenwick tree
	void make_fenwick_tree(void)
	{
		for (uint32_t i = 0; i < N_SYMBOLS; i++) 
		{
			uint32_t j = i | (i + 1);
			
			if (j < N_SYMBOLS)
				arr[j] += arr[i];
		}
	}

	// Convert Fenwick tree in arr into array of counts
	void unmake_fenwick_tree(void)
	{
		for (uint32_t i = N_SYMBOLS; i-- > 0;) 
		{
			uint32_t j = i | (i + 1);
			
			if (j < N_SYMBOLS)
				arr[j] -= arr[i];
		}
	}

	uint32_t value(uint32_t symbol)
	{
		uint32_t sum = arr[symbol];
		uint32_t j = symbol & (symbol + 1);

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];
		symbol &= symbol - 1;

		if (symbol <= j)		return sum;
		sum -= arr[symbol - 1];

		return sum;
	}

	// *******************
	void rescale()
	{
		unmake_fenwick_tree();

		while (total >= MAX_TOTAL)
		{
			total = 0;
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			{
				arr[i] = (arr[i] + 1) / 2;
				total += arr[i];
			}
		}

		make_fenwick_tree();
	}

public:
	CFenwickTreeModelFixedSize() : total(0)
	{
		static_assert(N_SYMBOLS <= 256);
	};

	~CFenwickTreeModelFixedSize()
	{
	};

	CFenwickTreeModelFixedSize(const CFenwickTreeModelFixedSize& c)
	{
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			arr[i] = c.arr[i];
		total = c.total;
	}

	CFenwickTreeModelFixedSize& operator=(const CFenwickTreeModelFixedSize& c)
	{
		if (this != &c)
		{
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				arr[i] = c.arr[i];
			total = c.total;
		}

		return *this;
	}

	void Init(const uint32_t* _init_stats)
	{
		if (_init_stats)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				arr[i] = _init_stats[i];
		else
			std::fill_n(arr, N_SYMBOLS, 1);

		total = std::accumulate(arr, arr + N_SYMBOLS, 0u);
//		rescale();

		make_fenwick_tree();
	}

	void Init(const CFenwickTreeModelFixedSize& c)
	{
		std::copy_n(c.arr, N_SYMBOLS, arr);
		total = c.total;
	}

	void GetFreq(uint32_t symbol, uint32_t& sym_freq, uint32_t& left_freq, uint32_t& totf)
	{
		left_freq = prefix_sum(symbol);
		sym_freq = value(symbol);
		totf = total;
	}

	void Update(uint32_t symbol)
	{
		increment(symbol);
		total += ADDER;

		if (total >= MAX_TOTAL)
			rescale();
	}

	uint32_t GetSym(uint32_t left_freq)
	{
		uint32_t i = 0;
		uint32_t j = N_SYMBOLS;

		j = LSB[j];

		for (; j > 0; j >>= 1)
		{
			if (i + j <= N_SYMBOLS && arr[i + j - 1] <= left_freq)
			{
				left_freq -= arr[i + j - 1];
				i += j;
			}
		}

		return i;
	}

/*	void GetSymFreqAndUpdate(int& left_freq, int& sym, int& sym_freq)
	{
	}*/

	uint32_t GetTotal()
	{
		return total;
	}
};


// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM, unsigned MAX_TOTAL, unsigned ADDER> class CRangeCoderModel
{
	union {
		CRangeEncoder<T_IO_STREAM>* rce;
		CRangeDecoder<T_IO_STREAM>* rcd;
	} rc;

	CSimpleModel<MAX_TOTAL, ADDER> simple_model;

	uint32_t no_symbols;

public:
	CRangeCoderModel(CBasicRangeCoder<T_IO_STREAM>* rcb, uint32_t _no_symbols, uint32_t* _init, bool compress) :
		no_symbols(_no_symbols)
	{
		simple_model.Init(no_symbols, _init);

		if (compress)
			rc.rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
		else
			rc.rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
	}

	CRangeCoderModel(const CRangeCoderModel& c)
	{
		simple_model.Init(c.simple_model);
		rc = c.rc;

		no_symbols = c.no_symbols;
//		totf = c.totf;
//		rescale = c.rescale;
	}

	~CRangeCoderModel()
	{
	}

	void Encode(const uint32_t x)
	{
		uint32_t syfreq, ltfreq, totf;
		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rc.rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	void EncodeExcluding(const uint32_t x, const uint32_t exc)
	{
		uint32_t syfreq, ltfreq, totf;
		simple_model.GetFreqExc(x, syfreq, ltfreq, totf, exc);
		rc.rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	uint32_t Decode()
	{
		uint32_t syfreq, ltfreq, totf;

		totf = simple_model.GetTotal();
		ltfreq = static_cast<uint32_t>(rc.rcd->GetCumulativeFreq(totf));

		uint32_t x = simple_model.GetSym(ltfreq);

		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rc.rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

	uint32_t DecodeExcluding(const uint32_t exc)
	{
		uint32_t syfreq, ltfreq, totf;

		totf = simple_model.GetTotalExc(exc);
		ltfreq = static_cast<uint32_t>(rc.rcd->GetCumulativeFreq(totf));

		uint32_t x = simple_model.GetSym(ltfreq);

		simple_model.GetFreqExc(x, syfreq, ltfreq, totf, exc);
		rc.rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

/*	CSimpleModel<MAX_TOTAL, ADDER>* GetSimpleModel()
	{
		return &simple_model;
	}*/

	void Init(uint32_t* init)
	{
		simple_model.Init(no_symbols, init);
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM, unsigned N_SYMBOLS, unsigned MAX_TOTAL, unsigned ADDER> class CRangeCoderModelFixedSize
{
	union {
		CRangeEncoder<T_IO_STREAM>* rce;
		CRangeDecoder<T_IO_STREAM>* rcd;
	} rc;

	CSimpleModelFixedSize<N_SYMBOLS, MAX_TOTAL, ADDER> simple_model;

public:
	CRangeCoderModelFixedSize(CBasicRangeCoder<T_IO_STREAM>* rcb, uint32_t* _init, bool compress)
	{
		simple_model.Init(_init);

		if (compress)
			rc.rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
		else
			rc.rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
	}

	CRangeCoderModelFixedSize(const CRangeCoderModelFixedSize& c)
	{
		simple_model.Init(c.simple_model);
		rc = c.rc;
	}

	~CRangeCoderModelFixedSize()
	{
	}

	void Encode(const uint32_t x)
	{
		uint32_t syfreq, ltfreq, totf;

		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rc.rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	void EncodeExcluding(const uint32_t x, const uint32_t exc)
	{
		uint32_t syfreq, ltfreq, totf;

		simple_model.GetFreqExc(x, syfreq, ltfreq, totf, exc);
		rc.rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	void EncodeExcluding(const uint32_t x, const uint32_t exc1, const uint32_t exc2)
	{
		uint32_t syfreq, ltfreq, totf;

		simple_model.GetFreqExc(x, syfreq, ltfreq, totf, exc1, exc2);
		rc.rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	uint32_t Decode()
	{
		uint32_t syfreq, ltfreq;
		uint32_t totf = simple_model.GetTotal();
		uint32_t x = 0;

		ltfreq = static_cast<uint32_t>(rc.rcd->GetCumulativeFreq(totf));

//		int x = simple_model.GetSym(ltfreq);
//		simple_model.GetFreq(x, syfreq, ltfreq, totf);

		simple_model.GetSymFreqAndUpdate(ltfreq, x, syfreq);
		rc.rcd->UpdateFrequency(syfreq, ltfreq, totf);
//		simple_model.Update(x);

		return x;
	}

	uint32_t DecodeExcluding(const uint32_t exc)
	{
		uint32_t syfreq;
		uint32_t totf = simple_model.GetTotalExc(exc);
		uint32_t ltfreq = static_cast<uint32_t>(rc.rcd->GetCumulativeFreq(totf));

		uint32_t x = simple_model.GetSymExc(ltfreq, exc);

		simple_model.GetFreqExc(x, syfreq, ltfreq, totf, exc);
		rc.rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

	uint32_t DecodeExcluding(const uint32_t exc1, const uint32_t exc2)
	{
		uint32_t syfreq;
		uint32_t totf = simple_model.GetTotalExc(exc1, exc2);
		uint32_t ltfreq = static_cast<uint32_t>(rc.rcd->GetCumulativeFreq(totf));

		uint32_t x = simple_model.GetSymExc(ltfreq, exc1, exc2);

		simple_model.GetFreqExc(x, syfreq, ltfreq, totf, exc1, exc2);
		rc.rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

/*	CSimpleModelFixedSize<N_SYMBOLS>* GetSimpleModel()
	{
		return &simple_model;
	}*/
};


// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM, unsigned N_SYMBOLS, unsigned MAX_TOTAL, unsigned ADDER> class CRangeCoderFenwickTreeFixedSize
{
	union {
		CRangeEncoder<T_IO_STREAM>* rce;
		CRangeDecoder<T_IO_STREAM>* rcd;
	} rc;

	CFenwickTreeModelFixedSize<N_SYMBOLS, MAX_TOTAL, ADDER> fenwick_tree_model;

public:
	CRangeCoderFenwickTreeFixedSize(CBasicRangeCoder<T_IO_STREAM>* rcb, uint32_t* _init, bool compress)
	{
		fenwick_tree_model.Init(_init);

		if (compress)
			rc.rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
		else
			rc.rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
	}

	CRangeCoderFenwickTreeFixedSize(const CRangeCoderFenwickTreeFixedSize& c)
	{
		fenwick_tree_model.Init(c.fenwick_tree_model);
		rc = c.rc;
	}

	~CRangeCoderFenwickTreeFixedSize()
	{
	}

	void Encode(const uint32_t x)
	{
		uint32_t syfreq, ltfreq, totf;

		fenwick_tree_model.GetFreq(x, syfreq, ltfreq, totf);
		rc.rce->EncodeFrequency(syfreq, ltfreq, totf);

		fenwick_tree_model.Update(x);
	}

	uint32_t Decode()
	{
		uint32_t syfreq;
		uint32_t totf = fenwick_tree_model.GetTotal();
		uint32_t ltfreq = static_cast<uint32_t>(rc.rcd->GetCumulativeFreq(totf));

		uint32_t x = fenwick_tree_model.GetSym(ltfreq);

		fenwick_tree_model.GetFreq(x, syfreq, ltfreq, totf);
		rc.rcd->UpdateFrequency(syfreq, ltfreq, totf);
		fenwick_tree_model.Update(x);

		return x;
	}
};

#undef CASE
} // namespace entropy_coder
// EOF
