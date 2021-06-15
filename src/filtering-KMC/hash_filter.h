#pragma once
#include "../common/libs/libdivide/libdivide.h"
#include <cinttypes>

class CHashModuloFilter
{
	uint32_t modulo;
	uint64_t hash_mm(uint64_t x) const
	{
		x ^= x >> 33;
		x *= 0xff51afd7ed558ccduLL;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53uLL;
		x ^= x >> 33;
		return x;
	}

	libdivide::divider<uint64_t> div;

public:
	CHashModuloFilter(uint32_t modulo) :
		modulo(modulo),
		div(modulo)
	{

	}
	
	bool checkModuloHash(uint64_t kmer) const
	{	
		switch (modulo)
		{
			case 1:
				return true;
			case 2:
				return hash_mm(kmer) % 2 == 0;
			case 3:
				return hash_mm(kmer) % 3 == 0;
			case 4:
				return hash_mm(kmer) % 4 == 0;
			case 5:
				return hash_mm(kmer) % 5 == 0;
			case 6:
				return hash_mm(kmer) % 6 == 0;
			case 7:
				return hash_mm(kmer) % 7 == 0;
			case 8:
				return hash_mm(kmer) % 8 == 0;
			case 9:
				return hash_mm(kmer) % 9 == 0;
			case 10:
				return hash_mm(kmer) % 10 == 0;
			case 11:
				return hash_mm(kmer) % 11 == 0;
			case 12:
				return hash_mm(kmer) % 12 == 0;
			case 13:
				return hash_mm(kmer) % 13 == 0;
			case 14:
				return hash_mm(kmer) % 14 == 0;
			case 15:
				return hash_mm(kmer) % 15 == 0;
			case 16:
				return hash_mm(kmer) % 16 == 0;
			case 17:
				return hash_mm(kmer) % 17 == 0;
			case 18:
				return hash_mm(kmer) % 18 == 0;
			case 19:
				return hash_mm(kmer) % 19 == 0;
			case 20:
				return hash_mm(kmer) % 20 == 0;
			default:
			{
				auto h = hash_mm(kmer);
				return h - modulo * (h / div) == 0;		
			}
		}
	}
};
