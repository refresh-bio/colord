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

#include <iostream>
#include <cmath>
#include <cstdint>

#include "defs.h"
#include <assert.h>

#define UNROLL_FREQUENCY_CODING

#define RC_64BIT

namespace entropy_coder
{
// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CBasicRangeCoder
{
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CRangeEncoder : public CBasicRangeCoder<T_IO_STREAM>
{
public:
#ifdef RC_64BIT
	typedef uint64_t Code;
	typedef uint64_t Freq;

	const Freq TopValue = 0x00ffffffffffffULL;
	const Code Mask = 0xff00000000000000ULL;
#else
	typedef uint32_t Code;
	typedef uint32_t Freq;

	const Freq TopValue = 0x00ffffULL;
	const Code Mask = 0xff000000ULL;
#endif

	const uint32_t MAX_SHIFT = 8 * sizeof(Code) - 8;
	const uint32_t CODE_BYTES = sizeof(Code);

	T_IO_STREAM& io_stream;

	Code low;
	Freq range;

	CRangeEncoder(T_IO_STREAM& _io_stream) : io_stream(_io_stream), low(0), range(0)
	{}

	void Start()
	{
		low = 0;
		range = Mask;
	}

/*	float EstimateCodeLen(Freq s, Freq t)
	{
		return -log2f((float)s / t);
	}*/

	void EncodeFrequency(Freq symFreq_, Freq cumFreq_, Freq totalFreqSum_)
	{
		assert(range > totalFreqSum_);
		range /= totalFreqSum_;
		low += range * cumFreq_;
		range *= symFreq_;

#ifndef UNROLL_FREQUENCY_CODING
		while (range <= TopValue)
		{
			assert(range != 0);
			if ((low ^ (low + range)) & Mask)
			{
				Freq r = (Freq)low;
				range = (r | TopValue) - r;
			}
			io_stream.PutByte(low >> MAX_SHIFT);
			low <<= 8, range <<= 8;
		}
#else
		if (range <= TopValue)
		{
			assert(range != 0);
			if ((low ^ (low + range)) & Mask)
			{
				Freq r = (Freq)low;
				range = (r | TopValue) - r;
			}
			io_stream.PutByte(low >> MAX_SHIFT);
			low <<= 8, range <<= 8;

			if (range <= TopValue)
			{
				assert(range != 0);
				if ((low ^ (low + range)) & Mask)
				{
					Freq r = (Freq)low;
					range = (r | TopValue) - r;
				}
				io_stream.PutByte(low >> MAX_SHIFT);
				low <<= 8, range <<= 8;

				if (range <= TopValue)
				{
					assert(range != 0);
					if ((low ^ (low + range)) & Mask)
					{
						Freq r = (Freq)low;
						range = (r | TopValue) - r;
					}
					io_stream.PutByte(low >> MAX_SHIFT);
					low <<= 8, range <<= 8;

					if (range <= TopValue)
					{
						assert(range != 0);
						if ((low ^ (low + range)) & Mask)
						{
							Freq r = (Freq)low;
							range = (r | TopValue) - r;
						}
						io_stream.PutByte(low >> MAX_SHIFT);
						low <<= 8, range <<= 8;

#ifdef RC_64BIT
						if (range <= TopValue)
						{
							assert(range != 0);
							if ((low ^ (low + range)) & Mask)
							{
								Freq r = (Freq)low;
								range = (r | TopValue) - r;
							}
							io_stream.PutByte(low >> MAX_SHIFT);
							low <<= 8, range <<= 8;

							if (range <= TopValue)
							{
								assert(range != 0);
								if ((low ^ (low + range)) & Mask)
								{
									Freq r = (Freq)low;
									range = (r | TopValue) - r;
								}
								io_stream.PutByte(low >> MAX_SHIFT);
								low <<= 8, range <<= 8;

								if (range <= TopValue)
								{
									assert(range != 0);
									if ((low ^ (low + range)) & Mask)
									{
										Freq r = (Freq)low;
										range = (r | TopValue) - r;
									}
									io_stream.PutByte(low >> MAX_SHIFT);
									low <<= 8, range <<= 8;

									if (range <= TopValue)
									{
										assert(range != 0);
										if ((low ^ (low + range)) & Mask)
										{
											Freq r = (Freq)low;
											range = (r | TopValue) - r;
										}
										io_stream.PutByte(low >> MAX_SHIFT);
										low <<= 8, range <<= 8;
									}
								}
							}
						}
#endif
					}
				}
			}
		}
#endif
	}

	void End()
	{
		for (uint32_t i = 0; i < CODE_BYTES; i++)
		{
			io_stream.PutByte(low >> MAX_SHIFT);
			low <<= 8;
		}
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CRangeDecoder : public CBasicRangeCoder<T_IO_STREAM>
{
public:
#ifdef RC_64BIT
	typedef uint64_t Code;
	typedef uint64_t Freq;

	const Freq TopValue = 0x00ffffffffffffULL;
	const Code Mask = 0xff00000000000000ULL;
#else
	typedef uint32_t Code;
	typedef uint32_t Freq;

	const Freq TopValue = 0x00ffffULL;
	const Code Mask = 0xff000000ULL;
#endif

	const uint32_t MAX_SHIFT = 8 * sizeof(Code) - 8;
	const uint32_t CODE_BYTES = sizeof(Code);
	const uint32_t CODE_BITS = 8 * sizeof(Code);

	T_IO_STREAM& io_stream;

	Code	low;
	Freq	range;

	CRangeDecoder(T_IO_STREAM& _io_stream) : io_stream(_io_stream), low(0), range(0)
	{
		buffer = 0;
	}

	void Start()
	{
		if (io_stream.Size() < CODE_BYTES)
			return;

		buffer = 0;
		for (uint32_t i = 1; i <= CODE_BYTES; ++i)
		{
			buffer |= (Code)io_stream.GetByte() << (CODE_BITS - i * 8);
		}

		low = 0;
		range = Mask;
	}

	Freq GetCumulativeFreq(Freq totalFreq_)
	{
		assert(totalFreq_ != 0);
		return (Freq)(buffer / (range /= totalFreq_));
	}

	void UpdateFrequency(Freq symFreq_, Freq lowEnd_, Freq /*totalFreq_*/)
	{
		Freq r = lowEnd_ * range;
		buffer -= r;
		low += r;
		range *= symFreq_;

#ifndef UNROLL_FREQUENCY_CODING
		while (range <= TopValue)
		{
			if ((low ^ (low + range)) & Mask)
			{
				Freq r = (Freq)low;
				range = (r | TopValue) - r;
			}

			buffer = (buffer << 8) + io_stream.GetByte();
			low <<= 8, range <<= 8;
		}
#else
		if (range <= TopValue)
		{
			if ((low ^ (low + range)) & Mask)
			{
				Freq r = (Freq)low;
				range = (r | TopValue) - r;
			}

			buffer = (buffer << 8) + io_stream.GetByte();
			low <<= 8, range <<= 8;

			if (range <= TopValue)
			{
				if ((low ^ (low + range)) & Mask)
				{
					Freq r = (Freq)low;
					range = (r | TopValue) - r;
				}

				buffer = (buffer << 8) + io_stream.GetByte();
				low <<= 8, range <<= 8;

				if (range <= TopValue)
				{
					if ((low ^ (low + range)) & Mask)
					{
						Freq r = (Freq)low;
						range = (r | TopValue) - r;
					}

					buffer = (buffer << 8) + io_stream.GetByte();
					low <<= 8, range <<= 8;

					if (range <= TopValue)
					{
						if ((low ^ (low + range)) & Mask)
						{
							Freq r = (Freq)low;
							range = (r | TopValue) - r;
						}

						buffer = (buffer << 8) + io_stream.GetByte();
						low <<= 8, range <<= 8;

#ifdef RC_64BIT
						if (range <= TopValue)
						{
							if ((low ^ (low + range)) & Mask)
							{
								Freq r = (Freq)low;
								range = (r | TopValue) - r;
							}

							buffer = (buffer << 8) + io_stream.GetByte();
							low <<= 8, range <<= 8;

							if (range <= TopValue)
							{
								if ((low ^ (low + range)) & Mask)
								{
									Freq r = (Freq)low;
									range = (r | TopValue) - r;
								}

								buffer = (buffer << 8) + io_stream.GetByte();
								low <<= 8, range <<= 8;

								if (range <= TopValue)
								{
									if ((low ^ (low + range)) & Mask)
									{
										Freq r = (Freq)low;
										range = (r | TopValue) - r;
									}

									buffer = (buffer << 8) + io_stream.GetByte();
									low <<= 8, range <<= 8;

									if (range <= TopValue)
									{
										if ((low ^ (low + range)) & Mask)
										{
											Freq r = (Freq)low;
											range = (r | TopValue) - r;
										}

										buffer = (buffer << 8) + io_stream.GetByte();
										low <<= 8, range <<= 8;
									}
								}
							}
						}
#endif
					}
				}
			}
		}
#endif
	}

	void End()
	{}

private:
	Code buffer;
};
} // namespace entropy_coder
// EOF
