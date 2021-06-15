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
#include <cinttypes>

template<unsigned N>
struct SwitchImpl
{};

template<>
struct SwitchImpl<10>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{
		case 9: left_freq += stats[8];
		case 8: left_freq += stats[7];
		case 7: left_freq += stats[6];
		case 6: left_freq += stats[5];
		case 5: left_freq += stats[4];
		case 4: left_freq += stats[3];
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};

template<>
struct SwitchImpl<9>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 8: left_freq += stats[7];
		case 7: left_freq += stats[6];
		case 6: left_freq += stats[5];
		case 5: left_freq += stats[4];
		case 4: left_freq += stats[3];
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};


template<>
struct SwitchImpl<8>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 7: left_freq += stats[6];
		case 6: left_freq += stats[5];
		case 5: left_freq += stats[4];
		case 4: left_freq += stats[3];
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};

template<>
struct SwitchImpl<7>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 6: left_freq += stats[5];
		case 5: left_freq += stats[4];
		case 4: left_freq += stats[3];
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};


template<>
struct SwitchImpl<6>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{
		case 5: left_freq += stats[4];
		case 4: left_freq += stats[3];
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};


template<>
struct SwitchImpl<5>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 4: left_freq += stats[3];
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};


template<>
struct SwitchImpl<4>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 3: left_freq += stats[2];
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};

template<>
struct SwitchImpl<3>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 2: left_freq += stats[1];
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};


template<>
struct SwitchImpl<2>
{
	static void get_freq(uint32_t symbol, uint32_t& left_freq, uint32_t* stats)
	{
		switch (symbol)
		{		
		case 1: left_freq += stats[0];
		case 0: break;
		}
	}
};


template<unsigned N>
struct IfImpl
{};

template<>
struct IfImpl<10>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		t += stats[3];
		if (left_freq < t)		return 3;
		t += stats[4];
		if (left_freq < t)		return 4;
		t += stats[5];
		if (left_freq < t)		return 5;
		t += stats[6];
		if (left_freq < t)		return 6;
		t += stats[7];
		if (left_freq < t)		return 7;
		t += stats[8];
		if (left_freq < t)		return 8;
		return 9;
	}

	static uint32_t get_sym_lf(uint32_t&left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		t += stats[3];
		if (left_freq < t) { left_freq = t - stats[3];  return 3; }
		t += stats[4];
		if (left_freq < t) { left_freq = t - stats[4];  return 4; }
		t += stats[5];
		if (left_freq < t) { left_freq = t - stats[5];  return 5; }
		t += stats[6];
		if (left_freq < t) { left_freq = t - stats[6];  return 6; }
		t += stats[7];
		if (left_freq < t) { left_freq = t - stats[7];  return 7; }
		t += stats[8];
		if (left_freq < t) { left_freq = t - stats[8];  return 8; }
		left_freq = t;
		return 9;
	}
};

template<>
struct IfImpl<9>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		t += stats[3];
		if (left_freq < t)		return 3;
		t += stats[4];
		if (left_freq < t)		return 4;
		t += stats[5];
		if (left_freq < t)		return 5;
		t += stats[6];
		if (left_freq < t)		return 6;
		t += stats[7];
		if (left_freq < t)		return 7;
		return 8;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		t += stats[3];
		if (left_freq < t) { left_freq = t - stats[3];  return 3; }
		t += stats[4];
		if (left_freq < t) { left_freq = t - stats[4];  return 4; }
		t += stats[5];
		if (left_freq < t) { left_freq = t - stats[5];  return 5; }
		t += stats[6];
		if (left_freq < t) { left_freq = t - stats[6];  return 6; }
		t += stats[7];
		if (left_freq < t) { left_freq = t - stats[7];  return 7; }
		left_freq = t;
		return 8;
	}
};

template<>
struct IfImpl<8>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		t += stats[3];
		if (left_freq < t)		return 3;
		t += stats[4];
		if (left_freq < t)		return 4;
		t += stats[5];
		if (left_freq < t)		return 5;
		t += stats[6];
		if (left_freq < t)		return 6;
		return 7;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		t += stats[3];
		if (left_freq < t) { left_freq = t - stats[3];  return 3; }
		t += stats[4];
		if (left_freq < t) { left_freq = t - stats[4];  return 4; }
		t += stats[5];
		if (left_freq < t) { left_freq = t - stats[5];  return 5; }
		t += stats[6];
		if (left_freq < t) { left_freq = t - stats[6];  return 6; }
		left_freq = t;
		return 7;
	}
};

template<>
struct IfImpl<7>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		t += stats[3];
		if (left_freq < t)		return 3;
		t += stats[4];
		if (left_freq < t)		return 4;
		t += stats[5];
		if (left_freq < t)		return 5;
		return 6;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		t += stats[3];
		if (left_freq < t) { left_freq = t - stats[3];  return 3; }
		t += stats[4];
		if (left_freq < t) { left_freq = t - stats[4];  return 4; }
		t += stats[5];
		if (left_freq < t) { left_freq = t - stats[5];  return 5; }
		left_freq = t;
		return 6;
	}
};

template<>
struct IfImpl<6>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		t += stats[3];
		if (left_freq < t)		return 3;
		t += stats[4];
		if (left_freq < t)		return 4;
		return 5;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		t += stats[3];
		if (left_freq < t) { left_freq = t - stats[3];  return 3; }
		t += stats[4];
		if (left_freq < t) { left_freq = t - stats[4];  return 4; }
		left_freq = t;
		return 5;
	}
};

template<>
struct IfImpl<5>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		t += stats[3];
		if (left_freq < t)		return 3;
		return 4;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		t += stats[3];
		if (left_freq < t) { left_freq = t - stats[3];  return 3; }
		left_freq = t;
		return 4;
	}
};

template<>
struct IfImpl<4>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		t += stats[2];
		if (left_freq < t)		return 2;
		return 3;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		t += stats[2];
		if (left_freq < t) { left_freq = t - stats[2];  return 2; }
		left_freq = t;
		return 3;
	}
};

template<>
struct IfImpl<3>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t)		return 0;
		t += stats[1];
		if (left_freq < t)		return 1;
		return 2;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		t += stats[1];
		if (left_freq < t) { left_freq = t - stats[1];  return 1; }
		left_freq = t;
		return 2;
	}
};

template<>
struct IfImpl<2>
{
	static uint32_t get_sym(uint32_t left_freq, uint32_t* stats)
	{
		if (left_freq < stats[0])		return 0;
		return 1;
	}

	static uint32_t get_sym_lf(uint32_t& left_freq, uint32_t* stats)
	{
		uint32_t t = stats[0];
		if (left_freq < t) { left_freq = t - stats[0];  return 0; }
		left_freq = t;
		return 1;
	}
};



