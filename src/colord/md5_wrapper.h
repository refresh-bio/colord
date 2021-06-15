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
#include "libs/md5/md5.h"
#include <cinttypes>
#include <vector>

class CMD5
{
	md5_state_t md5_state;
public:
	CMD5()
	{
		md5_init(&md5_state);
	}
	void Update(const uint8_t* data, size_t size)
	{
		md5_append(&md5_state, reinterpret_cast<const md5_byte_t*>(data), size);
	}

	std::vector<uint8_t> Get()
	{
		std::vector<uint8_t> res(16);
		md5_finish(&md5_state, reinterpret_cast<md5_byte_t*>(res.data()));
		return res;
	}
};
