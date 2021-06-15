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
