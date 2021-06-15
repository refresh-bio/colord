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

#include <mmintrin.h>
#include <cstdint>
#include <xmmintrin.h>
#include <iostream> 
#include <cstddef>

#include "defs.h"
#include "rc.h"

//#define FIND_EXT_SUPPORT

namespace entropy_coder 
{
template<typename MODEL> class CContextHM {
public:
	typedef struct {
		context_t ctx;
		MODEL* rcm;
#ifdef FIND_EXT_SUPPORT
		size_t counter;
#endif
	} item_t;

	typedef context_t key_type;
	typedef MODEL* value_type;
	typedef MODEL model_type;
	typedef size_t aux_type;

private:
	double max_fill_factor;

	size_t size;
	item_t* data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;

	size_t ht_memory;
	size_t ht_total;
	size_t ht_match;

	void restruct(void)
	{
		item_t* old_data = data;
		size_t old_allocated = allocated;

		allocated *= 2;
		size = 0;

		allocated_mask = allocated - 1ull;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		ht_memory += allocated * sizeof(item_t);

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].rcm != nullptr)
#ifdef FIND_EXT_SUPPORT
				insert(old_data[i].ctx, old_data[i].rcm, old_data[i].counter);
#else
				insert(old_data[i].ctx, old_data[i].rcm);
#endif

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}

	// Based on murmur64
	size_t hash(context_t ctx)
	{
		auto h = ctx;

		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;

		return h & allocated_mask;
	}

public:
	CContextHM()
	{
		ht_memory = 0;
		ht_total = 0;
		ht_match = 0;

		allocated = 1u << 20;
		allocated_mask = allocated - 1;

		size = 0;
		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		max_fill_factor = 0.4;

		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)(allocated * max_fill_factor);
	}

	~CContextHM()
	{
		if (data == nullptr)
			return;

		for (size_t i = 0; i < allocated; ++i)
			if (data[i].rcm)
				delete data[i].rcm;
		delete[] data;
	}

	size_t get_bytes() const {
		return ht_memory;
	}

#ifdef FIND_EXT_SUPPORT
	void debug_list(std::vector<CContextHM<MODEL>::item_t>& v_ctx)
	{
		v_ctx.clear();

		for (size_t i = 0; i < allocated; ++i)
			if (data[i].rcm)
				v_ctx.push_back(data[i]);

		sort(v_ctx.begin(), v_ctx.end(), [](auto& x, auto& y) {return x.counter > y.counter; });
	}
#endif

#ifdef FIND_EXT_SUPPORT
	bool insert(const context_t ctx, MODEL* rcm, size_t counter = 0)
#else
	bool insert(const context_t ctx, MODEL* rcm)
#endif
	{
		if (size >= size_when_restruct)
			restruct();

		size_t h = hash(ctx);

		if (data[h].rcm != nullptr)
		{
			do
			{
				h = (h + 1) & allocated_mask;
			} while (data[h].rcm != nullptr);
		}

		++size;

		data[h].ctx = ctx;
		data[h].rcm = rcm;
#ifdef FIND_EXT_SUPPORT
		data[h].counter = counter;
#endif

		return true;
	}

	MODEL* find(const context_t ctx)
	{
		size_t h = hash(ctx);

		if (data[h].rcm == nullptr)
			return nullptr;

		if (data[h].ctx == ctx)
			return data[h].rcm;

		h = (h + 1) & allocated_mask;

		while (data[h].rcm != nullptr)
		{
			if (data[h].ctx == ctx)
				return data[h].rcm;
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}

#ifdef FIND_EXT_SUPPORT
	MODEL* find_ext(const context_t ctx, size_t*& p_counter)
	{
		size_t h = hash(ctx);

		if (data[h].rcm == nullptr)
			return nullptr;

		if (data[h].ctx == ctx)
		{
			p_counter = &data[h].counter;
			return data[h].rcm;
		}

		h = (h + 1) & allocated_mask;

		while (data[h].rcm != nullptr)
		{
			if (data[h].ctx == ctx)
			{
				p_counter = &data[h].counter;
				return data[h].rcm;
			}
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}
#endif

	void prefetch(const context_t ctx)
	{
		size_t h = hash(ctx);

#ifdef _WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}

	size_t get_size(void) const
	{
		return size;
	}
};

template<typename MODEL> class CContextVec {
public:
	typedef MODEL* item_t;

	typedef context_t key_type;
	typedef MODEL* value_type;
	typedef size_t aux_type;

private:
	size_t size;
	item_t* data;

	void allocate()
	{
		if (size)
			data = new item_t[size];
		else
			data = nullptr;

		for (uint64_t i = 0; i < size; ++i)
			data[i] = nullptr;
	}

	void deallocate()
	{
		if (data == nullptr)
			return;

		for (size_t i = 0; i < size; ++i)
			if (data[i])
				delete data[i];
		delete[] data;
	}

	void reallocate(uint64_t _size)
	{
		auto old_data = data;
		auto old_size = size;

		size = _size;
		data = new item_t[size];
		for (uint64_t i = 0; i < old_size; ++i)
			data[i] = old_data[i];

		for (uint64_t i = old_size; i < size; ++i)
			data[i] = nullptr;

		delete[] old_data;
	}

public:
	CContextVec(uint64_t _size = 0) : size(_size)
	{
		allocate();
	}

	~CContextVec()
	{
		deallocate();
	}

/*	void resize(uint64_t _size)
	{
		deallocate();
		size = _size;
		allocate();
	}*/

	bool insert(const context_t ctx, MODEL* rcm)
	{
		if (ctx >= size)
			reallocate((uint64_t) (ctx * 1.1) + 1);

		data[ctx] = rcm;

		return true;
	}

	MODEL* find(const context_t ctx)
	{
		if (ctx >= size)
			return nullptr;

		return data[ctx];
	}

	void prefetch(const context_t ctx)
	{
#ifdef _WIN32
		_mm_prefetch((const char*)(data + ctx), _MM_HINT_T0);
#else
		__builtin_prefetch(data + ctx);
#endif
	}

	size_t get_size(void) const
	{
		return size;
	}
};

template<typename MODEL> class CContextVecEmb {
public:
	typedef MODEL item_t;

	typedef context_t key_type;
	typedef MODEL* value_type;
	typedef size_t aux_type;

private:
	size_t size;
	vector<item_t> data;

	void reallocate(uint64_t _size, MODEL &item)
	{
		size = _size;
		data.resize(size, item);
	}

public:
	CContextVecEmb(uint64_t _size = 0, MODEL *item = nullptr) : size(_size)
	{
		if(size)
			reallocate(size, *item);
	}

	~CContextVecEmb()
	{
	}

	MODEL* insert(const context_t ctx, MODEL &rcm)
	{
		if (ctx >= size)
			reallocate((uint64_t)(ctx * 1.1) + 1, rcm);

		return data.data() + ctx;
	}

	MODEL* find(const context_t ctx)
	{
		if (ctx >= size)
			return nullptr;

		return data.data() + ctx;
	}

	void prefetch(const context_t ctx)
	{
		#ifdef _WIN32
				_mm_prefetch((const char*)(data.data() + ctx), _MM_HINT_T0);
		#else
				__builtin_prefetch(data.data() + ctx);
		#endif
	}

	size_t get_size(void) const
	{
		return size;
	}
};


} // namespace entropy_coder
// EOF
