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
#include <cstddef>

#include <algorithm>
#include <vector>
#include <utility>
#include <iterator>

// ************************************************************************
// *** Global const iterator
template<typename CIntHashMapLP>
class const_cint_hash_map_lp_iterator
{
	friend CIntHashMapLP;

public:
	using pair_type = typename CIntHashMapLP::pair_type;
	using key_type = typename CIntHashMapLP::key_type;
	using value_type = typename CIntHashMapLP::value_type;
	using mapped_type = typename CIntHashMapLP::mapped_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
//	using pointer = value_type*;
//	using reference = value_type&;
	using pointer = pair_type*;
	using reference = pair_type&;
	using vector_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;

	const_cint_hash_map_lp_iterator() = default;

	const_cint_hash_map_lp_iterator(CIntHashMapLP* _p_hm)
	{
		p_hm = _p_hm;
		empty_key = _p_hm->empty_key;
		iter = _p_hm->data.end();
		iter_begin = _p_hm->data.begin();
		iter_end = _p_hm->data.end();

		pp.first = empty_key;
	}

	const_cint_hash_map_lp_iterator(CIntHashMapLP* _p_hm, vector_iterator_type _iter)
	{
		p_hm = _p_hm;
		iter = _iter;
		iter_begin = _p_hm->data.begin();
		iter_end = _p_hm->data.end();
		empty_key = _p_hm->empty_key;

		if (p_hm->is_key(*iter))
		{
			pp.first = p_hm->get_key(*iter);
			++iter;
			if (iter != iter_end)
				pp.second = p_hm->get_value(*iter);
		}
		else if(p_hm->is_empty(*iter))
		{
			increment();
		}
		else if (iter == iter_begin)	// && must be value here
		{
			for (auto p = iter_end; p != iter_begin; --p)
				if (p_hm->is_key(*(p - 1)))
				{
					pp.first = p_hm->get_key(*(p - 1));
					break;
				}
			pp.second = p_hm->get_value(*iter);
		}
	}

	const pair_type& operator*() const
	{
		return pp;
	}

	const pair_type* operator->() const
	{
		return &pp;
	}

	const_cint_hash_map_lp_iterator<CIntHashMapLP>& operator++()
	{
		increment();
		return *this;
	}

	const_cint_hash_map_lp_iterator<CIntHashMapLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	bool operator==(const const_cint_hash_map_lp_iterator<CIntHashMapLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	bool operator!=(const const_cint_hash_map_lp_iterator<CIntHashMapLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	CIntHashMapLP* p_hm;
	vector_iterator_type iter;
	vector_iterator_type iter_begin;
	vector_iterator_type iter_end;
	key_type empty_key;
	pair_type pp;

	void increment()
	{
		for (++iter; iter != iter_end; ++iter)
		{
			if (p_hm->is_key(*iter))
			{
				pp.first = p_hm->get_key(*iter);
				++iter;
				if (iter != iter_end)
					pp.second = p_hm->get_value(*iter);
				break;
			}
			else if (p_hm->is_value(*iter))
			{
				pp.second = p_hm->get_value(*iter);
				break;
			}
		}
	}
};

// ************************************************************************
// *** Global iterator
template <typename CIntHashMapLP>
class cint_hash_map_lp_iterator : public const_cint_hash_map_lp_iterator<CIntHashMapLP>
{
	friend CIntHashMapLP;

public:
	using pair_type = typename CIntHashMapLP::pair_type;
	using key_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::key_type;
	using value_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::value_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename CIntHashMapLP::VectorType::iterator;

	cint_hash_map_lp_iterator() = default;
	cint_hash_map_lp_iterator(CIntHashMapLP* _p_hm) : const_cint_hash_map_lp_iterator<CIntHashMapLP>(_p_hm)
	{
	}

	cint_hash_map_lp_iterator(CIntHashMapLP* _p_hm, vector_iterator_type _iter) : const_cint_hash_map_lp_iterator<CIntHashMapLP>(_p_hm, _iter)
	{
	}

	pair_type& operator*()
	{
		return const_cast<pair_type&>(*this->iter);
	}

	pair_type* operator->()
	{
		return const_cast<pair_type*>(&(*this->iter));
	}
};

// ************************************************************************
// *** Local const iterator
template<typename CIntHashMapLP>
class const_cint_hash_map_lp_local_iterator
{
	friend CIntHashMapLP;

public:
	using pair_type = typename CIntHashMapLP::pair_type;
	using key_type = typename CIntHashMapLP::key_type;
	using value_type = typename CIntHashMapLP::value_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;
	using vector_const_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;

	const_cint_hash_map_lp_local_iterator() = default;

	const_cint_hash_map_lp_local_iterator(CIntHashMapLP* _p_hm)
	{
		//		local_key = _p_hm->empty_key;
		//		empty_key = _p_hm->empty_key;
		iter = _p_hm->data.end();
		//		iter_begin = _p_hm->data.begin();
				iter_end = _p_hm->data.end();
		//		compare = _p_hm->compare;
	}

	const_cint_hash_map_lp_local_iterator(CIntHashMapLP* _p_hm, vector_iterator_type _iter)
	{
		iter = _iter;
		iter_begin = _p_hm->data.begin();
		iter_end = _p_hm->data.end();

		if (iter != iter_end)
		{
			pp.first = get_key(*iter);
			++iter;
			if (iter == iter_end)
				iter = iter_begin;
			pp.second = get_value(*iter);
		}
	}

	const pair_type& operator*() const
	{
		return pp;
	}

	const pair_type* operator->() const
	{
		return &pp;
	}

	const_cint_hash_map_lp_local_iterator<CIntHashMapLP>& operator++()
	{
		increment();
		return *this;
	}

	const_cint_hash_map_lp_local_iterator<CIntHashMapLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	bool operator==(const const_cint_hash_map_lp_local_iterator<CIntHashMapLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	bool operator!=(const const_cint_hash_map_lp_local_iterator<CIntHashMapLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	pair_type pp;
	vector_iterator_type iter;
	vector_iterator_type iter_begin;
	vector_iterator_type iter_end;

	const uint32_t value_marker_shift = sizeof(value_type) * 8 - 2;
	const uint32_t key_marker_shift = sizeof(value_type) * 8 - 1;
	const value_type empty_key = (~((value_type)0)) >> 1;
	//		const value_type marker_mask = ((value_type)3) << general_marker_shift;
	const value_type marker_key = ((value_type)1) << key_marker_shift;
	//		const value_type marker_value = (value_type)0;
	const value_type key_mask = ~(((value_type)~0) << key_marker_shift);
	const value_type value_mask = ~(((value_type)~0) << value_marker_shift);

	bool is_key(value_type x) const
	{
		return (x >> key_marker_shift) == 1;
	}

	bool is_empty(value_type x) const
	{
		return x == empty_key;
	}

	bool is_value(value_type x) const
	{
		return (x >> value_marker_shift) == 0;
	}

	value_type get_key(value_type x) const
	{
		return x & key_mask;
		//			return x;
	}

	value_type get_value(value_type x) const
	{
		//			return x & value_mask;
		return x;
	}

	void increment()
	{
		++iter;

		if (iter == iter_end)
			iter = iter_begin;

		if (!is_value(*iter))
			iter = iter_end;
		else
			pp.second = get_value(*iter);

		return;		
	}
};

// ************************************************************************
// *** Local const iterator
template<typename CIntHashMapLP>
class const_cint_hash_map_lp_local_value_iterator
{
	friend CIntHashMapLP;

public:
	using key_type = typename CIntHashMapLP::key_type;
	using value_type = typename CIntHashMapLP::mapped_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;
	using vector_const_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;
	using key_equal = typename CIntHashMapLP::key_equal;

	const_cint_hash_map_lp_local_value_iterator() = default;

	const_cint_hash_map_lp_local_value_iterator(CIntHashMapLP* _p_hm)
	{
		iter = _p_hm->data.end();
		iter_end = _p_hm->data.end();
	}

	const_cint_hash_map_lp_local_value_iterator(CIntHashMapLP* _p_hm, vector_iterator_type _iter)
	{
		iter = _iter;
		iter_begin = _p_hm->data.begin();
		iter_end = _p_hm->data.end();

		if (iter != iter_end)
		{
			++iter;
			if (iter == iter_end)
				iter = iter_begin;
		}
	}

	const value_type& operator*() const
	{
		return *iter;
	}

	const value_type* operator->() const
	{
		return &(*iter);
	}

	const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>& operator++()
	{
		increment();
		return *this;
	}

	const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	bool operator==(const const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	bool operator!=(const const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	vector_iterator_type iter;
	vector_iterator_type iter_begin;
	vector_iterator_type iter_end;

	void increment()
	{
		++iter;

		if (iter == iter_end)
			iter = iter_begin;

		if (!is_value(*iter))
			iter = iter_end;

		return;
	}
};


// ************************************************************************
// *** Local iterator
template <typename CIntHashMapLP>
class cint_hash_map_lp_local_iterator : public const_cint_hash_map_lp_local_iterator<CIntHashMapLP>
{
	friend CIntHashMapLP;

public:
	using pair_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::pair_type;
	using key_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::key_type;
	using value_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::value_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename CIntHashMapLP::VectorType::iterator;
	using vector_const_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;

	cint_hash_map_lp_local_iterator() = default;
	cint_hash_map_lp_local_iterator(CIntHashMapLP* _p_hm) : const_cint_hash_map_lp_local_iterator<CIntHashMapLP>(_p_hm)
	{
	}

	cint_hash_map_lp_local_iterator(CIntHashMapLP* _p_hm, vector_iterator_type _iter) : const_cint_hash_map_lp_local_iterator<CIntHashMapLP>(_p_hm, _iter)
	{
	}

	cint_hash_map_lp_local_iterator(CIntHashMapLP* _p_hm, vector_const_iterator_type _iter) : const_cint_hash_map_lp_local_iterator<CIntHashMapLP>(_p_hm, _iter)
	{
	}

	pair_type& operator*()
	{
//		return const_cast<pair_type&>(*this->iter);
		return const_cast<pair_type&>(this->pp);
	}

	pair_type* operator->()
	{
		return const_cast<pair_type*>(&(this->pp));
	}
};

// ************************************************************************
// *** Local iterator
template <typename CIntHashMapLP>
class cint_hash_map_lp_local_value_iterator : public const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>
{
	friend CIntHashMapLP;

public:
	using key_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::key_type;
	using value_type = typename const_cint_hash_map_lp_iterator<CIntHashMapLP>::mapped_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename CIntHashMapLP::VectorType::iterator;
	using vector_const_iterator_type = typename CIntHashMapLP::VectorType::const_iterator;

	cint_hash_map_lp_local_value_iterator() = default;
	cint_hash_map_lp_local_value_iterator(CIntHashMapLP* _p_hm) : const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>(_p_hm)
	{
	}

	cint_hash_map_lp_local_value_iterator(CIntHashMapLP* _p_hm, vector_iterator_type _iter) : const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>(_p_hm, _iter)
	{
	}

	cint_hash_map_lp_local_value_iterator(CIntHashMapLP* _p_hm, vector_const_iterator_type _iter) : const_cint_hash_map_lp_local_value_iterator<CIntHashMapLP>(_p_hm, _iter)
	{
	}

	value_type& operator*()
	{
		return const_cast<value_type&>(this->iter->second);
	}

	value_type* operator->()
	{
		return const_cast<value_type*>(&(this->iter->second));
	}
};

// ************************************************************************
// *** Hash map with linear probing (multikey)
template<typename Key_Value_t,
	typename Hash_t = std::hash<Key_Value_t>>
	class cint_hash_map_lp {
	public:
		using pair_type = std::pair<Key_Value_t, Key_Value_t>;
		using key_type = Key_Value_t;
		using mapped_type = Key_Value_t;
		using value_type = Key_Value_t;
		using hasher = Hash_t;
		using reference = Key_Value_t&;
		using const_reference = const Key_Value_t&;
		using size_type = size_t;
		using difference_type = ptrdiff_t;
		using cint_hash_map_lp_type = cint_hash_map_lp<Key_Value_t, Hash_t>;
		using iterator = cint_hash_map_lp_iterator<cint_hash_map_lp_type>;
		using const_iterator = const_cint_hash_map_lp_iterator<cint_hash_map_lp_type>;
		using local_iterator = cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>;
		using const_local_iterator = const_cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>;
		using local_value_iterator = cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>;
		using const_local_value_iterator = const_cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>;

	private:
		using VectorType = std::vector<value_type>;

	private:
		Hash_t hash;

		double max_fill_factor;
		size_t no_elements;
		size_t no_elements_unique;
		std::vector<value_type> data;
		size_t allocated;
		size_t size_when_restruct;
		size_t allocated_mask;

		const uint32_t value_marker_shift = sizeof(value_type) * 8 - 2;
		const uint32_t key_marker_shift = sizeof(value_type) * 8 - 1;
		const value_type empty_key = (~((value_type)0)) >> 1;
//		const value_type marker_mask = ((value_type)3) << general_marker_shift;
		const value_type marker_key = ((value_type)1) << key_marker_shift;
//		const value_type marker_value = (value_type)0;
		const value_type key_mask = ~(((value_type)~0) << key_marker_shift);
		const value_type value_mask = ~(((value_type)~0) << value_marker_shift);

		bool is_key(Key_Value_t x) const
		{
			return (x >> key_marker_shift) == 1;
		}

		bool is_empty(Key_Value_t x) const
		{
			return x == empty_key;
		}

		bool is_value(Key_Value_t x) const
		{
			return (x >> value_marker_shift) == 0;
		}

		Key_Value_t get_key(Key_Value_t x) const
		{
			return x & key_mask;
//			return x;
		}

		Key_Value_t get_value(Key_Value_t x) const
		{
//			return x & value_mask;
			return x;
		}

	public:
		friend class cint_hash_map_lp_iterator<cint_hash_map_lp_type>;
		friend class const_cint_hash_map_lp_iterator<cint_hash_map_lp_type>;
		friend class cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>;
		friend class const_cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>;
		friend class cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>;
		friend class const_cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>;

		virtual ~cint_hash_map_lp() = default;

		explicit cint_hash_map_lp(size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Hash_t _hash = Hash_t()) :
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);
		}

		cint_hash_map_lp(const cint_hash_map_lp<Key_Value_t, Hash_t>& src) = default;

		cint_hash_map_lp(cint_hash_map_lp<Key_Value_t, Hash_t>&& src) = default;

		template <typename InputIterator>
		cint_hash_map_lp(InputIterator first, InputIterator last,
			size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Hash_t _hash = Hash_t()) :
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);

			for (auto p = first; p != last; ++p)
				insert_fast(*p);
		}

		explicit cint_hash_map_lp(std::initializer_list<value_type> il,
			size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Hash_t _hash = Hash_t()) :
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);

			for (auto& x : il)
				insert_fast(x);
		}


		cint_hash_map_lp_type& operator=(
			const cint_hash_map_lp<Key_Value_t, Hash_t>& rhs)
		{
			if (this != &rhs)
			{
				data.clear();
				data.shrink_to_fit();

				hash = rhs.hash;
				data = rhs.data;
				no_elements = rhs.no_elements;
				no_elements_unique = rhs.no_elements_unique;
				allocated = rhs.allocated;
				size_when_restruct = rhs.size_when_restruct;
				allocated_mask = rhs.allocated_mask;
			}

			return *this;
		}

		cint_hash_map_lp_type& operator=(
			cint_hash_map_lp<Key_Value_t, Hash_t>&& rhs)
		{
			if (this != &rhs)
			{
				hash = rhs.hash;
				data = move(rhs.data);
				no_elements = rhs.no_elements;
				no_elements_unique = rhs.no_elements_unique;
				allocated = rhs.allocated;
				size_when_restruct = rhs.size_when_restruct;
				allocated_mask = rhs.allocated_mask;
			}

			return *this;
		}

		//	hash_map_type& operator=(std::initializer_list<value_type> il);

		void swap(cint_hash_map_lp& x)
		{
			swap(hash, x.hash);

			swap(max_fill_factor, x.max_fill_factor);
			swap(no_elements, x.no_elements);
			swap(no_elements_unique, x.no_elements_unique);
			swap(data, x.data);
			swap(allocated, x.allocated);
			swap(size_when_restruct, x.size_when_restruct);
			swap(allocated_mask, x.allocated_mask);
		}

		// *************************************
		iterator begin()
		{
			if (!no_elements)
				return end();

			auto p = cint_hash_map_lp_iterator<cint_hash_map_lp_type>(this, data.begin());
/*			if (p->first == empty_key)
				++p;*/

			return p;
		}

		iterator end()
		{
			return cint_hash_map_lp_iterator<cint_hash_map_lp_type>(this);
		}

		const_iterator begin() const
		{
			return const_cast<cint_hash_map_lp_type*>(this)->begin();
		}

		const_iterator end() const
		{
			return const_cast<cint_hash_map_lp_type*>(this)->end();
		}

		const_iterator cbegin() const
		{
			return begin();
		}

		const_iterator cend() const
		{
			return end();
		}

		// *************************************
		local_iterator local_end()
		{
			return cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>(this);
		}

		const_local_iterator local_end() const
		{
			return const_cast<cint_hash_map_lp_type*>(this)->local_end();
		}

		/*		const_local_iterator clocal_end() const
				{
					return local_end();
				}*/

		local_iterator local_begin(iterator x)
		{
			return cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>(this, x.iter);
		}

		local_iterator local_begin(const_iterator x)
		{
			return cint_hash_map_lp_local_iterator<cint_hash_map_lp_type>(this, x.iter);
		}

		// *************************************
		local_value_iterator local_value_end()
		{
			return cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>(this);
		}

		const_local_value_iterator local_value_end() const
		{
			return const_cast<cint_hash_map_lp_type*>(this)->local_value_end();
		}

		/*		const_local_iterator clocal_end() const
				{
					return local_end();
				}*/

		local_value_iterator local_value_begin(iterator x)
		{
			return cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>(this, x.iter);
		}

		local_value_iterator local_value_begin(const_iterator x)
		{
			return cint_hash_map_lp_local_value_iterator<cint_hash_map_lp_type>(this, x.iter);
		}

		// *************************************
		bool empty()
		{
			return no_elements == 0;
		}

		size_type allocated_size()
		{
			return allocated;
		}

		size_type size()
		{
			return no_elements;
		}

		size_type size_unique()
		{
			return no_elements_unique;
		}

		size_type max_size()
		{
			return data.max_size();
		}

		size_t reserve(size_t _requested_reserve)
		{
			if (_requested_reserve < allocated)
				return allocated;

			restruct(_requested_reserve);
		}

		//	T& operator[](const key_type& k);

		template <typename InputIterator>
		void insert(InputIterator first, InputIterator last)
		{
			for (auto p = first; p != last; ++p)
				insert_fast(*p);
		}

		void insert(std::initializer_list<value_type> il)
		{
			for (auto& x : il)
				insert_fast(x);
		}

		std::pair<iterator, bool> insert(const pair_type& x)
		{
			if (no_elements + no_elements_unique >= size_when_restruct)
				restruct(allocated * 2);

			size_t h = hash(x.first) & allocated_mask;
			bool key_observed = false;

			while (!is_empty(data[h]))
			{
				if (is_key(data[h]) && get_key(data[h]) == x.first)
				{
					key_observed = true;
					break;
				}

				h = (h + 1) & allocated_mask;
			}

			if (!key_observed)
			{
				assert(is_empty(data[h]));

				data[h] = marker_key + x.first;
				h = (h + 1) & allocated_mask;

				if (!is_empty(data[h]))
					push_forward(h);
//				data[h] = marker_value + x.second;
				data[h] = x.second;

				++no_elements_unique;
			}
			else
			{
				assert(is_key(data[h]));

				do {
					h = (h + 1) & allocated_mask;
				} while (is_value(data[h]));

				if (!is_empty(data[h]))
				{
					assert(is_key(data[h]));
					push_forward(h);
				}
//				data[h] = marker_value + x.second;
				data[h] = x.second;
			}

			++no_elements;

			return std::make_pair(iterator(this, data.begin() + h), true);
		}
		
		bool insert_fast(const pair_type& x)
		{
			if (no_elements + no_elements_unique >= size_when_restruct)
				restruct(allocated * 2);

			size_t h = hash(x.first) & allocated_mask;
			bool key_observed = false;

			while (!is_empty(data[h]))
			{
				if (is_key(data[h]) && get_key(data[h]) == x.first)
				{
					key_observed = true;
					break;
				}

				h = (h + 1) & allocated_mask;
			}

			if (!key_observed)
			{
				data[h] = marker_key + x.first;
				h = (h + 1) & allocated_mask;

				if (!is_empty(data[h]))
					push_forward(h);
//				data[h] = marker_value + x.second;
				data[h] = x.second;

				++no_elements_unique;
			}
			else
			{
				do {
					h = (h + 1) & allocated_mask;
				} while (is_value(data[h]));

				if (!is_empty(data[h]))
					push_forward(h);
//				data[h] = marker_value + x.second;
				data[h] = x.second;
			}

			++no_elements;

			return true;
		}
		
		local_iterator find(const key_type& key)
		{
			// !!! For unknown reasons turning off inline here speeds things up significantly, but only in "find", it does not affect "check"
			size_t pos = _find_noinline(key, hash(key) & allocated_mask);

			if (!is_empty(data[pos]))
				return local_iterator(this, data.begin() + pos);
			else
				return local_iterator(this);
		}

		bool check(const key_type& key)
		{
			return _check_noinline(key, hash(key) & allocated_mask);

/*			size_t pos = _find(key, hash(key) & allocated_mask);

			return compare(data[pos].first, key);*/
		}

		size_type count(const key_type& key) const
		{
			size_t pos = _find(key, hash(key) & allocated_mask);

			if (!is_key(data[pos]) || get_key(data[pos]) != key)
				return 0;

			size_t r = 0;

			for (pos = (pos + 1) & allocated_mask; is_value(data[pos]); pos = (pos + 1) & allocated_mask)
				++r;

			return r;
		}

		hasher hash_function() const
		{
			return hash;
		}

		void prefetch(const Key_Value_t& key)
		{
			size_t h = hash(key) & allocated_mask;

#ifdef _WIN32
			_mm_prefetch((const char*)(data.data() + h), _MM_HINT_T0);
#else
			__builtin_prefetch(data.data() + h);
#endif
		}

	private:
		void construct(size_t _init_reserved, double _max_fill_factor)
		{
			max_fill_factor = _max_fill_factor;

			if (max_fill_factor > 0.99)
				max_fill_factor = 0.99;
			else if (max_fill_factor < 0.01)
				max_fill_factor = 0.1;

			_reserve(_init_reserved);

			size_when_restruct = (size_t)(allocated * max_fill_factor);
			if (size_when_restruct == allocated)
				--size_when_restruct;
		}

		void _reserve(size_t requested_allocated)
		{
			allocated = requested_allocated;

			if (allocated < 8)
				allocated = 8;

			// Round up to the power of 2
			if ((allocated & (allocated - 1)))
			{
				while ((allocated & (allocated - 1)))
					allocated &= allocated - 1;
				allocated *= 2;
			}

			allocated_mask = allocated - 1ull;
			size_when_restruct = (size_t)(allocated * max_fill_factor);
			if (size_when_restruct == allocated)
				--size_when_restruct;

			data.resize(allocated, empty_key);

			no_elements = 0;
			no_elements_unique = 0;
		}

		void restruct(size_t new_allocated)
		{
			std::vector<value_type> old_data;
			old_data = move(data);
			size_t old_allocated = allocated;

			_reserve(new_allocated);

			no_elements = 0;
			no_elements_unique = 0;

			size_t pos;

			for (pos = 0; is_value(old_data[pos]) && pos < old_allocated; ++pos)
				;

			Key_Value_t key{};

			for (size_t i = pos; i < old_allocated; ++i)
			{
				if (is_key(old_data[i]))
					key = get_key(old_data[i]);
				else if (is_value(old_data[i]))
					insert_fast(std::make_pair(key, get_value(old_data[i])));
			}

			for (size_t i = 0; i < pos; ++i)
			{
				if (is_key(old_data[i]))
					key = get_key(old_data[i]);
				else if (is_value(old_data[i]))
					insert_fast(std::make_pair(key, get_value(old_data[i])));
			}
		}

		void push_forward(size_t h)
		{
			Key_Value_t x = data[h];
			h = (h + 1) & allocated_mask;

			while (!is_empty(data[h]))
			{
				std::swap(x, data[h]);
				h = (h + 1) & allocated_mask;
			}

			data[h] = x;
		}


		size_t _find(const Key_Value_t key, size_t pos)
		{
			if (is_empty(data[pos]) || (is_key(data[pos]) && get_key(data[pos]) == key))
				return pos;

			for (++pos; pos < allocated; ++pos)
				if(is_empty(data[pos]) || (is_key(data[pos]) && get_key(data[pos]) == key))
					return pos;

			for (pos = 0; pos < allocated; ++pos)
				if (is_empty(data[pos]) || (is_key(data[pos]) && get_key(data[pos]) == key))
					return pos;
					
			return pos;
		}

		NO_INLINE size_t _find_noinline(const Key_Value_t key, size_t pos)
		{
			if (is_empty(data[pos]) || (is_key(data[pos]) && get_key(data[pos]) == key))
				return pos;

			for (++pos; pos < allocated; ++pos)
				if (is_empty(data[pos]) || (is_key(data[pos]) && get_key(data[pos]) == key))
					return pos;

			for (pos = 0; pos < allocated; ++pos)
				if (is_empty(data[pos]) || (is_key(data[pos]) && get_key(data[pos]) == key))
					return pos;
					
			return pos;
		}
		
		NO_INLINE bool _check_noinline(const Key_Value_t key, size_t pos)
		{
			if (is_empty(data[pos]))
				return false;
			if (is_key(data[pos]) && get_key(data[pos]) == key)
				return true;

			for (++pos; pos < allocated; ++pos)
			{
				if (is_empty(data[pos]))
					return false;
				if (is_key(data[pos]) && get_key(data[pos]) == key)
					return true;
			}

			for (pos = 0; pos < allocated; ++pos)
			{
				if (is_empty(data[pos]))
					return false;
				if (is_key(data[pos]) && get_key(data[pos]) == key)
				return true;
			}
			
			return false;
		}

		bool _check(const Key_Value_t key, size_t pos)
		{
			if (is_empty(data[pos]))
				return false;
			if (is_key(data[pos]) && get_key(data[pos]) == key)
			return true;

			for (++pos; pos < allocated; ++pos)
			{
				if (is_empty(data[pos]))
					return false;
				if (is_key(data[pos]) && get_key(data[pos]) == key)
				return true;
			}

			for (pos = 0; pos < allocated; ++pos)
			{
				if (is_empty(data[pos]))
					return false;
				if (is_key(data[pos]) && get_key(data[pos]) == key)
				return true;
			}

			return false;
		}
};

// EOF
