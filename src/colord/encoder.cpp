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
#include "encoder.h"
#include "parallel_queue.h"
#include "queues_data.h"
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <set>
#include <unordered_map>

class CMmersSTL
{
	using hash_map_type = std::unordered_map<anchor_type, std::vector<uint32_t>>;
	hash_map_type mmers;

public:
	//unique
	uint32_t GetNMmers() { return static_cast<uint32_t>(mmers.size()); }

	explicit CMmersSTL(const read_t& read, uint32_t m)
	{
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;
			//auto [it, new_inserted] = mmers.emplace(mmer, std::vector<uint32_t>{});
			auto it = mmers.emplace(mmer, std::vector<uint32_t>{}).first;
			it->second.push_back(pos - m + 1);
		}
	}
	explicit CMmersSTL(const read_t& read, uint32_t m, const CMmersSTL& allowOnly)
	{
		const hash_map_type& include = allowOnly.mmers;
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;
			if (include.find(mmer) != include.end())
			{
				//auto [it, new_inserted] = mmers.emplace(mmer, std::vector<uint32_t>{});
				auto it = mmers.emplace(mmer, std::vector<uint32_t>{}).first;
				it->second.push_back(pos - m + 1);
			}
		}
	}
	void GetIntersection(const CMmersSTL& inParam,
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& outThis,
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& outParam
	)
	{
		outThis.clear();
		outParam.clear();
		const hash_map_type& in_mmers_this = this->mmers;
		const hash_map_type& in_mmers_param = inParam.mmers;

		const hash_map_type* _smaller_input = &in_mmers_this;
		const hash_map_type* _bigger_input = &in_mmers_param;

		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>* _smaller_output = &outThis;
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>* _bigger_output = &outParam;
		if (in_mmers_this.size() > in_mmers_param.size())
		{
			std::swap(_smaller_input, _bigger_input);
			std::swap(_smaller_output, _bigger_output);
		}
		const hash_map_type& smaller_input = *_smaller_input;
		const hash_map_type& bigger_input = *_bigger_input;
		std::vector<std::pair<uint64_t, std::vector<uint32_t>>>& smaller_output = *_smaller_output;
		std::vector<std::pair<uint64_t, std::vector<uint32_t>>>& bigger_output = *_bigger_output;

		for (const auto& [key, val] : smaller_input)
		{
			if (auto it = bigger_input.find(key); it != bigger_input.end())
			{
				smaller_output.emplace_back(key, val);
				std::vector<uint32_t> pos = it->second;
				bigger_output.emplace_back(key, pos);
			}
		}
	}
};

class CMmersHashMapLP
{
	using hash_map_type = hash_map_lp<uint64_t, uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;
	using hash_set_type = hash_set_lp<uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;	
	hash_map_type mmers;

public:
	uint32_t GetNMmers() { return static_cast<uint32_t>(mmers.size_unique()); }

	explicit CMmersHashMapLP(const read_t& read, uint32_t m) :
		mmers(~0ull, static_cast<size_t>((read_len(read) - m + 1) / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{})
	{
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;
//			mmers.insert(std::make_pair(mmer, pos - m + 1));
			mmers.insert_fast(std::make_pair(mmer, pos - m + 1));
		}
	}

	explicit CMmersHashMapLP(const read_t& read, uint32_t m, CMmersHashMapLP& allowOnly, CBloomFilter &bloom_mmers) :
		mmers(~0ull, 16, 0.8, std::equal_to<uint64_t>{}, MurMur64Hash{})
	{
		hash_map_type& include = allowOnly.mmers;
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		std::vector<std::pair<uint64_t, uint64_t>> v_mmers;
		v_mmers.reserve(read_len(read) - m + 1);

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;

			if (bloom_mmers.test(mmer))
				if (include.check(mmer))
					v_mmers.emplace_back(mmer, pos - m + 1);
		}

		mmers = hash_map_type(~0ull, static_cast<uint64_t>(v_mmers.size() / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});
		for (auto& x : v_mmers)
			mmers.insert_fast(x);
	}

	void GetIntersection(CMmersHashMapLP& inParam,
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& outThis,
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& outParam
	)
	{
		outThis.clear();
		outParam.clear();
		hash_map_type& in_mmers_this = this->mmers;
		hash_map_type& in_mmers_param = inParam.mmers;

		hash_map_type* _smaller_input = &in_mmers_this;
		hash_map_type* _bigger_input = &in_mmers_param;

		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>* _smaller_output = &outThis;
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>* _bigger_output = &outParam;
		if (in_mmers_this.size() > in_mmers_param.size())
		{
			std::swap(_smaller_input, _bigger_input);
			std::swap(_smaller_output, _bigger_output);
		}
		hash_map_type& smaller_input = *_smaller_input;
		hash_map_type& bigger_input = *_bigger_input;
		std::vector<std::pair<uint64_t, std::vector<uint32_t>>>& smaller_output = *_smaller_output;
		std::vector<std::pair<uint64_t, std::vector<uint32_t>>>& bigger_output = *_bigger_output;

		smaller_output.reserve(smaller_input.size());
		bigger_output.reserve(smaller_input.size());

		if (smaller_input.size() == smaller_input.size_unique()) //only unique values, it may be done simpler
		{
			for (auto outerIt = smaller_input.begin(); outerIt != smaller_input.end(); ++outerIt)
			{
				auto key = outerIt->first;
				if (auto it = bigger_input.find(key); it != bigger_input.local_end())
				{
					//auto smaller_it = smaller_input.find(key); //must be successfull, but it seems to be unnecessary, maybe there is a better way
/*					auto smaller_it = smaller_input.local_begin(outerIt); //it is guranted that local_begin work well in a case of unique records
					smaller_output.emplace_back(key, std::vector<uint32_t>{});
					for (; smaller_it != smaller_input.local_end(); ++smaller_it)
						smaller_output.back().second.push_back(static_cast<uint32_t>(smaller_it->second));*/
					auto smaller_it = smaller_input.local_value_begin(outerIt); //it is guranted that local_begin work well in a case of unique records
					smaller_output.emplace_back(key, std::vector<uint32_t>(smaller_it, smaller_input.local_value_end()));

					bigger_output.emplace_back(key, std::vector<uint32_t>{});
					for (; it != bigger_input.local_end(); ++it)
						bigger_output.back().second.push_back(static_cast<uint32_t>(it->second));
				}
			}
		}
		else
		{
			//hash_map_type done(~0ull, (smaller_input.size()) / 0.4, 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});

			hash_set_type done(~0ull, static_cast<size_t>(smaller_input.size() / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});
			for (auto outerIt = smaller_input.begin(); outerIt != smaller_input.end(); ++outerIt)
			{
				auto key = outerIt->first;
				if (done.check(key))
					continue;
				else
//					done.insert(key);
					done.insert_fast(key);
				if (auto it = bigger_input.find(key); it != bigger_input.local_end())
				{
					auto smaller_it = smaller_input.find(key); //in case of non unique values local_begin may produce erroneous results

					smaller_output.emplace_back(key, std::vector<uint32_t>{});
					for (; smaller_it != smaller_input.local_end(); ++smaller_it)
						smaller_output.back().second.push_back(static_cast<uint32_t>(smaller_it->second));

					bigger_output.emplace_back(key, std::vector<uint32_t>{});
					for (; it != bigger_input.local_end(); ++it)
						bigger_output.back().second.push_back(static_cast<uint32_t>(it->second));
				}
			}
		}
	}
};

//the idea of this class is to keep m-mers in hash table where values are positions
//but in some cases there is a lot of duplicated m-mers which are all stored as separate entries in HT
//for this reason in this class there is a limit of a numebr of stored duplicates
//if there is more duplicates in the input data the remaining positions are stored in a separate data structure (just a vector)
//this is possible with insert_up_to_n_duplicates method added to hash table which has some limitation (for example there may be no HT restruct to keep the results correct)
class CMmersHashMapDuplicateOptimizedLP
{
	using hash_map_type = hash_map_lp<uint64_t, uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;
	using hash_set_type = hash_set_lp<uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;
	hash_map_type mmers;

	const size_t insert_up_to = 20;

	std::vector<std::vector<uint64_t>> vec_for_highly_duplicated_elems;

	void insert_elem(anchor_type key, uint64_t val) {
		uint64_t x = vec_for_highly_duplicated_elems.size();
		auto y = mmers.insert_up_to_n_duplicates(key, val, x, insert_up_to);

		if (y == insert_up_to) {
			if (x == vec_for_highly_duplicated_elems.size())
				vec_for_highly_duplicated_elems.emplace_back();
			vec_for_highly_duplicated_elems[x].push_back(val);
		}
	}

public:
	uint32_t GetNMmers() {
		return static_cast<uint32_t>(mmers.size_unique());
	}

	size_t Size() {
		size_t x = 0;
		for (auto p : vec_for_highly_duplicated_elems)
			x += p.size();

		// - vec_for_highly_duplicated_elems.size() because this is the numebr of entries in HT that are just to store indexes of vectors
		return mmers.size() - vec_for_highly_duplicated_elems.size() + x;
	}

	explicit CMmersHashMapDuplicateOptimizedLP(const read_t& read, uint32_t m) :
		//+ 1 in initial size is crucial because we use insert_up_to_n_duplicates which does not allow restruct of HT
		mmers(~0ull, static_cast<size_t>((read_len(read) - m + 1) / 0.4) + 1, 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{})
	{
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;
			insert_elem(mmer, pos - m + 1);
		}
	}

	explicit CMmersHashMapDuplicateOptimizedLP(const read_t& read, uint32_t m, CMmersHashMapDuplicateOptimizedLP& allowOnly, CBloomFilter& bloom_mmers) :
		mmers(~0ull, 16, 0.8, std::equal_to<uint64_t>{}, MurMur64Hash{})
	{
		hash_map_type& include = allowOnly.mmers;
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		std::vector<std::pair<uint64_t, uint64_t>> v_mmers;
		v_mmers.reserve(read_len(read) - m + 1);

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;

			if (bloom_mmers.test(mmer))
				if (include.check(mmer))
					v_mmers.emplace_back(mmer, pos - m + 1);
		}
		//+ 1 in initial size is crucial because we use insert_up_to_n_duplicates which does not allow restruct of HT
		mmers = hash_map_type(~0ull, static_cast<uint64_t>(v_mmers.size() / 0.4) + 1, 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});
		for (auto& x : v_mmers)
			insert_elem(x.first, x.second);
	}

	void GetIntersection(CMmersHashMapDuplicateOptimizedLP& inParam,
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& outThis,
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& outParam
	)
	{
		outThis.clear();
		outParam.clear();
		hash_map_type& in_mmers_this = this->mmers;
		hash_map_type& in_mmers_param = inParam.mmers;

		std::vector<std::vector<uint64_t>>& vec_for_highly_duplicated_elems_this = this->vec_for_highly_duplicated_elems;
		std::vector<std::vector<uint64_t>>& vec_for_highly_duplicated_elems_param = inParam.vec_for_highly_duplicated_elems;

		hash_map_type* _smaller_input = &in_mmers_this;
		hash_map_type* _bigger_input = &in_mmers_param;

		std::vector<std::vector<uint64_t>>* _vec_for_highly_duplicated_elems_smaller = &vec_for_highly_duplicated_elems_this;
		std::vector<std::vector<uint64_t>>* _vec_for_highly_duplicated_elems_bigger = &vec_for_highly_duplicated_elems_param;

		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>* _smaller_output = &outThis;
		std::vector<std::pair<anchor_type, std::vector<uint32_t>>>* _bigger_output = &outParam;
		if (in_mmers_this.size() > in_mmers_param.size())
		{
			std::swap(_smaller_input, _bigger_input);
			std::swap(_smaller_output, _bigger_output);
			std::swap(_vec_for_highly_duplicated_elems_smaller, _vec_for_highly_duplicated_elems_bigger);
		}
		hash_map_type& smaller_input = *_smaller_input;
		hash_map_type& bigger_input = *_bigger_input;
		std::vector<std::vector<uint64_t>>& vec_for_highly_duplicated_elems_smaller = *_vec_for_highly_duplicated_elems_smaller;
		std::vector<std::vector<uint64_t>>& vec_for_highly_duplicated_elems_bigger = *_vec_for_highly_duplicated_elems_bigger;
		std::vector<std::pair<uint64_t, std::vector<uint32_t>>>& smaller_output = *_smaller_output;
		std::vector<std::pair<uint64_t, std::vector<uint32_t>>>& bigger_output = *_bigger_output;

		smaller_output.reserve(smaller_input.size());
		bigger_output.reserve(smaller_input.size());

		if (smaller_input.size() == smaller_input.size_unique()) //only unique values, it may be done simpler
		{
			for (auto outerIt = smaller_input.begin(); outerIt != smaller_input.end(); ++outerIt)
			{
				auto key = outerIt->first;
				if (auto it = bigger_input.find(key); it != bigger_input.local_end())
				{
					//auto smaller_it = smaller_input.find(key); //must be successfull, but it seems to be unnecessary, maybe there is a better way
/*					auto smaller_it = smaller_input.local_begin(outerIt); //it is guranted that local_begin work well in a case of unique records
					smaller_output.emplace_back(key, std::vector<uint32_t>{});
					for (; smaller_it != smaller_input.local_end(); ++smaller_it)
						smaller_output.back().second.push_back(static_cast<uint32_t>(smaller_it->second));*/
					auto smaller_it = smaller_input.local_value_begin(outerIt); //it is guranted that local_begin work well in a case of unique records
					smaller_output.emplace_back(key, std::vector<uint32_t>(smaller_it, smaller_input.local_value_end()));

					bigger_output.emplace_back(key, std::vector<uint32_t>{});
					for (size_t i = 0; it != bigger_input.local_end(); ++it, ++i) {
						assert(i <= insert_up_to - 1);
						if (i == insert_up_to - 1)
							for (auto x : vec_for_highly_duplicated_elems_bigger[it->second])
								bigger_output.back().second.push_back(static_cast<uint32_t>(x));
						else
							bigger_output.back().second.push_back(static_cast<uint32_t>(it->second));
					}
				}
			}
		}
		else
		{
			hash_set_type done(~0ull, static_cast<size_t>(smaller_input.size() / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});
			for (auto outerIt = smaller_input.begin(); outerIt != smaller_input.end(); ++outerIt)
			{
				auto key = outerIt->first;
				if (done.check(key))
					continue;
				else
					done.insert_fast(key);
				if (auto it = bigger_input.find(key); it != bigger_input.local_end())
				{
					auto smaller_it = smaller_input.find(key); //in case of non unique values local_begin may produce erroneous results

					smaller_output.emplace_back(key, std::vector<uint32_t>{});
					for (size_t i = 0; smaller_it != smaller_input.local_end(); ++smaller_it, ++i) {
						assert(i <= insert_up_to - 1);
						if (i == insert_up_to - 1)
							for (auto x : vec_for_highly_duplicated_elems_smaller[smaller_it->second])
								smaller_output.back().second.push_back(static_cast<uint32_t>(x));
						else
							smaller_output.back().second.push_back(static_cast<uint32_t>(smaller_it->second));
					}

					bigger_output.emplace_back(key, std::vector<uint32_t>{});
					for (size_t i = 0; it != bigger_input.local_end(); ++it, ++i) {
						assert(i <= insert_up_to - 1);
						if (i == insert_up_to - 1)
							for (auto x : vec_for_highly_duplicated_elems_bigger[it->second])
								bigger_output.back().second.push_back(static_cast<uint32_t>(x));
						else
							bigger_output.back().second.push_back(static_cast<uint32_t>(it->second));
					}
				}
			}
		}
	}
};

class CKmersHashSetLP
{
	using hash_set_type = hash_set_lp<uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;
	hash_set_type kmers;
public:
	CKmersHashSetLP(const std::vector<kmer_type>& kmers) :
		kmers(~0ull, static_cast<uint64_t>(kmers.size() / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{}) //TODO: 0.8?
	{
		for (const kmer_type kmer : kmers)
			this->kmers.insert(kmer);
	}

	size_t Size()
	{
		return kmers.size();
	}
	bool Check(kmer_type kmer)
	{
		return kmers.check(kmer);
	}

	const hash_set_type& GetKmers() const { return kmers; }
};

class CKmersHashMapLP
{
	using hash_map_type = hash_map_lp<uint64_t, uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;
	using hash_set_type = hash_set_lp<uint64_t, std::equal_to<uint64_t>, MurMur64Hash>;
	hash_map_type kmers;

public:
	uint32_t GetNMmers() { return static_cast<uint32_t>(kmers.size_unique()); }

	explicit CKmersHashMapLP(const read_t& read, uint32_t m) :
		kmers(~0ull, static_cast<size_t>((read_len(read) - m + 1) / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{})
	{
		if (read_len(read) < m)
			return;

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
		}

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;
			//			mmers.insert(std::make_pair(mmer, pos - m + 1));
			kmers.insert_fast(std::make_pair(mmer, pos - m + 1));
		}
	}

	explicit CKmersHashMapLP(const read_t& read, uint32_t m, CKmersHashSetLP& includeCanonical, CBloomFilter &bloomKmers, uint32_t modulo) :
		kmers(~0ull, static_cast<size_t>(includeCanonical.Size() / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{})
	{
		if (read_len(read) < m)
			return;

		MurMur64Hash mm;
		libdivide::divider<uint64_t> div(modulo);

		anchor_type mask = (1ull << (2 * m)) - 1;

		anchor_type mmer{}, rev{};
		uint32_t pos = 0;
		for (; pos < m - 1; ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];

			rev >>= 2;
			rev += ((uint64_t)(3 - read[pos])) << (2 * (m - 1));
		}

		for (; pos < read_len(read); ++pos)
		{
			assert(read[pos] < 4); //only ACGT
			mmer <<= 2;
			mmer += read[pos];
			mmer &= mask;

			rev >>= 2;
			rev += ((uint64_t)(3 - read[pos])) << (2 * (m - 1));
			auto can = mmer < rev ? mmer : rev;

			uint64_t h = mm(can);

			if(h - modulo * (h / div) == 0 && bloomKmers.test(can))
				if (includeCanonical.Check(can))
				{
					kmers.insert_fast(std::make_pair(mmer, pos - m + 1));
				}
		}
	}

	//-1 if not exists or is not unique, value in _map in the other case
	int64_t ExistsAndIsUnique(kmer_type kmer)
	{
		auto local_it = kmers.find(kmer);
		auto local_end = kmers.local_end();
		if (local_it == local_end)
			return -1;
		int64_t val = local_it->second;
		++local_it;
		if (local_it != local_end)
			return -1;
		return val;
	}
};



std::vector<std::tuple<anchor_type, uint32_t, uint32_t>> CEncoder::get_aligned_mmers_LIS(
	const std::vector<std::pair<uint32_t, anchor_type>>& sorted_mmers_enc_vec,
	const std::vector<std::pair<uint32_t, anchor_type>>& sorted_mmers_ref_vec,
	const std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& mmers_ref_read_vec)
{
	std::vector<int> lis_vec;
	lis_vec.reserve(sorted_mmers_enc_vec.size() * 2);		// rough estimation

	//for (auto [p, c] : sorted_mmers_enc_vec)
	for (const auto& elem : sorted_mmers_enc_vec)
	{
		auto c = elem.second;
		auto it = std::lower_bound(std::begin(mmers_ref_read_vec), std::end(mmers_ref_read_vec),
			c, [](const auto& e1, const auto& e2) { return e1.first < e2; }); //binary search, if it turn out to be not efficient enough consider using LUT
		assert(it != mmers_ref_read_vec.end());  //must be found, because sorted_mmers_enc_vec and sorted_mmers_ref_vec contains common m-mers
		assert(it->first == c);
		auto& positions = it->second;

		if (positions.size() == 1)
			lis_vec.emplace_back(positions.front());
		else
			for (auto vit = positions.rbegin(); vit != positions.rend(); ++vit)
				lis_vec.emplace_back(*vit);
	}
	std::vector<int> lis_res;
	LIS(lis_vec, lis_res);

	std::vector<std::tuple<anchor_type, uint32_t, uint32_t>> aligned_mmers_lis;
	aligned_mmers_lis.reserve(lis_res.size());

	uint32_t enc_pos{}, ref_pos{};
	for (uint32_t i = 0; i < lis_res.size(); ++i)
	{
		uint32_t pos_in_ref_read = lis_res[i];
		while (sorted_mmers_ref_vec[ref_pos++].first != pos_in_ref_read)
			;

		auto mmer = sorted_mmers_ref_vec[ref_pos - 1].second;
		while (sorted_mmers_enc_vec[enc_pos++].second != mmer)
			;
		aligned_mmers_lis.emplace_back(mmer, sorted_mmers_enc_vec[enc_pos - 1].first, sorted_mmers_ref_vec[ref_pos - 1].first);
	}

	return aligned_mmers_lis;
}

void CEncoder::AddPlainRead(const read_t& read)
{
	size_t read_l = read_len(read);
	stats.LogPlainRead(read_l);
	current_encoded_reads.emplace_back();
	current_encoded_reads.back().append(tuple_types::start_plain, 0, 0);
	for (size_t i = 0; i < read_l; ++i)	
		current_encoded_reads.back().append(tuple_types::plain, read[i], 1);	
}

void CEncoder::AddPlainReadWithN(const read_t& read)
{
	size_t read_l = read_len(read);
	stats.LogPlainReadWithN(read_l);
	current_encoded_reads.emplace_back();
	current_encoded_reads.back().append(tuple_types::start_plain_with_Ns, 0, 0);
	for (size_t i = 0; i < read_l; ++i)
		current_encoded_reads.back().append(tuple_types::plain, read[i], 1);
}


uint64_t CEncoder::get_number_of_matches(std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& first, std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& second)
{
	uint64_t res{};
	for (auto it = first.begin(); it != first.end(); ++it)
	{
		auto it2 = std::lower_bound(std::begin(second), std::end(second), it->first, [](const auto& e1, const auto& e2) {return e1.first < e2; });
		assert(it2 != second.end());
		assert(it2->first == it->first);
		res += it->second.size() * it2->second.size();
	}
	return res;
};

std::vector<std::pair<uint32_t, anchor_type>> CEncoder::Convert(const std::vector<std::pair<anchor_type, std::vector<uint32_t>>>& in, size_t len)
{
	std::vector<std::pair<uint32_t, anchor_type>> res;

	uint32_t res_size = 0;

	for (const auto& [mmer, vpos] : in)
		res_size += static_cast<uint32_t>(vpos.size());

	res.reserve(res_size);

	if (res_size < 0.05 * len)
	{
		for (const auto& [mmer, vpos] : in)
			for (auto pos : vpos)
				res.emplace_back(pos, mmer);

		std::sort(std::begin(res), std::end(res)); //pairs sorted by first element by default
	}
	else
	{
		std::vector<anchor_type> tmp(len, ~0ull);
		for (const auto& [mmer, vpos] : in)
			for (auto pos : vpos)
				tmp[pos] = mmer;

		for (size_t i = 0; i < len; ++i)
			if (tmp[i] != ~0ull)
				res.emplace_back(i, tmp[i]);
	}

	return res;
}

uint32_t CEncoder::MergeAnchors(const std::vector<std::tuple<anchor_type, uint32_t, uint32_t>>& inputAnchors, std::vector<Anchor>& res)
{
	res.clear();
	uint32_t start_index = 0;
	uint32_t tot_anch_len{};
	for (uint32_t i = 1; i < inputAnchors.size(); ++i)
	{
		//auto& [prev_mmer, prev_p1, prev_p2] = inputAnchors[i - 1];
		auto& prev_elem = inputAnchors[i - 1];
		auto& prev_p1 = std::get<1>(prev_elem);
		auto& prev_p2 = std::get<2>(prev_elem);

		//auto& [mmer, p1, p2] = inputAnchors[i];
		auto& cur_elem = inputAnchors[i];
		auto& p1 = std::get<1>(cur_elem);
		auto& p2 = std::get<2>(cur_elem);

		if (prev_p1 != p1 - 1 || prev_p2 != p2 - 1)
		{
			//const auto& [_, s1, s2] = inputAnchors[start_index];
			const auto& elem = inputAnchors[start_index];
			const auto& s1 = std::get<1>(elem);
			const auto& s2 = std::get<2>(elem);

			auto n_mmers = i - start_index;

			auto len = n_mmers + anchor_len - 1;
			res.emplace_back(len, s1, s2 );
			tot_anch_len += len;
			start_index = i;
		}
	}
	uint32_t i = static_cast<uint32_t>(inputAnchors.size());
	//const auto& [_, s1, s2] = inputAnchors[start_index];
	const auto& tmp = inputAnchors[start_index];
	const auto& s1 = std::get<1>(tmp);
	const auto& s2 = std::get<2>(tmp);

	auto n_mmers = i - start_index;

	auto len = n_mmers + anchor_len - 1;
	res.emplace_back(len, s1, s2);
	tot_anch_len += len;

	return tot_anch_len;
}

uint32_t CEncoder::AdjustAnchors(std::vector<Anchor>& anchors, uint32_t new_enc_start, uint32_t new_enc_end)
{
	uint32_t tot_len{};

	//find first anchor which ends after new start begin
	auto guard = std::numeric_limits<uint32_t>::max();
	uint32_t first_pos = guard;

	//replace with binary search is perfromance is issue here
	for (uint32_t i = 0; i < anchors.size(); ++i)
	{
		if (anchors[i].pos_in_enc + anchors[i].len > new_enc_start)
		{
			first_pos = i;
			break;
		}
	}
	if (first_pos == guard)
	{
		anchors.clear();
		return tot_len;
	}

	//if first anchor intersects with beginning of new part and its length would be smaller than anchor_len then we refuse such anchor (should we?)
	if (anchors[first_pos].pos_in_enc < new_enc_start && (anchors[first_pos].pos_in_enc + anchors[first_pos].len) - new_enc_start < anchor_len)
		++first_pos;

	//similar but for end
	uint32_t last_pos = guard;

	//replace with binary search is perfromance is issue here
	for (int32_t i = static_cast<int32_t>(anchors.size()) - 1; i >= 0; --i)
	{
		if (anchors[i].pos_in_enc < new_enc_end)
		{
			last_pos = i;
			break;
		}
	}
	if (last_pos == guard)
	{
		anchors.clear();
		return tot_len;
	}

	//if last anchor intersects with part end and its lenght is lower than anchor_len then we refuse such anchor (should we?)
	if (anchors[last_pos].pos_in_enc + anchors[last_pos].len > new_enc_end && new_enc_end - anchors[last_pos].pos_in_enc < anchor_len)
	{
		if (last_pos == 0)
		{
			anchors.clear();
			return tot_len;
		}
		--last_pos;
	}
	if (first_pos > last_pos)
	{
		anchors.clear();
		return tot_len;
	}

	anchors.erase(anchors.begin() + last_pos + 1, anchors.end());
	anchors.erase(anchors.begin(), anchors.begin() + first_pos);

	if (anchors.empty())
		return tot_len;

	//fix last end if needed
	if (anchors.back().pos_in_enc + anchors.back().len > new_enc_end)
	{
		anchors.back().len -= (anchors.back().pos_in_enc + anchors.back().len - new_enc_end);
	}

	for (uint32_t i = 0; i < anchors.size(); ++i)
	{
		if (i == 0 && anchors[0].pos_in_enc < new_enc_start)
		{
			auto diff = new_enc_start - anchors[0].pos_in_enc;
			anchors[0].len -= diff;
			anchors[0].pos_in_enc = 0;
			anchors[0].pos_in_ref += diff;
		}
		else
		{
			assert(anchors[i].pos_in_enc >= new_enc_start);
			anchors[i].pos_in_enc -= new_enc_start;
		}
		tot_len += anchors[i].len;
	}
	return tot_len;
}

AnalyseRefReadWithKmersRes CEncoder::AnalyseRefReadWithKmers(CKmers& enc_kmers, CBloomFilter& bloom_kmers, const read_t& enc_read, const read_t& ref_read, Candidate& candidate, CKmersHashSetLP& common_kmers)
{
	CKmersHashMapLP ref_kmers(ref_read, kmerLen, common_kmers, bloom_kmers, modulo);

	std::vector<Anchor> kmers_anchors;

	for (kmer_type kmer : common_kmers.GetKmers())
	{
		int64_t inEnc = enc_kmers.ExistsAndIsUnique(kmer);
		int64_t inRef = ref_kmers.ExistsAndIsUnique(kmer);
		if (inEnc == -1 || inRef == -1) //does not exists in one -> try rev. compl		
		{
			kmer = rev_compl(kmer, kmerLen);
			inEnc = enc_kmers.ExistsAndIsUnique(kmer);
			inRef = ref_kmers.ExistsAndIsUnique(kmer);
		}
		if (inEnc != -1 && inRef != -1)
			kmers_anchors.emplace_back(kmerLen, inEnc, inRef);
	}
	if (kmers_anchors.size() == 0) // no k-mer based anchors was found
		return AnalyseRefReadWithKmersRes::NoAnchors;


	std::sort(kmers_anchors.begin(), kmers_anchors.end(), [](const Anchor& lhs, const Anchor& rhs)
	{
		return lhs.pos_in_enc < rhs.pos_in_enc;
	});

	//if sorted k-mers positions in enc are not sorted in ref -> refuse
	//maybe it would be better to remove those k-mers that brokes sorted order

	if (!std::is_sorted(kmers_anchors.begin(), kmers_anchors.end(), [](const Anchor& lhs, const Anchor& rhs)
	{
		return lhs.pos_in_ref < rhs.pos_in_ref;
	}))
		return AnalyseRefReadWithKmersRes::CorespondingKmersIncompatibile; 

	//first remove overlapping k-mers, not very effective implementation (erases of a single element in the middle of a vector)
	//consider better implementation
	for (uint64_t i = 0; i < kmers_anchors.size() - 1; ++i)
	{
		if (kmers_anchors[i].pos_in_enc + kmers_anchors[i].len > kmers_anchors[i + 1].pos_in_enc ||
			kmers_anchors[i].pos_in_ref + kmers_anchors[i].len > kmers_anchors[i + 1].pos_in_ref)
		{
			kmers_anchors.erase(kmers_anchors.begin() + i + 1);
			--i;
			//			tot_anch_len -= kmerLen;
		}
	}

	//expand first kmer anchor to the left
	Anchor& first_anchor = kmers_anchors[0];

	while (first_anchor.pos_in_enc > 0 && first_anchor.pos_in_ref > 0)
	{
		if (enc_read[first_anchor.pos_in_enc - 1] == ref_read[first_anchor.pos_in_ref - 1])
		{
			--first_anchor.pos_in_enc;
			--first_anchor.pos_in_ref;
			++first_anchor.len;
			//			++tot_anch_len;
		}
		else
			break;
	}

	//In the best case all symbols between two k-mer based anchors are match, so it will join into a single anchor (again, currently erase in the middle of a vector in the implementation)
	for (uint64_t i = 0; i < kmers_anchors.size(); ++i)
	{
		if (i > 0)//expand to the left, except the first one 
		{
			auto prev_anch_end_enc = kmers_anchors[i - 1].pos_in_enc + kmers_anchors[i - 1].len;
			auto prev_anch_end_ref = kmers_anchors[i - 1].pos_in_ref + kmers_anchors[i - 1].len;

			while (true)
			{
				bool prev_reached_enc = kmers_anchors[i].pos_in_enc == prev_anch_end_enc;
				bool prev_reached_ref = kmers_anchors[i].pos_in_ref == prev_anch_end_ref;
				if (prev_reached_enc && prev_reached_ref)//merge anchors
				{
					kmers_anchors[i].len += kmers_anchors[i - 1].len;
					kmers_anchors.erase(kmers_anchors.begin() + i - 1);
					break;
				}
				if (prev_reached_enc || prev_reached_ref)
					break;
				if (enc_read[kmers_anchors[i].pos_in_enc - 1] != ref_read[kmers_anchors[i].pos_in_ref - 1])
					break;
				kmers_anchors[i].len++;
				kmers_anchors[i].pos_in_enc--;
				kmers_anchors[i].pos_in_ref--;
			}
		}
		//expand to the right
		if (i != kmers_anchors.size() - 1) //the last one is a special case of right expanding
		{
			auto next_anch_start_enc = kmers_anchors[i + 1].pos_in_enc;
			auto next_anch_start_ref = kmers_anchors[i + 1].pos_in_ref;

			auto pos_enc = kmers_anchors[i].pos_in_enc + kmers_anchors[i].len;
			auto pos_ref = kmers_anchors[i].pos_in_ref + kmers_anchors[i].len;

			while (true)
			{
				bool next_reached_enc = pos_enc == next_anch_start_enc;
				bool next_reached_ref = pos_ref == next_anch_start_ref;
				if (next_reached_enc && next_reached_ref) //merge anchors
				{
					kmers_anchors[i].len += kmers_anchors[i + 1].len;
					kmers_anchors.erase(kmers_anchors.begin() + i + 1);
					--i;
					break;
				}
				else if (next_reached_enc || next_reached_ref)
					break;
				if (enc_read[pos_enc] != ref_read[pos_ref])
					break;
				++pos_enc;
				++pos_ref;
				++kmers_anchors[i].len;
			}
		}
	}

	//expand after last anchor
	auto& last_anchor = kmers_anchors[kmers_anchors.size() - 1];
	auto pos_enc = last_anchor.pos_in_enc + last_anchor.len;
	auto pos_ref = last_anchor.pos_in_ref + last_anchor.len;
	while (pos_enc < read_len(enc_read) && pos_ref < read_len(ref_read))
	{
		if (enc_read[pos_enc] != ref_read[pos_ref])
			break;
		++pos_enc;
		++pos_ref;
		++last_anchor.len;
	}

	candidate.anchors = std::move(kmers_anchors);
	candidate.tot_anchor_len = 0;
	for (auto& a : candidate.anchors)
		candidate.tot_anchor_len += a.len; //maybe it could be calculated during steps above

	return AnalyseRefReadWithKmersRes::Accept;
}


AnalyseRefReadRes CEncoder::AnalyseRefRead(CMmers& encode_mmers, CBloomFilter &bloom_mmers, const read_t& enc_read, const read_t& ref_read, Candidate& candidate, int decision)
{
	CMmers candidate_mmers(ref_read, anchor_len, encode_mmers, bloom_mmers);

	std::vector<std::pair<anchor_type, std::vector<uint32_t>>> inter_mmers_enc_read_vec;
	std::vector<std::pair<anchor_type, std::vector<uint32_t>>> inter_mmers_ref_read_vec;
	encode_mmers.GetIntersection(candidate_mmers, inter_mmers_enc_read_vec, inter_mmers_ref_read_vec);

	std::sort(inter_mmers_ref_read_vec.begin(), inter_mmers_ref_read_vec.end());

	for (auto& x : inter_mmers_ref_read_vec)
		if(x.second.size() != 1)
			std::sort(x.second.begin(), x.second.end());

	//empty intersecion, in current implementation it should not be possible I think
	if (inter_mmers_enc_read_vec.size() == 0)
		return AnalyseRefReadRes::EmptyIntersection;

	if (decision != 0)
	{
		auto r = get_number_of_matches(inter_mmers_enc_read_vec, inter_mmers_ref_read_vec);
		decision = r > this->maxMatchesMultiplier * enc_read.size();
	}
	if (decision == 1) //too many matches
	{
		return AnalyseRefReadRes::TooManyMatches;
	}

	auto sorted_mmers_enc_vec = Convert(inter_mmers_enc_read_vec, read_len(enc_read));
	auto sorted_mmers_ref_vec = Convert(inter_mmers_ref_read_vec, read_len(ref_read));

	std::vector<std::tuple<anchor_type, uint32_t, uint32_t>> aligned_mmers_vec = get_aligned_mmers_LIS(sorted_mmers_enc_vec, sorted_mmers_ref_vec, inter_mmers_ref_read_vec);

	candidate.tot_anchor_len = MergeAnchors(aligned_mmers_vec, candidate.anchors);

	if (candidate.anchors.size() < minAnchors)
	{
		return AnalyseRefReadRes::TooLowAnchors;
	}
	return AnalyseRefReadRes::Accept;
}

EncodeCandidates CEncoder::prepareEncodeCandidates(const read_t& enc_read, const std::vector<uint32_t>& neighbours)
{
	EncodeCandidates res;
	//res.encode.enc_read_id = enc_id;
	//const auto& enc_read = reads.GetRead(enc_id);

	CMmers encode_mmers(enc_read, anchor_len);
	CBloomFilter bloom_mmers(enc_read, anchor_len);

	int decision = -1; //0 - encode, 1 - refuse, other - not know yet

	if (encode_mmers.GetNMmers() > minFractionOfMmersInEncodeToAlwaysEncode * read_len(enc_read))
		decision = 0;// encode
	else if (encode_mmers.GetNMmers() < minFractionOfMmersInEncode * read_len(enc_read))
		decision = 1;// refuse

	if (decision == 1)
	{
		stats.LogRefuse_not_enough_unique_mmers_in_enc_read();
		return res;
	}

	bool refused_too_many_matches = false;
	bool refused_too_low_anchors = false;

	for (const auto ref_read_id : neighbours)
	{
		Candidate candidate;
		candidate.ref_read_id = ref_read_id;
		auto ref_read = reference_reads.GetRefRead(ref_read_id);

		Candidate rev_compl_candidate;
		rev_compl_candidate.ref_read_id = ref_read_id;
		rev_compl_candidate.shouldReverse = true;
//		auto rev_compl_ref_read = get_rev_compl(ref_read);
		auto rev_compl_ref_read = reference_reads.GetRefRead(ref_read_id, true);

		MmerBasedAnchors(encode_mmers, bloom_mmers, enc_read, ref_read, rev_compl_ref_read, candidate, rev_compl_candidate, decision, res, refused_too_many_matches, refused_too_low_anchors);
	}

	if (res.candidates.empty())
	{
		if(refused_too_many_matches)
			stats.LogRefuse_to_many_matches();
		if (refused_too_low_anchors)
			stats.LogRefuse_too_low_anchors();
	}

	std::sort(std::begin(res.candidates), std::end(res.candidates), [](const auto& e1, const auto& e2) {
		return e1.tot_anchor_len > e2.tot_anchor_len;
	});

	return res;
}

bool CEncoder::KmerBasedAnchors(CKmersHashMapLP& enc_kmers, CBloomFilter& bloom_kmers, 
	const read_t& enc_read, const read_t& ref_read,
	const read_t& rev_compl_ref_read, Candidate& candidate, Candidate& rev_compl_candidate,  const std::vector<kmer_type>& common_kmers, EncodeCandidates& res)
{
	CKmersHashSetLP com_kmers(common_kmers);
	AnalyseRefReadWithKmersRes rev_cmpl_analyse_with_kmers_res = AnalyseRefReadWithKmers(enc_kmers, bloom_kmers, enc_read, rev_compl_ref_read, rev_compl_candidate, com_kmers);
	AnalyseRefReadWithKmersRes analyse_with_kmers_res = AnalyseRefReadWithKmers(enc_kmers, bloom_kmers, enc_read, ref_read, candidate, com_kmers);

	if (rev_cmpl_analyse_with_kmers_res == AnalyseRefReadWithKmersRes::Accept && analyse_with_kmers_res == AnalyseRefReadWithKmersRes::Accept)
	{
		if (candidate.tot_anchor_len > rev_compl_candidate.tot_anchor_len)
		{
			stats.LogNonRevChoosen();
			res.candidates.emplace_back(std::move(candidate));
		}
		else
		{
			stats.LogRevChoosen();
			res.candidates.emplace_back(std::move(rev_compl_candidate));
		}
	}
	else if (rev_cmpl_analyse_with_kmers_res == AnalyseRefReadWithKmersRes::Accept)
	{
		stats.LogRevChoosen();
		res.candidates.emplace_back(std::move(rev_compl_candidate));
	}
	else if (analyse_with_kmers_res == AnalyseRefReadWithKmersRes::Accept)
	{
		stats.LogNonRevChoosen();
		res.candidates.emplace_back(std::move(candidate));
	}
	else
		return false;
	return true;
}

void CEncoder::MmerBasedAnchors(CMmers& encode_mmers, CBloomFilter& bloom_mmers,
	const read_t& enc_read, const read_t& ref_read, const read_t& rev_compl_ref_read,
	Candidate& candidate, Candidate& rev_compl_candidate, int decision, EncodeCandidates& res,
	bool& refused_too_many_matches, bool& refuled_too_low_anchors)
{
	auto rev_cmpl_analyse_res = AnalyseRefRead(encode_mmers, bloom_mmers, enc_read, rev_compl_ref_read, rev_compl_candidate, decision);
	auto analyse_res = AnalyseRefRead(encode_mmers, bloom_mmers, enc_read, ref_read, candidate, decision);

	if (rev_cmpl_analyse_res == AnalyseRefReadRes::Accept && analyse_res == AnalyseRefReadRes::Accept)
	{
		if (candidate.tot_anchor_len > rev_compl_candidate.tot_anchor_len)
		{
			stats.LogNonRevChoosen();
			res.candidates.emplace_back(std::move(candidate));
		}
		else
		{
			stats.LogRevChoosen();
			res.candidates.emplace_back(std::move(rev_compl_candidate));
		}
	}
	else if (rev_cmpl_analyse_res == AnalyseRefReadRes::Accept)
	{
		stats.LogRevChoosen();
		res.candidates.emplace_back(std::move(rev_compl_candidate));
	}
	else if (analyse_res == AnalyseRefReadRes::Accept)
	{
		stats.LogNonRevChoosen();
		res.candidates.emplace_back(std::move(candidate));
	}
	else // non is accepted
	{
		if (analyse_res == AnalyseRefReadRes::TooManyMatches)
			refused_too_many_matches = true;
		else if (analyse_res == AnalyseRefReadRes::TooLowAnchors)
			refuled_too_low_anchors = true;

		if (rev_cmpl_analyse_res == AnalyseRefReadRes::TooManyMatches)
			refused_too_many_matches = true;
		else if (rev_cmpl_analyse_res == AnalyseRefReadRes::TooLowAnchors)
			refuled_too_low_anchors = true;
	}
}

EncodeCandidates CEncoder::prepareEncodeCandidatesHiFi(const read_t& enc_read, const std::vector<uint32_t>& neighbours, const std::vector<std::vector<kmer_type>>& common_kmers)
{
	EncodeCandidates res;
	//res.encode.enc_read_id = enc_id;
	//const auto& enc_read = reads.GetRead(enc_id);

	CMmers encode_mmers(enc_read, anchor_len);
	CBloomFilter bloom_mmers(enc_read, anchor_len);

	int decision = -1; //0 - encode, 1 - refuse, other - not know yet

	if (encode_mmers.GetNMmers() > minFractionOfMmersInEncodeToAlwaysEncode * read_len(enc_read))
		decision = 0;// encode
	else if (encode_mmers.GetNMmers() < minFractionOfMmersInEncode * read_len(enc_read))
		decision = 1;// refuse

	if (decision == 1)
	{
		stats.LogRefuse_not_enough_unique_mmers_in_enc_read();
		return res;
	}

	bool refused_too_many_matches = false;
	bool refuled_too_low_anchors = false;
	
	CKmersHashMapLP enc_kmers(enc_read, kmerLen);
	CBloomFilter bloom_kmers(enc_read, kmerLen, true);

	for (uint64_t i = 0; i < neighbours.size(); ++i)
	{
		const auto ref_read_id = neighbours[i];

		Candidate candidate;
		candidate.ref_read_id = ref_read_id;
		auto ref_read = reference_reads.GetRefRead(ref_read_id);

		Candidate rev_compl_candidate;
		rev_compl_candidate.ref_read_id = ref_read_id;
		rev_compl_candidate.shouldReverse = true;
		auto rev_compl_ref_read = reference_reads.GetRefRead(ref_read_id, true);

		if (KmerBasedAnchors(enc_kmers, bloom_kmers, enc_read, ref_read, rev_compl_ref_read, candidate, rev_compl_candidate, common_kmers[i], res))
			continue;

		//if k-mer based anchors analysis fails
		MmerBasedAnchors(encode_mmers, bloom_mmers, enc_read, ref_read, rev_compl_ref_read, candidate, rev_compl_candidate, decision, res, refused_too_many_matches, refuled_too_low_anchors);
	}
	if (res.candidates.empty())
	{
		if(refused_too_many_matches)
			stats.LogRefuse_to_many_matches();
		if (refuled_too_low_anchors)
			stats.LogRefuse_too_low_anchors();
	}
	std::sort(std::begin(res.candidates), std::end(res.candidates), [](const auto& e1, const auto& e2) {
		return e1.tot_anchor_len > e2.tot_anchor_len;
	});

	return res;
}

EditDistRes CEncoder::GetEditDist(read_view refPart, read_view encPart, uint32_t frag_no, uint32_t n_fragments)
{
	if (refPart.empty() || encPart.empty())
		return get_edit_dist_on_seq_empty(refPart, encPart);

	EditDistRes ed;

	uint32_t max_symbols_for_flank = static_cast<uint32_t>(encPart.size() * 2);
	if (frag_no == 0)
	{
		uint32_t ref_offset;
		ed = find_edit_dist_with_edlib_ex_odwr_reverse(refPart, encPart, max_symbols_for_flank, ref_offset, EDLIB_MODE_SHW);
		refactor_edit_script(refPart.substr(ref_offset), encPart, ed.editScript);
		ed.editScript = std::string(ref_offset, 'D') + ed.editScript;
	}
	else if (frag_no == n_fragments - 1)
	{
		uint32_t tmp;
		ed = find_edit_dist_with_edlib_ex_odwr(refPart.substr(0, max_symbols_for_flank), encPart, tmp, EDLIB_MODE_SHW);
		refactor_edit_script(refPart, encPart, ed.editScript);
	}
	else
	{
		ed = find_edit_dist_with_edlib_ex(refPart, encPart);
		refactor_edit_script(refPart, encPart, ed.editScript);
	}

	return ed;
}



uint32_t CEncoder::CountDeletions(std::string_view es)
{
	uint32_t nD{};
	for (auto c : es)
		if (c == 'D')
			++nD;
		else
			break;
	return nD;
}


std::string_view CEncoder::GetEditScriptEntropyInput(std::string_view es)
{
	//count deletions at the beginning
	uint32_t min_D_at_start_to_skip_them_all = 10; //number of deletions at the beginning to skip them during entropy estimation
	uint32_t nD = CountDeletions(es);

	std::string_view edit_script_entropy_input = es;
	if (nD >= min_D_at_start_to_skip_them_all)
	{
		edit_script_entropy_input = edit_script_entropy_input.substr(nD);
	}
	return edit_script_entropy_input;
}



bool CEncoder::EncodeWithEditScript(const EditDistRes& ed, read_view refPart, read_view encPart, uint32_t frag_no, uint32_t n_fragments)
{
	std::string_view edit_script_entropy_input = GetEditScriptEntropyInput(ed.editScript);

//	auto entropy_edit_script = entropy(edit_script_entropy_input);
	auto entropy_edit_script = CEntropy::entropy_es(edit_script_entropy_input);

	auto plain_entropy_input = encPart;
//	auto entropy_plain = entropy(plain_entropy_input);
	auto entropy_plain = CEntropy::entropy_dna(plain_entropy_input);

	return entropy_edit_script * edit_script_entropy_input.length() * editScriptCostMultiplier < entropy_plain * plain_entropy_input.length();
}

bool CEncoder::EncodeWithAlternativeRead(const std::vector<Candidate>& candidates, uint32_t level, read_view encPart, uint32_t cur_pos_in_encode_read, uint32_t end_enc,
	std::vector<Candidate>& alt_candidates)
{
	if (candidates.size() <= level + 1 || encPart.length() < this->minPartLenToConsiderAltRead || level >= maxRecurence)
		return false;


	alt_candidates = candidates;
	//recalculate and cut anchors 
	for (uint32_t i = level + 1; i < candidates.size(); ++i)
	{
		alt_candidates[i].tot_anchor_len = AdjustAnchors(alt_candidates[i].anchors, cur_pos_in_encode_read, end_enc);		
	}
	std::sort(alt_candidates.begin() + level + 1, alt_candidates.end(), [](const auto& e1, const auto& e2) {
		return e1.tot_anchor_len > e2.tot_anchor_len;
	});
	return alt_candidates[level + 1].tot_anchor_len != 0;
}

void singleEditScriptSymbolStore(char symb, uint32_t no_rep, es_t& encoded_read)
{
	bool store_with_loop = true;
	tuple_types type;
	int code = 0;
	if (symb == 'M')
	{
		type = tuple_types::match;
		if (no_rep >= min_anchor_len)
		{
			type = tuple_types::anchor;
			store_with_loop = false;
		}
	}
	else if (symb == 'D')
	{
		type = tuple_types::deletion;
		if (no_rep > 16)
		{
			type = tuple_types::skip;
			store_with_loop = false;
		}
	}
	else if (CMissmatchCoder::is_missmatch_es_symb(symb))
	{
		type = tuple_types::substitution;
		code = CMissmatchCoder::GetMissmatchCode(symb);
	}
	else if (symb == 'A' || symb == 'C' || symb == 'G' || symb == 'T')
	{
		type = tuple_types::insertion;
		code = SymbToBinMap[(uint8_t)symb];
	}
	else
	{
		std::cerr << "Fatal error: unexpected symbol in edit script: " << symb << "\n";
		exit(1);
	}

	if (store_with_loop)
		for (uint32_t i = 0; i < no_rep; ++i)
			encoded_read.append(type, code, 1);
	else
		encoded_read.append(type, code, no_rep);
}

void bigEditScriptToTuples(es_t& encoded_read, const std::string& currentBigEditScirpt)
{
	assert(currentBigEditScirpt.size());
	auto symb = currentBigEditScirpt[0];
	uint32_t no_rep = 1;
	for (uint32_t i = 1; i < currentBigEditScirpt.size(); ++i)
	{
		if (currentBigEditScirpt[i] != symb)
		{			
			singleEditScriptSymbolStore(symb, no_rep, encoded_read);

			symb = currentBigEditScirpt[i];
			no_rep = 1;
		}
		else
			++no_rep;
	}
	singleEditScriptSymbolStore(symb, no_rep, encoded_read);
}

void CEncoder::StoreFrag(std::string& currentBigEditScirpt, uint32_t level, es_t& encoded_read, uint32_t reference_read_id, uint32_t main_ref_id,
	uint32_t& last_pos_in_ref, uint32_t cur_pos_in_ref_read, bool& first, bool shouldReverse)
{
	if (currentBigEditScirpt.size()) //enything to write
	{
		if (level == 0)
		{
			if (reference_read_id != main_ref_id)
				encoded_read.append(tuple_types::alt_id, reference_read_id, shouldReverse);
			else if(!first) //if this is firs edit script store, then we do not need to go back to main ref as we are in main ref
				encoded_read.append(tuple_types::main_ref, 0, shouldReverse);
			
			bigEditScriptToTuples(encoded_read, currentBigEditScirpt);			
		}
		else
		{
			if (reference_read_id != main_ref_id)
				encoded_read.append(tuple_types::alt_id, reference_read_id, shouldReverse);
			else
				encoded_read.append(tuple_types::main_ref, 0, shouldReverse);
			std::string dels = std::string(last_pos_in_ref, 'D');						

			bigEditScriptToTuples(encoded_read, dels + currentBigEditScirpt);

		}
		last_pos_in_ref = cur_pos_in_ref_read;
		first = false;
	}
	currentBigEditScirpt.clear();
}

void CEncoder::EncodePart(uint32_t level, read_view encode_read, uint32_t frag_no, uint32_t n_fragments,
	uint32_t end_enc, uint32_t end_ref, std::vector<Candidate>& candidates,
	uint32_t reference_read_id, 
	read_t &ref_read,
	uint32_t main_ref_id, uint32_t cur_pos_in_ref_read, uint32_t cur_pos_in_encode_read,
	std::string& currentBigEditScirpt, es_t& encoded_read,
	uint32_t& last_pos_in_ref, bool& first)
{
	read_view refPart = read_view(ref_read).substr(cur_pos_in_ref_read, end_ref - cur_pos_in_ref_read);
	read_view encPart = encode_read.substr(cur_pos_in_encode_read, end_enc - cur_pos_in_encode_read);

	EditDistRes ed = GetEditDist(refPart, encPart, frag_no, n_fragments);

	bool decision;

	if (encPart.length() < minPartLenToConsiderAltRead)
		decision = entropyEstimator.EncodeWithEditScript(ed.editScript, encPart, refPart.length());
	else
		decision = EncodeWithEditScript(ed, refPart, encPart, frag_no, n_fragments);

	//if (EncodeWithEditScript(ed, refPart, encPart, frag_no, n_fragments))
//	if(entropyEstimator.EncodeWithEditScript(ed.editScript, encPart, refPart.length()))
	if(decision)
	{
		currentBigEditScirpt += ed.editScript;

		stats.LogEditScript(level, ed.editScript, end_enc - cur_pos_in_encode_read);

	}
	else
	{
		std::vector<Candidate> alt_candidates;

		if (EncodeWithAlternativeRead(candidates, level, encPart, cur_pos_in_encode_read, end_enc, alt_candidates))
		{
			if (frag_no == n_fragments - 1)
				stats.ComprStats(level).n_alternative_right_flank++;
			else if (frag_no == 0)
				stats.ComprStats(level).n_alternative_left_flank++;
			else
				stats.ComprStats(level).n_alternative_in_between++;

			StoreFrag(currentBigEditScirpt, level, encoded_read, reference_read_id, main_ref_id, last_pos_in_ref, cur_pos_in_ref_read, first, candidates[level].shouldReverse);

			AddEncodedReadWithCandidates(encPart, alt_candidates, level + 1, encoded_read, main_ref_id, first);

			if (frag_no != n_fragments - 1) // for the last fragment it does not matter how many symbols needs to be skiped
			{
				currentBigEditScirpt.append(end_ref - cur_pos_in_ref_read, 'D');
			}
		}

		else
		{
			stats.ComprStats(level).n_plain_symbols += encPart.length();

			for (uint8_t s : encPart)
				currentBigEditScirpt.push_back("ACGT"[s]);
			if (frag_no != n_fragments - 1) // for the last fragment it does not matter how many symbols needs to be skiped
			{
				uint32_t len_ref = (end_ref - cur_pos_in_ref_read);
				currentBigEditScirpt.append(len_ref, 'D'); //store deletions, except for the last part
			}
		}

	}
}

void CEncoder::AddEncodedReadWithCandidates(const read_view encode_read, std::vector<Candidate>& candidates, uint32_t level, es_t& encoded_read, uint32_t main_ref_id, bool& first)
{
	stats.LogLevel(level);
	uint32_t reference_read_id = candidates[level].ref_read_id;
	const auto& anchors = candidates[level].anchors;

	uint32_t last_pos_in_ref = 0;

	if (level == 0)
	{
		assert(reference_read_id == main_ref_id);		
		encoded_read.append(tuple_types::start_es, main_ref_id, candidates[0].shouldReverse);
	}

	std::string currentBigEditScirpt;

	currentBigEditScirpt.reserve(encode_read.size() / 8);			// rough estimation

	std::map<uint32_t, read_t> buffered_ref_reads;

	uint32_t n_fragments = static_cast<uint32_t>(anchors.size() * 2 + 1);
	uint32_t anch_id{};
	uint32_t cur_pos_in_ref_read = 0;
	uint32_t cur_pos_in_encode_read = 0;

	stats.ComprStats(level).n_left_flank_symb += anchors[0].pos_in_enc;
	stats.ComprStats(level).n_right_flank_symb += encode_read.length() - (anchors.back().pos_in_enc + anchors.back().len);

	for (uint32_t i = 0; i < n_fragments; ++i)
	{
		if (i % 2 == 0)
		{
			uint32_t end_enc = i == n_fragments - 1 ? static_cast<uint32_t>(encode_read.length()) : anchors[anch_id].pos_in_enc;
//			uint32_t end_ref = i == n_fragments - 1 ? static_cast<uint32_t>(read_len(reference_reads.GetRefRead(reference_read_id))) : anchors[anch_id].pos_in_ref;
			uint32_t end_ref = i == n_fragments - 1 ? static_cast<uint32_t>(reference_reads.GetRefReadLen(reference_read_id)) : anchors[anch_id].pos_in_ref;

			if (buffered_ref_reads.count(reference_read_id) == 0)
				buffered_ref_reads.emplace(reference_read_id, reference_reads.GetRefRead(reference_read_id, candidates[level].shouldReverse));


			EncodePart(level, encode_read, i, n_fragments, end_enc, end_ref, candidates, 
				reference_read_id, buffered_ref_reads[reference_read_id], main_ref_id, cur_pos_in_ref_read, cur_pos_in_encode_read, currentBigEditScirpt,
				encoded_read, last_pos_in_ref,
				first);
		}
		else //encode anchor
		{
			auto pos_ref = anchors[anch_id].pos_in_ref;
			auto pos_enc = anchors[anch_id].pos_in_enc;
			auto len = anchors[anch_id].len;
			currentBigEditScirpt.append(len, 'M');

			++stats.ComprStats(level).n_anchors;
			stats.ComprStats(level).n_symb_anchors += len;

			cur_pos_in_ref_read = pos_ref + len;
			cur_pos_in_encode_read = pos_enc + len;
			++anch_id;
		}
	}

	StoreFrag(currentBigEditScirpt, level, encoded_read, reference_read_id, main_ref_id, last_pos_in_ref, cur_pos_in_ref_read, first, candidates[level].shouldReverse);
}

void CEncoder::fixOverlapingInEncodeRead(std::vector<Anchor>& anchors)
{
	uint32_t curID = 0;
	while (curID < anchors.size() - 1)
	{
		uint32_t nextID = curID + 1;
		auto curEnd = anchors[curID].pos_in_enc + anchors[curID].len;

		//something need to be fixed
		if (anchors[nextID].pos_in_enc < curEnd)
		{
			uint32_t roznica = curEnd - anchors[nextID].pos_in_enc;

			//two anchors may share at most m - 1 symbos, so always it is possible to reduce lenght of latter anchor, that it is not empty
			anchors[nextID].pos_in_enc += roznica;
			assert(anchors[nextID].len > roznica);
			anchors[nextID].len -= roznica;
			anchors[nextID].pos_in_ref += roznica;
		}
		++curID;
	}
}


void CEncoder::fixOverlapingInReferenceRead(std::vector<Anchor>& anchors)
{
	uint32_t curID = 0;
	while (curID < anchors.size() - 1)
	{
		uint32_t nextID = curID + 1;
		auto curEnd = anchors[curID].pos_in_ref + anchors[curID].len;

		//something need to be fixed
		if (anchors[nextID].pos_in_ref < curEnd)
		{
			uint32_t roznica = curEnd - anchors[nextID].pos_in_ref;

			//two anchors may share at most m - 1 symbos, so always it is possible to reduce lenght of latter anchor, that it is not empty
			anchors[nextID].pos_in_ref += roznica;
			assert(anchors[nextID].len > roznica); 
			anchors[nextID].len -= roznica;
			anchors[nextID].pos_in_enc += roznica;
		}
		++curID;
	}
}


void CEncoder::processComprElem(const CCompressElem& compr_elem)
{

	const auto& neighbours = compr_elem.ref_reads;
	if (compr_elem.hasN)
	{
		AddPlainReadWithN(compr_elem.read);
		return;
	}
	entropyEstimator.LogRead(compr_elem.read);
	if (neighbours.empty())
	{
		AddPlainRead(compr_elem.read);
		return;
	}
	EncodeCandidates encodeCandidates = dataSource == DataSource::PBHiFi ? 
		prepareEncodeCandidatesHiFi(compr_elem.read, neighbours, compr_elem.common_kmers)
	: 
		prepareEncodeCandidates(compr_elem.read, neighbours);

	if (encodeCandidates.candidates.empty())
	{
		AddPlainRead(compr_elem.read);
		return;
	}

	for (auto& c : encodeCandidates.candidates)
	{
		fixOverlapingInReferenceRead(c.anchors);
		fixOverlapingInEncodeRead(c.anchors);
	}

	bool first = true; //control if id should be saved, not very nice solution
	
	current_encoded_reads.emplace_back();
	AddEncodedReadWithCandidates(compr_elem.read, encodeCandidates.candidates, 0, current_encoded_reads.back(), encodeCandidates.candidates[0].ref_read_id, first);
}


bool CEncoder::pop(CCompressPack& compress_pack)
{
	auto start = std::chrono::high_resolution_clock::now();
	bool is = compress_queue.Pop(compress_pack);
	waitOnQueue += std::chrono::high_resolution_clock::now() - start;
	return is;
}

void CEncoder::Encode()
{
	CCompressPack compress_pack;
	while (pop(compress_pack))
	{
		entropyEstimator.Reset();
		for (const auto& comp_elem : compress_pack.data)
			processComprElem(comp_elem);

		if (is_fastq)
		{
			auto current_encoded_reads_copy = current_encoded_reads;
			edit_script_for_qual_queue.Push(compress_pack.id, std::move(current_encoded_reads_copy));
		}

		compressed_queue.Push(compress_pack.id, std::move(current_encoded_reads));
	}
	compressed_queue.MarkCompleted();
	edit_script_for_qual_queue.MarkCompleted();
}
