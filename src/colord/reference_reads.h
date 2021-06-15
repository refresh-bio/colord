#pragma once
#include "utils.h"
#include <vector>
#include <mutex>
#include <cstring>

#ifdef USE_COMPACTED_READS

class CReferenceReads
{
	std::vector<read_t> ref_reads;
	uint32_t cur_read_pos{};

	uint8_t decompact_dir_lookup[256][4];
	uint8_t decompact_rc_lookup[256][4];

	read_t compact(const read_t& read) const
	{
		auto in_len = read_len(read);
		size_t out_len = (in_len + 3) / 4 + 1; //one for number of symbols in last byte
		read_t res(out_len);

		auto full_bytes = in_len / 4; //number of fully filled bytes in output
		uint8_t symbols_in_last_byte = in_len % 4;

		uint32_t in_pos = 0;
		for (uint32_t i = 0; i < full_bytes; ++i)
		{
			res[i] = read[in_pos++] << 6;
			res[i] += read[in_pos++] << 4;
			res[i] += read[in_pos++] << 2;
			res[i] += read[in_pos++];
		}

		switch (symbols_in_last_byte)
		{
		case 3:
			res[out_len - 2] += read[in_pos++] << 6;
			res[out_len - 2] += read[in_pos++] << 4;
			res[out_len - 2] += read[in_pos++] << 2;
			break;
		case 2:
			res[out_len - 2] += read[in_pos++] << 6;
			res[out_len - 2] += read[in_pos++] << 4;
			break;
		case 1:
			res[out_len - 2] += read[in_pos++] << 6;
			break;
		}
		assert(in_pos == in_len);
		res[out_len - 1] = symbols_in_last_byte;

		return res;
	}

	read_t decompact_dir(const read_t& read) const
	{
		auto symbols_in_last_byte = read.back();
		auto full_bytes = read.size() - 2 + !symbols_in_last_byte;
		auto out_len = full_bytes * 4 + symbols_in_last_byte + 1;

		read_t res(out_len);
		uint64_t out_pos = 0;

		auto ptr = res.data();

		auto read_ptr = &read[0];
		auto read_end = read_ptr + full_bytes;

		switch (full_bytes % 4)
		{
		case 3:
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
		case 2:
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
		case 1:
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
		}

		while (read_ptr != read_end)
		{
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
			memcpy(ptr, decompact_dir_lookup[*read_ptr++], 4);
			ptr += 4;
		}

/*		for (uint32_t i = 0; i < full_bytes; ++i)
		{
			memcpy(ptr, decompact_dir_lookup[read[i]], 4);
			ptr += 4;
		}*/

		out_pos = full_bytes * 4;

		switch (symbols_in_last_byte)
		{
		case 3:
			res[out_pos++] += read[read.size() - 2] >> 6;
			res[out_pos++] += (read[read.size() - 2] >> 4) & 3;
			res[out_pos++] += (read[read.size() - 2] >> 2) & 3;
			break;
		case 2:
			res[out_pos++] += (read[read.size() - 2] >> 6) & 3;
			res[out_pos++] += (read[read.size() - 2] >> 4) & 3;
			break;
		case 1:
			res[out_pos++] += (read[read.size() - 2] >> 6) & 3;
			break;
		}

		res[out_pos++] = 255;
		assert(out_pos == out_len);
		return res;
	}

	read_t decompact_rc(const read_t& read) const
	{
		auto symbols_in_last_byte = read.back();
		auto full_bytes = read.size() - 2 + !symbols_in_last_byte;
		auto out_len = full_bytes * 4 + symbols_in_last_byte + 1;

		read_t res(out_len);
//		uint32_t out_pos = 0;

		auto ptr = res.data() + symbols_in_last_byte;

/*		for (uint32_t i = 0; i < full_bytes; ++i)
		{
			memcpy(ptr, decompact_rc_lookup[read[full_bytes - i - 1u]], 4);
			ptr += 4;
		}*/

		auto read_ptr = &read[0] + full_bytes;
		auto read_start = &read[0];

		switch (full_bytes % 4)
		{
		case 3:
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
		case 2:
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
		case 1:
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
		}

		while (read_ptr != read_start)
		{
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
			memcpy(ptr, decompact_rc_lookup[*--read_ptr], 4);
			ptr += 4;
		}


		switch (symbols_in_last_byte)
		{
		case 3:
			res[2] += 3u - (read[read.size() - 2] >> 6);
			res[1] += 3u - ((read[read.size() - 2] >> 4) & 3);
			res[0] += 3u - ((read[read.size() - 2] >> 2) & 3);
			break;
		case 2:
			res[1] += 3u - ((read[read.size() - 2] >> 6) & 3);
			res[0] += 3u - ((read[read.size() - 2] >> 4) & 3);
			break;
		case 1:
			res[0] += 3u - ((read[read.size() - 2] >> 6) & 3);
			break;
		}

		res[out_len - 1u] = 255;
//		assert(out_pos == out_len);
		return res;
	}

public:
	CReferenceReads(uint32_t n_reads) :
		ref_reads(n_reads)
	{
		for (uint32_t i = 0; i < 256; ++i)
		{
			decompact_dir_lookup[i][0] = i >> 6;
			decompact_dir_lookup[i][1] = (i >> 4) & 3;
			decompact_dir_lookup[i][2] = (i >> 2) & 3;
			decompact_dir_lookup[i][3] = i & 3;

			decompact_rc_lookup[i][3] = 3u - (i >> 6);
			decompact_rc_lookup[i][2] = 3u - ((i >> 4) & 3);
			decompact_rc_lookup[i][1] = 3u - ((i >> 2) & 3);
			decompact_rc_lookup[i][0] = 3u - (i & 3);
		}
	}

	void Add(const read_t& ref_read)
	{
		ref_reads[cur_read_pos++] = compact(ref_read);
	}
	read_t GetRefRead(uint32_t id) const
	{
		return decompact_dir(ref_reads[id]);
	}

	read_t GetRefRead(uint32_t id, bool rev_comp) const
	{
		if(rev_comp)
			return decompact_rc(ref_reads[id]);
		else
			return decompact_dir(ref_reads[id]);
	}

	uint32_t GetRefReadLen(uint32_t id) const
	{
		auto symbols_in_last_byte = ref_reads[id].back();
		auto full_bytes = ref_reads[id].size() - 2 + !symbols_in_last_byte;
		auto out_len = full_bytes * 4 + symbols_in_last_byte + 1;

		return static_cast<uint32_t>(out_len) - 1;
	}

	void PrintMemoryUsage() const
	{
		uint64_t tot_size_bytes{};
		for (auto& r : ref_reads)
			tot_size_bytes += r.size();
		tot_size_bytes += sizeof(read_t) * ref_reads.size();
		std::cerr << "reference reads memory usage: " << (tot_size_bytes / 1024 / 1024) << "MiB\n";
	}
};
#else
class CReferenceReads
{
	std::vector<read_t> ref_reads;
	uint32_t cur_read_pos{};	
public:
	CReferenceReads(uint32_t n_reads) :
		ref_reads(n_reads)
	{
		
	}

	void Add(const read_t& ref_read)
	{				
		ref_reads[cur_read_pos++] = ref_read;
	}
	const read_t& GetRefRead(uint32_t id) const
	{		
		return ref_reads[id];		
	}

	void PrintMemoryUsage() const
	{
		uint64_t tot_size_bytes{};
		for (auto& r : ref_reads)
			tot_size_bytes += r.size();
		tot_size_bytes += sizeof(read_t) * ref_reads.size();
		std::cerr << "reference reads memory usage: " << (tot_size_bytes / 1024 / 1024) << "MiB\n";
	}
};
#endif