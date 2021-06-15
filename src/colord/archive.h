#pragma once

#include <cstdio>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <thread>
#include <mutex>

using namespace std;

class CArchive
{
	bool input_mode;
	FILE* f;
	size_t f_offset;

	struct part_t{
		size_t offset;
		size_t size;

		part_t() : offset(0), size(0)
		{};

		part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
		{};
	};

	typedef struct {
		string stream_name;
		size_t cur_id;
		size_t raw_size;
		size_t packed_size;
		size_t packed_data_size;
		vector<part_t> parts;
	} stream_t;

	map<size_t, stream_t> m_streams;
	mutex mtx;

	bool serialize();
	bool deserialize();
	size_t write_fixed(size_t x);
	size_t write(size_t x);
	size_t write(string s);
	size_t read_fixed(size_t& x);
	size_t read(size_t& x);
	size_t read(string& s);

public:
	CArchive(bool _input_mode);
	~CArchive();

	bool Open(string file_name);
	bool Close();

	int RegisterStream(string stream_name);
	int GetStreamId(string stream_name);
	
	size_t GetStreamPackedSize(int stream_id)
	{
		lock_guard<mutex> lck(mtx);

		if (stream_id < 0 || stream_id >= static_cast<int>(m_streams.size()))
			return 0;

		return m_streams[stream_id].packed_size;
	}

	size_t GetStreamPackedDataSize(int stream_id)
	{
		lock_guard<mutex> lck(mtx);

		if (stream_id < 0 || stream_id >= static_cast<int>(m_streams.size()))
			return 0;

		return m_streams[stream_id].packed_data_size;
	}

	bool AddPart(int stream_id, vector<uint8_t> &v_data, size_t metadata = 0);
	int AddPartPrepare(int stream_id);
	bool AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data, size_t metadata = 0);

	bool GetPart(int stream_id, vector<uint8_t> &v_data, size_t &metadata);
	void SetRawSize(int stream_id, size_t raw_size);
	size_t GetRawSize(int stream_id);

	size_t GetNoStreams()
	{
		lock_guard<mutex> lck(mtx);

		return m_streams.size();
	}
};

// EOF
