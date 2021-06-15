#include "archive.h"

#include <iostream>

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#endif

// ******************************************************************************
CArchive::CArchive(bool _input_mode)
{
	f = nullptr;
	input_mode = _input_mode;
}

// ******************************************************************************
CArchive::~CArchive()
{
	if (f)
		Close();
}

// ******************************************************************************
bool CArchive::Open(string file_name)
{
	lock_guard<mutex> lck(mtx);

	if (f)
		fclose(f);

	f = fopen(file_name.c_str(), input_mode ? "rb" : "wb");

	setvbuf(f, nullptr, _IOFBF, 64 << 20);

	if (!f)
		return false;

	if (input_mode)
		deserialize();

	f_offset = 0;

	return true;
}

// ******************************************************************************
bool CArchive::Close()
{
	lock_guard<mutex> lck(mtx);

	if (!f)
		return false;

	if (input_mode)
	{
		fclose(f);
		f = nullptr;
	}
	else
	{
		serialize();
		fclose(f);
		f = nullptr;
	}

	return true;
}

// ******************************************************************************
size_t CArchive::write_fixed(size_t x)
{
	fwrite(&x, 1, 8, f);

	return 8;
}

// ******************************************************************************
size_t CArchive::write(size_t x)
{
	int no_bytes = 0;

	for (size_t tmp = x; tmp; tmp >>= 8)
		++no_bytes;
	
	putc(no_bytes, f);

	for (int i = no_bytes; i; --i)
		putc((x >> ((i - 1) * 8)) & 0xff, f);

	return no_bytes + 1;
}

// ******************************************************************************
size_t CArchive::write(string s)
{
	fwrite(s.c_str(), 1, s.size(), f);
	putc(0, f);

	return s.size() + 1;
}

// ******************************************************************************
size_t CArchive::read_fixed(size_t& x)
{
	if (fread(&x, 8, 1, f))
		return 8;
	
	return 0;
}

// ******************************************************************************
size_t CArchive::read(size_t& x)
{
	int no_bytes = getc(f);

	x = 0;

	for (int i = 0; i < no_bytes; ++i)
	{
		x <<= 8;
		x += (size_t)getc(f);
	}

	return no_bytes + 1;
}

// ******************************************************************************
size_t CArchive::read(string& s)
{
	s.clear();

	while (true)
	{
		int c = getc(f);
		if (c == EOF)
			return 0;

		if (c == 0)
			return s.size() + 1;

		s.push_back((char)c);
	}

	return 0;
}

// ******************************************************************************
bool CArchive::serialize()
{
	size_t footer_size = 0;

	// Zapisuje informacje o po³o¿eniu kawa³ków strumieni
	footer_size += write(m_streams.size());

	for (auto& stream : m_streams)
	{
		size_t p = footer_size;

		footer_size += write(stream.second.stream_name);
		footer_size += write(stream.second.parts.size());
		footer_size += write(stream.second.raw_size);

		for (auto& part : stream.second.parts)
		{
			footer_size += write(part.offset);
			footer_size += write(part.size);
//			stream.second.packed_data_size += part.size;
		}

		stream.second.packed_size += footer_size - p;
	}

	write_fixed(footer_size);

	return true;
}

// ******************************************************************************
bool CArchive::deserialize()
{
	size_t footer_size;

	my_fseek(f, -8, SEEK_END);
	read_fixed(footer_size);

	my_fseek(f, -(long)(8 + footer_size), SEEK_END);

	// Odczytuje informacje o po³o¿eniu kawa³ków strumieni
	size_t n_streams;
	read(n_streams);

	for (size_t i = 0; i < n_streams; ++i)
	{
		m_streams[i] = stream_t();
		auto& stream_second = m_streams[i];

		read(stream_second.stream_name);
		read(stream_second.cur_id);
		read(stream_second.raw_size);

		stream_second.parts.resize(stream_second.cur_id);
		for (size_t j = 0; j < stream_second.cur_id; ++j)
		{
			read(stream_second.parts[j].offset);
			read(stream_second.parts[j].size);
		}

		stream_second.cur_id = 0;
	}
	
	my_fseek(f, 0, SEEK_SET);

	return true;
}

// ******************************************************************************
int CArchive::RegisterStream(string stream_name)
{
	lock_guard<mutex> lck(mtx);

	int id = (int) m_streams.size();

	m_streams[id] = stream_t();
	m_streams[id].cur_id = 0;
	m_streams[id].stream_name = stream_name;
	m_streams[id].raw_size = 0;
	m_streams[id].packed_size = 0;
	m_streams[id].packed_data_size = 0;

	return id;
}

// ******************************************************************************
int CArchive::GetStreamId(string stream_name)
{
	lock_guard<mutex> lck(mtx);

	for (auto& x : m_streams)
		if (x.second.stream_name == stream_name)
			return static_cast<int>(x.first);

	return -1;
}

// ******************************************************************************
bool CArchive::AddPart(int stream_id, vector<uint8_t> &v_data, size_t metadata)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].parts.push_back(part_t(f_offset, v_data.size()));

	f_offset += write(metadata);
	fwrite(v_data.data(), 1, v_data.size(), f);

	f_offset += v_data.size();

	m_streams[stream_id].packed_size += f_offset - m_streams[stream_id].parts.back().offset;
	m_streams[stream_id].packed_data_size += v_data.size();

	return true;
}

// ******************************************************************************
int CArchive::AddPartPrepare(int stream_id)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].parts.push_back(part_t(0, 0));

	return static_cast<int>(m_streams[stream_id].parts.size()) - 1;
}

// ******************************************************************************
bool CArchive::AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data, size_t metadata)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());

	f_offset += write(metadata);
	fwrite(v_data.data(), 1, v_data.size(), f);

	f_offset += v_data.size();

	m_streams[stream_id].packed_size += f_offset - m_streams[stream_id].parts[part_id].offset;
	m_streams[stream_id].packed_data_size += v_data.size();

	return true;
}

// ******************************************************************************
void CArchive::SetRawSize(int stream_id, size_t raw_size)
{
	lock_guard<mutex> lck(mtx);
	
	m_streams[stream_id].raw_size = raw_size;
}

// ******************************************************************************
size_t CArchive::GetRawSize(int stream_id)
{
	lock_guard<mutex> lck(mtx);
	
	return m_streams[stream_id].raw_size;
}

// ******************************************************************************
bool CArchive::GetPart(int stream_id, vector<uint8_t> &v_data, size_t &metadata)
{
	lock_guard<mutex> lck(mtx);
	
	auto& p = m_streams[stream_id];

	if (p.cur_id >= p.parts.size())
		return false;

	v_data.resize(p.parts[p.cur_id].size);

	my_fseek(f, p.parts[p.cur_id].offset, SEEK_SET);

	if(p.parts[p.cur_id].size != 0)
		read(metadata);
	else
	{
		metadata = 0;
		p.cur_id++;
		return true;
	}

	auto r = fread(v_data.data(), 1, p.parts[p.cur_id].size, f);

	p.cur_id++;

	if (r != p.parts[p.cur_id-1].size)
		return false;

	return r == p.parts[p.cur_id-1].size;
}

// EOF
