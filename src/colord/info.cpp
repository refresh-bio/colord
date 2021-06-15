#include "info.h"
#include "utils.h"
#include "archive.h"
#include <ctime>

void runInfo(const CInfoParams& params)
{
	
	CArchive archive(true);
	if (!archive.Open(params.inputFilePath))
	{
		std::cerr << "Error: cannot open archive: " << params.inputFilePath << "\n";
		exit(1);
	}
	
	CInfo info;
	int s_info = archive.GetStreamId("info");
	size_t meta;
	std::vector<uint8_t> _info;
	archive.GetPart(s_info, _info, meta);
	info.Deserialize(_info);
	archive.Close();

	std::cerr << "version major: " << info.version_major << "\n";
	std::cerr << "version minor: " << info.version_minor << "\n";

	std::cerr << "total bytes: " << info.total_bytes << "\n";
	std::cerr << "total bases: " << info.total_bases << "\n";
	std::cerr << "total reads: " << info.total_reads << "\n";
	time_t time = info.time;
	std::cerr << "time: " << asctime(localtime(&time)) << "\n";
	std::cerr << "command: " << info.full_command_line<< "\n";
}

