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
#include "colord_api.h"
#include <iostream>

void usage(char* prog_name)
{
	std::cerr << "Usage: " << prog_name << " <colord-archive> [reference_genome]\n";
}
int main(int argc, char** argv)
{
	if (argc < 2)
	{
		usage(argv[0]);
		return 1;
	}
	std::string inputPath;
	std::string refGenomePath;
	inputPath = argv[1];
	if (argc > 2)
		refGenomePath = argv[2];
	try
	{
		colord::DecompressionStream stream(inputPath, refGenomePath);
		auto info = stream.GetInfo();
		std::cerr << "Database info:\n\n";
		info.ToOstream(std::cerr);

		while (auto x = stream.NextRecord())
		{
			if (info.isFastq)
			{
				std::cout << "@" << x.ReadHeader() << "\n";
				std::cout << x.Read() << "\n";
				std::cout << "+" << x.QualHeader() << "\n";
				std::cout << x.Qual() << "\n";
			}
			else
			{
				std::cout << ">" << x.ReadHeader() << "\n";
				std::cout << x.Read() << "\n";
			}
		}
	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
		return 1;
	}	
	return 0;
}