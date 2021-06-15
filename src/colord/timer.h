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
#ifndef __APPLE__
//mkokot TODO: consider valid implementation for mac os, currently do not use this classes on mac os
#include <memory>


class CThreadWatchImpl;
class CThreadWatch
{		
	std::unique_ptr<CThreadWatchImpl> pimpl;
public:
	CThreadWatch();
	void startTimer();
	void stopTimer();
	double getElapsedTime();
	~CThreadWatch();
};
#else
//Fake impl for mac os
class CThreadWatch
{
public:
	CThreadWatch() {}
	void startTimer() {}
	void stopTimer() {}
	double getElapsedTime() { return 0.0; }
	~CThreadWatch() {}
};
#endif
