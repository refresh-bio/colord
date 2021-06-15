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
