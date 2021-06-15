#ifndef __APPLE__
#include <cstdio> // NULL
#include "timer.h"

#ifdef _WIN32
#include <windows.h>

typedef struct {
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stopWatch;

typedef struct
{
	ULARGE_INTEGER start;
	ULARGE_INTEGER stop;
} thread_watch_t;

class CThreadWatchImpl
{
	thread_watch_t timer_kernel, timer_user;
	LARGE_INTEGER frequency;
	double LIToSecs(LARGE_INTEGER & L);

public:
	CThreadWatchImpl();
	void startTimer();
	void stopTimer();
	double getElapsedTime();
};

// **********************************************************
double CThreadWatchImpl::LIToSecs(LARGE_INTEGER& L)
{
	return ((double)L.QuadPart / (double)frequency.QuadPart);
}

// **********************************************************
CThreadWatchImpl::CThreadWatchImpl()
{
	timer_kernel.start.QuadPart = 0;
	timer_kernel.stop.QuadPart = 0;
	timer_user.start.QuadPart = 0;
	timer_user.stop.QuadPart = 0;
	//	QueryPerformanceFrequency( &frequency );
	//	frequency = 1;		// 100ns
}

// **********************************************************
void CThreadWatchImpl::startTimer()
{
	FILETIME CreationTime, ExitTime, KernelTime, UserTime;
	GetThreadTimes(GetCurrentThread(), &CreationTime, &ExitTime, &KernelTime, &UserTime);

	timer_kernel.start.LowPart = KernelTime.dwLowDateTime;
	timer_kernel.start.HighPart = KernelTime.dwHighDateTime;
	timer_user.start.LowPart = UserTime.dwLowDateTime;
	timer_user.start.HighPart = UserTime.dwHighDateTime;
}

// **********************************************************
void CThreadWatchImpl::stopTimer()
{
	FILETIME CreationTime, ExitTime, KernelTime, UserTime;
	GetThreadTimes(GetCurrentThread(), &CreationTime, &ExitTime, &KernelTime, &UserTime);
	//    QueryPerformanceCounter(&timer.stop);
	timer_kernel.stop.LowPart = KernelTime.dwLowDateTime;
	timer_kernel.stop.HighPart = KernelTime.dwHighDateTime;
	timer_user.stop.LowPart = UserTime.dwLowDateTime;
	timer_user.stop.HighPart = UserTime.dwHighDateTime;
}

// **********************************************************
double CThreadWatchImpl::getElapsedTime()
{
	/*	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
	return LIToSecs( time) ;*/
	LARGE_INTEGER time;

	time.QuadPart = (timer_kernel.stop.QuadPart - timer_kernel.start.QuadPart);
	time.QuadPart += (timer_user.stop.QuadPart - timer_user.start.QuadPart);

	return (double)time.QuadPart / 1e7;			// 100ns clock
}

#else
#include <sys/time.h>
#include <sys/resource.h>


typedef timeval thread_watch_t;

// **********************************************************
class CThreadWatchImpl
{
	thread_watch_t start_kernel, start_user;
	thread_watch_t stop_kernel, stop_user;

public:
	CThreadWatchImpl();
	void startTimer();
	void stopTimer();
	double getElapsedTime();
};
#include <cstring>



// **********************************************************
CThreadWatchImpl::CThreadWatchImpl()
{
}

// **********************************************************
void CThreadWatchImpl::startTimer()
{
	rusage usage;
	getrusage(RUSAGE_THREAD, &usage);
	start_user = usage.ru_utime;
	start_kernel = usage.ru_stime;
}

// **********************************************************
void CThreadWatchImpl::stopTimer()
{
	rusage usage;
	getrusage(RUSAGE_THREAD, &usage);
	stop_user = usage.ru_utime;
	stop_kernel = usage.ru_stime;
}

// **********************************************************
double CThreadWatchImpl::getElapsedTime()
{
	double ret = 0.0;

	ret += stop_user.tv_sec + stop_user.tv_usec / 1000000.0;
	ret += stop_kernel.tv_sec + stop_kernel.tv_usec / 1000000.0;
	ret -= start_user.tv_sec + start_user.tv_usec / 1000000.0;
	ret -= start_kernel.tv_sec + start_kernel.tv_usec / 1000000.0;

	return ret;
}


#endif


// **********************************************************
CThreadWatch::CThreadWatch()
{
	pimpl = std::make_unique<CThreadWatchImpl>();
}

// **********************************************************
void CThreadWatch::startTimer()
{
	pimpl->startTimer();
}

// **********************************************************
void CThreadWatch::stopTimer()
{
	pimpl->stopTimer();
}

// **********************************************************
double CThreadWatch::getElapsedTime()
{
	return pimpl->getElapsedTime();
}
// **********************************************************
CThreadWatch::~CThreadWatch() = default; //must be here due to unique_ptr based pimpl


#endif
// ***** EOF
