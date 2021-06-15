/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

The source codes are based on codes written by Dennis and published:
http://allmybrain.com/2008/06/10/timing-cc-code-on-linux/

Version: 3.1.1
Date   : 2019-05-19
*/

#ifndef _TIMER_H
#define _TIMER_H

#ifdef WIN32
#include <windows.h>

typedef struct {
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stopWatch;

class CStopWatch {

private:
	stopWatch timer;
	LARGE_INTEGER frequency;
	double LIToSecs(LARGE_INTEGER & L);
public:
	CStopWatch();
	void startTimer();
	void stopTimer();
	double getElapsedTime();
};

typedef struct
{
	ULARGE_INTEGER start;
	ULARGE_INTEGER stop;
} thread_watch_t;


#else
#include <sys/time.h>
#include <sys/resource.h>
typedef struct {
	timeval start;
	timeval stop;
} stopWatch;

class CStopWatch {

private:
	stopWatch timer;
public:
	CStopWatch();
	void startTimer();
	void stopTimer();
	double getElapsedTime();
};


typedef timeval thread_watch_t;
#ifndef __APPLE__
// **********************************************************
#endif

#endif

#endif
// ***** EOF
