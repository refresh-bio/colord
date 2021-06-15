#pragma once
#include "defs.h"
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <string>
#include <ostream>
#include <map>
#include <ctime>

#ifdef MONITOR_QUEUES

class CQueueMonitor
{
	std::ostream& ostr;
	bool single_line;
	bool reporting;
	std::mutex mtx;

	struct desc_t
	{
		std::string name;
		size_t max_size;
		size_t cur_size;
	};

	std::map<uint32_t, desc_t> monitored_queues;

	void report()
	{
		if (!reporting)
			return;

		std::string sep = single_line ? "\t" : "\n";
		std::string rep;

		if (!single_line)
			rep = "**********\n";

		rep += std::to_string(std::time(nullptr)) + sep;

		for (auto& x : monitored_queues)
		{
			rep += std::to_string(x.first) + " ";
			rep += x.second.name + " ";
			rep += std::to_string(x.second.cur_size) + " ";
			rep += std::to_string(x.second.max_size) + sep;
		}

		if (!rep.empty())
			rep.back() = '\n';

		if (!single_line)
			rep += "**********\n";

		ostr << rep;
	}

public:
	CQueueMonitor(std::ostream &ostr, bool single_line = false, bool _reporting = true) : ostr(ostr), single_line(single_line), reporting(_reporting)
	{}

	void register_queue(uint32_t id, std::string name, size_t max_size)
	{
		std::lock_guard<std::mutex> lck(mtx);

		monitored_queues[id] = desc_t{ name, max_size, 0u };
	}

	void inc(uint32_t id, size_t val = 1)
	{
		std::lock_guard<std::mutex> lck(mtx);

		monitored_queues[id].cur_size += val;
		report();
	}

	void dec(uint32_t id, size_t val = 1)
	{
		std::lock_guard<std::mutex> lck(mtx);

		monitored_queues[id].cur_size -= val;
		report();
	}
};

#else
class CQueueMonitor
{
	void report()
	{		
	}

public:
	CQueueMonitor(std::ostream& ostr, bool single_line = false, bool _reporting = true) 
	{}

	void register_queue(uint32_t id, std::string name, size_t max_size)
	{
		
	}

	void inc(uint32_t id, size_t val = 1)
	{
		
	}

	void dec(uint32_t id, size_t val = 1)
	{
		
	}
};
#endif


//implements thread safe circular queue
template<typename T>
class CParallelQueue
{
	std::vector<T> data;
	bool full = false;
	bool is_completed = false;
	uint32_t n_writers;
	uint32_t start = 0;
	uint32_t end = 0;
	std::mutex mtx;
	std::condition_variable cv_push;
	std::condition_variable cv_pop;

	CQueueMonitor *qm;
	uint32_t q_id;

public:	
	CParallelQueue(uint32_t size, uint32_t n_writers = 1, CQueueMonitor* qm = nullptr, uint32_t q_id = 0) :
		data(size),
		n_writers(n_writers),
		qm(qm),
		q_id(q_id)
	{
	}

	void Push(T&& elem)
	{
		std::unique_lock lck(mtx);

		cv_push.wait(lck, [this] {return !full; });
		bool was_empty = start == end;
		data[end] = std::move(elem);
		end = (end + 1ul) % data.size();
		
		full = end == start;
		if (was_empty)
//			cv_pop.notify_all();
			cv_pop.notify_one();

		if (qm)
			qm->inc(q_id);
	}

	bool Pop(T& elem)
	{
		std::unique_lock lck(mtx);
		cv_pop.wait(lck, [this] {
			return start != end || full || is_completed;
		});
		if (is_completed && !full && start == end)
			return false;

		bool was_full = full;
		elem = std::move(data[start]);
		start = (start + 1ul) % data.size();
		full = false;
		if (was_full)
//			cv_push.notify_all();
			cv_push.notify_one();

		if(qm)
			qm->dec(q_id);

		return true;		
	}

	void MarkCompleted()
	{
		std::lock_guard lck(mtx);
		--n_writers;
		if (!n_writers)
		{
			is_completed = true;
			cv_pop.notify_all();
		}
	}
};




//implements thread safe circular queue
/*
It is assumed that there is exactly one producer who can steal idle threads
*/
template<typename T>
class CParallelQueuePopWaiting
{
	std::vector<T> data;
	bool full = false;
	bool is_completed = false;	
	uint32_t start = 0;
	uint32_t end = 0;
	std::mutex mtx;
	std::condition_variable cv_push;
	std::condition_variable cv_pop;
	uint32_t n_pop_waiting{};

	CQueueMonitor* qm;
	uint32_t q_id;

public:
	CParallelQueuePopWaiting(uint32_t size, CQueueMonitor *qm = nullptr, uint32_t q_id = 0) :
		data(size),
		qm(qm),
		q_id(q_id)
	{

	}
	uint32_t GetNWaitingOnPop() const
	{
		return n_pop_waiting;
	}
	void Push(T&& elem)
	{
		std::unique_lock lck(mtx);

		cv_push.wait(lck, [this] {return !full; });
		bool was_empty = start == end;
		data[end] = std::move(elem);
		end = (end + 1) % data.size();

		full = end == start;
		if (was_empty)
//			cv_pop.notify_all();
			cv_pop.notify_one();

		if (qm)
			qm->inc(q_id);
	}

	bool Pop(T& elem)
	{		
		std::unique_lock lck(mtx);
		bool first = true;
		cv_pop.wait(lck, [this, &first] {
			if (first)
				++n_pop_waiting;
			first = false;
			return start != end || full || is_completed;
		});
		--n_pop_waiting;
		if (is_completed && !full && start == end)
			return false;

		bool was_full = full;
		elem = std::move(data[start]);
		start = (start + 1) % data.size();
		full = false;
		if (was_full)
//			cv_push.notify_all();
			cv_push.notify_one();

		if (qm)
			qm->dec(q_id);

		return true;
	}

	void MarkCompleted()
	{
		std::lock_guard lck(mtx);		
		is_completed = true;
		cv_pop.notify_all();		
	}
};

//implements thread safe circular queue
template<typename T>
class CParallelPriorityQueue
{
/*	template<typename U>
	struct Elem
	{
		uint32_t priority;
		U data;
		bool operator<(const Elem& rhs) const
		{
			return priority > rhs.priority;
		}
	};*/
//	std::priority_queue<Elem<T>> data;

	std::map<uint32_t, T> map_data;

	uint32_t max_size;
	bool is_completed = false;
	uint32_t n_writers;
	
	std::mutex mtx;
	std::condition_variable cv_push;
	std::condition_variable cv_pop;
	uint32_t current_priority{};

	CQueueMonitor* qm;
	uint32_t q_id;

public:
	CParallelPriorityQueue(uint32_t size, uint32_t n_writers = 1, CQueueMonitor *qm = nullptr, uint32_t q_id = 0) :
		max_size(size),
		n_writers(n_writers),
		qm(qm),
		q_id(q_id)
	{

	}
	void Push(uint32_t priority, T&& elem)
	{
		std::unique_lock lck(mtx);		
		cv_push.wait(lck, [this, &priority] {
			return map_data.size() < max_size - 1 ||
				priority == current_priority;
		});
		//bool was_empty = data.empty();
#ifdef ENABLE_QUEUE_LOGING
		std::cerr << "Push pack with priority: " << priority << "\n";
#endif
//		data.emplace(Elem<T>{priority, std::move(elem)});
		map_data.emplace(priority, std::move(elem));

		if(priority == current_priority)		
			cv_pop.notify_all();

		cv_push.notify_all();

		if (qm)
			qm->inc(q_id);
	}

	bool Pop(T& elem)
	{
		//std::cerr << "Try Pop, expected priority: " << current_priority << "\n";
		std::unique_lock lck(mtx);
		cv_pop.wait(lck, [this] {
//			return is_completed || (!data.empty() && data.top().priority == current_priority);
			return is_completed || (!map_data.empty() && map_data.begin()->first == current_priority);
		});
//		if (is_completed && data.empty())
		if (is_completed && map_data.empty())
			return false;

		//bool was_full = data.size() == max_size;
//		elem = std::move(data.top().data);
//		data.pop();

//		elem = std::move(set_data.begin()->data);
		elem.swap(map_data.begin()->second);
//		elem = move(map_data.begin()->second);
		map_data.erase(map_data.begin());
#ifdef ENABLE_QUEUE_LOGING
		std::cerr << "Pop pack with priority: " << current_priority << "\n";
#endif
		++current_priority;

		cv_push.notify_all();
		
		if (qm)
			qm->dec(q_id);

		return true;
	}

	void MarkCompleted()
	{
		std::lock_guard lck(mtx);
		--n_writers;
		if (!n_writers)
		{
			is_completed = true;
			cv_pop.notify_all();
		}
	}
};
