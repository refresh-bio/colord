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

#include <functional>
namespace pooled_threads
{
	class ThreadPoolTask;

	class thread
	{
		bool _joinable = false;
		ThreadPoolTask* task;
		void Create(std::function<void()>&& f);
	public:
		thread(const thread&) = delete;
		thread& operator=(const thread&) = delete;

		thread(thread&& rhs) :
			_joinable(rhs._joinable), task(rhs.task)
		{
			rhs._joinable = false;
			rhs.task = nullptr;
		}
		thread& operator=(thread&& rhs)
		{
			_joinable = rhs._joinable;
			task = rhs.task;
			rhs._joinable = false;
			rhs.task = nullptr;
			return *this;
		}

		template<typename _Callable, typename... _Args>
		explicit thread(_Callable&& __f, _Args&&... __args) :
			_joinable(true)
		{
			Create(std::bind(std::forward<_Callable>(__f), std::forward<_Args>(__args)...));
		}

		bool joinable()
		{
			return _joinable;
		}

		void join();

		~thread();
	};
}