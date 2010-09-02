/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__FIFO_H)
#define __FIFO_H

#include <mitsuba/mitsuba.h>
#include <pthread.h>
#include <deque>

#if defined(__OSX__)
#include <mach/semaphore.h>
#include <mach/mach.h>
#include <mach/task.h>
#else
#include <semaphore.h>
#endif

MTS_NAMESPACE_BEGIN

/**
 * Simple synchronized FIFO class with an underlying 
 * deque container
 */
template <class T> class FIFO {
public:
	FIFO() {
#if defined(__OSX__)
		if (semaphore_create(::mach_task_self(), &m_available,
			SYNC_POLICY_FIFO, 0) != KERN_SUCCESS)
			SLog(EError, "Could not create a semaphore!");
#else
		if (sem_init(&m_available, PTHREAD_PROCESS_PRIVATE, 0))
			SLog(EError, "Could not create a semaphore!");
#endif
		pthread_mutex_init(&m_mutex, NULL);
	}

	void put(T value) {
		pthread_mutex_lock(&m_mutex);
		m_data.push_back(value);
		pthread_mutex_unlock(&m_mutex);
#if defined(__OSX__)
		semaphore_signal(m_available);
#else
		sem_post(&m_available);
#endif
	}

	T get() {
		T result;
#if defined(__OSX__)
		semaphore_wait(m_available);
#else
		while (sem_wait(&m_available) != 0)
			;
#endif
		pthread_mutex_lock(&m_mutex);
		result = m_data.front();
		m_data.pop_front();
		pthread_mutex_unlock(&m_mutex);
		return result;
	}

	int getSize() const {
		int result;
		pthread_mutex_lock(&m_mutex);
		result = m_data.size();
		pthread_mutex_unlock(&m_mutex);
		return result;
	}

	~FIFO() {
#if defined(__OSX__)
		semaphore_destroy(::mach_task_self(), m_available);
#else
		sem_destroy(&m_available);
#endif
		pthread_mutex_destroy(&m_mutex);
	}
private:
	std::deque<T> m_data;
	mutable pthread_mutex_t m_mutex;
#if defined(__OSX__)
	semaphore_t m_available;
#else
	sem_t m_available;
#endif
};

MTS_NAMESPACE_END

#endif /* __FIFO_H */
