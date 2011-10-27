/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/lock.h>
#include <errno.h>

#if !defined(WIN32)
#include <sys/time.h>
#else
#include <sys/timeb.h>

static int gettimeofday(struct timeval *tp, void *) {
	struct _timeb timebuf;
	_ftime(&timebuf);
	tp->tv_sec = (long) timebuf.time;
	tp->tv_usec = (long) timebuf.millitm * 1000;
	return 0;
}
#endif

MTS_NAMESPACE_BEGIN

Mutex::Mutex() {
	pthread_mutexattr_t mta;
	pthread_mutexattr_init(&mta);
	pthread_mutexattr_settype(&mta, PTHREAD_MUTEX_RECURSIVE);
	if (pthread_mutex_init(&m_mutex, &mta) != 0) 
		Log(EError, "Could not create a mutex");
	pthread_mutexattr_destroy(&mta);
}
	
Mutex::~Mutex() {
	int rv = pthread_mutex_destroy(&m_mutex);
	if (rv != 0) {
		switch (rv) {
			case EBUSY: Log(EWarn, "Could not destroy a mutex: held by another thread"); break;
			case EINVAL: Log(EError, "Could not destroy a mutex: data corruption"); break;
			default: Log(EError, "Could not destroy a mutex: unknown reason (%i)", rv);
		}
	}
}

ConditionVariable::ConditionVariable(Mutex *mutex) {
	m_mutex = mutex != NULL ? mutex : new Mutex();
	pthread_cond_init(&m_cond, NULL);
}

ConditionVariable::~ConditionVariable() {
	pthread_cond_destroy(&m_cond);
}

bool ConditionVariable::wait(int ms) {
	if (ms == -1) {
		wait();
		return true;
	} else {
		struct timeval tv;
		gettimeofday(&tv, NULL);

		struct timespec ts;
		ts.tv_sec = tv.tv_sec + ms/1000;
		ts.tv_nsec = tv.tv_usec * 1000 + (ms % 1000) * 1000000;
		if (ts.tv_nsec >= 1000000000) {
			ts.tv_nsec -= 1000000000;
			ts.tv_sec++;
		}
		int retval = pthread_cond_timedwait(&m_cond, &m_mutex->m_mutex, &ts);
#if defined(WIN32)
		/* Should be outside of the switch statement (depending
		   on the used environment, the constants ETIMEDOUT and WSAETIMEDOUT
		   are potentially identical */   
		if (retval == WSAETIMEDOUT)
			return false;
#endif
		switch (retval) {
			case 0: return true;
			case ETIMEDOUT: return false;
			case EINVAL: Log(EError, "Invalid condition variable/mutex or time interval (%i)!", ms);
			case EPERM: Log(EError, "Mutex is not owned by the current thread");
			default: Log(EError, "Unknown return value (%i) in ConditionVariable::wait(%i)!", retval, ms);
				return false;
		}
	}
}

WaitFlag::WaitFlag(bool flag) 
	: m_flag(flag) {
	pthread_mutex_init(&m_mutex, NULL);
	pthread_cond_init(&m_cond, NULL);
}

WaitFlag::~WaitFlag() {
	pthread_mutex_destroy(&m_mutex);
	pthread_cond_destroy(&m_cond);
}

void WaitFlag::set(bool value) {
	pthread_mutex_lock(&m_mutex);
	m_flag = value;
	if (m_flag) 
		pthread_cond_broadcast(&m_cond);
	pthread_mutex_unlock(&m_mutex);
}

void WaitFlag::wait() {
	pthread_mutex_lock(&m_mutex);
	/* Wait for a signal from the CV and
		release the mutex while waiting */
	while (!m_flag) 
		pthread_cond_wait(&m_cond, &m_mutex);
	pthread_mutex_unlock(&m_mutex);
}

bool WaitFlag::wait(int ms) {
	pthread_mutex_lock(&m_mutex);
	if (!m_flag) {
		struct timespec ts;
		struct timeval tv;
		gettimeofday(&tv, NULL);
		ts.tv_sec = tv.tv_sec + ms/1000;
		ts.tv_nsec = (tv.tv_usec + ms % 1000) * 1000;
		if (pthread_cond_timedwait(&m_cond, &m_mutex, &ts) != 0) {
			pthread_mutex_unlock(&m_mutex);
			return false;
		}
	}
	pthread_mutex_unlock(&m_mutex);
	return true;
}

MTS_IMPLEMENT_CLASS(ConditionVariable, false, Object)
MTS_IMPLEMENT_CLASS(WaitFlag, false, Object)
MTS_IMPLEMENT_CLASS(Mutex, false, Object)
MTS_NAMESPACE_END
