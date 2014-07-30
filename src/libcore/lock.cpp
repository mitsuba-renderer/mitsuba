/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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

#include <boost/thread/thread_time.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/condition_variable.hpp>

MTS_NAMESPACE_BEGIN

struct Mutex::MutexPrivate {
	boost::recursive_timed_mutex mutex;
};

Mutex::Mutex() : d(new MutexPrivate) {
}

Mutex::~Mutex() {
}

void Mutex::lock() {
	d->mutex.lock();
}

void Mutex::unlock() {
	d->mutex.unlock();
}

struct ConditionVariable::ConditionVariablePrivate {
	ref<Mutex> mutex;
	boost::condition_variable_any cond;

	ConditionVariablePrivate(Mutex * m) : mutex(m) { }

	// Helper to get the actual mutex implementation. This works only if
	// the mutex is not null

	inline boost::recursive_timed_mutex& mutexImpl() {
		return mutex->d->mutex;
	}
};

ConditionVariable::ConditionVariable(Mutex *mutex)
	: d(new ConditionVariablePrivate(mutex != NULL ? mutex : new Mutex())) {

}

ConditionVariable::~ConditionVariable() {
}

void ConditionVariable::signal() {
	d->cond.notify_one();
}

void ConditionVariable::broadcast() {
	d->cond.notify_all();
}

void ConditionVariable::wait() {
	d->cond.wait(d->mutexImpl());
}

bool ConditionVariable::wait(int ms) {
	if (ms == -1) {
		wait();
		return true;
	} else {
		const boost::posix_time::ptime timeout =
			boost::get_system_time() + boost::posix_time::milliseconds(ms);
		return d->cond.timed_wait(d->mutexImpl(), timeout);
	}
}

struct WaitFlag::WaitFlagPrivate {
	bool flag;
	boost::timed_mutex mutex;
	boost::condition_variable_any cond;

	WaitFlagPrivate(bool f) : flag(f) {}
};

WaitFlag::WaitFlag(bool flag)
	: d(new WaitFlagPrivate(flag)) {
}

WaitFlag::~WaitFlag() {
}

const bool & WaitFlag::get() const {
	return d->flag;
}

void WaitFlag::set(bool value) {
	boost::timed_mutex::scoped_lock lock(d->mutex);
	d->flag = value;
	if (d->flag)
		d->cond.notify_all();
}

void WaitFlag::wait() {
	boost::timed_mutex::scoped_lock lock(d->mutex);
	// Wait for a signal from the CV and release the mutex while waiting
	while (!d->flag)
		d->cond.wait(d->mutex);
}

bool WaitFlag::wait(int ms) {
	boost::timed_mutex::scoped_lock lock(d->mutex);
	if (!d->flag) {
		const boost::posix_time::ptime timeout =
			boost::get_system_time() + boost::posix_time::milliseconds(ms);
		return d->cond.timed_wait(lock, timeout);
	}
	return true;
}

MTS_IMPLEMENT_CLASS(ConditionVariable, false, Object)
MTS_IMPLEMENT_CLASS(WaitFlag, false, Object)
MTS_IMPLEMENT_CLASS(Mutex, false, Object)
MTS_NAMESPACE_END
