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

#if !defined(__LOCK_H)
#define __LOCK_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Thin wrapper around the recursive pthreads lock
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE Mutex : public Object {
	friend class ConditionVariable;
public:
	/// Create a new mutex object
	Mutex();

	/// Lock the mutex
	inline void lock() {
		pthread_mutex_lock(&m_mutex);
	}

	/// Unlock the mutex
	inline void unlock() {
		pthread_mutex_unlock(&m_mutex);
	}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Mutex();
private:
	pthread_mutex_t m_mutex;
};

/**
 * \brief Wait flag synchronization primitive. Can be used to
 * wait for a certain event to occur.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE WaitFlag : public Object {
public:
	/**
	 * \brief Create a new wait flag
	 * \param flag
	 *    Initial state of the flag. If set to true, \ref wait()
	 *    will immediately return.
	 */
	WaitFlag(bool flag = false);

	/// Return the current flag value
	inline const bool &get() const { return m_flag; }

	/// Set the value of the flag
	void set(bool value);

	/// Wait for the flag to be set to true
	void wait();

	/**
	 * \brief Temporarily wait for the flag to be set to true
	 * 
	 * Similar to \ref wait(), but also uses a time value given
	 * in milliseconds. A return value of \c false signals
	 * that a timeout has occurred.
	 *
	 * \param ms Maximum waiting time in milliseconds
	 */
	bool wait(int ms);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~WaitFlag();
private:
	bool m_flag;
	pthread_mutex_t m_mutex;
	pthread_cond_t m_cond;
};

/**
 * \brief Condition variable synchronization primitive. Can be used to
 * wait for a condition to become true in a safe way.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE ConditionVariable : public Object {
public:
	/**
	 * \brief Create a new condition variable. Also takes a 
	 * mutex, which is later used by wait(). If none is specified,
	 * a new mutex instance will be created.
	 */
	ConditionVariable(Mutex *mutex = NULL);

	/**
	 * \brief Send a signal, which wakes up at least one of 
	 * the waiting threads. 
	 * 
	 * The calling thread does not have to hold the lock, 
	 * but more predictable scheduling will occur if this is the 
	 * case.
	 */
	inline void signal() {
		pthread_cond_signal(&m_cond);
	}

	/**
	 * \brief Send a signal, which wakes up any waiting threads. 
	 *
	 * The calling thread does not have to hold the lock, but more
	 * predictable scheduling will occur if this is the case.
	 */
	inline void broadcast() {
		pthread_cond_broadcast(&m_cond);
	}

	/** 
	 * \brief Wait for a signal and release the lock in the meanwhile.
	 *
	 * Assumes that the lock specified in the constructor has
	 * previously been acquired. After returning, the lock is
	 * held again.
	 */
	inline void wait() {
		pthread_cond_wait(&m_cond, &m_mutex->m_mutex);
	}

	/**
	 * \brief Temporarily wait for a signal and release the lock in the meanwhile.
	 *
	 * Similar to wait(), but also uses a time value given
	 * in milliseconds. A return value of \c false signals
	 * that a timeout has occurred. The lock is held after
	 * returning in either case.
	 *
	 * \param ms Maximum waiting time in milliseconds
	 */
	bool wait(int ms);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~ConditionVariable();
private:
	bool m_flag;
	ref<Mutex> m_mutex;
	pthread_cond_t m_cond;
};

MTS_NAMESPACE_END

#endif /* __LOCK_H */

