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

#pragma once
#if !defined(__MITSUBA_CORE_LOCK_H_)
#define __MITSUBA_CORE_LOCK_H_

#include <mitsuba/mitsuba.h>

#include <boost/scoped_ptr.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Thin wrapper around the recursive boost thread lock
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE Mutex : public Object {
    friend class ConditionVariable;
public:
    /// Create a new mutex object
    Mutex();

    /// Lock the mutex
    void lock();

    /// Unlock the mutex
    void unlock();

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Mutex();
private:
    struct MutexPrivate;
    boost::scoped_ptr<MutexPrivate> d;
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
    const bool &get() const;

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
    struct WaitFlagPrivate;
    boost::scoped_ptr<WaitFlagPrivate> d;
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
    void signal();

    /**
     * \brief Send a signal, which wakes up any waiting threads.
     *
     * The calling thread does not have to hold the lock, but more
     * predictable scheduling will occur if this is the case.
     */
    void broadcast();

    /**
     * \brief Wait for a signal and release the lock in the meanwhile.
     *
     * Assumes that the lock specified in the constructor has
     * previously been acquired. After returning, the lock is
     * held again.
     */
    void wait();

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
    struct ConditionVariablePrivate;
    boost::scoped_ptr<ConditionVariablePrivate> d;
};

/**
 * \brief Simple RAII-style locking of a Mutex. On construction it locks the
 * mutex and unlocks it on destruction. Based on boost::lock_guard,
 * assumes the Mutex will outlive the lock.
 *
 * \ingroup libcore
 */
class LockGuard {
public:
    explicit LockGuard(Mutex * m_) : m(m_) {
        m->lock();
    }

    ~LockGuard() {
        m->unlock();
    }
private:
    Mutex* const m;

    explicit LockGuard(LockGuard&);
    LockGuard& operator=(LockGuard&);
};

/**
 * \brief In addition to providing RAII-style locking, UniqueLock also allows
 * for deferred locking until lock() is called explicitly. unlock() is only
 * called by the destructor if the object has locked the mutex.
 * Based on boost::unique_lock, assumes the Mutex will outlive the lock.
 */
class UniqueLock {
public:
    explicit UniqueLock(Mutex * mutex, bool acquire_lock = true)
    : m(mutex), is_locked(false) {
        if (acquire_lock)
            lock();
    }

    ~UniqueLock() {
        if (ownsLock())
            m->unlock();
    }

    void lock() {
        SAssert(!ownsLock() && m != NULL);
        m->lock();
        is_locked = true;
    }

    void unlock() {
        SAssert(ownsLock() && m != NULL);
        m->unlock();
        is_locked = false;
    }

    Mutex * release() {
        Mutex * const mutex = m;
        m = static_cast<Mutex*>(NULL);
        is_locked = false;
        return mutex;
    }

    inline bool operator!() const {
        return !ownsLock();
    }

    inline bool ownsLock() const {
        return is_locked;
    }

private:
    Mutex* m;
    bool is_locked;

    explicit UniqueLock(UniqueLock&);
    UniqueLock& operator=(UniqueLock&);
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_LOCK_H_ */
