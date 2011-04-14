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

#if !defined(__TLS_H)
#define __TLS_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/tls.h mitsuba/mitsuba.h
 * \brief Thin wrapper around posix thread local storage. 
 *     Stores references to Object instances.
 *
 * This class implements a simple reference-counting thread-local
 * pointer storage.
 * \sa PrimitiveThreadLocal
 * \ingroup libcore
 */
template <typename T> class ThreadLocal {
public:
	/// Create a new (empty) thread-local reference
	inline ThreadLocal() {
		pthread_key_create(&m_key, ThreadLocal::deletePtr);
	}

	/**
	 * \brief Destroy the thread-local reference
	 *
	 * Beware: all threads that ever used this TLS (other than
	 * the caller of the destructor) must either either died or 
	 * called \ref cleanup() prior to the destructor invocation -- 
	 * otherwise, the reference counts will be incorrect and
	 * memory leaks will ensue..
	 */
	inline ~ThreadLocal() {
		cleanup();
		pthread_key_delete(m_key);
	}

	/**
	 * \brief Make the thread-local reference point to the given 
	 * instance and update the reference counts
	 */
	inline void set(T *value) {
		T *previous = get();
		if (previous != NULL)
			previous->decRef();
		if (value != NULL)
			value->incRef();
		pthread_setspecific(m_key, value);
	}

	/// Return the currently stored reference 
	inline T *get() {
		return static_cast<T *>(pthread_getspecific(m_key));
	}
	
	/// Return the currently stored reference (const version)
	inline const T *get() const {
		return static_cast<const T *>(pthread_getspecific(m_key));
	}

	/// Set the reference to <tt>NULL</tt> and update the reference counts.
	void cleanup() {
		set(NULL);
	}
private:
	static void deletePtr(void *ptr) {
		static_cast<T *>(ptr)->decRef();
	}
	
	pthread_key_t m_key;
};

/**
 * \headerfile mitsuba/core/tls.h mitsuba/mitsuba.h
 * \brief Thin wrapper around posix thread local storage. 
 *     Stores heap-allocated data other than Object instances.
 * \sa ThreadLocal
 * \ingroup libcore
 */
template <typename T> class PrimitiveThreadLocal {
public:
	/// Create a new thread local storage object
	inline PrimitiveThreadLocal() {
		pthread_key_create(&m_key, PrimitiveThreadLocal::deletePtr);
	}

	/**
	 * \brief Destroy the thread local storage object
	 *
	 * Beware: all threads that ever used this TLS (other than
	 * the caller of the destructor) must either either died or 
	 * called \ref cleanup() prior to the destructor invocation -- 
	 * otherwise, memory leaks will ensue..
	 */
	inline ~PrimitiveThreadLocal() {
		cleanup();
		pthread_key_delete(m_key);
	}

	/// Update the data associated with the current thread
	inline void set(T &value) {
		get() = value;
	}

	/// Return a reference to the data associated with the current thread
	inline T &get() {
		T *value = static_cast<T *>(pthread_getspecific(m_key));
		if (value == NULL) {
			value = new T();
			pthread_setspecific(m_key, value);
		}

		return *value;
	}

	/**
	 * \brief Return a reference to the data associated with the 
	 * current thread (const version)
	 */
	inline const T &get() const {
		T *value = static_cast<T *>(pthread_getspecific(m_key));
		if (value == NULL) {
			value = new T();
			pthread_setspecific(m_key, value);
		}

		return *value;
	}

	/// Delete the data record associated with the current thread
	void cleanup() {
		T *value = static_cast<T *>(pthread_getspecific(m_key));
		if (value != NULL)
			delete value;
		pthread_setspecific(m_key, NULL);
	}
private:
	static void deletePtr(void *ptr) {
		delete static_cast<T *>(ptr);
	}

	pthread_key_t m_key;
};

MTS_NAMESPACE_END

#endif /* __TLS_H */
