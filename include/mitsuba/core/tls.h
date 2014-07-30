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
#if !defined(__MITSUBA_CORE_TLS_H_)
#define __MITSUBA_CORE_TLS_H_

#include <mitsuba/mitsuba.h>
#include <boost/scoped_ptr.hpp>

MTS_NAMESPACE_BEGIN

///! \cond
namespace detail {
	/* Some forward declarations pertaining to TLS management */

	extern MTS_EXPORT_CORE void initializeGlobalTLS();
	extern MTS_EXPORT_CORE void destroyGlobalTLS();
	extern MTS_EXPORT_CORE void initializeLocalTLS();
	extern MTS_EXPORT_CORE void destroyLocalTLS();

	class MTS_EXPORT_CORE ThreadLocalBase {
	public:
		/// Functor to allocate memory for a TLS object
		typedef void *(*ConstructFunctor)();
		/// Functor to release memory of a TLS object
		typedef void (*DestructFunctor)(void *);
		/// Construct a new thread local storage object
		ThreadLocalBase(const ConstructFunctor &constructFunctor,
			const DestructFunctor &destructFfunctor);
		/// Destroy the thread local storage object
		~ThreadLocalBase();
		/// Return the data value associated with the current thread
		void *get();
		/// Return the data value associated with the current thread (const version)
		const void *get() const;
		/// Like the other \c get(), but also returns whether the TLS object existed before
		void *get(bool &existed);
		/// Like the other \c get(), but also returns whether the TLS object existed before (const version)
		const void *get(bool &existed) const;
	protected:
		struct ThreadLocalPrivate;
		mutable boost::scoped_ptr<ThreadLocalPrivate> d;
	};
}; // namespace tls
///! \endcond

/**
 * \headerfile mitsuba/core/tls.h mitsuba/mitsuba.h
 * \brief Thin wrapper around boost thread local storage.
 *     Stores references to Object instances.
 * \sa PrimitiveThreadLocal
 * \ingroup libcore
 *
 * This class implements a reference counting thread local storage object which captures
 * references to subclasses of \ref Object. In comparison to an API like <tt>boost::thread_specific_ptr</tt>
 * it has a much nicer cleanup mechanism. Held references are destroyed when the owning thread dies \a or
 * when the \c ThreadLocal instance is freed, whichever occurs first.
 */
template <typename ValueType> class ThreadLocal {
public:
	/// Construct a new thread local storage object
	ThreadLocal() : m_base(&ThreadLocal::construct, &ThreadLocal::destruct) { }

	/// Update the data associated with the current thread
	inline void set(ValueType *ptr) { ((ref<ValueType> *) m_base.get())->operator=(ptr); }

	/// Return a reference to the data associated with the current thread
	inline ValueType *get() { return ((ref<ValueType> *) m_base.get())->get(); }

	/**
	 * \brief Return a reference to the data associated with the
	 * current thread (const version)
	 */
	inline const ValueType *get() const { return ((const ref<ValueType> *) m_base.get())->get(); }
protected:
	inline static void *construct() {
		return new ref<ValueType>();
	}

	inline static void destruct(void *data) {
		delete static_cast<ref<ValueType> *>(data);
	}
protected:
	detail::ThreadLocalBase m_base;
};

/**
 * \headerfile mitsuba/core/tls.h mitsuba/mitsuba.h
 * \brief Thin wrapper around posix thread local storage.
 *     Stores heap-allocated data other than Object instances.
 * \sa ThreadLocal
 * \ingroup libcore
 *
 * This class implements a thread local storage object for POD-style data structures.
 * In comparison to an API like <tt>boost::thread_specific_ptr</tt> it has a much nicer
 * cleanup mechanism. Held references are destroyed when the owning thread dies \a or
 * when the \c PrimitiveThreadLocal instance is freed, whichever occurs first.
 */
template <typename ValueType> class PrimitiveThreadLocal {
public:
	/// Construct a new thread local storage object
	PrimitiveThreadLocal() : m_base(&PrimitiveThreadLocal::construct,
		&PrimitiveThreadLocal::destruct) { }

	/// Update the data associated with the current thread
	inline void set(ValueType &value) {
		get() = value;
	}

	/// Return a reference to the data associated with the current thread
	inline ValueType &get() {
		return *((ValueType *) m_base.get());
	}

	/**
	 * \brief Return a reference to the data associated with the
	 * current thread (const version)
	 */
	inline const ValueType &get() const {
		return *((const ValueType *) m_base.get());
	}
protected:
	inline static void *construct() {
		return new ValueType();
	}

	inline static void destruct(void *data) {
		if (data)
			delete static_cast<ValueType *>(data);
	}
protected:
	detail::ThreadLocalBase m_base;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_TLS_H_ */
