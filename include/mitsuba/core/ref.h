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
#if !defined(__MITSUBA_CORE_REF_H_)
#define __MITSUBA_CORE_REF_H_

#include "util.h"
#include <set>

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/ref.h mitsuba/mitsuba.h
 * \brief Reference counting helper
 *
 * The \a ref refeference template is a simple wrapper to store a
 * pointer to an object. It takes care of increasing and decreasing
 * the reference count of the object. When the last reference goes
 * out of scope, the associated object will be deallocated.
 *
 * \author Wenzel Jakob
 * \ingroup libcore
 */
template <typename T> class ref {
public:
	/// Create a NULL reference
	ref() : m_ptr(NULL) { }

	/// Construct a reference from a pointer
	ref(T *ptr) : m_ptr(ptr) { if (m_ptr) ((Object *) m_ptr)->incRef(); }

	/// Copy-constructor
	ref(const ref &pRef) : m_ptr(pRef.m_ptr) { if (m_ptr) ((Object *) m_ptr)->incRef(); }

	/// Destroy this reference
	~ref() { if (m_ptr) ((Object *) m_ptr)->decRef(); }

	/// Overwrite this reference with another reference
	inline ref& operator= (const ref& r) {
		if (m_ptr == r.m_ptr)
			return *this;
		if (m_ptr)
			((Object *) m_ptr)->decRef();
		m_ptr = r.m_ptr;
		if (m_ptr)
			((Object *) m_ptr)->incRef();
		return *this;
	}

	/// Overwrite this reference with a pointer to another object
	inline ref& operator= (T *ptr) {
		if (m_ptr == ptr)
			return *this;
		if (m_ptr)
			((Object *) m_ptr)->decRef();
		m_ptr = ptr;
		if (m_ptr)
			((Object *) m_ptr)->incRef();
		return *this;
	}

	/// Compare this reference with another reference
	inline bool operator== (const ref &pref) const { return (m_ptr == pref.m_ptr); }

	/// Compare this reference with another reference
	inline bool operator!= (const ref &pref) const { return (m_ptr != pref.m_ptr); }

	/// Compare this reference with a pointer
	inline bool operator== (const T* ptr) const { return (m_ptr == ptr); }

	/// Compare this reference with a pointer
	inline bool operator!= (const T* ptr) const { return (m_ptr != ptr); }

	/// Check whether this is a NULL reference
	inline bool operator!() const { return (m_ptr == NULL); }

	/// Access the object referenced by this reference
	inline T* operator-> () { return m_ptr; }

	/// Access the object referenced by this reference
	inline const T* operator-> () const { return m_ptr; }

	/// Return a C++ reference to the referenced object
	inline T& operator*() { return *m_ptr; }

	/// Return a C++ reference to the referenced object
	inline const T& operator*() const { return *m_ptr; }

	/// Return a pointer to the referenced object
	inline operator T* () { return m_ptr; }

	/// Return a pointer to the referenced object
	inline T* get() { return m_ptr; }

	/// Return a pointer to the referenced object
	inline const T* get() const { return m_ptr; }

	/**
	 * Return a string representation of this reference
	 */
	inline std::string toString() const {
		if (m_ptr == NULL)
			return formatString("ref<%s>[null]",
				T::m_theClass->getName().c_str());
		else
			return formatString("ref<%s>[ref=%i, ptr=%s]",
					m_ptr->getClass()->getName().c_str(),
					m_ptr->getRefCount(), m_ptr->toString().c_str());
	}
private:
	T *m_ptr;
};

/// \cond

/// Comparison operator for references
template <typename T> struct ref_comparator {
	bool operator() (const ref<T>& lhs, const ref<T>& rhs) const {
		return lhs.get() < rhs.get();
	}
};

/// \endcond

/**
 * \brief Simple reference-counted vector container based on \c std::vector and \ref ref
 */
template <typename T> class ref_vector : public std::vector< ref<T> > {
public:
	typedef std::vector< ref<T> > parent_type;

	ref_vector() : parent_type() {}
	ref_vector(size_t size) : parent_type(size) {}
	ref_vector(const ref_vector &vec) : parent_type(vec) {}

	/// Remove all duplicates without changing the order
	inline void ensureUnique() {
		std::set<T *> seen;

		typename parent_type::iterator it1 = this->begin(), it2 = this->begin();
		for (; it1 < this->end(); ++it1) {
			if (seen.find(it1->get()) != seen.end())
				continue;
			seen.insert(it1->get());
			if (it1 != it2)
				*it2++ = *it1;
			else
				it2++;
		}
		this->erase(it2, this->end());
	}

	/// Check if a certain pointer is contained in the vector
	inline bool contains(const T *ptr) const {
		for (typename parent_type::const_iterator it = this->begin();
				it != this->end(); ++it)
			if (it->get() == ptr)
				return true;
		return false;
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_REF_H_ */
