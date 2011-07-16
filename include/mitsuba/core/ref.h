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

#if !defined(__REFERENCE_H)
#define __REFERENCE_H

#include "util.h"

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
	inline ref& operator= (const ref& pref) {
		if (m_ptr == pref.m_ptr)
			return *this;
		T* tmp = m_ptr;
		m_ptr = pref.m_ptr;
		if (m_ptr)
			((Object *) m_ptr)->incRef();
		if (tmp)
			((Object *) tmp)->decRef();
		return *this;
	}
	
	/// Overwrite this reference with a pointer to another object
	inline ref& operator= (T *ptr) {
		if (m_ptr == ptr)
			return *this;
		T* tmp = m_ptr;
		m_ptr = ptr;
		if (m_ptr)
			((Object *) m_ptr)->incRef();
		if (tmp)
			((Object *) tmp)->decRef();
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

MTS_NAMESPACE_END

#endif /* __REFERENCE_H */
