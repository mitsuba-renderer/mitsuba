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

#if !defined(__OBJECT_H)
#define __OBJECT_H

#include <mitsuba/core/class.h>
#include <pthread.h>

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/object.h mitsuba/mitsuba.h
 * \brief Parent of all Mitsuba classes.
 *
 * Contains functions relevant to every object such as reference counting, 
 * limited type introspection and lifetime management.
 *
 * \sa ref, Class
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Object {
public:
	/// Construct a new object 
	Object();

	/// Return the current reference count
	inline int getRefCount() const;

	/** \brief Increase the reference count of the
	 * object by one.
	 */
	void incRef() const;

	/** \brief Decrease the reference count of 
	 * the object and possibly deallocate it. 
	 *
	 * The object will automatically be deallocated once 
	 * the reference count reaches zero.
	 */
	void decRef() const;

	/// Retrieve this object's class
	virtual const Class *getClass() const;

	/**
	 * \brief Return a human-readable string representation
	 * of the object's contents.
	 *
	 * This function is mainly useful for debugging purposes
	 * and should ideally be implemented by all subclasses.
	 * The default implementation simply returns <tt>MyObject[unknown]</tt>,
	 * where <tt>MyObject</tt> is the name of the subclass.
	 */
	virtual std::string toString() const;
protected:
	/** \brief Virtual private deconstructor.
	 * (Will only be called by \ref ref)
	 */
	virtual ~Object();
public:
	static Class *m_theClass; ///< Pointer to the object's class descriptor
private:
#ifndef WIN32
	volatile mutable int m_refCount;
#else
	volatile mutable LONG m_refCount;
#endif
};

inline int Object::getRefCount() const {
	return m_refCount;
}

MTS_NAMESPACE_END

#endif /* __OBJECT_H */
