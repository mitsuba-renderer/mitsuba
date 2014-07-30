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
#if !defined(__MITSUBA_CORE_OBJECT_H_)
#define __MITSUBA_CORE_OBJECT_H_

#include <mitsuba/core/class.h>

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
	void decRef(bool autoDeallocate = true) const;

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

	/** \brief Initializes the built-in reference count
	 * debugger (if enabled)
	 */
	static void staticInitialization();

	/// Free the memory taken by staticInitialization()
	static void staticShutdown();
protected:
	/** \brief Virtual private deconstructor.
	 * (Will only be called by \ref ref)
	 */
	virtual ~Object();
public:
	static Class *m_theClass; ///< Pointer to the object's class descriptor
private:
#if !defined(_MSC_VER)
	volatile mutable int m_refCount;
#else
	volatile mutable long m_refCount;
#endif
};

inline int Object::getRefCount() const {
	return m_refCount;
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_OBJECT_H_ */
