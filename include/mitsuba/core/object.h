/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__OBJECT_H)
#define __OBJECT_H

#include <mitsuba/core/class.h>
#include <pthread.h>

MTS_NAMESPACE_BEGIN

/** \brief Base class of all Mitsuba classes.
 *
 * Contains functions relevant to every object
 * such as reference counting, limited type
 * introspection and lifetime management.
 *
 * @see ref, Class
 */
class MTS_EXPORT_CORE Object {
public:
	/// Object constructor
	Object();

	/// Get the reference count
	inline int getRefCount() const;

	/** \brief Increase the reference count of the
	 * object.
	 */
	void incRef() const;

	/** \brief Decrease the reference count object
	 * of the object. It will automatically be de-
	 * allocated once the reference count reaches
	 * zero.
	 */
	void decRef() const;

	/// Retrieve this object's class
	virtual const Class *getClass() const;

	/// Convert to string
	virtual std::string toString() const;
protected:
	/** \brief Virtual private deconstructor.
	 * (Will only be called by 'ref')
	 */
	virtual ~Object();
private:
	mutable int m_refCount;
#if defined(WIN32)
	mutable CRITICAL_SECTION m_refLock;
#else
	mutable pthread_mutex_t m_refLock;
#endif
public:
	static Class *m_theClass;
};

inline int Object::getRefCount() const {
	return m_refCount;
}

MTS_NAMESPACE_END

#endif /* __OBJECT_H */
