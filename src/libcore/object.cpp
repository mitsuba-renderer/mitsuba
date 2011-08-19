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

#include <mitsuba/mitsuba.h>

#define DEBUG_REFCOUNTS 0

MTS_NAMESPACE_BEGIN

Class *MTS_CLASS(Object) = new Class("Object", false, "");

Object::Object()
 : m_refCount(0) {
}

void Object::incRef() const {
#if DEBUG_REFCOUNTS == 1
	if (Class::rttiIsInitialized())
		cout << this << ": Increasing reference count (" << getClass()->getName() << ") -> "
			<< m_refCount + 1 << endl;
#endif
#if defined(_WIN32)
	InterlockedIncrement(&m_refCount);
#else
	__sync_fetch_and_add(&m_refCount, 1);
#endif
}

void Object::decRef() const {
#if DEBUG_REFCOUNTS == 1
	if (Class::rttiIsInitialized()) {
		cout << this << ": Decreasing reference count (" << getClass()->getName() << ") -> "
			<< m_refCount - 1 << endl;
	}
#endif
#if defined(_WIN32)
	int count = InterlockedDecrement(&m_refCount);
#else
	int count = __sync_sub_and_fetch(&m_refCount, 1);
#endif
	AssertEx(count >= 0, "Reference count is below zero!");
	if (count == 0) {
#if DEBUG_REFCOUNTS == 1
		if (Class::rttiIsInitialized())
			cout << this << ": Deleting an instance of " << 
				getClass()->getName() << endl;
#endif
		delete this;
	}
}

const Class *Object::getClass() const {
    return m_theClass;
}

std::string Object::toString() const {
	std::ostringstream oss;
	oss << getClass()->getName();
	oss << "[unknown]";
	return oss.str();
}

Object::~Object() {
	if (m_refCount > 0) {
		Log(EWarn, "Deleting %s with reference count %i!", 
			toString().c_str(), m_refCount);
	}
}

MTS_NAMESPACE_END
