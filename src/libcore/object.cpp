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

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

Class *Object::m_theClass = new Class("Object", false, "");

Object::Object()
 : m_refCount(0) {
}

void Object::incRef() const {
#if defined(_WIN32)
	InterlockedIncrement(&m_refCount);
#else
	__sync_fetch_and_add(&m_refCount, 1);
#endif
}

void Object::decRef() const {
#if defined(_WIN32)
	int count = InterlockedDecrement(&m_refCount);
#else
	int count = __sync_sub_and_fetch(&m_refCount, 1);
#endif
	AssertEx(count >= 0, "Reference count is below zero!");
	if (count == 0) 
		delete this;
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
