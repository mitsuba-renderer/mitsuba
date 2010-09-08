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
#if defined(WIN32)
	InitializeCriticalSection(&m_refLock);
#else
	pthread_mutex_init(&m_refLock, NULL);
#endif
}

void Object::incRef() const {
#if defined(WIN32)
	EnterCriticalSection(&m_refLock);
	m_refCount++;
	LeaveCriticalSection(&m_refLock);
#else
	pthread_mutex_lock(&m_refLock);
	m_refCount++;
	pthread_mutex_unlock(&m_refLock);
#endif
}

void Object::decRef() const {
#if defined(WIN32)
	EnterCriticalSection(&m_refLock);
	int count = --m_refCount;
	LeaveCriticalSection(&m_refLock);
#else
	pthread_mutex_lock(&m_refLock);
	int count = --m_refCount;
	pthread_mutex_unlock(&m_refLock);
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
#if defined(WIN32)
	DeleteCriticalSection(&m_refLock);
#else
	pthread_mutex_destroy(&m_refLock);
#endif
}

MTS_NAMESPACE_END
