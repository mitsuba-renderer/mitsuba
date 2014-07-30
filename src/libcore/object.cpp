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

#include <mitsuba/mitsuba.h>
#include <set>

#if defined(_MSC_VER)
# include <intrin.h>
#endif

#define DEBUG_REFCOUNTS 0

#if DEBUG_REFCOUNTS == 1
# include <boost/thread/mutex.hpp>
#endif

MTS_NAMESPACE_BEGIN

#if DEBUG_REFCOUNTS == 1

namespace
{
class RefCountTracker {
	typedef boost::lock_guard<boost::mutex> Lock;
public:
	RefCountTracker() { }

	inline void add(const Object *obj) {
		Lock lock(m_mutex);
		m_objects.insert(obj);
	}

	inline void remove(const Object *obj) {
		Lock lock(m_mutex);
		m_objects.erase(obj);
	}

	~RefCountTracker() {
		Lock lock(m_mutex);
		if (m_objects.size() == 0) {
			cout << "Reference count debugger -- all is good! (no leaked objects)" << endl;
		} else {
			cout << endl;
			cout << "Reference count debugger: leaked " << m_objects.size() << " objects!" << endl;
			cout << "A detailed summary follows:" << endl << endl;

			for (std::set<const Object *>::iterator it = m_objects.begin();
					it != m_objects.end(); ++it) {
				cout << "=======================================================" << endl;
				cout << " Leaked object at address " << (*it) << ":" << endl << endl;
				cout << " Class type: " << (*it)->getClass()->getName() << endl;
				cout << "=======================================================" << endl;
			}
		}
	}
private:
	std::set<const Object *> m_objects;
	boost::mutex m_mutex;
};

static RefCountTracker *__ref_tracker = NULL;

} // namespace
#endif

Class *MTS_CLASS(Object) = new Class("Object", false, "");

Object::Object()
 : m_refCount(0) {
#if DEBUG_REFCOUNTS == 1
	if (__ref_tracker)
		__ref_tracker->add(this);
#endif
}

void Object::incRef() const {
#if DEBUG_REFCOUNTS == 1
	if (Class::rttiIsInitialized())
		cout << this << ": Increasing reference count (" << getClass()->getName() << ") -> "
			<< (int) (m_refCount + 1) << endl;
#endif
#if defined(_MSC_VER)
	_InterlockedIncrement(&m_refCount);
#else
	__sync_fetch_and_add(&m_refCount, 1);
#endif
}

void Object::decRef(bool autoDeallocate) const {
#if DEBUG_REFCOUNTS == 1
	if (Class::rttiIsInitialized()) {
		cout << this << ": Decreasing reference count (" << getClass()->getName() << ") -> "
			<< std::dec << (int) (m_refCount - 1) << endl;
	}
#endif
#if defined(_MSC_VER)
	int count = _InterlockedDecrement(&m_refCount);
#else
	int count = __sync_sub_and_fetch(&m_refCount, 1);
#endif
	AssertEx(count >= 0, "Reference count is below zero!");
	if (count == 0 && autoDeallocate) {
#if DEBUG_REFCOUNTS == 1
		if (Class::rttiIsInitialized())
			cout << this << ": Deleting an instance of " <<
				getClass()->getName() << endl;
		if (__ref_tracker)
			__ref_tracker->remove(this);
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
			toString().c_str(), (int) m_refCount);
	}
}

void Object::staticInitialization() {
#if DEBUG_REFCOUNTS == 1
	if (!__ref_tracker)
		__ref_tracker = new RefCountTracker();
#endif
}

void Object::staticShutdown() {
#if DEBUG_REFCOUNTS == 1
	if (__ref_tracker) {
		delete __ref_tracker;
		__ref_tracker = NULL;
	}
#endif
}

MTS_NAMESPACE_END
