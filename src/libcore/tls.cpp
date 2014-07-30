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

#include <mitsuba/core/tls.h>

#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/unordered_set.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>

#if defined(__OSX__)
# include <pthread.h>
#endif

MTS_NAMESPACE_BEGIN

namespace mi = boost::multi_index;

/* The native TLS classes on Linux/MacOS/Windows only support a limited number
   of dynamically allocated entries (usually 1024 or 1088). Furthermore, they
   do not provide appropriate cleanup semantics when the TLS object or one of
   the assocated threads dies. The custom TLS code provided in Mitsuba has no
   such limits (caching in various subsystems of Mitsuba may create a huge amount,
   so this is a big deal) as well as nice cleanup semantics. The implementation
   is designed to make the \c get() operation as fast as as possible at the cost
   of more involved locking when creating or destroying threads and TLS objects */
namespace detail {

/// A single TLS entry + cleanup hook
struct TLSEntry {
	void *data;
	void (*destructFunctor)(void *);

	inline TLSEntry() : data(NULL), destructFunctor(NULL) { }
};

/// boost multi-index element to act as replacement of map<Key,T>
template<typename T1, typename T2>
struct mutable_pair {
	mutable_pair(const T1 &f, const T2 &s) : first(f), second(s) { }

	T1 first;
	mutable T2 second;
};

/// Per-thread TLS entry map
struct PerThreadData {
	typedef mutable_pair<void *, TLSEntry> MapData;
	typedef mi::member<MapData, void *, &MapData::first> key_member;
	struct seq_tag {};
	struct key_tag {};

	typedef mi::multi_index_container<MapData,
		mi::indexed_by<
			mi::hashed_unique<mi::tag<key_tag>, key_member>,
			mi::sequenced<mi::tag<seq_tag> >
		>
	> Map;
	typedef mi::index<Map, key_tag>::type::iterator         key_iterator;
	typedef mi::index<Map, seq_tag>::type::reverse_iterator reverse_iterator;

	Map map;
	boost::recursive_mutex mutex;
};

/// List of all PerThreadData data structures (one for each thread)
boost::unordered_set<PerThreadData *> ptdGlobal;
/// Lock to protect ptdGlobal
boost::mutex ptdGlobalLock;

#if defined(__WINDOWS__)
__declspec(thread) PerThreadData *ptdLocal = NULL;
#elif defined(__LINUX__)
__thread PerThreadData *ptdLocal = NULL;
#elif defined(__OSX__)
pthread_key_t ptdLocal;
#endif

struct ThreadLocalBase::ThreadLocalPrivate {
	ConstructFunctor constructFunctor;
	DestructFunctor destructFunctor;

	ThreadLocalPrivate(const ConstructFunctor &constructFunctor,
			const DestructFunctor &destructFunctor) : constructFunctor(constructFunctor),
			destructFunctor(destructFunctor) { }

	~ThreadLocalPrivate() {
		/* The TLS object was destroyed. Walk through all threads
		   and clean up where necessary */
		boost::lock_guard<boost::mutex> guard(ptdGlobalLock);

		for (boost::unordered_set<PerThreadData *>::iterator it = ptdGlobal.begin();
				it != ptdGlobal.end(); ++it) {
			PerThreadData *ptd = *it;
			boost::unique_lock<boost::recursive_mutex> lock(ptd->mutex);

			PerThreadData::Map::iterator it2 = ptd->map.find(this);
			TLSEntry entry;

			if (it2 != ptd->map.end()) {
				entry = it2->second;
				ptd->map.erase(it2);
			}

			lock.unlock();

			if (entry.data)
				destructFunctor(entry.data);
		}
	}

	/// Look up a TLS entry. The goal is to make this operation very fast!
	std::pair<void *, bool> get() {
		bool existed = true;
		void *data;

#if defined(__OSX__)
		PerThreadData *ptd = (PerThreadData *) pthread_getspecific(ptdLocal);
#else
		PerThreadData *ptd = ptdLocal;
#endif
		if (EXPECT_NOT_TAKEN(!ptd))
			throw std::runtime_error("Internal error: call to ThreadLocalPrivate::get() "
				" precedes the construction of thread-specific data structures!");

		/* This is an uncontended thread-local lock (i.e. not to worry) */
		boost::lock_guard<boost::recursive_mutex> guard(ptd->mutex);
		PerThreadData::key_iterator it = ptd->map.find(this);
		if (EXPECT_TAKEN(it != ptd->map.end())) {
			data = it->second.data;
		} else {
			/* This is the first access from this thread */
			TLSEntry entry;
			entry.data = data = constructFunctor();
			entry.destructFunctor = destructFunctor;
			ptd->map.insert(PerThreadData::MapData(this, entry));
			existed = false;
		}

		return std::make_pair(data, existed);
	}
};

ThreadLocalBase::ThreadLocalBase(
		const ConstructFunctor &constructFunctor, const DestructFunctor &destructFunctor)
		: d(new ThreadLocalPrivate(constructFunctor, destructFunctor)) { }

ThreadLocalBase::~ThreadLocalBase() { }

void *ThreadLocalBase::get() {
	return d->get().first;
}

const void *ThreadLocalBase::get() const {
	return d->get().first;
}

void *ThreadLocalBase::get(bool &existed) {
	std::pair<void *, bool> result = d->get();
	existed = result.second;
	return result.first;
}

const void *ThreadLocalBase::get(bool &existed) const {
	std::pair<void *, bool> result = d->get();
	existed = result.second;
	return result.first;
}

void initializeGlobalTLS() {
#if defined(__OSX__)
	pthread_key_create(&ptdLocal, NULL);
#endif
}

void destroyGlobalTLS() {
#if defined(__OSX__)
	pthread_key_delete(ptdLocal);
	memset(&ptdLocal, 0, sizeof(pthread_key_t));
#endif
}

/// A new thread was started -- set up TLS data structures
void initializeLocalTLS() {
	boost::lock_guard<boost::mutex> guard(ptdGlobalLock);
#if defined(__OSX__)
	PerThreadData *ptd = (PerThreadData *) pthread_getspecific(ptdLocal);
	if (!ptd) {
		ptd = new PerThreadData();
		ptdGlobal.insert(ptd);
		pthread_setspecific(ptdLocal, ptd);
	}
#else
	if (!ptdLocal) {
		ptdLocal = new PerThreadData();
		ptdGlobal.insert(ptdLocal);
	}
#endif
}

/// A thread has died -- destroy any remaining TLS entries associated with it
void destroyLocalTLS() {
	boost::lock_guard<boost::mutex> guard(ptdGlobalLock);

#if defined(__OSX__)
	PerThreadData *ptd = (PerThreadData *) pthread_getspecific(ptdLocal);
#else
	PerThreadData *ptd = ptdLocal;
#endif

	boost::unique_lock<boost::recursive_mutex> lock(ptd->mutex);

	// Destroy the data in reverse order of creation
	for (PerThreadData::reverse_iterator it = mi::get<PerThreadData::seq_tag>(ptd->map).rbegin();
		it != mi::get<PerThreadData::seq_tag>(ptd->map).rend(); ++it) {
		TLSEntry &entry = it->second;
		entry.destructFunctor(entry.data);
	}

	lock.unlock();
	ptdGlobal.erase(ptd);
	delete ptd;
#if defined(__OSX__)
	pthread_setspecific(ptdLocal, NULL);
#else
	ptdLocal = NULL;
#endif
}

} /* namespace detail */


MTS_NAMESPACE_END
