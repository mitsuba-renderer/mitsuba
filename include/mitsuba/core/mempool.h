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
#if !defined(__MITSUBA_CORE_MEMPOOL_H_)
#define __MITSUBA_CORE_MEMPOOL_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/// Create a new memory pool with an initial set of 128 entries
#define MTS_MEMPOOL_GRANULARITY 128

/// Set this to one to track down memory pool-related issues
#define MTS_DEBUG_MEMPOOL 0

/**
 * \brief Basic memory pool for efficient allocation and deallocation
 * of objects of the same type.
 *
 * This class attempts to keep most instances contiguous in memory, while
 * having only minimal interaction with the underlying allocator.
 *
 * \ingroup libcore
 */
template <typename T> class BasicMemoryPool {
public:
	/// Create a new memory pool with an initial set of 128 entries
	BasicMemoryPool(size_t nEntries = MTS_MEMPOOL_GRANULARITY) : m_size(0) {
		increaseCapacity(nEntries);
	}

	/// Destruct the memory pool and release all entries
	~BasicMemoryPool() {
		for (size_t i=0; i<m_cleanup.size(); ++i)
			freeAligned(m_cleanup[i]);
	}

	/// Acquire an entry
	inline T *alloc() {
		if (EXPECT_NOT_TAKEN(m_free.empty()))
			increaseCapacity();
		T *result = m_free.back();
		m_free.pop_back();
		return result;
	}

	void assertNotContained(T *ptr) {
#if MTS_DEBUG_MEMPOOL == 1
		if (std::find(m_free.begin(), m_free.end(), ptr) != m_free.end())
			SLog(EError, "BasicMemoryPool:assertNotContained(): Memory pool inconsistency!");
#endif
	}

	/// Release an entry
	inline void release(T *ptr) {
#if MTS_DEBUG_MEMPOOL == 1
		if (std::find(m_free.begin(), m_free.end(), ptr) != m_free.end())
			SLog(EError, "BasicMemoryPool::release(): Memory pool "
				"inconsistency. Tried to release %s", ptr->toString().c_str());
#endif
		m_free.push_back(ptr);
	}

	/// Return the total size of the memory pool
	inline size_t size() const {
		return m_size;
	}

	/// Check if every entry has been released
	bool unused() const {
		return m_free.size() == m_size;
	}

	/// Return a human-readable description
	std::string toString() const {
		std::ostringstream oss;
		oss << "BasicMemoryPool[size=" << m_size << ", free=" << m_free.size() << "]";
		return oss.str();
	}
private:
	void increaseCapacity(size_t nEntries = MTS_MEMPOOL_GRANULARITY) {
		T *ptr = static_cast<T *>(allocAligned(sizeof(T) * nEntries));
		for (size_t i=0; i<nEntries; ++i)
			m_free.push_back(&ptr[i]);
		m_cleanup.push_back(ptr);
		m_size += nEntries;
	}
private:
	std::vector<T *> m_free;
	std::vector<T *> m_cleanup;
	size_t m_size;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_MEMPOOL_H_ */
