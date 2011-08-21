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

#if !defined(__DEBUG_H)
#define __DEBUG_H

#include <mitsuba/core/timer.h>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Statistics collection
// -----------------------------------------------------------------------

/// Size (in characters) of the console-based progress message
#define PROGRESS_MSG_SIZE 56

/**
 * Specifies the number of internal counters associated with each 
 * \ref StatsCounter instance.
 *
 * This is needed for SMP/ccNUMA systems where different processors might 
 * be contending for a cache line containing a counter. The solution used
 * here tries to ensure that every processor has  its own local counter.
 */
#define NUM_COUNTERS       128   // Must be a power of 2

/// Bitmask for \ref NUM_COUNTERS
#define NUM_COUNTERS_MASK (NUM_COUNTERS-1)

/// Determines the multiples (e.g. 1000, 1024) and units of a \ref StatsCounter
enum EStatsType {
	ENumberValue = 0, ///< Simple unitless number, e.g. # of rays
	EByteCount,       ///< Number of read/written/transferred bytes
	EPercentage,      ///< Percentage with respect to a base counter
	EAverage          ///< Average value with respect to a base counter
};

/**
 * \brief Counter data structure, which is suitable for ccNUMA/SMP machines
 *
 * This counter takes up at least one cache line to reduce false sharing.
 */
struct CacheLineCounter { 
#if (defined(WIN32) && !defined(WIN64)) || (defined(__POWERPC__) && !defined(_LP64))
	// WIN32 & Darwin (PPC/32) doesn't support atomic 64 bit increment operations
	// -> restrict counters to 32bit :(
	uint32_t value;
	uint32_t unused2;
#else
	uint64_t value;
#endif
	char unused[120];
};

/** \brief General-purpose statistics counter
 * 
 * This class implements a simple counter, which can be used to track various
 * quantities within Mitsuba. At various points during the execution, it is
 * possible to then call \ref Statistics::printStats() to get a human-readable 
 * report of their values.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE StatsCounter {
public:
	/** 
	 * \brief Create a new statistics counter 
	 * 
	 * \param category Category of the counter when shown in the statistics summary
	 * \param name     Name of the counter when shown in the statistics summary
	 * \param type     Characterization of the quantity that will be measured
	 * \param initial  Initial value of the counter
	 * \param base     Initial value of the base counter (only for <tt>type == EPercentage</tt> and <tt>EAverage</tt>)
	 */
	StatsCounter(const std::string &category, const std::string &name,
		EStatsType type = ENumberValue, uint64_t initial = 0L, uint64_t base = 0L);

	/// Free all storage used by the counter
	~StatsCounter();

	/// Increment the counter value by one
#if defined(MTS_NO_STATISTICS)
	inline uint64_t operator++() { return 0L; }
#elif defined(WIN64)
	inline uint64_t operator++() {
		const int offset = Thread::getID() & NUM_COUNTERS_MASK;
		InterlockedExchangeAdd64(reinterpret_cast<LONG64 *>(&m_value[offset].value), 1);
		return m_value[offset].value;
}
#elif defined(WIN32)
	inline uint32_t operator++() {
		const int offset = Thread::getID() & NUM_COUNTERS_MASK;
		InterlockedExchangeAdd(reinterpret_cast<LONG *>(&m_value[offset].value), 1);
		return m_value[offset].value;
}
#elif defined(__POWERPC__) && !defined(_LP64)
	inline uint64_t operator++() {
		return (uint64_t) __sync_fetch_and_add(&m_value[Thread::getID() & NUM_COUNTERS_MASK].value, 1);
	}
#else
	inline uint64_t operator++() {
		return __sync_fetch_and_add(&m_value[Thread::getID() & NUM_COUNTERS_MASK].value, 1);
	}
#endif

	/// Increment the counter by the specified amount
#ifdef MTS_NO_STATISTICS
	inline void operator+=(size_t amount) { }
#elif defined(WIN64)
	inline void operator+=(size_t amount) {
		InterlockedExchangeAdd64(reinterpret_cast<LONG64 *>(&m_value[Thread::getID() & NUM_COUNTERS_MASK].value), amount);
	}
#elif defined(WIN32)
	inline void operator+=(size_t amount) {
		InterlockedExchangeAdd(reinterpret_cast<LONG *>(&m_value[Thread::getID() & NUM_COUNTERS_MASK].value), amount);
	}
#else
	inline void operator+=(size_t amount) {
		__sync_fetch_and_add(&m_value[Thread::getID() & NUM_COUNTERS_MASK].value, amount);
	}
#endif
	
	/// Increment the base counter by the specified amount (only for use with EPercentage/EAverage)
#ifdef MTS_NO_STATISTICS
	inline void incrementBase(size_t amount = 1) { }
#elif defined(WIN64)
	inline void incrementBase(size_t amount = 1) {
		InterlockedExchangeAdd64(reinterpret_cast<LONG64 *>(&m_base[Thread::getID() & NUM_COUNTERS_MASK].value), amount);
	}
#elif defined(WIN32)
	inline void incrementBase(size_t amount = 1) {
		InterlockedExchangeAdd(reinterpret_cast<LONG *>(&m_base[Thread::getID() & NUM_COUNTERS_MASK].value), amount);
	}
#else
	inline void incrementBase(size_t amount = 1) {
		__sync_fetch_and_add(&m_base[Thread::getID() & NUM_COUNTERS_MASK].value, amount);
	}
#endif

	/// Return the name of this counter
	inline const std::string &getName() const { return m_name; }

	/// Return the category of this counter
	inline const std::string &getCategory() const { return m_category; }

	/// Return the type of this counter
	inline EStatsType getType() const { return m_type; }

	/// Return the value of this counter as 64-bit unsigned integer
#ifdef MTS_NO_STATISTICS
	inline uint64_t getValue() const { return 0L; }
#else
	inline uint64_t getValue() const {
		uint64_t result = 0;
		for (int i=0; i<NUM_COUNTERS; ++i)
			result += m_value[i].value;
		return result;
	}
#endif

	/// Get the reference number (only used with the EPercentage/EAverage counter type)
#ifdef MTS_NO_STATISTICS
	inline uint64_t getBase() const { return 0L; }
#else
	inline uint64_t getBase() const {
		uint64_t result = 0;
		for (int i=0; i<NUM_COUNTERS; ++i)
			result += m_base[i].value;
		return result;
	}
#endif

	/// Reset the stored counter values
	inline void reset() {
		for (int i=0; i<NUM_COUNTERS; ++i) {
			m_value[i].value = m_base[i].value = 0;
		}
	}

	/// Sorting by name (for the statistics)
	bool operator<(const StatsCounter &v) const;
private:
	std::string m_category;
	std::string m_name;
	EStatsType m_type;
	CacheLineCounter *m_value;
	CacheLineCounter *m_base;
};

/** \brief General-purpose progress reporter
 *
 * This class is used to track the progress of various operations that might
 * take longer than, say, a second. It provides interactive feedback on both
 * the console binaries and within Mitsuba's Qt-based GUI.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE ProgressReporter {
public:
	/**
	 * Construct a new progress reporter. 
	 * 'ptr' is a custom pointer payload to be submitted with progress messages
	 */
	ProgressReporter(const std::string &title, long long total, const void *ptr);

	/// Reset the progress reporter to its initial state
	void reset();

	/// Update the progress display
	void update(long long value);

	/// Return the current value
	inline long long getValue() const { return m_value; }

	/// Finish
	inline void finish() { 
		if (m_value < m_total)
			update(m_total);
	}

	/// Enable/disable progress bars
	static void setEnabled(bool enabled);

	/// Check whether progress bars are enabled
	static inline bool isEnabled() { return m_enabled; }
private:
	/// Convert a time value to a human-readable format
	void printTime(Float time, std::ostream &os) const;
private:
	static bool m_enabled;
	std::string m_title;
	long long m_total, m_value;
	unsigned int m_lastMs;
	int m_percentage;
	int m_fillSize, m_fillPos;
	char m_string[PROGRESS_MSG_SIZE];
	ref<Timer> m_timer;
	const void *m_ptr;
};

/** \brief Collects various rendering statistics and presents them
 * in a human-readable form.
 *
 * \remark Only the \ref getInstance(), \ref getStats(), and 
 * \ref printStats() functions are implemented in the Python bindings.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Statistics : public Object {
public:
	/// Return the global stats collector instance
	inline static Statistics *getInstance() { return m_instance; }

	/// Register a counter with the statistics collector
	void registerCounter(const StatsCounter *ctr);

	/// Record that a plugin has been loaded
	void logPlugin(const std::string &pname, const std::string &descr);

	/// Print a summary of gathered statistics
	void printStats();

	/// Return a string containing gathered statistics
	std::string getStats();

	/// Initialize the global statistics collector
	static void staticInitialization();

	/// Free the memory taken by staticInitialization()
	static void staticShutdown();

	MTS_DECLARE_CLASS()
protected:
	/// Create a statistics instance
	Statistics();
	/// Virtual destructor
	virtual ~Statistics() { }
private:
	struct compareCategory {
		bool operator()(const StatsCounter *c1, const StatsCounter *c2) {
			if (c1->getCategory() == c2->getCategory())
				return c1->getName() <= c2->getName();
			return c1->getCategory() < c2->getCategory();
		}
	};
	 
	static ref<Statistics> m_instance;
	std::vector<const StatsCounter *> m_counters;
	std::vector<std::pair<std::string, std::string> > m_plugins;
	ref<Mutex> m_mutex;
};

MTS_NAMESPACE_END

#endif /* __DEBUG_H */
