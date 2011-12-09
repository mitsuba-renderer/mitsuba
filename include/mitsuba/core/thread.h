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

#if !defined(__THREAD_H)
#define __THREAD_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/thread.h mitsuba/mitsuba.h
 * \brief Cross-platform thread implementation
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Thread : public Object {
public:
	/// Possible priority values for \ref Thread::setPriority()
	enum EThreadPriority {
		EIdlePriority = 0,
		ELowestPriority,
		ELowPriority,
		ENormalPriority,
		EHighPriority,
		EHighestPriority,
		ERealtimePriority
	};

	/**
	 * \brief Create a new thread object
	 * \param name An identifying name of this thread 
	 *   (will be shown in debug messages)
	 * \param stackSize Initial stack size of the thread
	 *   (0 = default)
	 * \remark Note that it is currently not possible to
	 *         construct Thread instances from Python
	 */
	Thread(const std::string &name, 
		unsigned int stackSize = 0);

	/**
	 * \brief Set the thread priority
	 *
	 * This does not always work -- for instance, Linux 
	 * requires root privileges for this operation.
	 *
	 * \return \c true upon success.
	 */
	bool setPriority(EThreadPriority priority);

	/// Return the thread priority
	inline EThreadPriority getPriority() const { return m_priority; }

	/**
	 * \brief Specify whether or not this thread is critical
	 * 
	 * When an thread marked critical crashes from an uncaught 
	 * exception, the whole process is brought down. 
	 * The default is \c false.
	 */
	inline void setCritical(bool critical) { m_critical = critical; }

	/// Return the value of the critical flag
	inline bool getCritical() const { return m_critical; }

	/// Return the thread's stack size
	inline int getStackSize() const { return m_stackSize; }

	/// Return the thread ID
#if defined(WIN32)
	inline static int getID() { return (int) GetCurrentThreadId(); }
#elif MTS_USE_ELF_TLS == 1
	inline static int getID() { return m_id; }
#else
	inline static int getID() { 
		return getThread()->m_id;
	}
#endif

	/// Return the name of this thread
	inline const std::string &getName() const { return m_name; }

	/// Set the name of this thread
	inline void setName(const std::string &name) { m_name = name; }

	/// Return the parent thread
	inline Thread *getParent() { return m_parent; }

	/// Return the parent thread (const version)
	inline const Thread *getParent() const { return m_parent.get(); }

	/// Set the logger instance used to process log messages from this thread
	inline void setLogger(Logger *logger) { m_logger = logger; }

	/// Return the thread's logger instance
	inline Logger *getLogger() { return m_logger; }

	/// Set the thread's file resolver
	inline void setFileResolver(FileResolver *fresolver) { m_fresolver = fresolver; }

	/// Return the thread's file resolver
	inline FileResolver *getFileResolver() { return m_fresolver; }

	/// Return the current thread
	inline static Thread *getThread() { return m_self->get(); }

	/// Is this thread still running?
	inline bool isRunning() const { return m_running; }

	/// Start the thread
	void start();

	/**
	 * \brief Detach the thread and release resources
	 * 
	 * After a call to this function, \ref join()
	 * cannot be used anymore. This releases resources, which
	 * would otherwise be held until a call to \ref join().
	 */
	void detach();

	/// Wait until the thread finishes
	void join();

	/// Return a string representation
	virtual std::string toString() const;

	/// Sleep for a certain amount of time
	static void sleep(unsigned int ms);

	/// Initialize the threading system
	static void staticInitialization();

	/// Shut down the threading system
	static void staticShutdown();

	/// Initialize Mitsuba's threading system for simultaneous use of OpenMP
	static void initializeOpenMP(size_t threadCount);

	/**
	 * \brief Register an unmanaged thread with Mitsuba (i.e. one that 
	 * doesn't derive from \c mitsuba::Thread)
	 *
	 * Should be called from the thread in question. The function returns
	 * a Mitsuba handle to the thread
	 */
	static Thread *registerUnmanagedThread(const std::string &name);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Thread();

	/// Thread dispatch function
	static void *dispatch(void *par);

	/**
	 * Exit the thread, should be called from
	 * inside the thread
	 */
	void exit();

	/// Yield to another processor
	void yield();

	/// The thread's run method
	virtual void run() = 0;
private:
	ref<Thread> m_parent;
	ref<Logger> m_logger;
	ref<FileResolver> m_fresolver;
	ref<Mutex> m_joinMutex;
	std::string m_name;
	unsigned int m_stackSize;
	bool m_running, m_joined;
	EThreadPriority m_priority;
	pthread_t m_thread;
	static ThreadLocal<Thread> *m_self;
	bool m_critical;
#if defined(__LINUX__) || defined(__OSX__)
	static int m_idCounter;
	static ref<Mutex> m_idMutex;
#if MTS_USE_ELF_TLS == 1
	static __thread int m_id
		__attribute__((tls_model("global-dynamic")));
#else
	int m_id;
#endif
#endif
};

/// Like omp_get_thread_num(), but works on all platforms
extern MTS_EXPORT_CORE int mts_get_thread_num();
/// Like omp_get_max_threads(), but works on all platforms
extern MTS_EXPORT_CORE int mts_get_max_threads();

MTS_NAMESPACE_END

#endif /* __THREAD_H */
