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

#if !defined(__THREAD_H)
#define __THREAD_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_CORE Thread : public Object {
public:
	/// Available thread priorities
	enum EThreadPriority {
		EIdlePriority = 0,
		ELowestPriority,
		ELowPriority,
		ENormalPriority,
		EHighPriority,
		EHighestPriority,
		ERealtimePriority
	};

	/// Create a new thread object
	Thread(const std::string &name, 
		unsigned int stackSize = 0);

	/// Set the thread priority
	void setPriority(EThreadPriority priority);

	/**
	 * Set the critical flag. When an thread marked critical crashes
	 * from an uncaught exception, the whole process is terminated
	 * (default: false).
	 */
	inline void setCritical(bool critical) { m_critical = critical; }

	/// Return the value of the critical flag
	inline bool getCritical() const { return m_critical; }

	/// Return the thread priority
	inline EThreadPriority getPriority() const { return m_priority; }

	/// Return the thread's stack size
	inline int getStackSize() const { return m_stackSize; }

	/// Return the thread's ID
#if defined(__OSX__)
	inline static int getID() { return getThread()->m_id; }
#elif defined(WIN32)
	inline static int getID() { return (int) GetCurrentThreadId(); }
#else
	inline static int getID() { return m_id; }
#endif

	/// Return the thread's name
	inline const std::string &getName() const { return m_name; }

	/// Set the thread's name
	inline void setName(const std::string &name) { m_name = name; }

	/// Return the parent thread
	inline Thread *getParent() { return m_parent; }

	/// Return the parent thread
	inline const Thread *getParent() const { return m_parent.get(); }

	/// Set the thread's logger
	inline void setLogger(Logger *logger) { m_logger = logger; }

	/// Return the thread's logger
	inline Logger *getLogger() { return m_logger; }

	/// Return the current thread
	inline static Thread *getThread() {
		return m_self->get();
	}

	/// Is this thread still running?
	inline bool isRunning() const { return m_running; }

	/// Start the thread
	void start();

	/**
	 * Detach the thread - after this, <tt>join()</tt>
	 * cannot be used anymore. Releases resources,
	 * would otherwise be held until a call to <tt>join().</tt>
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
#endif
#if defined(__OSX__)
	int m_id;
#elif defined(__LINUX__)
	static int __thread m_id;
#endif
};

MTS_NAMESPACE_END

#endif /* __THREAD_H */
