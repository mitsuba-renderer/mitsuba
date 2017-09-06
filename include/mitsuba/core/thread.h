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
#if !defined(__MITSUBA_CORE_THREAD_H_)
#define __MITSUBA_CORE_THREAD_H_

#include <mitsuba/mitsuba.h>
#include <boost/scoped_ptr.hpp>

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
     * \remark Note that it is currently not possible to
     *         construct Thread instances from Python
     */
    Thread(const std::string &name);

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
    EThreadPriority getPriority() const;

    /**
     * \brief Set the core affinity
     *
     * This function provides a hint to the operating system
     * scheduler that the thread should preferably run
     * on the specified processor core. By default, the parameter
     * is set to -1, which means that there is no affinity.
     */
    void setCoreAffinity(int core);

    /// Return the core affinity
    int getCoreAffinity() const;

    /**
     * \brief Specify whether or not this thread is critical
     *
     * When an thread marked critical crashes from an uncaught
     * exception, the whole process is brought down.
     * The default is \c false.
     */
    void setCritical(bool critical);

    /// Return the value of the critical flag
    bool getCritical() const;

    /// Return the thread ID
    static int getID();

    /// Return the name of this thread
    const std::string &getName() const;

    /// Set the name of this thread
    void setName(const std::string &name);

    /// Return the parent thread
    Thread *getParent();

    /// Return the parent thread (const version)
    const Thread *getParent() const;

    /// Set the logger instance used to process log messages from this thread
    void setLogger(Logger *logger);

    /// Return the thread's logger instance
    Logger *getLogger();

    /// Set the thread's file resolver
    void setFileResolver(FileResolver *fresolver);

    /// Return the thread's file resolver
    FileResolver *getFileResolver();

    /// Return the current thread
    static Thread *getThread();

    /// Is this thread still running?
    bool isRunning() const;

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

    /**
     * \brief Register a thread crash handler
     *
     * A crash handler is called whenever a thread fails with an uncaught
     * exception. This can be used to implement more useful error messages
     * in certain circumstances
     */
    static void registerCrashHandler(bool (*handler)(void));

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Thread();

    /// Thread dispatch function
    static void dispatch(Thread *thread);

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
    struct ThreadPrivate;
    boost::scoped_ptr<ThreadPrivate> d;
};

#if defined(MTS_OPENMP)
#if defined(__OSX__)
/// Variant of \c omp_get_max_threads that works on OSX
extern MTS_EXPORT_CORE int mts_omp_get_max_threads();

/// Variant of \c omp_get_thread_num that works on OSX
extern MTS_EXPORT_CORE int mts_omp_get_thread_num();
#else
#define mts_omp_get_max_threads omp_get_max_threads
#define mts_omp_get_thread_num omp_get_thread_num
#endif
#else
#define mts_omp_get_max_threads() 1
#define mts_omp_get_thread_num() 0
#endif

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_THREAD_H_ */
