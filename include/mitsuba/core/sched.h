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
#if !defined(__MITSUBA_CORE_SCHED_H_)
#define __MITSUBA_CORE_SCHED_H_

#include <mitsuba/core/serialization.h>
#include <mitsuba/core/lock.h>
#include <deque>

/**
 * Uncomment this to enable scheduling debug messages
 */
//#define DEBUG_SCHED 1

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract work unit -- represents a small amount of information
 * that encodes part of a larger processing task.
 *
 * Instead of the usual serialization function and unserialization
 * constructor, implementations of this class supply \ref load()
 * and \ref save() methods that can be used for essentially the
 * same purpose, but without requiring any memory allocations.
 *
 * \sa WorkProcessor
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE WorkUnit : public Object {
public:
    /// Copy the content of another work unit of the same type
    virtual void set(const WorkUnit *workUnit) = 0;

    /// Fill the work unit with content acquired from a binary data stream
    virtual void load(Stream *stream) = 0;

    /// Serialize a work unit to a binary data stream
    virtual void save(Stream *stream) const = 0;

    /// Return a string representation
    virtual std::string toString() const = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~WorkUnit() { }
};

/**
 * \brief Abstract work result -- represents the result of a
 * processed \ref WorkUnit instance.
 *
 * Instead of the usual serialization function and unserialization
 * constructor, implementations of this class supply \ref load()
 * and \ref save() methods that can be used for essentially the
 * same purpose, but without requiring any memory allocations.
 *
 * \sa WorkProcessor
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE WorkResult : public Object {
public:
    /// Fill the work result with content acquired from a binary data stream
    virtual void load(Stream *stream) = 0;

    /// Serialize a work result to a binary data stream
    virtual void save(Stream *stream) const = 0;

    /// Return a string representation
    virtual std::string toString() const = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~WorkResult() { }
};

/**
 * \brief Abstract work processor -- takes work units and turns them into
 *\ref WorkResult instances.
 *
 * When executing a parallel task using Mitsuba's scheduling system, the
 * actual work is done in an implementation of this interface.
 *
 * The class is serializable so that it can be sent over the network if
 * required. It is possible to keep local state in \ref WorkProcessor
 * instances (e.g. scratch space for computations), though anything not
 * returned in the form of a \ref WorkResult will eventually be lost.
 * Each \ref Worker (both locally and remotely) has its own \ref WorkProcessor,
 * and therefore no form of locking is required within instances of this class.
 *
 * \sa WorkUnit
 * \sa WorkResult
 * \sa ParallelProcess
 * \sa Scheduler
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE WorkProcessor : public SerializableObject {
    friend class Scheduler;
public:
    /// Create a work unit of the proper type and size.
    virtual ref<WorkUnit> createWorkUnit() const = 0;

    /// Create a work result of the proper type and size
    virtual ref<WorkResult> createWorkResult() const = 0;

    /**
     * \brief Create a copy of this work processor instance.
     *
     * \remark In practice, before the cloned work processor
     * is actually used, its \ref prepare() method will be called.
     * Therefore, any state that is initialized in \ref prepeare()
     * does not have to be copied.
     */
    virtual ref<WorkProcessor> clone() const = 0;

    /**
     * \brief Called once before processing starts.
     *
     * This is useful for allocating scratch space or resolving references
     * to resource objects. Lengthy computations should be performed in
     * process() instead of here, since this this method will be called
     * while the central scheduler lock is held. A thrown exception will
     * lead to the termination of the parallel process.
     */
    virtual void prepare() = 0;

    /**
     * \brief Process a work unit and store the computed results.
     *
     * The <tt>active</tt> parameter can be used to signal a premature
     * stop of the execution flow. In this case, the work result is allowed
     * to be undefined (it will simply be ignored). A thrown exception will
     * lead to the termination of the parallel process.
     */
    virtual void process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop) = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~WorkProcessor() { }
    /// Protected constructors
    inline WorkProcessor() { }
    inline WorkProcessor(Stream *stream, InstanceManager *manager)
        : SerializableObject(stream, manager) { }

    /**
     * \brief Look up a named resource, which has been bound to
     * the associated parallel process.
     *
     * Throws an exception if the resource is not known / bound.
     */
    SerializableObject *getResource(const std::string &name);
protected:
    std::map<std::string, SerializableObject *> m_resources;
};

/**
 * \brief Abstract parallelizable task.
 *
 * Instances of this class model a larger piece of work that can be split
 * into independent `units' and subsequently farmed out over a cluster or
 * processed locally. After the work units have been completed, the results
 * are pieced back together to a solution of the original large-scale problem.
 *
 * This class implements the core logic running on the central scheduling
 * server, i.e. the part that is responsible for generating work units and
 * accepting their results. The module that performs the actual computation
 * is an instance of \ref WorkProcessor, which is also specified here.
 * Finally, the this class references `resources', which denote
 * chunks of globally shared read-only data required during execution.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE ParallelProcess : public Object {
    friend class Scheduler;
public:
    /// Binding from local resource names to global resource IDs
    typedef std::map<std::string, int> ResourceBindings;

    /// Return codes used by generateWork() and getReturnStatus()
    enum EStatus {
        EUnknown, ///< Unknown return status
        EPause,   ///< Temporarily, no work units can be created
        ESuccess, ///< The process finished / a piece of work was generated
        EFailure  ///< The process failed / no more work is available
    };

    /**
     * \brief Generate a piece of work.
     *
     * Takes a pre-allocated \ref WorkUnit instance of
     * the appropriate sub-type and size (as specified by
     * \ref ParallelProcess::getWorkUnitName()) and
     * fills it with the appropriate content. Returns ESuccess
     * on success and EFailure or EPause when no more work is
     * left -- in that case, the work unit will be ignored and
     * the process completed (\ref EFailure) or temporarily
     * paused (\ref EPause). When \ref EPause was used,
     * resubmission via \ref Scheduler::schedule() will
     * be required once more work is available. In some cases, it
     * is useful to distribute 'nearby' pieces of work to the same
     * processor -- the \c worker parameter can be used to
     * implement this.
     * This function should run as quickly as possible, since it
     * will be executed while the scheduler mutex is held. A
     * thrown exception will lead to the termination of the
     * parallel process.
     *
     * \param unit   Work unit data structure to be filled
     * \param worker ID of the worker executing this function
     */
    virtual EStatus generateWork(WorkUnit *unit, int worker) = 0;

    /**
     * \brief Called whenever a work unit has been completed.
     *
     * Note that this function may concurrently be executed by
     * multiple threads. Also, processing of work results will
     * generally be out of order with respect to the creation
     * in \ref generateWork().
     *
     * When a work unit is only partially completed due to
     * a call to \ref Scheduler::cancel(), the second
     * parameter is set to true.
     * A thrown exception will lead to the termination of
     * the parallel process.
     *
     * \param result Work result to be processed
     * \param cancelled Was the associated work unit not fully completed
     */
    virtual void processResult(const WorkResult *result,
        bool cancelled) = 0;

    /**
     * \brief Called when the parallel process is canceled by
     * \ref Scheduler::cancel().
     *
     * The default implementation does nothing.
     */
    virtual void handleCancellation();

    /**
     * \brief Query the return status of a process after its
     * execution has finished.
     *
     * Returns one of \ref Success, \ref Failure or \ref Unknown.
     * (\ref EUnknown means that the process is either still running
     * or has never been scheduled).
     */
    inline EStatus getReturnStatus() const { return m_returnStatus; }

    /**
     * \brief Create an instance of the algorithm responsible
     * for executing the work units of this parallel process.
     */
    virtual ref<WorkProcessor> createWorkProcessor() const = 0;

    /**
     * \brief Bind a resource to this parallel process.
     *
     * Takes a resource ID as given by the scheduler and associates it
     * with a name. This name can later be used by the work processor
     * to access the resource data.
     *
     * \param name Process-specific name of the resource
     * \param id Resource ID as returned by \ref Scheduler::registerResource()
     * \sa WorkProcessor::getResource
     */
    virtual void bindResource(const std::string &name, int id);

    /**
     * \brief Is this process strictly local?
     *
     * If a process is marked as local, it shouldn't be distributed
     * to remote processing nodes. The default implementation
     * returns false.
     */
    virtual bool isLocal() const;

    /**
     * \brief Return the log level for events associated with this process.
     *
     * By default, this is set to EDebug
     */
    inline ELogLevel getLogLevel() const { return m_logLevel; }

    /**
     * \brief Return a list of all bound resources
     */
    inline const ResourceBindings &getResourceBindings() const { return m_bindings; }

    /**
     * \brief Return a list of plugins required by this parallel process.
     *
     * This is required so that remote machines can load the plugins before
     * they accept work from this process. The default implementation just
     * returns all plugins that are loaded in the current application.
     */
    virtual std::vector<std::string> getRequiredPlugins();

    MTS_DECLARE_CLASS()
protected:
    /// Protected constructor
    inline ParallelProcess() : m_returnStatus(EUnknown),
        m_logLevel(EDebug) { }
    /// Virtual destructor
    virtual ~ParallelProcess() { }
protected:
    ResourceBindings m_bindings;
    EStatus m_returnStatus;
    ELogLevel m_logLevel;
};

class Worker;

/**
 * \brief Centralized task scheduler implementation.
 *
 * Accepts parallelizable jobs and distributes their computational load
 * both locally and remotely. This is done by associating different types
 * of \ref Worker instances with the scheduler. These try to acquire work
 * units from the scheduler, which are then executed on the current machine
 * or sent to remote nodes over a network connection.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Scheduler : public Object {
    friend class Worker;
public:
    /**
     * \brief Schedule a parallelizable process for execution.
     *
     * If the scheduler is currently running and idle, its execution
     * will begin immediately. Returns \c false if the process
     * is already scheduled and has not yet terminated and \c true
     * in any other case.
     */
    bool schedule(ParallelProcess *process);

    /**
     * \brief Block until the process has successfully been completed
     * or canceled prematurely.
     *
     * Returns false if the process does not exist or has already
     * finished by the time \ref wait() is invoked.
     */
    bool wait(const ParallelProcess *process);

    /**
     * \brief Cancel the execution of a parallelizable process.
     *
     * Upon return, no more work from this process is running.
     * Returns false if the process does not exist (anymore).
     */
    inline bool cancel(ParallelProcess *proc) {
        return cancel(proc, false);
    }

    /**
     * \brief Register a serializable resource with the scheduler.
     *
     * A resource should be thought of as a constant state that is shared
     * amongst all processing nodes. Resources can be reused by
     * subsequent parallel processes, and consequently do not have to be
     * re-transmitted over the network. Returns a resource ID, which can be
     * used to reference the associated data.
     */
    int registerResource(SerializableObject *resource);

    /**
     * \brief Register a \a multiple resource with the scheduler.
     *
     * \a Multi means that in comparison to the previous method, a separate
     * instance is provided  for every core. An example where this is useful
     * is to distribute random generator state when performing parallel
     * Monte Carlo simulations. \c resources must be a vector whose
     * length is equal to \ref getCoreCount().
     */
    int registerMultiResource(std::vector<SerializableObject *> &resources);

    /**
     * \brief Increase the reference count of a previously registered resource.
     *
     * The resource must be unregistered an additional time after calling
     * this function.
     *
     * \sa unregisterResource
     */
    void retainResource(int resourceID);

    /**
     * \brief Unregister a resource from the scheduler
     *
     * Note that the resource's won't be removed until all processes using
     * it have terminated)
     *
     * \return \c false if the resource could not be found
     */
    bool unregisterResource(int id);

    /**
     * \brief Return the ID of a registered resource
     *
     * Throws an exception if the resource cannot be found.
     */
    int getResourceID(const SerializableObject *resource) const;

    /// Register a worker with the scheduler
    void registerWorker(Worker *processor);

    /// Unregister a worker from the scheduler
    void unregisterWorker(Worker *processor);

    /// Get the number of workers
    size_t getWorkerCount() const;

    /// Get the number of local workers
    size_t getLocalWorkerCount() const;

    /// Retrieve one of the workers by index
    Worker *getWorker(int index);

    /// Start all workers and begin execution of any scheduled processes
    void start();

    /**
     * \brief Pause the distribution of work units and shut down all
     * running workers.
     *
     * Any currently scheduled work units are still completed.
     * Processing can be resumed via \ref start().
     */
    void pause();

    /// Cancel all running processes and free memory used by resources
    void stop();

    /// Return the total number of cores exposed through this scheduler
    size_t getCoreCount() const;

    /// Does the scheduler have one or more local workers?
    bool hasLocalWorkers() const;

    /// Does the scheduler have one or more remote workers?
    bool hasRemoteWorkers() const;

    /// Return a pointer to the scheduler of this process
    inline static Scheduler *getInstance() { return m_scheduler; }

    /// Has the scheduler been started?
    inline bool isRunning() const { return m_running; }

    /// Is the scheduler currently executing work?
    bool isBusy() const;

    /// Initialize the scheduler of this process -- called once in main()
    static void staticInitialization();

    /// Free the memory taken by staticInitialization()
    static void staticShutdown();

    MTS_DECLARE_CLASS()
public:
    // Public, but shouldn't be part of the documentation
    /// \cond
    struct ProcessRecord {
        /* Unique ID value assigned to this process */
        int id;
        /* Current number of in-flight work units */
        int inflight;
        /* Is the parallel process still generating work */
        bool morework;
        /* Was the process cancelled using \c cancel()?*/
        bool cancelled;
        /* Is the process currently in the queue? */
        bool active;
        /* Signaled every time a work unit arrives */
        ref<ConditionVariable> cond;
        /* Set when the process is done/canceled */
        ref<WaitFlag> done;
        /* Log level for events associated with this process */
        ELogLevel logLevel;

        inline ProcessRecord(int id, ELogLevel logLevel, Mutex *mutex)
         : id(id), inflight(0), morework(true), cancelled(false),
            active(true), logLevel(logLevel) {
            cond = new ConditionVariable(mutex);
            done = new WaitFlag();
        }
    };

    /**
     * Data structure, which contains a piece of work
     * as well as the information required to either execute
     * it locally or submit it to a processing node over the
     * network.
     */
    struct Item {
        int id;
        int workerIndex;
        int coreOffset;
        ParallelProcess *proc;
        ProcessRecord *rec;
        ref<WorkProcessor> wp;
        ref<WorkUnit> workUnit;
        ref<WorkResult> workResult;
        bool stop;

        inline Item() : id(-1), workerIndex(-1), coreOffset(-1),
            proc(NULL), rec(NULL), stop(false) {
        }

        std::string toString() const;
    };

    struct ResourceRecord {
        std::vector<SerializableObject *> resources;
        ref<MemoryStream> stream;
        int refCount;
        bool multi;

        inline ResourceRecord(SerializableObject *resource)
         : resources(1), refCount(1), multi(false) {
            resources[0] = resource;
        }

        inline ResourceRecord(std::vector<SerializableObject *> resources)
         : resources(resources), refCount(1), multi(true) {
        }
    };

    /// A list of status codes returned by acquireWork()
    enum EStatus {
        /// Sucessfully acquired a work unit
        EOK,
        /// There is currently no work (and onlyTry was set to true)
        ENone,
        /// The scheduler is shutting down
        EStop
    };
    /// \endcond

    /// Look up a resource by ID & core index
    SerializableObject *getResource(int id, int coreIndex = -1);

    /// Return a resource in the form of a binary data stream
    const MemoryStream *getResourceStream(int id);

    /**
     * \brief Test whether this is a multi-resource,
     * i.e. different for every core.
     */
    bool isMultiResource(int id) const;
protected:
    /// Protected constructor
    Scheduler();

    /// Virtual destructor
    virtual ~Scheduler();

    /**
     * Acquire a piece of work from the scheduler -- internally
     * used by the different worker implementations.
     */
    EStatus acquireWork(Item &item, bool local, bool onlyTry, bool keepLock);

    /// Release the main scheduler lock -- internally used by the remote worker
    inline void releaseLock() { m_mutex->unlock(); }

    inline void releaseWork(Item &item) {
        ProcessRecord *rec = item.rec;
        try {
            item.proc->processResult(item.workResult, item.stop);
        } catch (const std::exception &ex) {
            Log(EWarn, "Caught an exception - canceling process %i: %s",
                item.id, ex.what());
            cancel(item.proc, true);
            return;
        }
        LockGuard lock(m_mutex);
        --rec->inflight;
        rec->cond->signal();
        if (rec->inflight == 0 && !rec->morework && !item.stop)
            signalProcessTermination(item.proc, item.rec);
    }

    /**
     * Cancel the execution of a parallelizable process. Upon
     * return, no more work from this process is running. When
     * the second parameter is set to true, the number of in-flight
     * work units for this process is reduced by one.
     */
    bool cancel(ParallelProcess *proc, bool reduceInflight);

    /**
     * Internally used to prepare a Scheduler::Item structure
     * when only the process ID is known.
     */
    inline void setProcessByID(Item &item, int id) {
        LockGuard lock(m_mutex);
        ParallelProcess *proc = m_idToProcess[id];
        if (proc == NULL) {
            Log(EError, "Process %i is not locally known!", id);
        };
        item.proc = proc;
        item.id = id;
        item.rec = m_processes[proc];
        item.wp = proc->createWorkProcessor();
        const ParallelProcess::ResourceBindings &bindings = item.proc->getResourceBindings();
        for (ParallelProcess::ResourceBindings::const_iterator it = bindings.begin();
            it != bindings.end(); ++it)
            item.wp->m_resources[(*it).first] = m_scheduler->getResource((*it).second, item.coreOffset);
        try {
            item.wp->prepare();
            item.workUnit = item.wp->createWorkUnit();
            item.workResult = item.wp->createWorkResult();
        } catch (std::exception &) {
            throw;
        }
    }

    /// Announces the termination of a process
    void signalProcessTermination(ParallelProcess *proc, ProcessRecord *rec);
private:
    /// Global scheduler instance
    static ref<Scheduler> m_scheduler;
    /// Mutex, which protects local data structures
    mutable ref<Mutex> m_mutex;
    /// CV used to signal the availability of work
    ref<ConditionVariable> m_workAvailable;
    /// Scheduled processes in FIFO order
    std::deque<int> m_localQueue, m_remoteQueue;
    /// Set of all currently scheduled processes
    std::map<const ParallelProcess *, ProcessRecord *> m_processes;
    /// Maps process IDs to processes
    std::map<int, ParallelProcess *> m_idToProcess;
    /// List of shared resources
    std::map<int, ResourceRecord *> m_resources;
    /// List of all active workers
    std::vector<Worker *> m_workers;
    int m_resourceCounter, m_processCounter;
    bool m_running;
};

/**
 * \brief Base class of all worker implementations
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Worker : public Thread {
    friend class Scheduler;
public:
    /// Return the number of cores exposed by this worker
    inline size_t getCoreCount() const { return m_coreCount; }
    /// Is this a remote worker?
    inline bool isRemoteWorker() const { return m_isRemote; };

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Worker() { }

    /// Protected constructor
    Worker(const std::string &name);

    /* Decrement reference counts to any referenced objects */
    virtual void clear();

    /// Used internally by the scheduler
    virtual void start(Scheduler *scheduler,
        int workerIndex, int coreOffset);

    /**
     * \brief Called to inform a worker that a resource is no longer in use.
     *
     * The remote worker uses this to notify the machine on the other end that
     * the memory used by this resource can now be released.
     */
    virtual void signalResourceExpiration(int id) = 0;

    /**
     * \brief Called to inform a worker that a process has been cancelled
     *
     * Guaranteed to be called while the Scheduler's main lock is held.
     */
    virtual void signalProcessCancellation(int id) = 0;

    /**
     * \brief Called to inform a worker that a process has successfully been
     * completed and any associated resources can be freed.
     */
    virtual void signalProcessTermination(int id) = 0;

    /* Inline functions to access protected members of Scheduler */
    inline Scheduler::EStatus acquireWork(bool local,
            bool onlyTry = false, bool keepLock = false) {
        return m_scheduler->acquireWork(m_schedItem, local, onlyTry, keepLock);
    }

    void releaseSchedulerLock() {
        return m_scheduler->releaseLock();
    }

    /// Release a processed work unit
    inline void releaseWork(Scheduler::Item &item) {
        m_scheduler->releaseWork(item);
    }

    /// Initialize the m_schedItem data structure when only the process ID is known
    void setProcessByID(Scheduler::Item &item, int id) {
        return m_scheduler->setProcessByID(item, id);
    }

    /**
     * Cancel the currently scheduled parallel process and possibly
     * reduce the number of in-flight work units
     * Returns false if the process does not exist (anymore).
     */
    inline void cancel(bool reduceInflight) {
        m_scheduler->cancel(m_schedItem.proc, reduceInflight);
    }
protected:
    Scheduler *m_scheduler;
    Scheduler::Item m_schedItem;
    size_t m_coreCount;
    bool m_isRemote;
};

/**
 * \brief Acquires work from the scheduler and executes
 * it locally.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE LocalWorker : public Worker {
public:
    /**
     * \brief Create a new local worker thread
     *
     * \param coreID
     *   When an CPU core ID (>=0) is specified here, the worker
     *   thread will attempt to register core affinity with the
     *   operating system. Passing -1 disables this.
     *
     * \param name
     *   An identifying string for this thread
     *
     * \param priority
     *   The desired thread priority (not supported on some
     *   operating systems)
     */
    LocalWorker(int coreID, const std::string &name,
        Thread::EThreadPriority priority = Thread::ENormalPriority);

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~LocalWorker();
    /* Worker implementation */
    virtual void run();
    virtual void signalResourceExpiration(int id);
    virtual void signalProcessCancellation(int id);
    virtual void signalProcessTermination(int id);
};

/**
 * \brief Dummy work unit without contents
 * \ingroup libcore
 */
class MTS_EXPORT_CORE DummyWorkUnit : public WorkUnit {
public:
    void set(const WorkUnit *workUnit);
    void load(Stream *stream);
    void save(Stream *stream) const;
    std::string toString() const;

    MTS_DECLARE_CLASS()
protected:
    virtual ~DummyWorkUnit() { }
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SCHED_H_ */
