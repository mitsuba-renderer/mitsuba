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

#include <mitsuba/core/sched.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/mstream.h>

#include <boost/thread/thread.hpp>

MTS_NAMESPACE_BEGIN

SerializableObject *WorkProcessor::getResource(const std::string &name) {
    if (m_resources.find(name) == m_resources.end())
        Log(EError, "Could not find a resource named \"%s\"!", name.c_str());
    return m_resources[name];
}

void ParallelProcess::bindResource(const std::string &name, int id) {
    m_bindings[name] = id;
}

std::vector<std::string> ParallelProcess::getRequiredPlugins() {
    return PluginManager::getInstance()->getLoadedPlugins();
}

void ParallelProcess::handleCancellation() {
}

bool ParallelProcess::isLocal() const {
    return false;
}

/* ==================================================================== */
/*                              Scheduler                               */
/* ==================================================================== */

ref<Scheduler> Scheduler::m_scheduler;

Scheduler::Scheduler() {
    m_mutex = new Mutex();
    m_workAvailable = new ConditionVariable(m_mutex);
    m_resourceCounter = 0;
    m_processCounter = 0;
    m_running = false;
}

Scheduler::~Scheduler() {
    for (size_t i=0; i<m_workers.size(); ++i)
        m_workers[i]->decRef();
}

void Scheduler::registerWorker(Worker *worker) {
    LockGuard lock(m_mutex);
    m_workers.push_back(worker);
    worker->incRef();
}

void Scheduler::unregisterWorker(Worker *worker) {
    LockGuard lock(m_mutex);
    m_workers.erase(std::remove(m_workers.begin(), m_workers.end(), worker),
        m_workers.end());
    worker->decRef();
}

Worker *Scheduler::getWorker(int index) {
    Worker *result = NULL;
    LockGuard lock(m_mutex);
    if (index < (int) m_workers.size()) {
        result = m_workers[index];
    } else {
        Log(EError, "Scheduler::getWorker() - out of bounds");
    }
    return result;
}

size_t Scheduler::getWorkerCount() const {
    size_t count;
    LockGuard lock(m_mutex); // make valgrind/helgrind happy
    count = m_workers.size();
    return count;
}

size_t Scheduler::getLocalWorkerCount() const {
    size_t count = 0;
    LockGuard lock(m_mutex);
    for (size_t i=0; i<m_workers.size(); ++i) {
        if (m_workers[i]->getClass() == MTS_CLASS(LocalWorker))
            count++;
    }
    return count;
}

bool Scheduler::isBusy() const {
    bool result;
    LockGuard lock(m_mutex); // make valgrind/helgrind happy
    result = m_processes.size() > 0;
    return result;
}

int Scheduler::registerResource(SerializableObject *object) {
    LockGuard lock(m_mutex);
    int resourceID = m_resourceCounter++;
    ResourceRecord *rec = new ResourceRecord(object);
    if (hasRemoteWorkers()) {
        ref<InstanceManager> manager = new InstanceManager();
        rec->stream = new MemoryStream();
        rec->stream->setByteOrder(Stream::ENetworkByteOrder);
        manager->serialize(rec->stream, rec->resources[0]);
    }
    m_resources[resourceID] = rec;
    object->incRef();
#if defined(DEBUG_SCHED)
    if (rec->stream.get())
        Log(EDebug, "Registered resource %i: %s (%i KB)", resourceID, object->getClass()->getName().c_str(),
            rec->stream->getPos() / 1024);
    else
        Log(EDebug, "Registered resource %i: %s", resourceID, object->getClass()->getName().c_str());
#endif
    return resourceID;
}

int Scheduler::registerMultiResource(std::vector<SerializableObject *> &objects) {
    if (objects.size() != getCoreCount())
        Log(EError, "registerMultiResource() : resource vector does not have the right size!");
    LockGuard lock(m_mutex);
    int resourceID = m_resourceCounter++;
    ResourceRecord *rec = new ResourceRecord(objects);
    m_resources[resourceID] = rec;
    for (size_t i=0; i<objects.size(); ++i)
        objects[i]->incRef();
#if defined(DEBUG_SCHED)
    Log(EDebug, "Registered multi resource %i: %s", resourceID, objects[0]->getClass()->getName().c_str());
#endif
    return resourceID;
}

void Scheduler::retainResource(int id) {
    LockGuard lock(m_mutex);
    if (m_resources.find(id) == m_resources.end()) {
        Log(EError, "retainResource(): could not find the resource with ID %i!", id);
    }
    ResourceRecord *rec = m_resources[id];
    rec->refCount++;
}

bool Scheduler::unregisterResource(int id) {
    LockGuard lock(m_mutex);
    if (m_resources.find(id) == m_resources.end()) {
        Log(EWarn, "unregisterResource(): could not find the resource with ID %i!", id);
        return false;
    }
    ResourceRecord *rec = m_resources[id];
    if (--rec->refCount == 0) {
#if defined(DEBUG_SCHED)
        Log(EDebug, "Freeing resource %i", id);
#endif
        for (size_t i=0; i<rec->resources.size(); ++i)
            rec->resources[i]->decRef();
        m_resources.erase(id);
        delete rec;
        for (size_t i=0; i<m_workers.size(); ++i)
            m_workers[i]->signalResourceExpiration(id);
    }
    return true;
}

SerializableObject *Scheduler::getResource(int id, int coreIndex) {
    SerializableObject *result = NULL;

    LockGuard lock(m_mutex);
    std::map<int, ResourceRecord *>::iterator it = m_resources.find(id);
    if (it == m_resources.end()) {
        Log(EError, "getResource(): could not find the resource with ID %i!", id);
    }
    ResourceRecord *rec = (*it).second;
    if (rec->multi) {
        if (coreIndex == -1) {
            Log(EError, "getResource(): tried to look up multi resource %i without specifying a core index!", id);
        }
        result = rec->resources.at(coreIndex);
    } else {
        result = rec->resources[0];
    }
    return result;
}

bool Scheduler::isMultiResource(int id) const {
    LockGuard lock(m_mutex);
    std::map<int, ResourceRecord *>::const_iterator it = m_resources.find(id);
    if (it == m_resources.end()) {
        Log(EError, "getResourceStream(): could not find the resource with ID %i!", id);
    }
    bool result = (*it).second->multi;
    return result;
}

const MemoryStream *Scheduler::getResourceStream(int id) {
    LockGuard lock(m_mutex);
    std::map<int, ResourceRecord *>::iterator it = m_resources.find(id);
    if (it == m_resources.end()) {
        Log(EError, "getResourceStream(): could not find the resource with ID %i!", id);
    }
    ResourceRecord *rec = (*it).second;
    if ((*it).second->multi) {
        Log(EError, "getResourceStream(): only standard resource lookups are permitted!");
    }

    if (!rec->stream) {
        ref<InstanceManager> manager = new InstanceManager();
        rec->stream = new MemoryStream();
        rec->stream->setByteOrder(Stream::ENetworkByteOrder);
        manager->serialize(rec->stream, rec->resources[0]);
    }
    return rec->stream;
}

int Scheduler::getResourceID(const SerializableObject *obj) const {
    LockGuard lock(m_mutex);
    std::map<int, ResourceRecord *>::const_iterator it = m_resources.begin();
    for (; it!=m_resources.end(); ++it) {
        ResourceRecord *rec = (*it).second;
        for (size_t j=0; j<rec->resources.size(); ++j) {
            if (rec->resources[j] == obj) {
                int id = (*it).first;
                return id;
            }
        }
    }
    Log(EError, "Resource could not be found!");
    return -1; // Never reached
}

bool Scheduler::schedule(ParallelProcess *process) {
    LockGuard lock(m_mutex);

    if (process->isLocal() && !hasLocalWorkers()) {
        Log(EError, "Cannot schedule a local process when "
            "there are no local workers!");
    }

    if (m_processes.find(process) != m_processes.end()) {
        ProcessRecord *rec = m_processes[process];
        if (rec->morework && !rec->active) {
            /* Paused process - reactivate */
#if defined(DEBUG_SCHED)
            Log(rec->logLevel, "Waking inactive process %i..", rec->id);
#endif
            rec->active = true;
            m_localQueue.push_back(rec->id);
            if (!process->isLocal())
                m_remoteQueue.push_back(rec->id);
            m_workAvailable->broadcast();
            return true;
        }
        /* The process is still active */
        return false;
    }

    /* First, check that all resources are available and increase
       their reference count */
    const ParallelProcess::ResourceBindings &bindings = process->getResourceBindings();
    for (ParallelProcess::ResourceBindings::const_iterator it = bindings.begin();
        it != bindings.end(); ++it) {
        if (m_resources.find((*it).second) == m_resources.end()) {
            Log(EError, "Unable to find resource %i (%s) referenced by %s",
                (*it).second, (*it).first.c_str(), process->toString().c_str());
        }
        m_resources[(*it).second]->refCount++;
    }
    ProcessRecord *rec = new ProcessRecord(m_processCounter++,
        process->getLogLevel(), m_mutex);
    m_processes[process] = rec;
#if defined(DEBUG_SCHED)
    Log(rec->logLevel, "Scheduling process %i: %s..", rec->id, process->toString().c_str());
#endif
    process->m_returnStatus = ParallelProcess::EUnknown;
    m_idToProcess[rec->id] = process;
    m_localQueue.push_back(rec->id);
    if (!process->isLocal())
        m_remoteQueue.push_back(rec->id);
    process->incRef();
    m_workAvailable->broadcast();
    return true;
}

bool Scheduler::hasRemoteWorkers() const {
    bool hasRemoteWorkers = false;
    LockGuard lock(m_mutex);
    for (size_t i=0; i<m_workers.size(); ++i)
        hasRemoteWorkers |= m_workers[i]->isRemoteWorker();
    return hasRemoteWorkers;
}

bool Scheduler::hasLocalWorkers() const {
    bool hasLocalWorkers = false;
    LockGuard lock(m_mutex);
    for (size_t i=0; i<m_workers.size(); ++i)
        hasLocalWorkers |= !m_workers[i]->isRemoteWorker();
    return hasLocalWorkers;
}

bool Scheduler::wait(const ParallelProcess *process) {
    UniqueLock lock(m_mutex);

    std::map<const ParallelProcess *, ProcessRecord *>::iterator it =
        m_processes.find(process);
    if (it == m_processes.end()) {
        /* The process is not known */
        return false;
    }

    ProcessRecord *rec = (*it).second;
    /* Increase the WaitFlag reference count as otherwise,
       it might be deleted before having a chance to verify
       that the flag has really changed */

#if defined(DEBUG_SCHED)
    Log(rec->logLevel, "Waiting for process %i", rec->id);
#endif

    WaitFlag *flag = rec->done;
    flag->incRef();
    lock.unlock();
    flag->wait();

    lock.lock();
    flag->decRef();
    lock.unlock();
    return true;
}

bool Scheduler::cancel(ParallelProcess *process, bool reduceInflight) {
    UniqueLock lock(m_mutex);
    std::map<const ParallelProcess *, ProcessRecord *>::iterator it =
        m_processes.find(process);
    if (it == m_processes.end()) {
#if defined(DEBUG_SCHED)
        Log(EDebug, "Scheduler::cancel() - the process is not currently running");
#endif
        return false;
    }

    ProcessRecord *rec = (*it).second;
    if (reduceInflight) {
        --rec->inflight;
        rec->cond->signal();
    }

    if (rec->cancelled) {
#if defined(DEBUG_SCHED)
        Log(rec->logLevel, "Scheduler::cancel() - the process is already being cancelled. "
            "Waiting until this has happened..");
#endif
        lock.unlock();
        wait(process);
        return true;
    }

#if defined(DEBUG_SCHED)
    Log(rec->logLevel, "Cancelling process %i (%i work units in flight)..", rec->id, rec->inflight);
#endif

    for (size_t i=0; i<m_workers.size(); ++i)
        m_workers[i]->signalProcessCancellation(rec->id);

    /* Ensure that this process won't be scheduled again */
    m_localQueue.erase(std::remove(m_localQueue.begin(), m_localQueue.end(), rec->id),
        m_localQueue.end());
    m_remoteQueue.erase(std::remove(m_remoteQueue.begin(), m_remoteQueue.end(), rec->id),
        m_remoteQueue.end());

    /* Ensure that the process won't be considered 'done' when the
       last in-flight work unit is returned */
    rec->morework = true;
    rec->cancelled = true;

    /* Now wait until no more work from this process circulates and release
       the lock while waiting. */
    while (rec->inflight != 0)
        rec->cond->wait();

    /* Decrease the reference count of all bound resources */
    const ParallelProcess::ResourceBindings &bindings = process->getResourceBindings();
    for (ParallelProcess::ResourceBindings::const_iterator it = bindings.begin();
        it != bindings.end(); ++it) {
        unregisterResource((*it).second);
    }

    m_processes.erase(process);
    m_idToProcess.erase(rec->id);
    process->m_returnStatus = ParallelProcess::EFailure;

    try {
        process->handleCancellation();
    } catch (const std::exception &ex) {
        Log(EWarn, "Process %i's cancellation handler threw an exception.", ex.what());
    }

    /* Wake up any threads waiting on this process */
    rec->done->set(true);
    process->decRef();

#if defined(DEBUG_SCHED)
    Log(rec->logLevel, "Process %i was cancelled.", rec->id);
#endif

    delete rec;

    return true;
}

Scheduler::EStatus Scheduler::acquireWork(Item &item,
        bool local, bool onlyTry, bool keepLock) {
    UniqueLock lock(m_mutex);
    std::deque<int> &queue = local ? m_localQueue : m_remoteQueue;
    while (true) {
        if (onlyTry && queue.size() == 0) {
            return ENone;
        }

        /* Wait until work is available and return false
           if stop() is called */
        while (queue.size() == 0 && m_running)
            m_workAvailable->wait();

        if (!m_running) {
            return EStop;
        }

        /* Try to create a work unit from the parallel
           process currently on top of the queue */
        ParallelProcess::EStatus wStatus;
        try {
            int id = queue.front();
            if (item.id != id) {
                /* First work unit from this parallel process - establish
                   connections to referenced resources and prepare the
                   work processor */
                setProcessByID(item, id);
            }

            wStatus = item.proc->generateWork(item.workUnit, item.workerIndex);
        } catch (const std::exception &ex) {
            Log(EWarn, "Caught an exception - canceling process %i: %s",
                item.id, ex.what());
            cancel(item.proc);
            continue;
        }

        if (wStatus == ParallelProcess::ESuccess) {
            break;
        } else if (wStatus == ParallelProcess::EFailure) {
#if defined(DEBUG_SCHED)
            if (item.rec->morework)
                Log(item.rec->logLevel, "Process %i has finished generating work", item.rec->id);
#endif
            item.rec->morework = false;
            item.rec->active = false;
            queue.pop_front();
            if (item.rec->inflight == 0)
                signalProcessTermination(item.proc, item.rec);
        } else if (wStatus == ParallelProcess::EPause) {
#if defined(DEBUG_SCHED)
            Log(item.rec->logLevel, "Pausing process %i", item.rec->id);
#endif
            item.rec->active = false;
            queue.pop_front();
        }
    }

    item.rec->inflight++;
    item.stop = false;

    if (!keepLock)
        lock.unlock();
    else
        lock.release(); /* Avoid the automatic unlocking upon destruction */

    boost::this_thread::yield();
    return EOK;
}

void Scheduler::signalProcessTermination(ParallelProcess *proc, ProcessRecord *rec) {
#if defined(DEBUG_SCHED)
    Log(rec->logLevel, "Process %i is complete.", rec->id);
#endif
    for (size_t i=0; i<m_workers.size(); ++i)
        m_workers[i]->signalProcessTermination(rec->id);
    /* The parallel process has been completed. Decrease the reference count
        of all used resources */
    const ParallelProcess::ResourceBindings &bindings = proc->getResourceBindings();
    for (ParallelProcess::ResourceBindings::const_iterator it = bindings.begin();
        it != bindings.end(); ++it) {
        unregisterResource((*it).second);
    }
    rec->done->set(true);
    m_processes.erase(proc);
    m_localQueue.erase(std::remove(m_localQueue.begin(), m_localQueue.end(), rec->id),
        m_localQueue.end());
    m_remoteQueue.erase(std::remove(m_remoteQueue.begin(), m_remoteQueue.end(), rec->id),
        m_remoteQueue.end());
    proc->m_returnStatus = ParallelProcess::ESuccess;
    m_idToProcess.erase(rec->id);
    delete rec;
    proc->decRef();
}

void Scheduler::start() {
    Assert(!m_running);
#if defined(DEBUG_SCHED)
    Log(EDebug, "Starting ..");
#endif
    m_running = true;
    if (m_workers.size() == 0)
        Log(EError, "Cannot start the scheduler - there are no registered workers!");

    int coreIndex = 0;
    for (size_t i=0; i<m_workers.size(); ++i) {
        m_workers[i]->start(this, (int) i, coreIndex);
        coreIndex += (int) m_workers[i]->getCoreCount();
    }
}

void Scheduler::pause() {
    Assert(m_running);
#if defined(DEBUG_SCHED)
    Log(EDebug, "Pausing ..");
#endif
    UniqueLock lock(m_mutex);
    m_running = false;
    /* Wake up any workers waiting for work units */
    m_workAvailable->broadcast();
    lock.unlock();
    /* Return when all of them have finished */
    for (size_t i=0; i<m_workers.size(); ++i)
        m_workers[i]->join();
    /* Decrement reference counts to any referenced objects */
    for (size_t i=0; i<m_workers.size(); ++i)
        m_workers[i]->clear();
}

void Scheduler::stop() {
    if (m_running)
        pause();
#if defined(DEBUG_SCHED)
    Log(EDebug, "Stopping ..");
#endif
    LockGuard lock(m_mutex);
    for (std::map<const ParallelProcess *, ProcessRecord *>::iterator
            it = m_processes.begin(); it != m_processes.end(); ++it) {
        (*it).first->decRef();
        (*it).second->done->set(true);
        delete (*it).second;
    }
    m_processes.clear();
    m_idToProcess.clear();
    m_localQueue.clear();
    m_remoteQueue.clear();
    for (std::map<int, ResourceRecord *>::iterator
        it = m_resources.begin(); it != m_resources.end(); ++it) {
        ResourceRecord *rec = (*it).second;
        for (size_t i=0; i<rec->resources.size(); ++i)
            rec->resources[i]->decRef();
        delete rec;
    }
    m_resources.clear();
}

size_t Scheduler::getCoreCount() const {
    size_t coreCount = 0;
    LockGuard lock(m_mutex);
    for (size_t i=0; i<m_workers.size(); ++i)
        coreCount += m_workers[i]->getCoreCount();
    return coreCount;
}

std::string Scheduler::Item::toString() const {
    std::ostringstream oss;
    oss << "Scheduler::Item[" << endl
        << "  id=" << rec->id << "," << endl
        << "  coreOffset=" << coreOffset << "," << endl
        << "  proc=" << (proc == NULL ? "null" : indent(proc->toString()).c_str()) << "," << endl
        << "  wp=" << (wp == NULL ? "null" : indent(wp->toString()).c_str()) << "," << endl
        << "  workUnit=" << (wp == NULL ? "null": indent(workUnit->toString()).c_str()) << endl
        << "]";
    return oss.str();
}

void Scheduler::staticInitialization() {
    m_scheduler = new Scheduler();
}

void Scheduler::staticShutdown() {
    m_scheduler->stop();
    m_scheduler = NULL;
}

/* ==================================================================== */
/*                         Worker implementations                       */
/* ==================================================================== */

Worker::Worker(const std::string &name) : Thread(name), m_coreCount(0), m_isRemote(false) {
}

void Worker::clear() {
    m_schedItem.wp = NULL;
    m_schedItem.workUnit = NULL;
    m_schedItem.workResult = NULL;
    m_schedItem.id = -1;
}

void Worker::start(Scheduler *scheduler, int workerIndex, int coreOffset) {
    m_schedItem.workerIndex = workerIndex;
    m_schedItem.coreOffset = coreOffset;
    m_scheduler = scheduler;
    Thread::start();
}

LocalWorker::LocalWorker(int coreID, const std::string &name,
        Thread::EThreadPriority priority) : Worker(name) {
    if (coreID >= 0)
        setCoreAffinity(coreID);
    m_coreCount = 1;
#if !defined(__LINUX__)
    /* Don't set thead priority on Linux, since it uses
       dynamic priorities */
    setPriority(priority);
#endif
}

LocalWorker::~LocalWorker() {
}

void LocalWorker::run() {
    while (acquireWork(true) != Scheduler::EStop) {
        try {
            m_schedItem.wp->process(m_schedItem.workUnit, m_schedItem.workResult, m_schedItem.stop);
        } catch (const std::exception &ex) {
            m_schedItem.stop = true;
            releaseWork(m_schedItem);
            ELogLevel warnLogLevel = Thread::getThread()->getLogger()->getErrorLevel() == EError
                ? EWarn : EInfo;
            Log(warnLogLevel, "Caught an exception - canceling process %i: %s",
                m_schedItem.id, ex.what());
            cancel(false);
            continue;
        }
        releaseWork(m_schedItem);
    }
}

void LocalWorker::signalResourceExpiration(int id) {
    /* No-op for local workers */
}

void LocalWorker::signalProcessTermination(int id) {
    /* No-op for local workers */
}

void LocalWorker::signalProcessCancellation(int id) {
    if (m_schedItem.id == id)
        m_schedItem.stop = true;
}


/* ==================================================================== */
/*                        Work unit implementations                    */
/* ==================================================================== */

void DummyWorkUnit::set(const WorkUnit *workUnit) { }
void DummyWorkUnit::load(Stream *stream) { }
void DummyWorkUnit::save(Stream *stream) const { }
std::string DummyWorkUnit::toString() const {
    return "DummyWorkUnit[]";
}

MTS_IMPLEMENT_CLASS(Worker, true, Thread)
MTS_IMPLEMENT_CLASS(WorkUnit, true, Object)
MTS_IMPLEMENT_CLASS(DummyWorkUnit, false, WorkUnit)
MTS_IMPLEMENT_CLASS(WorkResult, true, Object)
MTS_IMPLEMENT_CLASS(LocalWorker, false, Worker)
MTS_IMPLEMENT_CLASS(WorkProcessor, true, Object)
MTS_IMPLEMENT_CLASS(Scheduler, false, Object)
MTS_IMPLEMENT_CLASS(ParallelProcess, true, Object)
MTS_NAMESPACE_END
