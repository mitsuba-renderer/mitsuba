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

#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/version.h>

MTS_NAMESPACE_BEGIN

class CancelThread : public Thread {
public:
    CancelThread(ParallelProcess *proc) : Thread("cthr"), m_proc(proc) { }

    void run() {
        Scheduler::getInstance()->cancel(m_proc);
        m_proc = NULL;
    }
protected:
    virtual ~CancelThread() { }
private:
    ref<ParallelProcess> m_proc;
};

RemoteWorker::RemoteWorker(const std::string &name, Stream *stream) : Worker(name), m_stream(stream) {
    const size_t dataLength = strlen(MTS_VERSION)+3;
    char *data = (char *) alloca(dataLength);
    strncpy(data, MTS_VERSION, strlen(MTS_VERSION)+1);
    data[dataLength-2] = SPECTRUM_SAMPLES;
#ifdef DOUBLE_PRECISION
    data[dataLength-1] = 1;
#else
    data[dataLength-1] = 0;
#endif
    m_stream->writeShort(StreamBackend::EHello);
    m_stream->write(data, dataLength);
    m_stream->flush();

    int msg = m_stream->readShort();
    if (msg == StreamBackend::EIncompatible)
        Log(EError, "The server reported a version or configuration mismatch -- unable to connect!");
    else if (msg != StreamBackend::EHello)
        Log(EError, "Received an invalid response!");
    m_coreCount = m_stream->readShort();
    m_nodeName = m_stream->readString();
    m_mutex = new Mutex();
    m_finishCond = new ConditionVariable(m_mutex);
    m_memStream = new MemoryStream();
    m_memStream->setByteOrder(Stream::ENetworkByteOrder);
    m_reader = new RemoteWorkerReader(this);
    m_reader->start();
    m_inFlight = 0;
    m_isRemote = true;
    Log(EDebug, "Connection to \"%s\" established (%i cores).",
        m_nodeName.c_str(), m_coreCount);
}

RemoteWorker::~RemoteWorker() {
    Log(EDebug, "Shutting down");
    if (!m_reader || !m_mutex || !m_memStream)
        return;

    LockGuard lock(m_mutex);
    m_reader->shutdown();
    m_memStream->writeShort(StreamBackend::EQuit);
    try {
        flush();
    } catch (std::runtime_error &e) {
        Log(EWarn, "Could not flush buffer: %s", e.what());
    }
    m_reader->join();
}

void RemoteWorker::start(Scheduler *scheduler, int workerIndex, int coreOffset) {
    Worker::start(scheduler, workerIndex, coreOffset);
    m_reader->m_schedItem.coreOffset = coreOffset;
}

void RemoteWorker::flush() {
    m_memStream->seek(0);
    m_memStream->copyTo(m_stream);
    m_memStream->reset();
    m_stream->flush();
}

void RemoteWorker::run() {
    Scheduler::EStatus status;

    while ((status = acquireWork(false, true, true)) != Scheduler::EStop) {
        if (status == Scheduler::ENone) {
            flush();
            if ((status = acquireWork(false, false, true)) == Scheduler::EStop)
                break;
        }
        /* Acquire the lock each iteration, release it at the end of each one */
        LockGuard lock(m_mutex);

        const int id = m_schedItem.rec->id;
        if (m_processes.find(id) == m_processes.end()) {
            /* The backend has not yet seen this process - submit
               all information required to receive and execute work
               units on the other side */
            std::vector<std::pair<int, const MemoryStream *> > resources;
            std::vector<std::pair<int, const SerializableObject *> > multiResources;

            /* First, look up all resources required by this process (the scheduler lock
               needs to be held for that, so do it quickly) */
            const ParallelProcess::ResourceBindings &bindings = m_schedItem.proc->getResourceBindings();
            for (ParallelProcess::ResourceBindings::const_iterator it = bindings.begin();
                it != bindings.end(); ++it) {
                int resID = (*it).second;

                if (m_resources.find(resID) == m_resources.end()) {
                    if (!m_scheduler->isMultiResource(resID)) {
                        resources.push_back(std::pair<int, const MemoryStream *>(resID,
                            m_scheduler->getResourceStream(resID)));
                    } else {
                        for (size_t i=0; i<m_coreCount; ++i)
                            multiResources.push_back(std::pair<int, const SerializableObject *>(resID,
                                m_scheduler->getResource(resID, (int) (m_schedItem.coreOffset + i))));
                    }
                }
                m_resources.insert(resID);
            }
            /* We can safely release the scheduler lock now. The local message buffer lock is still
               held and thus, there is no danger of sending a cancellation message for a process which
               the remote side has not even seen yet. */
            releaseSchedulerLock();

            std::vector<std::string> plugins = m_schedItem.proc->getRequiredPlugins();
            for (size_t i=0; i<plugins.size(); ++i) {
                if (m_plugins.find(plugins[i]) == m_plugins.end()) {
                    /* Ask the remote side to load any plugins, which might be required first */
                    m_memStream->writeShort(StreamBackend::EEnsurePluginLoaded);
                    m_memStream->writeString(plugins[i]);
                    m_plugins.insert(plugins[i]);
                }
            }

            m_memStream->writeShort(StreamBackend::ENewProcess);
            m_memStream->writeInt(id);
            m_memStream->writeInt(m_schedItem.proc->getLogLevel());

            ref<InstanceManager> manager = new InstanceManager();
            manager->serialize(m_memStream, m_schedItem.wp);
            m_processes.insert(id);

            for (size_t i=0; i<resources.size(); ++i) {
                int resID = resources[i].first;
                const MemoryStream *resStream = resources[i].second;
                Log(EDebug, "Sending resource %i to \"%s\" (%i KB)", resID, m_nodeName.c_str(),
                    resStream->getPos() / 1024);
                m_memStream->writeShort(StreamBackend::ENewResource);
                m_memStream->writeInt(resID);
                m_memStream->writeSize(resStream->getPos());
                m_memStream->write(resStream->getData(), resStream->getPos());
            }

            for (size_t i=0; i<multiResources.size(); i += m_coreCount) {
                int resID = multiResources[i].first;
                ref<MemoryStream> resStream = new MemoryStream();
                ref<InstanceManager> manager = new InstanceManager();
                resStream->setByteOrder(Stream::ENetworkByteOrder);
                for (size_t j=0; j<m_coreCount; ++j)
                    manager->serialize(resStream, multiResources[i+j].second);
                Log(EDebug, "Sending multi resource %i to \"%s\" (%i KB)", resID, m_nodeName.c_str(),
                    resStream->getPos() / 1024);
                m_memStream->writeShort(StreamBackend::ENewMultiResource);
                m_memStream->writeInt(resID);
                m_memStream->writeSize(resStream->getPos());
                m_memStream->write(resStream->getData(), resStream->getPos());
            }

            for (ParallelProcess::ResourceBindings::const_iterator it = bindings.begin();
                it != bindings.end(); ++it) {
                m_memStream->writeShort(StreamBackend::EBindResource);
                m_memStream->writeInt(id);
                m_memStream->writeString((*it).first);
                m_memStream->writeInt((*it).second);
            }
        } else {
            releaseSchedulerLock();
        }

        m_memStream->writeShort(StreamBackend::EWorkUnit);
        m_memStream->writeInt(id);
        m_schedItem.workUnit->save(m_memStream);

        if (++m_inFlight == MTS_BACKLOG_FACTOR * m_coreCount) {
            flush();
            /* There are now too many packets in transit. Wait
               until this clears up a bit before attempting to
               send more work */
            while (m_inFlight > MTS_CONTINUE_FACTOR * m_coreCount)
                m_finishCond->wait();
        }
    }
    LockGuard lock(m_mutex);
    flush();
}

void RemoteWorker::signalResourceExpiration(int id) {
    LockGuard lock(m_mutex);
    if (m_resources.find(id) == m_resources.end()) {
        return;
    }
    m_memStream->writeShort(StreamBackend::EResourceExpired);
    m_memStream->writeInt(id);
    flush();
    m_resources.erase(id);
}

void RemoteWorker::signalProcessCancellation(int id) {
    LockGuard lock(m_mutex);
    if (m_processes.find(id) == m_processes.end()) {
        return;
    }
    m_memStream->writeShort(StreamBackend::EProcessCancelled);
    m_memStream->writeInt(id);
    flush();
    m_processes.erase(id);
}

void RemoteWorker::signalProcessTermination(int id) {
    LockGuard lock(m_mutex);
    if (m_processes.find(id) == m_processes.end()) {
        return;
    }
    m_memStream->writeShort(StreamBackend::EProcessTerminated);
    m_memStream->writeInt(id);
    flush();
    m_processes.erase(id);
}

void RemoteWorker::clear() {
    Worker::clear();
    m_reader->m_schedItem.wp = NULL;
    m_reader->m_schedItem.workUnit = NULL;
    m_reader->m_schedItem.workResult = NULL;
    m_reader->m_schedItem.id = -1;
}


RemoteWorkerReader::RemoteWorkerReader(RemoteWorker *worker)
 : Thread(formatString("%s_r", worker->getName().c_str())),
    m_parent(worker), m_shutdown(false), m_currentID(-1) {
    m_stream = m_parent->m_stream;
    setCritical(true);
}

void RemoteWorkerReader::run() {
    int id=-1; short msg=-1;

    while (true) {
        try {
            msg = m_stream->readShort();
            id = m_stream->readInt();

            if (id != m_currentID) {
                m_parent->setProcessByID(m_schedItem, id);
                m_currentID = id;
            }

            switch (msg) {
                case StreamBackend::EWorkResult:
                    m_schedItem.workResult->load(m_stream);
                    m_schedItem.stop = false;
                    m_parent->releaseWork(m_schedItem);
                    m_parent->signalCompletion();
                    break;
                case StreamBackend::ECancelledWorkResult:
                    m_schedItem.stop = true;
                    m_parent->releaseWork(m_schedItem);
                    m_parent->signalCompletion();
                    break;
                case StreamBackend::EProcessCancelled: {
                        Log(EWarn, "Process %i encountered a problem on node \"%s\"."
                            " - Cancelling the process..", id, m_parent->getNodeName().c_str());
                        /* We can't block here waiting for the process to terminate, since
                        we need to listen for canceled in-flight work units. Handle
                        the cancellation notification in a separate thread */
                        CancelThread *thr = new CancelThread(m_schedItem.proc);
                        thr->incRef();
                        thr->start();
                        m_joinThreads.push_back(thr);
                    }
                    break;
                default:
                    Log(EError, "Received an unknown message (type %i)", id);
            };
        } catch (std::runtime_error &e) {
            if (!m_shutdown)
                throw e;
            break;
        }
    }
    for (size_t i=0; i<m_joinThreads.size(); ++i) {
        m_joinThreads[i]->join();
        m_joinThreads[i]->decRef();
    }
}

/* ==================================================================== */
/*                         Stream server backend                        */
/* ==================================================================== */

StreamBackend::StreamBackend(const std::string &thrName, Scheduler *scheduler,
        const std::string &nodeName, Stream *stream, bool detach) : Thread(thrName),
        m_scheduler(scheduler), m_nodeName(nodeName), m_stream(stream), m_detach(detach) {
    m_sendMutex = new Mutex();
    m_memStream = new MemoryStream();
    m_memStream->setByteOrder(Stream::ENetworkByteOrder);
}

StreamBackend::~StreamBackend() { }

void StreamBackend::run() {
    if (m_detach)
        detach();

    if (m_stream->getClass()->derivesFrom(MTS_CLASS(SocketStream))) {
        SocketStream *sstream = static_cast<SocketStream *>(m_stream.get());
        Log(EInfo, "Incoming connection from %s", sstream->getPeer().c_str());
    }

    short msg = m_stream->readShort();
    if (msg != EHello) {
        Log(EWarn, "Received invalid data -- dropping the connection!");
        return;
    }

    const size_t dataLength = strlen(MTS_VERSION)+3;
    char *data    = (char *) alloca(dataLength),
         *refData = (char *) alloca(dataLength);
    strncpy(refData, MTS_VERSION, strlen(MTS_VERSION)+1);
    refData[dataLength-2] = SPECTRUM_SAMPLES;
#ifdef DOUBLE_PRECISION
    refData[dataLength-1] = 1;
#else
    refData[dataLength-1] = 0;
#endif
    m_stream->read(data, dataLength);

    if (memcmp(data, refData, dataLength) != 0) {
        m_stream->writeShort(EIncompatible);
        m_stream->flush();
        Log(EWarn, "The client either the wrong version, or it is compiled "
            "using different configuration flags -- dropping the connection!");
        return;
    }

    Log(EDebug, "Program versions match.");
    m_memStream->writeShort(EHello);
    m_memStream->writeShort((short) m_scheduler->getCoreCount());
    m_memStream->writeString(m_nodeName);
    m_memStream->seek(0);
    m_memStream->copyTo(m_stream);
    m_stream->flush();
    bool running = true;

    try {
        while (running) {
            msg = m_stream->readShort();
            switch (msg) {
                case ENewProcess: {
                        int id = m_stream->readInt();
                        ELogLevel logLevel = (ELogLevel) m_stream->readInt();
                        ref<InstanceManager> manager = new InstanceManager();
                        ref<WorkProcessor> wp = static_cast<WorkProcessor *>(manager->getInstance(m_stream));
                        RemoteProcess *rp = new RemoteProcess(id, logLevel, this, wp);
                        rp->incRef();
                        m_processes[id] = rp;
                    }
                    break;
                case ENewResource: {
                        int id = m_stream->readInt();
                        size_t size = m_stream->readSize();
                        ref<InstanceManager> manager = new InstanceManager();
                        ref<MemoryStream> mstream = new MemoryStream(size);
                        mstream->setByteOrder(Stream::ENetworkByteOrder);
                        m_stream->copyTo(mstream, size);
                        mstream->seek(0);
                        ref<SerializableObject> res = static_cast<SerializableObject *>(manager->getInstance(mstream));
                        m_resources[id] = m_scheduler->registerResource(res);
                    }
                    break;
                case ENewMultiResource: {
                        int id = m_stream->readInt();
                        size_t size = m_stream->readSize();
                        ref<InstanceManager> manager = new InstanceManager();
                        ref<MemoryStream> mstream = new MemoryStream(size);
                        mstream->setByteOrder(Stream::ENetworkByteOrder);
                        m_stream->copyTo(mstream, size);
                        mstream->seek(0);
                        size_t coreCount = m_scheduler->getCoreCount();
                        std::vector<SerializableObject *> objects(coreCount);
                        for (size_t i=0; i<coreCount; ++i)
                            objects[i] = static_cast<SerializableObject *>(manager->getInstance(mstream));
                        m_resources[id] = m_scheduler->registerMultiResource(objects);
                    }
                    break;
                case EEnsurePluginLoaded: {
                        std::string name = m_stream->readString();
                        PluginManager::getInstance()->ensurePluginLoaded(name);
                    }
                    break;
                case EBindResource: {
                        int procID = m_stream->readInt();
                        std::string resName = m_stream->readString();
                        int resID = m_stream->readInt();
                        RemoteProcess *rp = m_processes[procID];
                        rp->bindResource(resName, m_resources[resID]);
                    }
                    break;
                case EWorkUnit : {
                        int id = m_stream->readInt();
                        RemoteProcess *rp = m_processes[id];
                        WorkUnit *wu = rp->getEmptyWorkUnit();
                        wu->load(m_stream);
                        rp->putFullWorkUnit(wu);
                        m_scheduler->schedule(rp);
                    }
                    break;
                case EProcessTerminated : {
                        int id = m_stream->readInt();
                        RemoteProcess *rp = m_processes[id];
                        rp->setDone();
                        rp->decRef();
                        m_scheduler->schedule(rp);
                        m_processes.erase(id);
                    }
                    break;
                case EProcessCancelled: {
                        int id = m_stream->readInt();
                        RemoteProcess *rp = m_processes[id];
                        m_scheduler->cancel(rp);
                        m_processes.erase(id);
                        rp->decRef();
                    }
                    break;
                case EResourceExpired: {
                        int id = m_stream->readInt();
                        int localID = m_resources[id];
                        m_scheduler->unregisterResource(localID);
                        m_resources.erase(id);
                    }
                    break;
                case EQuit: running = false; break;
                default: Log(EError, "Received an unknown message type: %i", msg);
            }
        }
    } catch (std::exception &e) {
        Log(EWarn, "Uncaught exception \"%s\" - cleaning up!", e.what());
    } catch (...) {
        Log(EWarn, "Uncaught exception (unknown type) - cleaning up!");
    }

    for (std::map<int, RemoteProcess *>::const_iterator it = m_processes.begin();
        it != m_processes.end(); ++it) {

        Log(EWarn, "Cancelling stray process %i", (*it).first);
        m_scheduler->cancel((*it).second);
        (*it).second->decRef();
    }

    for (std::map<int, int>::const_iterator it = m_resources.begin();
        it != m_resources.end(); ++it) {

        Log(EWarn, "Removing stray resource %i", (*it).first);
        m_scheduler->unregisterResource((*it).second);
    }

    if (m_stream->getClass()->derivesFrom(MTS_CLASS(SocketStream))) {
        SocketStream *sstream = static_cast<SocketStream *>(m_stream.get());
        Log(EInfo, "Closing connection to %s - received %i KB / sent %i KB",
            sstream->getPeer().c_str(), (int) (sstream->getReceivedBytes() / 1024),
            (int) (sstream->getSentBytes() / 1024));
    }
}

void StreamBackend::sendCancellation(int id, int numLost) {
    Log(EInfo, "Notifying the remote side about the cancellation of process %i", id);

    LockGuard lock(m_sendMutex);
    m_memStream->reset();
    m_memStream->writeShort(EProcessCancelled);
    m_memStream->writeInt(id);
    for (int i=0; i<numLost; ++i) {
        m_memStream->writeShort(ECancelledWorkResult);
        m_memStream->writeInt(id);
    }
    try {
        m_memStream->seek(0);
        m_memStream->copyTo(m_stream);
        m_stream->flush();
    } catch (std::exception &) {
        Log(EWarn, "Connection error - could not submit cancellation notification");
        /* A connection failure occurred - this will eventually be
           caught and handled in run() and is therefore ignored for now */
    }
}

void StreamBackend::sendWorkResult(int id, const WorkResult *result, bool cancelled) {
    LockGuard lock(m_sendMutex);
    m_memStream->reset();
    m_memStream->writeShort(cancelled ? ECancelledWorkResult : EWorkResult);
    m_memStream->writeInt(id);
    if (!cancelled)
        result->save(m_memStream);
    try {
        m_memStream->seek(0);
        m_memStream->copyTo(m_stream);
        m_stream->flush();
    } catch (std::exception &) {
        Log(EWarn, "Connection error - could not submit work result");
        /* A connection failure occurred - this will eventually be
           caught and handled in run() and is therefore ignored for now */
    }
}

/* ==================================================================== */
/*                            Remote process                            */
/* ==================================================================== */

RemoteProcess::RemoteProcess(int id, ELogLevel logLevel,
        StreamBackend *backend, WorkProcessor *wp)
        : m_id(id), m_backend(backend), m_wp(wp) {
    m_mutex = new Mutex();
    m_logLevel = logLevel;
    m_done = false;
}

RemoteProcess::~RemoteProcess() {
#if defined(DEBUG_SCHED)
    Log(m_logLevel, "Destroying the remote process %i", m_id);
#endif
    for (size_t i=0; i<m_full.size(); ++i)
        m_full[i]->decRef();
    for (size_t i=0; i<m_empty.size(); ++i)
        m_empty[i]->decRef();
}

ParallelProcess::EStatus RemoteProcess::generateWork(WorkUnit *unit, int worker) {
    EStatus status;

    LockGuard lock(m_mutex);
    if (m_full.size() > 0) {
        unit->set(m_full.front());
        m_empty.push_back(m_full.front());
        m_full.pop_front();
        status = ESuccess;
    } else {
        status = m_done ? EFailure : EPause;
    }
    return status;
}

void RemoteProcess::processResult(const WorkResult *result, bool cancelled) {
    m_backend->sendWorkResult(m_id, result, cancelled);
}

ref<WorkProcessor> RemoteProcess::createWorkProcessor() const {
    return m_wp->clone();
}

/* Executed while the main scheduler lock is held. */
void RemoteProcess::handleCancellation() {
    /* Also acquire the local queue mutex, purge all queued
       work units and inform the remote side how many were lost */
    LockGuard lock(m_mutex);
    m_backend->sendCancellation(m_id, (int) m_full.size());
    m_empty.insert(m_empty.end(), m_full.begin(), m_full.end());
    m_full.clear();
}

MTS_IMPLEMENT_CLASS(RemoteWorker, false, Worker)
MTS_IMPLEMENT_CLASS(RemoteWorkerReader, false, Thread)
MTS_IMPLEMENT_CLASS(StreamBackend, false, Thread)
MTS_IMPLEMENT_CLASS(RemoteProcess, false, ParallelProcess)
MTS_NAMESPACE_END
