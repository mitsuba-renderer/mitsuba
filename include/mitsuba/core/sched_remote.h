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
#if !defined(__MITSUBA_CORE_SCHED_REMOTE_H_)
#define __MITSUBA_CORE_SCHED_REMOTE_H_

#include <mitsuba/core/sched.h>
#include <set>

/// Default port of <tt>mtssrv</tt>
#define MTS_DEFAULT_PORT 7554

/** How many work units should be sent to a remote worker
   at a time? This is a multiple of the worker's core count */
#define MTS_BACKLOG_FACTOR 3

/** Once the back log factor drops below this value (also a
   multiple of the core size), the stream processor will
   continue sending batches of work units */
#define MTS_CONTINUE_FACTOR 2

MTS_NAMESPACE_BEGIN

class RemoteWorkerReader;
class StreamBackend;

/**
 * \brief Acquires work from the scheduler and forwards
 * it to a processing node reachable through a \ref Stream.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE RemoteWorker : public Worker {
	friend class RemoteWorkerReader;
public:
	/**
	 * \brief Construct a new remote worker with the given name and
	 * communication stream
	 */
	RemoteWorker(const std::string &name, Stream *stream);

	/// Return the name of the node on the other side
	inline const std::string &getNodeName() const { return m_nodeName; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~RemoteWorker();
	/* Worker implementation */
	virtual void run();
	virtual void clear();
	virtual void signalResourceExpiration(int id);
	virtual void signalProcessCancellation(int id);
	virtual void signalProcessTermination(int id);
	virtual void start(Scheduler *scheduler, int workerIndex, int coreOffset);
	void flush();

	inline void signalCompletion() {
		LockGuard lock(m_mutex);
		m_inFlight--;
		m_finishCond->signal();
	}
protected:
	ref<Mutex> m_mutex;
	ref<ConditionVariable> m_finishCond;
	ref<MemoryStream> m_memStream;
	ref<Stream> m_stream;
	ref<RemoteWorkerReader> m_reader;

	/* List of processes and resources that are
	   currently active at the remote node */
	std::set<int> m_resources;
	std::set<int> m_processes;
	std::set<std::string> m_plugins;
	std::string m_nodeName;
	size_t m_inFlight;
};

/**
 * \brief Communication helper thread required by \ref RemoteWorker.
 *
 * Constantly waits for finished work units sent by the processing node.
 */
class MTS_EXPORT_CORE RemoteWorkerReader : public Thread {
	friend class RemoteWorker;
public:
	RemoteWorkerReader(RemoteWorker *parent);

	inline void shutdown() { m_shutdown = true; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~RemoteWorkerReader() { }
	/// Thread body
	void run();
private:
	std::vector<Thread *> m_joinThreads;
	RemoteWorker *m_parent;
	ref<Stream> m_stream;
	bool m_shutdown;
	int m_currentID;
	Scheduler::Item m_schedItem;
};

/**
 * \brief Parallel process facade used to insert work units from a
 * remote scheduler into the local one.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE RemoteProcess : public ParallelProcess {
public:
	/**
	 * \brief Create a new remote process
	 *
	 * \param id        Identification number for this process
	 * \param logLevel  Log level for events associated with this process
	 * \param backend   The responsible server-side communication backend
	 * \param proc      Work processor instance for use with this process
	 */
	RemoteProcess(int id, ELogLevel logLevel,
		StreamBackend *backend, WorkProcessor *proc);

	/* ParallelProcess interface implementation */
	EStatus generateWork(WorkUnit *unit, int worker);
	void processResult(const WorkResult *result,
		bool cancelled);
	ref<WorkProcessor> createWorkProcessor() const;
	void handleCancellation();

	/// Get an empty work unit from the process (or create one)
	inline WorkUnit *getEmptyWorkUnit() {
		ref<WorkUnit> wu;
		LockGuard lock(m_mutex);
		if (m_empty.empty()) {
			wu = m_wp->createWorkUnit();
			wu->incRef();
		} else {
			wu = m_empty.back();
			m_empty.pop_back();
		}
		return wu;
	}

	/// Make a full work unit available to the process
	inline void putFullWorkUnit(WorkUnit *wu) {
		LockGuard lock(m_mutex);
		m_full.push_back(wu);
	}

	/// Mark the process as finished
	inline void setDone() {
		LockGuard lock(m_mutex);
		m_done = true;
	}

	MTS_DECLARE_CLASS()
protected:
	// Virtual destructor
	virtual ~RemoteProcess();
private:
	int m_id;
	ref<StreamBackend> m_backend;
	std::vector<WorkUnit *> m_empty;
	std::deque<WorkUnit *> m_full;
	ref<WorkProcessor> m_wp;
	ref<Mutex> m_mutex;
	bool m_done;
};

/**
 * \brief Network processing communication backend
 *
 * Attaches to the end of a stream, accepts work units and forwards
 * them to the local scheduler. Can be used to create network processing nodes.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE StreamBackend : public Thread {
	friend class RemoteProcess;
	friend class RemoteWorker;
	friend class RemoteWorkerReader;
public:
	/**
	 * \brief Create a new stream backend
	 *
	 * \param name
	 *    Name of the created thread
	 * \param scheduler
	 *    Scheduler instance used to process work units
	 * \param nodeName
	 *    Exposed name of this node
	 * \param stream
	 *    Stream used for communications
	 * \param detach
	 *    Should the associated thread be joinable or detach instead?
	 */
	StreamBackend(const std::string &name, Scheduler *scheduler,
		const std::string &nodeName, Stream *stream, bool detach);

	MTS_DECLARE_CLASS()
protected:
	enum EMessage {
		EUnknown = 0,
		ENewProcess,
		ENewResource,
		ENewMultiResource,
		EBindResource,
		EWorkUnit,
		EWorkResult,
		ECancelledWorkResult,
		EProcessTerminated,
		EProcessCancelled,
		EEnsurePluginLoaded,
		EResourceExpired,
		EQuit,
		EIncompatible,
		EHello = 0x1bcd
	};

	/// Virtual destructor
	virtual ~StreamBackend();
	virtual void run();
	void sendWorkResult(int id, const WorkResult *result, bool cancelled);
	void sendCancellation(int id, int numLost);
private:
	Scheduler *m_scheduler;
	std::string m_nodeName;
	ref<Stream> m_stream;
	ref<MemoryStream> m_memStream;
	std::map<int, RemoteProcess *> m_processes;
	std::map<int, int> m_resources;
	ref<Mutex> m_sendMutex;
	bool m_detach;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SCHED_REMOTE_H_ */
