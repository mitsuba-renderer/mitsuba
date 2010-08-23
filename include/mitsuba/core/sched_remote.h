#if !defined(__SCHED_REMOTE_H)
#define __SCHED_REMOTE_H

#include <mitsuba/core/sched.h>

/* How many work units should be sent to a remote worker
   at a time? This is a multiple of the worker's core count */
#define BACKLOG_FACTOR 3

/* Once the back log factor drops below this value (also a
   multiple of the core size), the stream processor will
   continue sending batches of work units */
#define CONTINUE_FACTOR 2

MTS_NAMESPACE_BEGIN

class RemoteWorkerReader;
class StreamBackend;
class FileResolver;

/**
 * Remote worker thread. Acquires work from the scheduler and forwards
 * it to a processing node reachable over some form of stream (usually
 * a <tt>SocketStream</tt>).
 */
class MTS_EXPORT_CORE RemoteWorker : public Worker {
	friend class RemoteWorkerReader;
public:
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
		m_mutex->lock();
		m_inFlight--;
		m_finishCond->signal();
		m_mutex->unlock();
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
 * Remote worker helper thread - constantly waits for finished
 * work units sent by the processing node.
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
 * 'Fake' parallel process used to insert work units from a
 * remote scheduler into the local one.
 */
class MTS_EXPORT_CORE RemoteProcess : public ParallelProcess {
public:
	RemoteProcess(int id, ELogLevel logLevel, 
		StreamBackend *backend, WorkProcessor *proc);

	/* ParallelProcess interface implementation */
	EStatus generateWork(WorkUnit *unit, int worker);
	void processResult(const WorkResult *result,	
		bool cancelled);
	ref<WorkProcessor> createWorkProcessor() const;
	void handleCancellation();

	inline WorkUnit *getEmptyWorkUnit() {
		ref<WorkUnit> wu;
		m_mutex->lock();
		if (m_empty.empty()) {
			wu = m_wp->createWorkUnit();
			wu->incRef();
		} else {
			wu = m_empty.back();
			m_empty.pop_back();
		}
		m_mutex->unlock();
		return wu;
	}

	inline void putFullWorkUnit(WorkUnit *wu) {
		m_mutex->lock();
		m_full.push_back(wu);
		m_mutex->unlock();
	}

	inline void setDone() {
		m_mutex->lock();
		m_done = true;
		m_mutex->unlock();
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
 * Stream backend - attaches to the end of a stream, accepts work units
 * and forwards them to the local scheduler. Can be used to create
 * remote processing nodes.
 */
class MTS_EXPORT_CORE StreamBackend : public Thread {
	friend class RemoteProcess;
public:
	enum EMessage {
		EUnknown = 0,
		ENewProcess,
		ENewResource,
		ENewManifoldResource,
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

	/**
	 * Create a new stream backend
	 *
	 * @param name
	 *    Name of the created thread
	 * @param scheduler
	 *    Scheduler instance used to process work units
	 * @param nodeName
	 *    Exposed name of this node
	 * @param stream
	 *    Stream used for communications
	 * @param detach
	 *    Should the associated thread be joinable or detach instead?
	 */
	StreamBackend(const std::string &name, Scheduler *scheduler, 
		const std::string &nodeName, Stream *stream, bool detach);

	MTS_DECLARE_CLASS()
protected:
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
	ref<FileResolver> m_resolver;
	bool m_detach;
};

MTS_NAMESPACE_END

#endif /* __SCHED_REMOTE_H */
