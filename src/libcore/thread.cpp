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

#include <mitsuba/core/lock.h>
#include <mitsuba/core/fresolver.h>
#include <errno.h>
#include <omp.h>

MTS_NAMESPACE_BEGIN


#if defined(_MSC_VER)
namespace
{
// Helper function to set a native thread name. MSDN:
//   http://msdn.microsoft.com/en-us/library/xcb2z8hs.aspx

const DWORD MS_VC_EXCEPTION=0x406D1388;

#pragma pack(push,8)
struct THREADNAME_INFO
{
	DWORD dwType;     // Must be 0x1000.
	LPCSTR szName;    // Pointer to name (in user addr space).
	DWORD dwThreadID; // Thread ID (-1=caller thread).
	DWORD dwFlags;    // Reserved for future use, must be zero.
};
#pragma pack(pop)

void SetThreadName(const char* threadName, DWORD dwThreadID = -1)
{
	THREADNAME_INFO info;
	info.dwType     = 0x1000;
	info.szName     = threadName;
	info.dwThreadID = dwThreadID;
	info.dwFlags    = 0;

	__try
	{
		RaiseException( MS_VC_EXCEPTION, 0, sizeof(info)/sizeof(ULONG_PTR),
			(ULONG_PTR*)&info );
	}
	__except(EXCEPTION_EXECUTE_HANDLER)
	{
	}
}


} // namespace
#endif // _MSC_VER


/**
 * Dummy class to associate a thread identity with the main thread
 */
class MainThread : public Thread {
public:
	MainThread() : Thread("main") {
	}

	virtual void run() {
		Log(EError, "The main thread is already running!");
	}

	MTS_DECLARE_CLASS()
protected:
	virtual ~MainThread() { }
};

class UnmanagedThread : public Thread {
public:
	UnmanagedThread(const std::string &name)
		: Thread(name) { }

	virtual void run() {
		Log(EError, "The unmanaged thread is already running!");
	}

	MTS_DECLARE_CLASS()
protected:
	virtual ~UnmanagedThread() { }
};


ThreadLocal<Thread> *Thread::m_self = NULL;

#if defined(__LINUX__) || defined(__OSX__)
int Thread::m_idCounter;
ref<Mutex> Thread::m_idMutex;
#endif

#if MTS_USE_ELF_TLS == 1
__thread int Thread::m_id 
	__attribute__((tls_model("global-dynamic")));
#endif

Thread::Thread(const std::string &name, unsigned int stackSize) 
 : m_name(name), m_stackSize(stackSize), m_running(false), m_joined(false),
   m_priority(ENormalPriority), m_critical(false) {
	m_joinMutex = new Mutex();
	memset(&m_thread, 0, sizeof(pthread_t));
}

void Thread::start() {
	if (m_running)
		Log(EError, "Thread is already running!");
	if (!m_self)
		Log(EError, "Threading has not been initialized!");

	Log(EDebug, "Spawning thread \"%s\"", m_name.c_str());

	/* Configure thread attributes */
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	if (m_stackSize != 0)
		pthread_attr_setstacksize(&attr, m_stackSize);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	m_parent = Thread::getThread();

	/* Inherit the parent thread's logger if none was set */
	if (!m_logger)
		m_logger = m_parent->getLogger();

	/* Inherit the parent thread's file resolver if none was set */
	if (!m_fresolver)
		m_fresolver = m_parent->getFileResolver();

	m_running = true;
	m_joined = false;

	incRef();
	if (pthread_create(&m_thread, &attr, &Thread::dispatch, this))
		Log(EError, "Could not create thread!");
}

bool Thread::setPriority(EThreadPriority priority) {
	m_priority = priority;
	if (!m_running)
		return true;

	Float factor;
	switch (priority) {
		case EIdlePriority: factor = 0.0f; break;
		case ELowestPriority: factor = 0.2f; break;
		case ELowPriority: factor = 0.4f; break;
		case EHighPriority: factor = 0.6f; break;
		case EHighestPriority: factor = 0.8f; break;
		case ERealtimePriority: factor = 1.0f; break;
		default: factor = 0.0f; break;
	}

	struct sched_param param;
	int policy;
	int retval = pthread_getschedparam(m_thread, &policy, &param);
	if (retval) {
		Log(EWarn, "Error in pthread_getschedparam(): %s!", strerror(retval));
		return false;
	}

	int min = sched_get_priority_min(policy);
	int max = sched_get_priority_max(policy);

	if (min == max) {
		Log(EWarn, "Could not adjust the thread priority -- valid range is zero!");
		return false;
	}

	param.sched_priority = (int) (min + (max-min)*factor);

	retval = pthread_setschedparam(m_thread, policy, &param);
	if (retval) {
		Log(EWarn, "Could not adjust the thread priority to %i: %s!", param.sched_priority, strerror(retval));
		return false;
	}
	return true;
}

void *Thread::dispatch(void *par) {
	Thread *thread = static_cast<Thread *>(par);
	Thread::m_self->set(thread);

	if (thread->getPriority() != ENormalPriority)
		thread->setPriority(thread->getPriority());

#if MTS_USE_ELF_TLS == 1
	m_idMutex->lock();
	m_id = ++m_idCounter;
	m_idMutex->unlock();
#elif defined(__LINUX__) || defined(__OSX__)
	m_idMutex->lock();
	thread->m_id = ++m_idCounter;
	m_idMutex->unlock();
#elif defined(_MSC_VER)
	if (!thread->getName().empty()) {
		const std::string threadName = "Mitsuba: " + thread->getName();
		SetThreadName(threadName.c_str());
	}
#endif

	try {
		thread->run();
	} catch (std::exception &e) {
		ELogLevel warnLogLevel = thread->getLogger()->getErrorLevel() == EError
			? EWarn : EInfo;
		Log(warnLogLevel, "Fatal error: uncaught exception: \"%s\"", e.what());
		if (thread->m_critical)
			_exit(-1);
	} catch (...) {
		ELogLevel warnLogLevel = thread->getLogger()->getErrorLevel() == EError
			? EWarn : EInfo;
		Log(warnLogLevel, "Fatal error - uncaught exception (unknown type)");
		if (thread->m_critical)
			_exit(-1);
	}

	thread->exit();

	return NULL;
}


void Thread::join() {
	/* Only one call to join() */
	m_joinMutex->lock();
	if (m_joined) {
		m_joinMutex->unlock();
		return;
	}
	m_joinMutex->unlock();
	int retval = pthread_join(m_thread, NULL);
	switch (retval) {
		case EINVAL: Log(EError, "Thread::join() - the thread is not joinable");
		case ESRCH:
			Log(EWarn, "Thread::join() - the thread does not exist");
			break;
		case EDEADLK: Log(EError, "Thread::join() - a deadlock was detected / "
			"thread tried to join on itself");
		default: break;
	}
	m_joined = true;
}

void Thread::detach() {
	int retval = pthread_detach(m_thread);
	switch (retval) {
		case EINVAL: Log(EError, "Thread::detach() - the thread is not joinable");
		case ESRCH: Log(EError, "Thread::detach() - the thread does not exist");
		default: break;
	}
}

void Thread::sleep(unsigned int ms) {
#if defined(WIN32)
	SleepEx(ms, FALSE);
#else
	struct timeval tv;
	tv.tv_sec = ms/1000;
	tv.tv_usec = (ms * 1000) % 1000000;
	select(0, NULL, NULL, NULL, &tv);
#endif
}

void Thread::yield() {
#if defined(WIN32)
	Yield();
#else
	sched_yield();
#endif
}

void Thread::exit() {
	Log(EDebug, "Thread \"%s\" has finished", m_name.c_str());
	m_running = false;
	decRef();
	m_self->set(NULL);
	pthread_exit(NULL);
}

std::string Thread::toString() const {
	std::ostringstream oss;
	oss << "Thread[" << endl
		<< "  name = \"" << m_name << "\"," << endl
		<< "  running = " << m_running << "," << endl
		<< "  joined = " << m_joined << "," << endl
		<< "  priority = " << m_priority << "," << endl
		<< "  critical = " << m_critical << "," << endl
		<< "  stackSize = " << m_stackSize << endl
		<< "]";
	return oss.str();
}

void Thread::staticInitialization() {
#if defined(__OSX__)
	__mts_autorelease_init();
#endif

	m_self = new ThreadLocal<Thread>();
	Thread *mainThread = new MainThread();
#if defined(__LINUX__) || defined(__OSX__)
	m_idMutex = new Mutex();
	m_idCounter = 0;

	#if MTS_USE_ELF_TLS == 1
		m_id = 0;
	#else
		mainThread->m_id = 0;
	#endif
#endif
	mainThread->m_running = true;
	mainThread->m_thread = pthread_self();
	mainThread->m_joinMutex = new Mutex();
	mainThread->m_joined = false;
	mainThread->m_fresolver = new FileResolver();
	m_self->set(mainThread);

}

static std::vector<Thread *> __unmanagedThreads;

Thread *Thread::registerUnmanagedThread(const std::string &name) {
	Thread *thread = getThread();
	if (!thread) {
		thread = new UnmanagedThread(name);
		thread->m_running = false;
		thread->m_thread = pthread_self();
		thread->m_joinMutex = new Mutex();
		thread->m_joined = false;
		thread->incRef();
		m_self->set(thread);
		#pragma omp critical
			__unmanagedThreads.push_back((UnmanagedThread *) thread);
	}
	return thread;
}

void Thread::staticShutdown() {
	for (size_t i=0; i<__unmanagedThreads.size(); ++i)
		__unmanagedThreads[i]->decRef();
	__unmanagedThreads.clear();
	getThread()->m_running = false;
	m_self->set(NULL);
	delete m_self;
	m_self = NULL;
#if defined(__LINUX__) || defined(__OSX__)
	m_idMutex = NULL;
#endif
#if defined(__OSX__)
	__mts_autorelease_shutdown();
#endif
}

#if defined(__OSX__)
PrimitiveThreadLocal<int> __threadID;
int __threadCount = 0;

int mts_get_thread_num() {
	return __threadID.get();
}
int mts_get_max_threads() {
	return __threadCount;
}
#else
int mts_get_thread_num() {
	return omp_get_thread_num();
}
int mts_get_max_threads() {
	return omp_get_max_threads();
}
#endif

void Thread::initializeOpenMP(size_t threadCount) {
	ref<Logger> logger = Thread::getThread()->getLogger();
	ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

#if defined(__OSX__)
	__threadCount = (int) threadCount;
#endif

	omp_set_num_threads((int) threadCount);
	int counter = 0;

	#pragma omp parallel
	{
		Thread *thread = Thread::getThread();
		if (!thread) {
			#pragma omp critical
			{
				thread = new UnmanagedThread(
					formatString("omp%i", counter));
				#if MTS_BROKEN_OPENMP == 1
				__threadID.set(counter);
				#endif
				counter++;
			}
			thread->m_running = false;
			thread->m_thread = pthread_self();
			thread->m_joinMutex = new Mutex();
			thread->m_joined = false;
			thread->m_fresolver = fResolver;
			thread->m_logger = logger;
			thread->incRef();
			m_self->set(thread);
			#pragma omp critical
			__unmanagedThreads.push_back((UnmanagedThread *) thread);
		}
	}
}

Thread::~Thread() {
	if (m_running)
		Log(EWarn, "Destructor called while Thread '%s' is still running", m_name.c_str());
}

MTS_IMPLEMENT_CLASS(Thread, true, Object)
MTS_IMPLEMENT_CLASS(MainThread, false, Thread)
MTS_IMPLEMENT_CLASS(UnmanagedThread, false, Thread)
MTS_NAMESPACE_END
