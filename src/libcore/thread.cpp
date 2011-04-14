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

class OpenMPThread : public Thread {
public:
	OpenMPThread(int threadIdx)
		: Thread(formatString("omp%i", threadIdx)) { }

	virtual void run() {
		Log(EError, "The OpenMP thread is already running!");
	}

	MTS_DECLARE_CLASS()
protected:
	virtual ~OpenMPThread() { }
};


ThreadLocal<Thread> *Thread::m_self = NULL;

#if defined(__LINUX__) || defined(__OSX__)
int Thread::m_idCounter;
ref<Mutex> Thread::m_idMutex;
#endif

#if defined(__LINUX__)
int __thread Thread::m_id;
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

#if defined(__OSX__)
	m_idMutex->lock();
	thread->m_id = ++m_idCounter;
	m_idMutex->unlock();
#elif defined(__LINUX__)
	m_idMutex->lock();
	m_id = ++m_idCounter;
	m_idMutex->unlock();
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
		<< "  name=\"" << m_name << "\"," << endl
		<< "  running=" << m_running << "," << endl
		<< "  joined=" << m_joined << "," << endl
		<< "  priority=" << m_priority << "," << endl
		<< "  critical=" << m_critical << "," << endl
		<< "  stackSize=" << m_stackSize << endl
		<< "]";
	return oss.str();
}

void Thread::staticInitialization() {
#if defined(__OSX__)
	__ubi_autorelease_init();
#endif

	m_self = new ThreadLocal<Thread>();
#if defined(__LINUX__) || defined(__OSX__)
	m_idMutex = new Mutex();
	m_idCounter = 0;
#endif
	Thread *mainThread = new MainThread();
	mainThread->m_running = true;
	mainThread->m_thread = pthread_self();
	mainThread->m_joinMutex = new Mutex();
	mainThread->m_joined = false;
	mainThread->m_fresolver = new FileResolver();
	m_self->set(mainThread);
#if defined(__OSX__)
	mainThread->m_id = 0;
#elif defined(__LINUX__)
	m_id = 0;
#endif
}

static std::vector<OpenMPThread *> __ompThreads;

void Thread::staticShutdown() {
	for (size_t i=0; i<__ompThreads.size(); ++i)
		__ompThreads[i]->decRef();
	__ompThreads.clear();
	getThread()->m_running = false;
	m_self->set(NULL);
	delete m_self;
	m_self = NULL;
#if defined(__LINUX__) || defined(__OSX__)
	m_idMutex = NULL;
#endif
#if defined(__OSX__)
	__ubi_autorelease_shutdown();
#endif
}

void Thread::initializeOpenMP(size_t threadCount) {
	ref<Logger> logger = Thread::getThread()->getLogger();
	ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

	omp_set_num_threads(threadCount);

	#pragma omp parallel
	{
		Thread *thread = Thread::getThread();
		if (!thread) {
			thread = new OpenMPThread(omp_get_thread_num());
			thread->m_running = true;
			thread->m_thread = pthread_self();
			thread->m_joinMutex = new Mutex();
			thread->m_joined = false;
			thread->m_fresolver = fResolver;
			thread->m_logger = logger;
			thread->incRef();
			m_self->set(thread);
			#pragma omp critical
			__ompThreads.push_back((OpenMPThread *) thread);
		}
	}
}

Thread::~Thread() {
	if (m_running)
		Log(EWarn, "Destructor called while Thread '%s' is still running", m_name.c_str());
}

MTS_IMPLEMENT_CLASS(Thread, true, Object)
MTS_IMPLEMENT_CLASS(MainThread, false, Thread)
MTS_IMPLEMENT_CLASS(OpenMPThread, false, Thread)
MTS_NAMESPACE_END
