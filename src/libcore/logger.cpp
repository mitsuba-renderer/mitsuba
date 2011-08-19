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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/lock.h>
#include <stdarg.h>

#if defined(__OSX__)
#include <sys/sysctl.h>
#endif

MTS_NAMESPACE_BEGIN

Logger::Logger(ELogLevel level)
 : m_logLevel(level), m_errorLevel(EError), m_warningCount(0) {
	m_mutex = new Mutex();
}

Logger::~Logger() {
	for (size_t i=0; i<m_appenders.size(); ++i)
		m_appenders[i]->decRef();
}

void Logger::setFormatter(Formatter *formatter) {
	m_mutex->lock();
	m_formatter = formatter;
	m_mutex->unlock();
}

void Logger::setLogLevel(ELogLevel level) {
	m_logLevel = level;
}

void Logger::setErrorLevel(ELogLevel level) {
	Assert(m_errorLevel <= EError);
	m_errorLevel = level;
}

void Logger::log(ELogLevel level, const Class *theClass, 
	const char *file, int line, const char *fmt, ...) {

	if (level < m_logLevel)
		return;

	char tmp[512], *msg = tmp;
	va_list iterator;

#if defined(WIN32)
	va_start(iterator, fmt);
	size_t size = _vscprintf(fmt, iterator) + 1;

	if (size >= sizeof(tmp)) 
		msg = new char[size];

	vsnprintf_s(msg, size, size-1, fmt, iterator);
	va_end(iterator);
#else
	va_start(iterator, fmt);
	size_t size = vsnprintf(tmp, sizeof(tmp), fmt, iterator);
	va_end(iterator);

	if (size >= sizeof(tmp)) {
		/* Overflow! -- dynamically allocate memory */
		msg = new char[size+1];
		va_start(iterator, fmt);
		vsnprintf(msg, size+1, fmt, iterator);
		va_end(iterator);
	}
#endif

	if (m_formatter == NULL) {
		std::cerr << "PANIC: Logging has not been properly initialized!" << std::endl;
		exit(-1);
	}

	std::string text = m_formatter->format(level, theClass, 
		Thread::getThread(), msg, file, line);

	if (msg != tmp)
		delete[] msg;

	if (level < m_errorLevel) {
		m_mutex->lock();
		if (level >= EWarn)
			m_warningCount++;
		for (size_t i=0; i<m_appenders.size(); ++i)
			m_appenders[i]->append(level, text);
		m_mutex->unlock();
	} else {
#if defined(__LINUX__)
		/* A critical error occurred: trap if we're running in a debugger */
		
		char exePath[PATH_MAX];
		pid_t ppid = getppid();
		memset(exePath, 0, PATH_MAX);
		if (readlink(formatString("/proc/%i/exe", ppid).c_str(), exePath, PATH_MAX) != -1) {
			if (!strcmp(exePath, "/usr/bin/gdb")) {
				__asm__ ("int $3");
			}
		}
#elif defined(__OSX__)
		int                 mib[4];
		struct kinfo_proc   info;
		size_t              size;
		info.kp_proc.p_flag = 0;
		mib[0] = CTL_KERN;
		mib[1] = KERN_PROC;
		mib[2] = KERN_PROC_PID;
		mib[3] = getpid();
		size = sizeof(info);
		sysctl(mib, sizeof(mib) / sizeof(*mib), &info, &size, NULL, 0);
		bool runningInDebugger = (info.kp_proc.p_flag & P_TRACED) != 0;

		if (runningInDebugger)
			__asm__ ("int $3");
#elif defined(WIN32)
		if (IsDebuggerPresent()) 
			__debugbreak();
#endif

		throw std::runtime_error(text);
	}
}

void Logger::logProgress(Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta, const void *ptr) {
	m_mutex->lock();
	for (size_t i=0; i<m_appenders.size(); ++i)
		m_appenders[i]->logProgress(
			progress, name, formatted, eta, ptr);
	m_mutex->unlock();
}

void Logger::addAppender(Appender *appender) {
	appender->incRef();
	m_mutex->lock();
	m_appenders.push_back(appender);
	m_mutex->unlock();
}

void Logger::removeAppender(Appender *appender) {
	m_mutex->lock();
	m_appenders.erase(std::remove(m_appenders.begin(), 
		m_appenders.end(), appender), m_appenders.end());
	m_mutex->unlock();
	appender->decRef();
}

void Logger::clearAppenders() {
	m_mutex->lock();
	for (size_t i=0; i<m_appenders.size(); ++i)
		m_appenders[i]->decRef();
	m_appenders.clear();
	m_mutex->unlock();
}

void Logger::staticInitialization() {
	Logger *logger = new Logger(EInfo);
	ref<Appender> appender = new StreamAppender(&std::cout);
	ref<Formatter> formatter = new DefaultFormatter();
	logger->addAppender(appender);
	logger->setFormatter(formatter);
	Thread::getThread()->setLogger(logger);
}

void Logger::staticShutdown() {
	Thread::getThread()->setLogger(NULL);
}

MTS_IMPLEMENT_CLASS(Logger, false, Object)
MTS_NAMESPACE_END
