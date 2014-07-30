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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/lock.h>
#include <stdarg.h>

#if defined(__OSX__)
# include <sys/sysctl.h>
#elif defined(__WINDOWS__)
# include <windows.h>
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
	LockGuard lock(m_mutex);
	m_formatter = formatter;
}

void Logger::setLogLevel(ELogLevel level) {
	LockGuard lock(m_mutex);
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

#if defined(__WINDOWS__)
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
		std::cerr << "PANIC: Logging has not been properly initialized!" << endl;
		exit(-1);
	}

	std::string text = m_formatter->format(level, theClass,
		Thread::getThread(), msg, file, line);

	if (msg != tmp)
		delete[] msg;

	if (level < m_errorLevel) {
		LockGuard lock(m_mutex);
		if (level >= EWarn)
			m_warningCount++;
		for (size_t i=0; i<m_appenders.size(); ++i)
			m_appenders[i]->append(level, text);
	} else {
#if defined(__LINUX__)
		/* A critical error occurred: trap if we're running in a debugger */

		char exePath[PATH_MAX];
		pid_t ppid = getppid();
		memset(exePath, 0, PATH_MAX);
		if (readlink(formatString("/proc/%i/exe", ppid).c_str(), exePath, PATH_MAX) != -1) {
			if (!strcmp(exePath, "/usr/bin/gdb")) {
#if defined(__i386__) || defined(__x86_64__)
				__asm__ ("int $3");
#else
				__builtin_trap();
#endif
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
#elif defined(__WINDOWS__)
		if (IsDebuggerPresent())
			__debugbreak();
#endif

		DefaultFormatter fmt;
		fmt.setHaveDate(false);
		fmt.setHaveLogLevel(false);
		text = fmt.format(level, theClass,
			Thread::getThread(), msg, file, line);
		throw std::runtime_error(text);
	}
}

void Logger::logProgress(Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta, const void *ptr) {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_appenders.size(); ++i)
		m_appenders[i]->logProgress(
			progress, name, formatted, eta, ptr);
}

void Logger::addAppender(Appender *appender) {
	appender->incRef();
	LockGuard lock(m_mutex);
	m_appenders.push_back(appender);
}

void Logger::removeAppender(Appender *appender) {
	LockGuard lock(m_mutex);
	m_appenders.erase(std::remove(m_appenders.begin(),
		m_appenders.end(), appender), m_appenders.end());
	appender->decRef();
}

bool Logger::readLog(std::string &target) {
	bool success = false;
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_appenders.size(); ++i) {
		Appender *appender = m_appenders[i];
		if (appender->getClass()->derivesFrom(MTS_CLASS(StreamAppender))) {
			StreamAppender *streamAppender =
				static_cast<StreamAppender *>(appender);
			if (streamAppender->logsToFile()) {
				streamAppender->readLog(target);
				success = true;
				break;
			}
		}
	}
	return success;
}

void Logger::clearAppenders() {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_appenders.size(); ++i)
		m_appenders[i]->decRef();
	m_appenders.clear();
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
