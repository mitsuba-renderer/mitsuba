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
#if !defined(__MITSUBA_CORE_LOGGER_H_)
#define __MITSUBA_CORE_LOGGER_H_

#include <mitsuba/core/formatter.h>

// -----------------------------------------------------------------------
//  Logging
// -----------------------------------------------------------------------

MTS_NAMESPACE_BEGIN

/*! \addtogroup libcore */
/*! @{ */

/// Write a Log message to the console (to be used within subclasses of <tt>Object</tt>)
#define Log(level, fmt, ...) do { \
		mitsuba::Thread *thread = mitsuba::Thread::getThread(); \
		if (EXPECT_NOT_TAKEN(thread == NULL)) \
			throw std::runtime_error("Null thread pointer"); \
		mitsuba::Logger *logger = thread->getLogger(); \
		if (logger != NULL && level >= logger->getLogLevel()) \
			logger->log(level, m_theClass, \
				__FILE__, __LINE__, fmt, ## __VA_ARGS__); \
	} while (0)

/**
 * \brief Write a Log message to the console (static version - to be used
 * outside of classes that derive from Object)
 */
#define SLog(level, fmt, ...) do { \
		mitsuba::Thread *thread = mitsuba::Thread::getThread(); \
		if (EXPECT_NOT_TAKEN(thread == NULL)) \
			throw std::runtime_error("Null thread pointer"); \
		mitsuba::Logger *logger = thread->getLogger(); \
		if (logger != NULL && level >= logger->getLogLevel()) \
			logger->log(level, NULL, \
				__FILE__, __LINE__, fmt, ## __VA_ARGS__); \
	} while (0)

/*! @} */

#ifdef MTS_NDEBUG
#define Assert(cond) ((void) 0)
#define AssertEx(cond, explanation) ((void) 0)
#define SAssert(cond) ((void) 0)
#define SAssertEx(cond, explanation) ((void) 0)
#else

/* Assertions */
/*! \addtogroup libcore */
/*! @{ */

/// Assert that a condition is true (to be used \a inside of classes that derive from \ref Object)
#define Assert(cond) do { \
		if (!(cond)) Log(EError, "Assertion \"%s\" failed in %s:%i", \
		#cond, __FILE__, __LINE__); \
	} while (0)

/// ``Static'' assertion (to be used \a outside of classes that derive from \ref Object)
#define SAssert(cond) do { \
		if (!(cond)) SLog(EError, "Assertion \"%s\" failed in %s:%i", \
		#cond, __FILE__, __LINE__); \
	} while (0)

/// Assertion with a customizable error explanation
#define AssertEx(cond, explanation) do { \
		if (!(cond)) Log(EError, "Assertion \"%s\" failed in %s:%i (" explanation ")", \
		#cond, __FILE__, __LINE__); \
	} while (0)

/// Static assertion with a customizable error explanation (see \ref SLog)
#define SAssertEx(cond, explanation) do { \
		if (!(cond)) SLog(EError, "Assertion \"%s\" failed in %s:%i (" explanation ")", \
		#cond, __FILE__, __LINE__); \
	} while (0)
#endif

/// Throw an exception reporting that the given function is not implemented
#define NotImplementedError(funcName) \
	throw std::runtime_error(formatString("%s::" funcName "(): Not implemented!", \
			getClass()->getName().c_str()));

/*! @} */

/**
 * \headerfile mitsuba/core/logger.h mitsuba/mitsuba.h
 * \brief Responsible for processing log messages
 *
 * Upon receiving a log message, the Logger class invokes
 * a Formatter to convert it into a human-readable form.
 * Following that, it sends this information to every
 * registered Appender.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Logger : public Object {
	friend class Statistics;
public:
	/// Construct a new logger with the given minimum log level
	Logger(ELogLevel logLevel = EDebug);

	/**
	 * \brief Process a log message
	 * \param level Log level of the message
	 * \param theClass Class descriptor of the message creator
	 * \param fileName Source file of the message creator
	 * \param lineNumber Source line number of the message creator
	 * \param fmt printf-style string formatter
	 * \note This function is not exposed in the Python bindings.
	 *       Instead, please use \cc mitsuba.core.Log
	 */
	void log(ELogLevel level, const Class *theClass,
		const char *fileName, int lineNumber,
		const char *fmt, ...);

	/**
	 * \brief Process a progress message
	 * \param progress Percentage value in [0,100]
	 * \param name Title of the progress message
	 * \param formatted Formatted string representation of the message
	 * \param eta Estimated time until 100% is reached.
	 * \param ptr Custom pointer payload. This is used to express the
	 *    context of a progress message. When rendering a scene, it
	 *    will usually contain a pointer to the associated \c RenderJob.
	 * \remark The \c ptr argument is missing in the Python bindings
	 */
	void logProgress(Float progress, const std::string &name,
		const std::string &formatted, const std::string &eta,
		const void *ptr);

	/// Set the log level (everything below will be ignored)
	void setLogLevel(ELogLevel level);

	/**
	 * \brief Set the error log level (this level and anything
	 * above will throw exceptions).
	 *
	 * The value provided here can be used for instance to turn
	 * warnings into errors. But \a level must always be
	 * less than \ref EError, i.e. it isn't possible to
	 * cause errors not to throw an exception.
	 */
	void setErrorLevel(ELogLevel level);

	/// Return the current log level
	inline ELogLevel getLogLevel() const { return m_logLevel; }

	/// Return the current error level
	inline ELogLevel getErrorLevel() const { return m_errorLevel; }

	/// Add an appender to this logger
	void addAppender(Appender *appender);

	/// Remove an appender from this logger
	void removeAppender(Appender *appender);

	/// Remove all appenders from this logger
	void clearAppenders();

	/// Return the number of registered appenders
	inline size_t getAppenderCount() const { return m_appenders.size(); }

	/// Return one of the appenders
	inline Appender *getAppender(size_t index) { return m_appenders[index]; }

	/// Return one of the appenders
	inline const Appender *getAppender(size_t index) const { return m_appenders[index]; }

	/// Set the logger's formatter implementation
	void setFormatter(Formatter *formatter);

	/// Return the logger's formatter implementation
	inline Formatter *getFormatter() { return m_formatter; }

	/**
	 * \brief Return the contents of the log file as a string (if it exists)
	 *
	 * \return \c true upon success
	 */
	bool readLog(std::string &target);

	/// Return the number of warnings reported so far
	inline size_t getWarningCount() const { return m_warningCount; }

	/// Initialize logging
	static void staticInitialization();

	/// Shutdown logging
	static void staticShutdown();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Logger();
private:
	ELogLevel m_logLevel;
	ELogLevel m_errorLevel;
	ref<Formatter> m_formatter;
	ref<Mutex> m_mutex;
	std::vector<Appender *> m_appenders;
	size_t m_warningCount;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_LOGGER_H_ */
