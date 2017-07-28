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
#if !defined(__MITSUBA_CORE_FORMATTER_H_)
#define __MITSUBA_CORE_FORMATTER_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/// Available Log message types
enum ELogLevel {
    ETrace = 0,   ///< Trace message, for extremely verbose debugging
    EDebug = 100, ///< Debug message, usually turned off
    EInfo = 200,  ///< More relevant debug / information message
    EWarn = 300,  ///< Warning message
    EError = 400  ///< Error message, causes an exception to be thrown
};

/** \brief Abstract interface for converting log information into
 * a human-readable format
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Formatter : public Object {
public:
    /**
     * \brief Turn a log message into a human-readable format
     * \param logLevel  The importance of the debug message
     * \param theClass  Originating class or NULL
     * \param thread    Thread, which is reponsible for creating the message
     * \param text      Text content associated with the log message
     * \param file      File, which is responsible for creating the message
     * \param line      Associated line within the source file
     */

    virtual std::string format(ELogLevel logLevel, const Class *theClass,
            const Thread *thread, const std::string &text,
            const char *file, int line) = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Formatter() { }
};

/** \brief The default formatter used to turn log messages into
 * a human-readable form
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE DefaultFormatter : public Formatter {
    friend class Logger;
public:
    /// Create a new default formatter
    DefaultFormatter();

    std::string format(ELogLevel logLevel, const Class *theClass,
            const Thread *thread, const std::string &text,
            const char *file, int line);

    /// Should date information be included? The default is yes.
    inline void setHaveDate(bool value) { m_haveDate = value; }

    /// Should thread information be included? The default is yes.
    inline void setHaveThread(bool value) { m_haveThread = value; }

    /// Should log level information be included? The default is yes.
    inline void setHaveLogLevel(bool value) { m_haveLogLevel = value; }

    /// Should class information be included? The default is yes.
    inline void setHaveClass(bool value) { m_haveClass = value; }

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~DefaultFormatter() { }
protected:
    bool m_haveDate;
    bool m_haveLogLevel;
    bool m_haveThread;
    bool m_haveClass;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_FORMATTER_H_ */
