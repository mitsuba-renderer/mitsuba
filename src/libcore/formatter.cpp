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
#include <mitsuba/core/thread.h>
#include <boost/filesystem.hpp>
#include <ctime>

MTS_NAMESPACE_BEGIN

DefaultFormatter::DefaultFormatter()
 : m_haveDate(true), m_haveLogLevel(true), m_haveThread(true), m_haveClass(true) {
}

std::string DefaultFormatter::format(ELogLevel logLevel, const Class *theClass,
        const Thread *thread, const std::string &text, const char *file, int line) {
    std::ostringstream oss;
    char buffer[128];

    /* Date/Time */
    if (m_haveDate) {
        time_t theTime = std::time(NULL);
        strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S ", std::localtime(&theTime));
        oss << buffer;
    }

    /* Log level */
    if (m_haveLogLevel) {
        switch (logLevel) {
            case ETrace: oss << "TRACE "; break;
            case EDebug: oss << "DEBUG "; break;
            case EInfo:  oss << "INFO  "; break;
            case EWarn:  oss << "WARN  "; break;
            case EError: oss << "ERROR "; break;
            default:     oss << "CUSTM "; break;
        }
    }

    /* Thread */
    if (thread && m_haveThread) {
        oss << thread->getName();

        for (int i=0; i<(5 - (int) thread->getName().size()); i++)
            oss << ' ';
    }

    /* Class */
    if (m_haveClass) {
        if (theClass)
            oss << "[" << theClass->getName() << "] ";
        else if (line != -1 && file)
            oss << "[" << fs::path(file).filename().string() << ":" << line << "] ";
    }

    /* Text */
    oss << text;

    return oss.str();
}

MTS_IMPLEMENT_CLASS(Formatter, true, Object)
MTS_IMPLEMENT_CLASS(DefaultFormatter, false, Formatter)
MTS_NAMESPACE_END
