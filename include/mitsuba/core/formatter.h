/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__FORMATTER_H)
#define __FORMATTER_H

MTS_NAMESPACE_BEGIN

/// Available Log message types
enum ELogLevel {
	/// Trace message, for extremely verbose debugging
	ETrace = 0,
	/// Debug message, usually turned off
	EDebug = 100,
	/// More relevant debug / information message
	EInfo = 200,
	/// Warning message
	EWarn = 300,
	/// Error message, causes an exception to be thrown
	EError = 400
};

class Thread;

/** \brief The Formatter class defines an interface for converting
 * log information into a human-readable format
 */
class MTS_EXPORT_CORE Formatter : public Object {
public:
	/// Format a line of text
	virtual std::string format(ELogLevel pLogLevel, const Class *pClass,
			const Thread *pThread, const std::string &pText, 
			const char *pFile, int pLine) = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Formatter() { }
};

/** \brief The default formatter implementation
 */
class MTS_EXPORT_CORE DefaultFormatter : public Formatter {
public:
	/// Create a new default formatter
	DefaultFormatter();

	/// Format a line of text
	std::string format(ELogLevel pLogLevel, const Class *pClass,
			const Thread *pThread, const std::string &pText, 
			const char *pFile, int pLine);

	/// Should the date be included?
	inline void setHaveDate(bool value) { m_haveDate = value; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~DefaultFormatter() { }
protected:
	bool m_haveDate;
};

MTS_NAMESPACE_END

#endif /* __FORMATTER_H */
