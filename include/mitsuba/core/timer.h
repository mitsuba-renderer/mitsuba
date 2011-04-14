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

#if !defined(__TIMER_H)
#define __TIMER_H

#include <mitsuba/mitsuba.h>

#ifndef WIN32
#include <sys/time.h>
#endif

MTS_NAMESPACE_BEGIN

/** \brief Platform independent milli/microsecond timer
 * \ingroup libcore
 */
class MTS_EXPORT_CORE Timer : public Object {
public:
	/// Create a new timer and reset it
	Timer();

	/// Reset the timer
	void reset();

	/// Return the milliseconds which have passed since the last reset
	unsigned int getMilliseconds() const;

	/// Return the microseconds which have passed since the last reset
	unsigned int getMicroseconds() const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Timer();
private:
#ifdef WIN32
	__int64 m_start;
	__int64 m_frequency;
#else
	struct timeval m_start;
#endif
};

MTS_NAMESPACE_END

#endif /* __TIMER_H */
