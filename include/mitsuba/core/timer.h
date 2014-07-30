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
#if !defined(__MITSUBA_CORE_TIMER_H_)
#define __MITSUBA_CORE_TIMER_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Platform independent milli/micro/nanosecond timer
 * \ingroup libcore
 *
 * This class implements a simple cross-platform timer with nanosecond resolution.
 * It operates similarly to a good stop watch: it can be started and stopped and
 * records both the time since the last \ref start() invocation, and the total time
 * collected in separate intervals.
 *
 * \author Edgar Velazquez-Armendariz
 */
class MTS_EXPORT_CORE Timer : public Object {
public:
	/**
	 * \brief Create a new timer and start it unless the optional
	 * \c start argument is set to \c false.
	 */
	Timer(bool start = true);

	/// Start the timer
	void start();

	/**
	 * \brief Reset the timer, including the total elapsed time across
	 * all intervals (and restart it by default)
	 */
	void reset(bool restart = true);

	/// Stop the timer and return the total elapsed time across all intervals in seconds
	Float stop();

	/// Return the number of nanoseconds that the timer has ticked so far (in total)
	uint64_t getNanoseconds() const;

	/// Return the number of microseconds that the timer has ticked so far (in total)
	unsigned int getMicroseconds() const;

	/// Return the number of milliseconds that the timer has ticked so far (in total)
	unsigned int getMilliseconds() const;

	/// Return the number of seconds that the timer has ticked so far (in total)
	Float getSeconds() const;

	/// Return the number of nanoseconds that have elapsed since the \ref start() invocation
	uint64_t getNanosecondsSinceStart() const;

	/// Return the number of microseconds that have elapsed since the last \ref start() invocation
	unsigned int getMicrosecondsSinceStart() const;

	/// Return the number of milliseconds that have elapsed since the last \ref start() invocation
	unsigned int getMillisecondsSinceStart() const;

	/// Return the number of seconds that have elapsed since the last \ref start() invocation
	Float getSecondsSinceStart() const;

	/**
	 * \brief "Lap"-style interface
	 *
	 * This function is the atomic equivalent to stopping the
	 * timer, recording the time passed since it was started,
	 * and restarting it. The resulting time value in seconds
	 * is returned.
	 */
	Float lap();

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Timer();

	/// Return the time in nanoseconds since the last \ref start() invocation
	double timeSinceStart() const;
private:
	/// Time since the last call to start() in nanoseconds
	double m_startTime;

	/// Total time the timer has been active in nanoseconds
	double m_elapsed;

	/// Is the timer currently active?
	bool m_active;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_TIMER_H_ */
