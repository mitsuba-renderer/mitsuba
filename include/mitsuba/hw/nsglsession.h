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

#if !defined(__MITSUBA_HW_NSGLSESSION_H_)
#define __MITSUBA_HW_NSGLSESSION_H_

#include <mitsuba/hw/session.h>

MTS_NAMESPACE_BEGIN

/** \brief A MacOS X (NSGL) windowing environment session
 * \ingroup libhw
 */
class MTS_EXPORT_HW NSGLSession : public Session {
public:
	/// Create a new session
	NSGLSession();

	/// Initialize the session
	void init();

	/// Shut the session down
	void shutdown();

	/// Process all events and call event callbacks
	void processEvents();

	/**
	 * \brief Process all events and call event callbacks.
	 *
	 * This function will run until the \c stop parameter is set
	 * to \c true from within an event callback.
	 */
	void processEventsBlocking(bool &stop);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NSGLSession();
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_NSGLSESSION_H_ */
