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

#if !defined(__WGLSESSION_H)
#define __WGLSESSION_H

#include <mitsuba/hw/session.h>

MTS_NAMESPACE_BEGIN

/** \brief Windows (WGL) windowing environment session
 * \ingroup libhw
 */
class MTS_EXPORT_HW WGLSession : public Session {
	friend class WGLDevice;
	friend class WGLRenderer;
public:
	/// Create a new session
	WGLSession();

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
	virtual ~WGLSession();
protected:
	HINSTANCE m_hinstance;
	std::string m_wndClassName;
};

MTS_NAMESPACE_END

#endif /* __WGLSESSION_H */
