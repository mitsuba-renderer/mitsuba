/*************************************************************************

	Apollyon, a set of C++ utilities including standard templates and
	a OpenGL-based rendering subsystem.

   					Copyright (C) 2005 Wenzel Jakob

	Apollyon is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation. This program
	is distributed WITHOUT ANY WARRANTY; without even the implied
	warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the GNU Lesser General Public License for more details.

 *************************************************************************/

#if !defined(__NSGLSESSION_H)
#define __NSGLSESSION_H

#include <mitsuba/hw/session.h>

MTS_NAMESPACE_BEGIN

/** \brief A MacOS X (NSGL) windowing environment session
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

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NSGLSession();
};

MTS_NAMESPACE_END

#endif /* __NSGLSESSION_H */
