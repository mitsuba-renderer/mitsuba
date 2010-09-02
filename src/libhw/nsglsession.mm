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

#include <mitsuba/hw/nsgldevice.h>
#include <Cocoa/Cocoa.h>

MTS_NAMESPACE_BEGIN

NSGLSession::NSGLSession()
 : Session() {
}

NSGLSession::~NSGLSession() {
    Log(EDebug, "Destroying NSGL session");
	if (m_initialized)
		shutdown();
}

void NSGLSession::init() {
	Session::init();
    Log(EDebug, "Initializing NSGL session");
	m_initialized = true;
}

void NSGLSession::shutdown() {
	Session::shutdown();
    Log(EDebug, "Shutting down NSGL session");
	m_initialized = false;
}

void NSGLSession::processEvents() {
	std::vector<Device *>::iterator it = m_devices.begin();

	/* Separate autorelease pool */
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	NSDate *distantPast = [NSDate distantPast];
	NSEvent *event;

	/* NSRunloop emulation */
	do {
		event = [NSApp nextEventMatchingMask: NSAnyEventMask untilDate: distantPast inMode: NSDefaultRunLoopMode dequeue: YES];
		if (event == nil)
			break;
		[NSApp sendEvent: event];
	} while (true);

	[pool release];

	/* After the events have been delivered to the devices,
	   process them */
	for (; it!=m_devices.end(); ++it) {
		NSGLDevice *device = static_cast<NSGLDevice *>(*it);
		device->processEvents();
	}
}

MTS_IMPLEMENT_CLASS(NSGLSession, false, Session)
MTS_NAMESPACE_END
