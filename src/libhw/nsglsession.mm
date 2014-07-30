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
	/* Separate autorelease pool */
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	NSEvent *event;

	/* NSRunloop emulation */
	do {
		event = [NSApp nextEventMatchingMask: NSAnyEventMask
			untilDate: nil inMode: NSDefaultRunLoopMode dequeue: YES];
		if (event == nil)
			break;
		[NSApp sendEvent: event];
	} while (true);

	[pool release];

	/* After the events have been delivered to the devices, process them */
	for (std::vector<Device *>::iterator it = m_devices.begin(); it!=m_devices.end(); ++it)
		static_cast<NSGLDevice *>(*it)->processEvents();
}

void NSGLSession::processEventsBlocking(bool &stop) {
	/* Separate autorelease pool */
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	NSEvent *event;

	/* NSRunloop emulation */
	while (true) {
		event = [NSApp nextEventMatchingMask: NSAnyEventMask
			untilDate: nil inMode: NSDefaultRunLoopMode dequeue: YES];

		if (event != nil) {
			[NSApp sendEvent: event];
		} else {
			/* No more event -- process them */
			for (std::vector<Device *>::iterator it = m_devices.begin(); it!=m_devices.end(); ++it)
				static_cast<NSGLDevice *>(*it)->processEvents();

			if (stop)
				break;

			/* Wait for the next event (blocking) */
			event = [NSApp nextEventMatchingMask: NSAnyEventMask
				untilDate: [NSDate distantFuture] inMode: NSDefaultRunLoopMode dequeue: YES];

			if (event != nil)
				[NSApp sendEvent: event];
		}
	}

	[pool release];

	/* After the events have been delivered to the devices, process them */
	for (std::vector<Device *>::iterator it = m_devices.begin(); it!=m_devices.end(); ++it)
		static_cast<NSGLDevice *>(*it)->processEvents();
}

MTS_IMPLEMENT_CLASS(NSGLSession, false, Session)
MTS_NAMESPACE_END
