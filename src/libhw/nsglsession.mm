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
