#include <mitsuba/hw/session.h>
#if  defined(WIN32)
#include <mitsuba/hw/wglsession.h>
#elif defined(__OSX__)
#include <mitsuba/hw/nsglsession.h>
#else
#include <mitsuba/hw/x11session.h>
#endif

MTS_NAMESPACE_BEGIN

Session::Session() {
	m_initialized = false;
}

Session *Session::create() {
#if defined(WIN32)
	return new WGLSession();
#elif defined(__OSX__)
	return new NSGLSession();
#else
	return new X11Session();
#endif
}

void Session::init() {
	Assert(!m_initialized);
}

void Session::shutdown() {
	Assert(m_initialized);
}

MTS_IMPLEMENT_CLASS(Session, true, Object)
MTS_NAMESPACE_END
