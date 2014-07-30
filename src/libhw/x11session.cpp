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

#include <mitsuba/hw/x11device.h>

MTS_NAMESPACE_BEGIN

X11Session::X11Session()
 : Session(), m_display(NULL), m_screen(0), m_hasVidMode(false), m_hasGLX(false) {
}

X11Session::~X11Session() {
	Log(EDebug, "Destroying X11 session");
	if (m_initialized)
		shutdown();
}

static int X11SessionErrorHandler(Display *pDisplay, XErrorEvent *pEvent) {
	SLog(EWarn, "Xlib error: Error code %d, request code %d",
			pEvent->error_code, pEvent->request_code);
	return 0;
}

void X11Session::setDisplayName(const std::string &pDisplayName) {
	m_displayName = pDisplayName;
}

void X11Session::init() {
	int temp, minor, major;

	Session::init();

	Log(EDebug, "Initializing X11 session");

	m_display = XOpenDisplay(m_displayName == "" ? NULL : m_displayName.c_str());
	m_hasGLX = m_hasVidMode = false;

	if (m_display == NULL) {
		Log(EError, "Cannot open the display");
	}

	XSetErrorHandler(X11SessionErrorHandler);

	if (!XF86VidModeQueryVersion(m_display, &major, &minor)) {
		Log(EWarn, "VidMode extension is not supported");
	} else {
		m_hasVidMode = true;
		Log(EDebug, "VidMode extension %i.%i found", major, minor);
	}

	if(!XQueryExtension(m_display, "GLX", &temp, &temp, &temp)) {
		Log(EError, "OpenGL is not supported");
	} else {
		if (glXQueryVersion(m_display, &major, &minor)) {
			if (major == 1 && minor < 1) {
				Log(EWarn, "GLX Version is too old (1.1 or higher is required)");
			} else {
				m_hasGLX = true;
				Log(EDebug, "GLX Version %i.%i found", major, minor);
			}
		} else {
			Log(EWarn, "Cannot query the GLX version");
		}
	}

	m_screen = DefaultScreen(m_display);
	m_root = RootWindow(m_display, m_screen);
	m_initialized = true;
}

void X11Session::shutdown() {
	Session::shutdown();

	Log(EDebug, "Shutting down X11 session");
	XCloseDisplay(m_display);
	m_initialized = false;
}

void X11Session::processEvents() {
	XEvent event;
	while (XPending(m_display)) {
		XNextEvent(m_display, &event);
		Window window = event.xany.window;

		std::vector<Device *>::iterator it = m_devices.begin();

		for (; it!=m_devices.end(); ++it) {
			X11Device *device = static_cast<X11Device *>(*it);
			if (device->getWindow() == window) {
				device->processEvent(event);
				break;
			}
		}
	}
}

void X11Session::processEventsBlocking(bool &stop) {
	XEvent event;
	while (true) {
		if (!XPending(m_display) && stop)
			break;
		XNextEvent(m_display, &event);
		Window window = event.xany.window;

		std::vector<Device *>::iterator it = m_devices.begin();

		for (; it!=m_devices.end(); ++it) {
			X11Device *device = static_cast<X11Device *>(*it);
			if (device->getWindow() == window) {
				device->processEvent(event);
				break;
			}
		}
	}
}


MTS_IMPLEMENT_CLASS(X11Session, false, Session)
MTS_NAMESPACE_END
