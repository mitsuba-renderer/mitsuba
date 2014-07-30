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

#include <mitsuba/hw/wglrenderer.h>

MTS_NAMESPACE_BEGIN

WGLRenderer::WGLRenderer(WGLSession *session)
 : GLRenderer(session) {
}

WGLRenderer::~WGLRenderer() {
	if (m_initialized)
		shutdown();
}

void WGLRenderer::init(Device *device, Renderer *other) {
	WGLDevice *wglDevice = static_cast<WGLDevice *>(device);

	if (m_session == NULL) {
		Log(EDebug, "Using an existing WGL context");
		m_context = wglGetCurrentContext();
		if (m_context == NULL)
			Log(EError, "Unable to retrieve the current WGL context!");
		m_borrowed = true;
	} else {
		/* Create a GL context */
		Log(EDebug, "Initializing WGL renderer");
		m_context = wglCreateContext(wglDevice->getDeviceContext());

		if (other != NULL) {
			Assert(other->getClass() == m_theClass);
			if (wglShareLists(m_context,
				static_cast<WGLRenderer *>(other)->m_context) != TRUE)
				Log(EError, "Unable to set up context sharing: %s",
					lastErrorText().c_str());
		}

		device->makeCurrent(this);
		m_borrowed = false;
	}

	GLRenderer::init(device);

	m_initialized = true;
}

void WGLRenderer::shutdown() {
	GLRenderer::shutdown();

	if (!m_borrowed) {
		Log(EDebug, "Shutting down WGL Renderer");
		wglDeleteContext(m_context);
	}

	m_initialized = false;
}

void *WGLRenderer::lookupExtension(const std::string &name) const {
	return (void *) wglGetProcAddress(name.c_str());
}

MTS_IMPLEMENT_CLASS(WGLRenderer, true, GLRenderer)
MTS_NAMESPACE_END
