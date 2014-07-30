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

#include <mitsuba/hw/nsglrenderer.h>
#include <dlfcn.h>

MTS_NAMESPACE_BEGIN

NSGLRenderer::NSGLRenderer(NSGLSession *session)
 : GLRenderer(session) {
}

NSGLRenderer::~NSGLRenderer() {
	Log(EDebug, "Destroying NSGL renderer");

	if (m_initialized)
		shutdown();
}

void NSGLRenderer::init(Device *device, Renderer *other) {
    NSGLDevice *nsglDevice = static_cast<NSGLDevice *>(device);

	if (m_session == NULL) {
		Log(EDebug, "Using an existing NSGL context");
		m_context = [NSOpenGLContext currentContext];
		m_borrowed = true;
	} else {
		NSOpenGLContext *otherContext = nil;
		NSOpenGLPixelFormat *format;
		Log(EDebug, "Initializing NSGL renderer");

		if (other != NULL) {
			Assert(other->getClass() == m_theClass);
			otherContext = static_cast<NSGLRenderer *>(other)->m_context;
		    NSOpenGLPixelFormatAttribute attribs[] = { 0 };
			format = [[NSOpenGLPixelFormat alloc] initWithAttributes:attribs];
		} else {
			format = ((NSOpenGLPixelFormat *) nsglDevice->getPixelFormat());
		}

		/* Create a GL context */
		m_context = [[NSOpenGLContext alloc]
			initWithFormat: format
			shareContext: otherContext];

	//	long vsync = 1;
	//	[m_context setValues: &vsync forParameter: NSOpenGLCPSwapInterval];
		if (other != NULL)
			[format release];

		if (m_context == nil)
			Log(EError, "Could not create NSGL rendering context");

		device->makeCurrent(this);
		m_borrowed = false;
	}

	GLRenderer::init(device);

	m_initialized = true;
}

void NSGLRenderer::shutdown() {
	GLRenderer::shutdown();
	if (!m_borrowed) {
		Log(EDebug, "Shutting down NSGL renderer");
		[m_context release];
	}
	m_initialized = false;
}

void *NSGLRenderer::lookupExtension(const std::string &name) const {
	std::string symName = "_" + name;
	return dlsym(RTLD_DEFAULT, symName.c_str());
}

MTS_IMPLEMENT_CLASS(NSGLRenderer, true, GLRenderer)
MTS_NAMESPACE_END
