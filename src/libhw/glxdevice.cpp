#include <mitsuba/hw/glxdevice.h>

#ifndef GLX_ARB_multisample
#define GLX_ARB_multisample
#define GLX_SAMPLE_BUFFERS_ARB		100000
#define GLX_SAMPLES_ARB				100001
#endif

MTS_NAMESPACE_BEGIN

GLXDevice::GLXDevice(X11Session *session)
 : X11Device(session) {
}

GLXDevice::~GLXDevice() {
	if (m_initialized)
		shutdown();
}

void GLXDevice::flip() {
	Assert(m_initialized);

	Device::flip();

	X11Session *session = static_cast<X11Session *>(getSession());
	glFinish();
	if (m_doubleBuffer) 
		glXSwapBuffers(session->m_display, m_window);
}

XVisualInfo *GLXDevice::createVisual() {
	int attribs[64], i=0;
	X11Session *session = static_cast<X11Session *>(getSession());

	if (!session->m_hasGLX)
        Log(EError, "GLX support is required for hardware rendering!");

	attribs[i++] = GLX_RGBA;
	attribs[i++] = GLX_RED_SIZE; attribs[i++] = m_redBits;
	attribs[i++] = GLX_GREEN_SIZE; attribs[i++] = m_greenBits;
	attribs[i++] = GLX_BLUE_SIZE; attribs[i++] = m_blueBits;
	attribs[i++] = GLX_ALPHA_SIZE; attribs[i++] = m_alphaBits;
	attribs[i++] = GLX_DEPTH_SIZE; attribs[i++] = m_depthBits;
	attribs[i++] = GLX_STENCIL_SIZE; attribs[i++] = m_stencilBits;

	if (m_doubleBuffer) {
		attribs[i++] = GLX_DOUBLEBUFFER;
	}

	if (m_fsaa > 1) {
		attribs[i++] = GLX_SAMPLE_BUFFERS_ARB; attribs[i++] = 1;
		attribs[i++] = GLX_SAMPLES_ARB; attribs[i++] = m_fsaa;
	}

	attribs[i] = None;

	XVisualInfo *visinfo = glXChooseVisual(session->m_display, session->m_screen, attribs);

    if (visinfo == NULL)
        Log(EError, "Could not find a matching visual!");

	return visinfo;
}

MTS_IMPLEMENT_CLASS(GLXDevice, false, X11Device)
MTS_NAMESPACE_END
