#include <mitsuba/hw/glxrenderer.h>

MTS_NAMESPACE_BEGIN

GLXRenderer::GLXRenderer(X11Session *session)
 : GLRenderer(session) {
}

GLXRenderer::~GLXRenderer() {
	if (m_initialized)
		shutdown();
}

void GLXRenderer::init(Device *device, Renderer *other) {
    X11Session *session = static_cast<X11Session *>(m_session.get());
    GLXDevice *glxDevice = static_cast<GLXDevice *>(device);

	if (session == NULL) {
		Log(EDebug, "Using an existing GLX context");
		m_context = glXGetCurrentContext();
		m_borrowed = true;
	} else {
		GLXContext otherContext = NULL;
		Log(EDebug, "Initializing GLX renderer");

		if (other != NULL) {
			Assert(other->getClass() == m_theClass);
			otherContext = static_cast<GLXRenderer *>(other)->m_context;
		}

		/* Create a GL context */
		m_context = glXCreateContext(session->m_display, glxDevice->getVisual(), otherContext, True);
		if (m_context == NULL)
			Log(EError, "Could not create GLX context: failed on the client side!");
		else if (m_context == (GLXContext) BadMatch)
			Log(EError, "Could not create GLX context: bad match with shared context!");
		else if (m_context == (GLXContext) BadValue)
			Log(EError, "Could not create GLX context: bad visual!");
		else if (m_context == (GLXContext) BadAlloc)
			Log(EError, "Could not create GLX context: not enough resources!");

		glxDevice->makeCurrent(this);
		m_borrowed = false;
	}

	GLRenderer::init(glxDevice);

	m_initialized = true;
}

void GLXRenderer::shutdown() {
	GLRenderer::shutdown();

	if (!m_borrowed) {
		Log(EDebug, "Shutting down GLX Renderer");
		X11Session *session = static_cast<X11Session *>(m_session.get());   
		glXDestroyContext(session->m_display, m_context);
	}

	m_initialized = false;
}

MTS_IMPLEMENT_CLASS(GLXRenderer, false, GLRenderer)
MTS_NAMESPACE_END
