#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/glsync.h>

#define MTS_SYNC_TIMEOUT 100000000  // 100ms

MTS_NAMESPACE_BEGIN

extern GLEWContext *glewGetContext();

GLSync::GLSync() : GPUSync(), m_sync(0) {
}

void GLSync::init() {
	if (m_sync != 0)
		cleanup();
	m_sync = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
	glFlush();
	if (m_sync == 0)
		Log(EError, "Unable to create a memory sync object!");
}

void GLSync::cleanup() {
	glDeleteSync(m_sync);
	m_sync = 0;
}

void GLSync::wait() {
	GLenum retval = glClientWaitSync(m_sync, 
		GL_SYNC_FLUSH_COMMANDS_BIT, MTS_SYNC_TIMEOUT);

	while (true) {
		switch (retval) {
			case GL_ALREADY_SIGNALED:
			case GL_CONDITION_SATISFIED:
				return;
			case GL_WAIT_FAILED:
				break;
			default:
				Log(EError, "glClientWaitSync: unexpected return value!");
		}

		retval = glClientWaitSync(m_sync, 0, MTS_SYNC_TIMEOUT);
	}
}

void GLSync::enqueueWait() {
	uint64_t timeout = 0xFFFFFFFFFFFFFFFFULL;
	glWaitSync(m_sync, 0, timeout);
	m_sync = 0;
}

GLSync::~GLSync() {
	if (m_sync != 0)
		cleanup();
}

MTS_IMPLEMENT_CLASS(GLSync, false, GPUSync)
MTS_NAMESPACE_END
