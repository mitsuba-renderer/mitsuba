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

#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/glsync.h>

#define MTS_SYNC_TIMEOUT 100000000  // 100ms

MTS_NAMESPACE_BEGIN

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
