/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#include <Carbon/Carbon.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gltexture.h>
#include <mitsuba/hw/glgeometry.h>
#include <mitsuba/hw/glprogram.h>
#include <mitsuba/hw/glsync.h>

MTS_NAMESPACE_BEGIN

PrimitiveThreadLocal<GLEWContext> glewContext;

GLEWContext *glewGetContext() {
	return &glewContext.get();
}

GLRenderer::GLRenderer(Session *session)
 : Renderer(session) {
}

GLRenderer::~GLRenderer() {
}

void GLRenderer::init(Device *device, Renderer *other) {
	Renderer::init(device, other);

	m_driverRenderer = (char *) glGetString(GL_RENDERER);
	m_driverVendor = (char *) glGetString(GL_VENDOR);
	m_driverVersion = (char *) glGetString(GL_VERSION);

	Log(m_logLevel, "OpenGL renderer : %s", m_driverRenderer.c_str());
	Log(m_logLevel, "OpenGL vendor   : %s", m_driverVendor.c_str());
	Log(m_logLevel, "OpenGL version  : %s", m_driverVersion.c_str());

	/* OpenGL extensions */
	GLenum err = glewInit();
	if (err != GLEW_OK) 
		Log(EError, "GLEW Error: %s\n", glewGetErrorString(err));

	if (glewIsSupported("GL_EXT_framebuffer_object")) {
		m_capabilities->setSupported(
			RendererCapabilities::ERenderToTexture, true);
		Log(m_logLevel, "Capabilities: Framebuffers objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Framebuffers objects are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_shading_language_100")) {
		m_capabilities->setSupported(
			RendererCapabilities::EShadingLanguage, true);
		Log(m_logLevel, "Capabilities: GLSL is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: GLSL is NOT supported!");
	}

	if (glewIsSupported("GL_ARB_texture_float")) {
		m_capabilities->setSupported(
			RendererCapabilities::EFloatingPointTextures, true);
		Log(m_logLevel, "Capabilities: Floating point textures are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Floating point textures are NOT supported!");
	}

	bool leopardWorkaround = false;
#if defined(__OSX__)
	/* Color render buffers sort-of work in Leopard/8600M or 9600M, but the extension is not reported */

	SInt32 MacVersion;
	if (Gestalt(gestaltSystemVersion, &MacVersion) == noErr) {
		if (MacVersion >= 0x1050 && MacVersion < 0x1060) {
			Log(EInfo, "Enabling Leopard floating point color buffer workaround");
			leopardWorkaround = true;
		}
	}
#endif

	if (glewIsSupported("GL_ARB_color_buffer_float") || leopardWorkaround) {
		m_capabilities->setSupported(
			RendererCapabilities::EFloatingPointBuffer, true);
		Log(m_logLevel, "Capabilities: Floating point color buffers are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Floating point color buffers are NOT supported!");
	}

	if (glewIsSupported("GL_EXT_framebuffer_blit")) {
		m_capabilities->setSupported(
			RendererCapabilities::EBufferBlit, true);
		Log(m_logLevel, "Capabilities: Fast buffer blitting is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Fast buffer blitting is NOT supported!");
	}

	if (glewIsSupported("GL_EXT_framebuffer_multisample") && 
		glewIsSupported("GL_EXT_framebuffer_blit") &&
		glewIsSupported("GL_ARB_texture_multisample")) {
		m_capabilities->setSupported(
			RendererCapabilities::EMultisampleRenderToTexture, true);
		Log(m_logLevel, "Capabilities: Multisample framebuffer objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Multisample framebuffer objects are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_vertex_buffer_object")) {
		m_capabilities->setSupported(
			RendererCapabilities::EVertexBufferObjects, true);
		Log(m_logLevel, "Capabilities: Vertex buffer objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Vertex buffer objects are NOT supported!");
	}

	if (glewIsSupported("GL_EXT_geometry_shader4")) {
		m_capabilities->setSupported(
			RendererCapabilities::EGeometryShaders, true);
		Log(m_logLevel, "Capabilities: Geometry shaders are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Geometry shaders are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_sync")) {
		m_capabilities->setSupported(
			RendererCapabilities::ESyncObjects, true);
		Log(m_logLevel, "Capabilities: Synchronization objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Synchronization objects are NOT supported!");
	}

	if (glewIsSupported("GL_NV_vertex_buffer_unified_memory")) {
		m_capabilities->setSupported(
			RendererCapabilities::EBindless, true);
		Log(m_logLevel, "Capabilities: Bindless rendering is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Bindless rendering is NOT supported!");
	}

	/* Hinting */
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	/* Disable color value clamping */
	if (m_capabilities->isSupported(
			RendererCapabilities::EFloatingPointBuffer) && !leopardWorkaround) {
		glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
	}

	/* Clip to viewport */
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	setBlendMode(EBlendNone);

	checkError();
}

void GLRenderer::shutdown() {
	Renderer::shutdown();
}

GPUTexture *GLRenderer::createGPUTexture(const std::string &name,
		Bitmap *bitmap) {
	return new GLTexture(name, bitmap);
}

GPUGeometry *GLRenderer::createGPUGeometry(const TriMesh *mesh) {
	return new GLGeometry(mesh);
}

GPUProgram *GLRenderer::createGPUProgram(const std::string &name) {
	return new GLProgram(name);
}
	
GPUSync *GLRenderer::createGPUSync() {
	return new GLSync();
}

void GLRenderer::clear() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	checkError();
}

void GLRenderer::checkError(bool onlyWarn) {
	int glError = glGetError();

	if (glError)
		Log(onlyWarn ? m_warnLogLevel : EError, "OpenGL Error : %s", gluErrorString(glError));
}

void GLRenderer::beginDrawingMeshes(bool transmitOnlyPositions) {
	m_transmitOnlyPositions = transmitOnlyPositions;
	glEnableClientState(GL_VERTEX_ARRAY);

	if (!transmitOnlyPositions) {
		glEnableClientState(GL_NORMAL_ARRAY);
		glClientActiveTexture(GL_TEXTURE0);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glClientActiveTexture(GL_TEXTURE1);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
#if defined(MTS_HAS_VERTEX_COLORS)
		glEnableClientState(GL_COLOR_ARRAY);
#endif
	}

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
#if defined(MTS_HAS_VERTEX_COLORS)
		const int stride = sizeof(GLfloat) * 14;
#else
		const int stride = sizeof(GLfloat) * 11;
#endif
		glEnableClientState(GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV);
		glEnableClientState(GL_ELEMENT_ARRAY_UNIFIED_NV);
		glVertexFormatNV(3, GL_FLOAT, stride);
		if (!transmitOnlyPositions) {
			glNormalFormatNV(GL_FLOAT, stride);
			glClientActiveTexture(GL_TEXTURE0);
			glTexCoordFormatNV(2, GL_FLOAT, stride);
			glClientActiveTexture(GL_TEXTURE1);
			glTexCoordFormatNV(3, GL_FLOAT, stride);
			glColorFormatNV(3, GL_FLOAT, stride);
		}
	}
}

void GLRenderer::drawTriMesh(const TriMesh *mesh) {
	std::map<const TriMesh *, GPUGeometry *>::iterator it = m_geometry.find(mesh);
	if (it != m_geometry.end()) {
		/* Draw using vertex buffer objects */
		GLGeometry *geometry = static_cast<GLGeometry *>((*it).second);
		if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
			glBufferAddressRangeNV(GL_VERTEX_ARRAY_ADDRESS_NV, 0, 
				geometry->m_vertexAddr, geometry->m_vertexSize);
			if (!m_transmitOnlyPositions) {
				glBufferAddressRangeNV(GL_NORMAL_ARRAY_ADDRESS_NV, 0, geometry->m_vertexAddr+3*sizeof(GLfloat), 
					geometry->m_vertexSize - 3*sizeof(GLfloat));
				glBufferAddressRangeNV(GL_TEXTURE_COORD_ARRAY_ADDRESS_NV, 0, geometry->m_vertexAddr+6*sizeof(GLfloat), 
					geometry->m_vertexSize - 6*sizeof(GLfloat));
				glBufferAddressRangeNV(GL_TEXTURE_COORD_ARRAY_ADDRESS_NV, 1, geometry->m_vertexAddr+8*sizeof(GLfloat), 
					geometry->m_vertexSize - 8*sizeof(GLfloat));
#if defined(MTS_HAS_VERTEX_COLORS)
				glBufferAddressRangeNV(GL_COLOR_ARRAY_ADDRESS_NV, 0, geometry->m_vertexAddr+11*sizeof(GLfloat),
					geometry->m_vertexSize - 11*sizeof(GLfloat));
#endif
			}
			glBufferAddressRangeNV(GL_ELEMENT_ARRAY_ADDRESS_NV, 0, 
				geometry->m_indexAddr, geometry->m_indexSize);
		} else {
#if defined(MTS_HAS_VERTEX_COLORS)
			const int stride = sizeof(GLfloat) * 14;
#else
			const int stride = sizeof(GLfloat) * 11;
#endif

			glBindBuffer(GL_ARRAY_BUFFER, geometry->m_vertexID);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->m_indexID);

			/* Set up the vertex/normal arrays */
			glVertexPointer(3, GL_FLOAT, stride, (GLfloat *) 0);

			if (!m_transmitOnlyPositions) {
				glNormalPointer(GL_FLOAT, stride, (GLfloat *) 0 + 3);

				glClientActiveTexture(GL_TEXTURE0);
				glTexCoordPointer(2, GL_FLOAT, stride, (GLfloat *) 0 + 6);

				/* Pass 'dpdu' as second set of texture coordinates */
				glClientActiveTexture(GL_TEXTURE1);
				glTexCoordPointer(3, GL_FLOAT, stride, (GLfloat *) 0 + 8);

#if defined(MTS_HAS_VERTEX_COLORS)
				glColorPointer(3, GL_FLOAT, stride, (GLfloat *) 0 + 11);
#endif
			}
		}

		/* Draw all triangles */
		glDrawElements(GL_TRIANGLES, (GLsizei) (mesh->getTriangleCount() * 3), 
			GL_UNSIGNED_INT, (GLvoid *) 0);
	
		if (!m_capabilities->isSupported(RendererCapabilities::EBindless)) {
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	} else {
		/* Draw the old-fashioned way without VBOs */
		const GLchar *vertices = (const GLchar *) mesh->getVertexBuffer();
		const GLchar *indices  = (const GLchar *) mesh->getTriangles();
		GLenum dataType = sizeof(Float) == 4 ? GL_FLOAT : GL_DOUBLE;

		glVertexPointer(3, dataType, sizeof(Vertex), vertices);

		if (!m_transmitOnlyPositions) {
			glNormalPointer(dataType, sizeof(Vertex), 
				vertices + sizeof(Float) * 3);

#if defined(MTS_HAS_VERTEX_COLORS)
			glColorPointer(3, GL_FLOAT, sizeof(Vertex),
				vertices + sizeof(Float) * 14);
#endif
			glClientActiveTexture(GL_TEXTURE0);
			glTexCoordPointer(2, dataType, sizeof(Vertex), 
				vertices + sizeof(Float) * 6);

			/* Pass 'dpdu' as second set of texture coordinates */
			glClientActiveTexture(GL_TEXTURE1);
			glTexCoordPointer(3, dataType, sizeof(Vertex), 
				vertices + sizeof(Float) * 8);
		}

		/* Draw all triangles */
		glDrawElements(GL_TRIANGLES, (GLsizei) (mesh->getTriangleCount()*3), 
			GL_UNSIGNED_INT, indices);
	}
}

void GLRenderer::endDrawingMeshes() {
	glDisableClientState(GL_VERTEX_ARRAY);
	if (!m_transmitOnlyPositions) {
#if defined(MTS_HAS_VERTEX_COLORS)
		glDisableClientState(GL_COLOR_ARRAY);
#endif
		glDisableClientState(GL_NORMAL_ARRAY);
		glClientActiveTexture(GL_TEXTURE1);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glClientActiveTexture(GL_TEXTURE0);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	}

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		glDisableClientState(GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV);
		glDisableClientState(GL_ELEMENT_ARRAY_UNIFIED_NV);
	} else {
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
}
	
void GLRenderer::drawAll() {
	GLRenderer::beginDrawingMeshes(true);
	std::map<const TriMesh *, GPUGeometry *>::iterator it;
	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		for (it = m_geometry.begin(); it != m_geometry.end(); ++it) {
			const TriMesh *mesh = static_cast<const TriMesh *>((*it).first);
			const GLGeometry *geometry = static_cast<const GLGeometry *>((*it).second);
			glBufferAddressRangeNV(GL_VERTEX_ARRAY_ADDRESS_NV, 0, 
				geometry->m_vertexAddr, geometry->m_vertexSize);
			glBufferAddressRangeNV(GL_ELEMENT_ARRAY_ADDRESS_NV, 0, 
				geometry->m_indexAddr, geometry->m_indexSize);
			glDrawElements(GL_TRIANGLES, (GLsizei) (mesh->getTriangleCount()*3), 
				GL_UNSIGNED_INT, 0);
		}
	} else {
		for (it = m_geometry.begin(); it != m_geometry.end(); ++it) {
			const TriMesh *mesh = static_cast<const TriMesh *>((*it).first);
			const GLGeometry *geometry = static_cast<const GLGeometry *>((*it).second);

			glBindBuffer(GL_ARRAY_BUFFER, geometry->m_vertexID);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->m_indexID);

			/* Set up the vertex/normal arrays */
			glVertexPointer(3, GL_FLOAT, sizeof(GLfloat) * 11, (GLfloat *) 0);
			glDrawElements(GL_TRIANGLES, (GLsizei) (mesh->getTriangleCount()*3), 
				GL_UNSIGNED_INT, (GLvoid *) 0);
		}
	}
	GLRenderer::endDrawingMeshes();
}

void GLRenderer::blitTexture(const GPUTexture *tex, bool flipVertically,
		bool centerHoriz, bool centerVert, const Vector2i &offset) {
	tex->bind();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	if (tex->getType() == GPUTexture::ETexture2D) {
		GLint viewport[4];	
		glGetIntegerv(GL_VIEWPORT, viewport);
		Vector2i scrSize = Vector2i(viewport[2], viewport[3]);
		Vector2i texSize = Vector2i(tex->getSize().x, tex->getSize().y);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glBegin(GL_QUADS);

		Vector2i upperLeft(0), lowerRight(0);
		if (centerHoriz)
			upperLeft.x = (scrSize.x - texSize.x)/2;
		if (centerVert)
			upperLeft.y = (scrSize.y - texSize.y)/2;
		upperLeft += offset;
		lowerRight = upperLeft + texSize;

		if (flipVertically)
			std::swap(upperLeft.y, lowerRight.y);

		const float zDepth = -1.0f; // just before the far plane
		glTexCoord2f(0.0f, 0.0f);
		glVertex3f(upperLeft.x, upperLeft.y, zDepth);
		glTexCoord2f(1.0f, 0.0f);
		glVertex3f(lowerRight.x, upperLeft.y, zDepth);
		glTexCoord2f(1.0f, 1.0f);
		glVertex3f(lowerRight.x, lowerRight.y, zDepth);
		glTexCoord2f(0.0f, 1.0f);
		glVertex3f(upperLeft.x, lowerRight.y, zDepth);
		glEnd();
	} else if (tex->getType() == GPUTexture::ETextureCubeMap) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		/* From OpenTK */
		glBegin(GL_QUADS);
		// 0 -x
		glTexCoord3f(-1.0f, +1.0f, -1.0f);
		glVertex2f(-1.0f, +0.333f);
		glTexCoord3f(-1.0f, +1.0f, +1.0f);
		glVertex2f(-0.5f, +0.333f);
		glTexCoord3f(-1.0f, -1.0f, +1.0f);
		glVertex2f(-0.5f, -0.333f);
		glTexCoord3f(-1.0f, -1.0f, -1.0f);
		glVertex2f(-1.0f, -0.333f);

		// 1 +z
		glTexCoord3f(-1.0f, +1.0f, +1.0f);
		glVertex2f(-0.5f, +0.333f);
		glTexCoord3f(+1.0f, +1.0f, +1.0f);
		glVertex2f(+0.0f, +0.333f);
		glTexCoord3f(+1.0f, -1.0f, +1.0f);
		glVertex2f(+0.0f, -0.333f);
		glTexCoord3f(-1.0f, -1.0f, +1.0f);
		glVertex2f(-0.5f, -0.333f);

		// 2 +x
		glTexCoord3f(+1.0f, +1.0f, +1.0f);
		glVertex2f(+0.0f, +0.333f);
		glTexCoord3f(+1.0f, +1.0f, -1.0f);
		glVertex2f(+0.5f, +0.333f);
		glTexCoord3f(+1.0f, -1.0f, -1.0f);
		glVertex2f(+0.5f, -0.333f);
		glTexCoord3f(+1.0f, -1.0f, +1.0f);
		glVertex2f(+0.0f, -0.333f);

		// 3 -z
		glTexCoord3f(+1.0f, +1.0f, -1.0f);
		glVertex2f(+0.5f, +0.333f);
		glTexCoord3f(-1.0f, +1.0f, -1.0f);
		glVertex2f(+1.0f, +0.333f);
		glTexCoord3f(-1.0f, -1.0f, -1.0f);
		glVertex2f(+1.0f, -0.333f);
		glTexCoord3f(+1.0f, -1.0f, -1.0f);
		glVertex2f(+0.5f, -0.333f);

		// 4 +y
		glTexCoord3f(-1.0f, +1.0f, -1.0f);
		glVertex2f(-0.5f, +1.0f);
		glTexCoord3f(+1.0f, +1.0, -1.0f);
		glVertex2f(+0.0f, +1.0);
		glTexCoord3f(+1.0f, +1.0f, +1.0f);
		glVertex2f(+0.0f, +0.333f);
		glTexCoord3f(-1.0f, +1.0f, +1.0f);
		glVertex2f(-0.5f, +0.333f);

		// 5 -y
		glTexCoord3f(-1.0f, -1.0f, +1.0f);
		glVertex2f(-0.5f, -0.333f);
		glTexCoord3f(+1.0f, -1.0, +1.0f);
		glVertex2f(+0.0f, -0.333f);
		glTexCoord3f(+1.0f, -1.0f, -1.0f);
		glVertex2f(+0.0f, -1.0f);
		glTexCoord3f(-1.0f, -1.0f, -1.0f);
		glVertex2f(-0.5f, -1.0f);
		glEnd();
	}
	tex->unbind();
}

void GLRenderer::blitQuad(bool flipVertically) {
	GLint viewport[4];	
	glGetIntegerv(GL_VIEWPORT, viewport);
	Vector2 scrSize(viewport[2], viewport[3]);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	const Float zDepth = -1.0f;
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, flipVertically ? 1.0f : 0.0f);
	glVertex3f(0.0f, 0.0f, zDepth);
	glTexCoord2f(1.0f, flipVertically ? 1.0f : 0.0f);
	glVertex3f(scrSize.x, 0.0f, zDepth);
	glTexCoord2f(1.0f, flipVertically ? 0.0f : 1.0f);
	glVertex3f(scrSize.x, scrSize.y, zDepth);
	glTexCoord2f(0.0f, flipVertically ? 0.0f : 1.0f);
	glVertex3f(0.0f, scrSize.y, zDepth);
	glEnd();
}

void GLRenderer::setCamera(const ProjectiveCamera *camera) {
	GLfloat temp1[16], temp2[16];
	Matrix4x4 *view = const_cast<Matrix4x4 *>(camera->getViewTransform().getMatrix());
	Matrix4x4 *proj = const_cast<Matrix4x4 *>(camera->getGLProjectionTransform().getMatrix());

	int pos=0;
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			temp1[pos]=(GLfloat) proj->m[i][j];
			temp2[pos++]=(GLfloat) view->m[i][j];
		}
	}

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(temp1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0f, 1.0f, -1.0f);
	glMultMatrixf(temp2);
}

void GLRenderer::setCamera(const ProjectiveCamera *camera, const Point2 &jitter) {
	GLfloat temp1[16], temp2[16];
	ref<Matrix4x4> view = const_cast<Matrix4x4 *>(camera->getViewTransform().getMatrix());
	ref<Matrix4x4> proj = const_cast<Matrix4x4 *>(camera->getGLProjectionTransform(jitter).getMatrix());

	int pos=0;
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			temp1[pos]=(GLfloat) proj->m[i][j];
			temp2[pos++]=(GLfloat) view->m[i][j];
		}
	}

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(temp1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0f, 1.0f, -1.0f);
	glMultMatrixf(temp2);
}

void GLRenderer::setCamera(const Matrix4x4 *proj, const Matrix4x4 *view) {
	GLfloat temp1[16], temp2[16];
	int pos=0;
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			temp1[pos]=(GLfloat) proj->m[i][j];
			temp2[pos++]=(GLfloat) view->m[i][j];
		}
	}

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(temp1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0f, 1.0f, -1.0f);
	glMultMatrixf(temp2);
}

void GLRenderer::setDepthOffset(Float value) {
	if (value == 0)
		glDisable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(4.0f, (GLfloat) value);
}

void GLRenderer::setDepthMask(bool value) {
	glDepthMask(value ? GL_TRUE : GL_FALSE);
}

void GLRenderer::setDepthTest(bool value) {
	if (value)
		glEnable(GL_DEPTH_TEST);
	else
		glDisable(GL_DEPTH_TEST);
}

void GLRenderer::setColorMask(bool value) {
	GLboolean flag = value ? GL_TRUE : GL_FALSE;
	glColorMask(flag, flag, flag, flag);
}

void GLRenderer::flush() {
	glFlush();
}

void GLRenderer::finish() {
	glFinish();
}
	
void GLRenderer::setColor(const Spectrum &spec) {
	Float r, g, b;
	spec.toLinearRGB(r, g, b);
	glColor4f(r, g, b, 1);
}

void GLRenderer::setBlendMode(EBlendMode mode) {
	switch (mode) {
		case EBlendNone:
			glDisable(GL_BLEND);
			break;
		case EBlendAdditive:
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
			break;
		case EBlendAlpha:
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			break;
		default:
			Log(EError, "Invalid blend mode!");
	}
}

void GLRenderer::setCullMode(ECullMode mode) {
	switch (mode) {
		case ECullNone:
			glDisable(GL_CULL_FACE);
			break;
		case ECullFront:
			glEnable(GL_CULL_FACE);
			glCullFace(GL_FRONT);
			break;
		case ECullBack:
			glEnable(GL_CULL_FACE);
			glCullFace(GL_BACK);
			break;
		default:
			Log(EError, "Invalid culling mode!");
	}
}

MTS_IMPLEMENT_CLASS(GLRenderer, true, Renderer)
MTS_NAMESPACE_END
