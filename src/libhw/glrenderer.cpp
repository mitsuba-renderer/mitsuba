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
#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gltexture.h>
#include <mitsuba/hw/glgeometry.h>
#include <mitsuba/hw/glprogram.h>
#include <mitsuba/hw/glsync.h>
#include <mitsuba/hw/font.h>
#include <boost/algorithm/string.hpp>

static mitsuba::PrimitiveThreadLocal<GLEWContextStruct> glewContext;

GLEWContextStruct *glewGetContext() {
	return &glewContext.get();
}

MTS_NAMESPACE_BEGIN

/* Helper functions */
namespace {
	FINLINE void loadMatrix(const Matrix4x4 &mat) {
		GLfloat temp[16];
		int pos = 0;

		for (int j=0; j<4; j++)
			for (int i=0; i<4; i++)
				temp[pos++] = (GLfloat) mat(i, j);

		glLoadMatrixf(temp);
	}

	FINLINE Matrix4x4 fetchMatrix(GLenum which) {
		GLfloat temp[16];
		Matrix4x4 mat;
		int pos = 0;

		glGetFloatv(which, temp);

		for (int j=0; j<4; j++)
			for (int i=0; i<4; i++)
				mat(i, j) = (Float) temp[pos++];

		return mat;
	}

	FINLINE void multMatrix(const Matrix4x4 &mat) {
		GLfloat temp[16];
		int pos = 0;

		for (int j=0; j<4; j++)
			for (int i=0; i<4; i++)
				temp[pos++] = (GLfloat) mat(i, j);

		glMultMatrixf(temp);
	}
}

GLRenderer::GLRenderer(Session *session)
 : Renderer(session) { }

GLRenderer::~GLRenderer() { }

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

	if (glewIsSupported("GL_ARB_color_buffer_float")) {
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
		Log((m_warnLogLevel == EWarn) ? EInfo : m_warnLogLevel,
			"Capabilities: Multisample framebuffer objects are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_vertex_buffer_object")) {
		m_capabilities->setSupported(
			RendererCapabilities::EVertexBufferObjects, true);
		Log(m_logLevel, "Capabilities: Vertex buffer objects are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Vertex buffer objects are NOT supported!");
	}

	if (glewIsSupported("GL_EXT_geometry_shader4") && glewIsSupported("GL_EXT_gpu_shader4")) {
		m_capabilities->setSupported(
			RendererCapabilities::EGeometryShaders, true);
		Log(m_logLevel, "Capabilities: Geometry shaders are supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Geometry shaders are NOT supported!");
	}

	if (glewIsSupported("GL_ARB_shader_texture_lod") || glewIsSupported("GL_EXT_gpu_shader4")) {
		m_capabilities->setSupported(
			RendererCapabilities::ECustomTextureFiltering, true);
		Log(m_logLevel, "Capabilities: Custom texture filtering is supported.");
	} else {
		Log(m_warnLogLevel, "Capabilities: Custom texture filtering is NOT supported.");
	}

	bool radeonOnOSX = false;

#if defined(__OSX__)
	/* Synchronization objects cause problem with ATI cards on OSX -- ignore
	   them even if the driver claims to support it */
	radeonOnOSX = boost::to_lower_copy(m_driverRenderer).find("radeon") != std::string::npos;
#endif

	if (glewIsSupported("GL_ARB_sync") && !radeonOnOSX) {
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
		Log((m_warnLogLevel == EWarn) ? EInfo : m_warnLogLevel,
			"Capabilities: Bindless rendering is NOT supported!");
	}

	/* Hinting */
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	/* Disable color value clamping */
	if (m_capabilities->isSupported(
			RendererCapabilities::EFloatingPointBuffer)) {
		glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);
		glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
	}

	/* Clip to viewport */
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	setBlendMode(EBlendNone);

	glEnable(GL_POINT_SMOOTH);

	m_normalsEnabled = false;
	m_texcoordsEnabled = false;
	m_tangentsEnabled = false;
	m_colorsEnabled = false;
	m_stride = -1;
	m_queuedTriangles = 0;
	m_transmitOnlyPositions = false;

	checkError();
}

void GLRenderer::shutdown() {
	Renderer::shutdown();
}

GPUTexture *GLRenderer::createGPUTexture(const std::string &name,
		Bitmap *bitmap) {
	return new GLTexture(name, bitmap);
}

GPUGeometry *GLRenderer::createGPUGeometry(const Shape *shape) {
	ref<TriMesh> mesh = const_cast<Shape *>(shape)->createTriMesh();
	if (!mesh)
		return NULL;
	return new GLGeometry(mesh);
}

GPUProgram *GLRenderer::createGPUProgram(const std::string &name) {
	return new GLProgram(name);
}

GPUSync *GLRenderer::createGPUSync() {
	return new GLSync();
}

void GLRenderer::reconfigure(const Device *device) {
	glViewport(0, 0, device->getSize().x, device->getSize().y);
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
	m_stride = -1;

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		glEnableClientState(GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV);
		glEnableClientState(GL_ELEMENT_ARRAY_UNIFIED_NV);
	}
}

void GLRenderer::drawMesh(const TriMesh *mesh) {
	std::map<const Shape *, GPUGeometry *>::iterator it = m_geometry.find(mesh);
	if (it != m_geometry.end()) {
		GLRenderer::drawMesh((*it).second);
	} else {
		/* This shape is not resident in GPU memory. Draw the slow way.. */
		const GLchar *positions = (const GLchar *) mesh->getVertexPositions();
		const GLchar *normals = (const GLchar *) mesh->getVertexNormals();
		const GLchar *texcoords = (const GLchar *) mesh->getVertexTexcoords();
		const GLchar *tangents = (const GLchar *) mesh->getUVTangents();
		const GLchar *colors = (const GLchar *) mesh->getVertexColors();
		const GLint *indices  = (const GLint *) mesh->getTriangles();
		GLenum dataType = sizeof(Float) == 4 ? GL_FLOAT : GL_DOUBLE;

		glVertexPointer(3, dataType, 0, positions);

		if (!m_transmitOnlyPositions) {
			if (mesh->hasVertexNormals()) {
				if (!m_normalsEnabled) {
					glEnableClientState(GL_NORMAL_ARRAY);
					m_normalsEnabled = true;
				}
				glNormalPointer(dataType, 0, normals);
			} else if (m_normalsEnabled) {
				glDisableClientState(GL_NORMAL_ARRAY);
				m_normalsEnabled = false;
			}

			glClientActiveTexture(GL_TEXTURE0);
			if (mesh->hasVertexTexcoords()) {
				if (!m_texcoordsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_texcoordsEnabled = true;
				}
				glTexCoordPointer(2, dataType, 0, texcoords);
			} else if (m_texcoordsEnabled) {
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_texcoordsEnabled = false;
			}

			/* Pass 'dpdu' as second set of texture coordinates */
			glClientActiveTexture(GL_TEXTURE1);
			if (mesh->hasUVTangents()) {
				if (!m_tangentsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_tangentsEnabled = true;
				}
				glTexCoordPointer(3, dataType, sizeof(Vector), tangents);
			} else if (m_tangentsEnabled) {
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_tangentsEnabled = false;
			}

			if (mesh->hasVertexColors()) {
				if (!m_colorsEnabled) {
					glEnableClientState(GL_COLOR_ARRAY);
					m_colorsEnabled = true;
				}
				glColorPointer(3, dataType, 0, colors);
			} else if (m_colorsEnabled) {
				glDisableClientState(GL_COLOR_ARRAY);
				m_colorsEnabled = false;
			}
		}

		size_t size = mesh->getTriangleCount();
		if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
			/* Draw all triangles */
			glDrawElements(GL_TRIANGLES, (GLsizei) (mesh->getTriangleCount()*3),
				GL_UNSIGNED_INT, indices);
			m_queuedTriangles += size;
		} else {
			/* Spoon-feed them (keeps the OS responsive) */
			size_t size = mesh->getTriangleCount(), cur = 0;
			while (cur < size) {
				size_t drawAmt = std::min(size - cur,
						MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
				if (drawAmt > 0)
					glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3),
						GL_UNSIGNED_INT, indices + cur * 3);
				m_queuedTriangles += drawAmt; cur += drawAmt;
				if (cur < size)
					finish();
			}
		}
	}
}

void GLRenderer::drawMesh(const GPUGeometry *_geo) {
	const GLGeometry *geo = static_cast<const GLGeometry *>(_geo);
	const TriMesh *mesh   = geo->getTriMesh();

	GLuint indexSize    = geo->m_size[GLGeometry::EIndexID];
	GLuint vertexSize   = geo->m_size[GLGeometry::EVertexID];

	/* Draw using vertex buffer objects (bindless if supported) */
	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		GLuint64 indexAddr  = geo->m_addr[GLGeometry::EIndexID];
		GLuint64 vertexAddr = geo->m_addr[GLGeometry::EVertexID];

		int stride = geo->m_stride;
		if (stride != m_stride) {
			glVertexFormatNV(3, GL_FLOAT, stride);
			glNormalFormatNV(GL_FLOAT, stride);
			glClientActiveTexture(GL_TEXTURE0);
			glTexCoordFormatNV(2, GL_FLOAT, stride);
			glClientActiveTexture(GL_TEXTURE1);
			glTexCoordFormatNV(3, GL_FLOAT, stride);
			glColorFormatNV(3, GL_FLOAT, stride);
			m_stride = stride;
		}

		glBufferAddressRangeNV(GL_VERTEX_ARRAY_ADDRESS_NV,
				0, vertexAddr, vertexSize);

		if (!m_transmitOnlyPositions) {
			int pos = 3 * sizeof(GLfloat);

			if (mesh->hasVertexNormals()) {
				if (!m_normalsEnabled) {
					glEnableClientState(GL_NORMAL_ARRAY);
					m_normalsEnabled = true;
				}
				glBufferAddressRangeNV(GL_NORMAL_ARRAY_ADDRESS_NV, 0,
					vertexAddr + pos, vertexSize - pos);

				pos += 3 * sizeof(GLfloat);
			} else if (m_normalsEnabled) {
				glDisableClientState(GL_NORMAL_ARRAY);
				m_normalsEnabled = false;
			}

			if (mesh->hasVertexTexcoords()) {
				glClientActiveTexture(GL_TEXTURE0);
				if (!m_texcoordsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_texcoordsEnabled = true;
				}
				glBufferAddressRangeNV(GL_TEXTURE_COORD_ARRAY_ADDRESS_NV, 0,
					vertexAddr + pos, vertexSize - pos);

				pos += 2 * sizeof(GLfloat);
			} else if (m_texcoordsEnabled) {
				glClientActiveTexture(GL_TEXTURE0);
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_texcoordsEnabled = false;
			}

			/* Pass 'dpdu' as second set of texture coordinates */
			if (mesh->hasUVTangents()) {
				glClientActiveTexture(GL_TEXTURE1);
				if (!m_tangentsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_tangentsEnabled = true;
				}

				glBufferAddressRangeNV(GL_TEXTURE_COORD_ARRAY_ADDRESS_NV, 1,
					vertexAddr + pos, vertexSize - pos);
				pos += 3 * sizeof(GLfloat);
			} else if (m_tangentsEnabled) {
				glClientActiveTexture(GL_TEXTURE1);
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_tangentsEnabled = false;
			}

			if (mesh->hasVertexColors()) {
				if (!m_colorsEnabled) {
					glEnableClientState(GL_COLOR_ARRAY);
					m_colorsEnabled = true;
				}

				glBufferAddressRangeNV(GL_COLOR_ARRAY_ADDRESS_NV, 0,
					vertexAddr + pos, vertexSize - pos);
			} else if (m_colorsEnabled) {
				glDisableClientState(GL_COLOR_ARRAY);
				m_colorsEnabled = false;
			}
		}
		glBufferAddressRangeNV(GL_ELEMENT_ARRAY_ADDRESS_NV, 0,
			indexAddr, indexSize);
	} else {
		glBindBuffer(GL_ARRAY_BUFFER, geo->m_id[GLGeometry::EVertexID]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geo->m_id[GLGeometry::EIndexID]);
		int stride = geo->m_stride;

		/* Set up the vertex/normal arrays */
		glVertexPointer(3, GL_FLOAT, stride, (GLfloat *) 0);

		if (!m_transmitOnlyPositions) {
			int pos = 3;
			if (mesh->hasVertexNormals()) {
				if (!m_normalsEnabled) {
					glEnableClientState(GL_NORMAL_ARRAY);
					m_normalsEnabled = true;
				}
				glNormalPointer(GL_FLOAT, stride, (GLfloat *) 0 + pos);
				pos += 3;
			} else if (m_normalsEnabled) {
				glDisableClientState(GL_NORMAL_ARRAY);
				m_normalsEnabled = false;
			}

			if (mesh->hasVertexTexcoords()) {
				glClientActiveTexture(GL_TEXTURE0);
				if (!m_texcoordsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_texcoordsEnabled = true;
				}
				glTexCoordPointer(2, GL_FLOAT, stride, (GLfloat *) 0 + pos);
				pos += 2;
			} else if (m_texcoordsEnabled) {
				glClientActiveTexture(GL_TEXTURE0);
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_texcoordsEnabled = false;
			}

			/* Pass 'dpdu' as second set of texture coordinates */
			if (mesh->hasUVTangents()) {
				glClientActiveTexture(GL_TEXTURE1);
				if (!m_tangentsEnabled) {
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					m_tangentsEnabled = true;
				}
				glTexCoordPointer(3, GL_FLOAT, stride, (GLfloat *) 0 + pos);
				pos += 3;
			} else if (m_tangentsEnabled) {
				glClientActiveTexture(GL_TEXTURE1);
				glDisableClientState(GL_TEXTURE_COORD_ARRAY);
				m_tangentsEnabled = false;
			}

			if (mesh->hasVertexColors()) {
				if (!m_colorsEnabled) {
					glEnableClientState(GL_COLOR_ARRAY);
					m_colorsEnabled = true;
				}
				glColorPointer(3, GL_FLOAT, stride, (GLfloat *) 0 + pos);
			} else if (m_colorsEnabled) {
				glDisableClientState(GL_COLOR_ARRAY);
				m_colorsEnabled = false;
			}
		}
	}

	size_t size = mesh->getTriangleCount();
	if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
		/* Draw all triangles */
		glDrawElements(GL_TRIANGLES, (GLsizei) (size * 3),
			GL_UNSIGNED_INT, (GLvoid *) 0);
		m_queuedTriangles += size;
	} else {
		/* Spoon-feed them (keeps the OS responsive) */
		size_t size = mesh->getTriangleCount(), cur = 0;
		while (cur < size) {
			size_t drawAmt = std::min(size - cur,
					MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
			if (drawAmt > 0)
				glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3),
					GL_UNSIGNED_INT, (GLuint *) 0 + cur * 3);
			m_queuedTriangles += drawAmt; cur += drawAmt;
			if (cur < size)
				finish();
		}
	}
}

void GLRenderer::endDrawingMeshes() {
	glDisableClientState(GL_VERTEX_ARRAY);
	if (m_normalsEnabled) {
		glDisableClientState(GL_NORMAL_ARRAY);
		m_normalsEnabled = false;
	}

	if (m_texcoordsEnabled) {
		glClientActiveTexture(GL_TEXTURE0);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		m_texcoordsEnabled = false;
	}

	if (m_tangentsEnabled) {
		glClientActiveTexture(GL_TEXTURE1);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		m_tangentsEnabled = false;
	}

	if (m_colorsEnabled) {
		glDisableClientState(GL_COLOR_ARRAY);
		m_colorsEnabled = false;
	}

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		glDisableClientState(GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV);
		glDisableClientState(GL_ELEMENT_ARRAY_UNIFIED_NV);
	} else {
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
}

void GLRenderer::drawAll(const std::vector<TransformedGPUGeometry> &allGeometry) {
	Matrix4x4 curObjTrafo;
	curObjTrafo.setIdentity();

	glMatrixMode(GL_MODELVIEW);
	Matrix4x4 backup = fetchMatrix(GL_MODELVIEW_MATRIX);

	GLRenderer::beginDrawingMeshes(true);

	if (m_capabilities->isSupported(RendererCapabilities::EBindless)) {
		for (std::vector<TransformedGPUGeometry>::const_iterator it = allGeometry.begin();
				it != allGeometry.end(); ++it) {
			const GLGeometry *geo  = static_cast<const GLGeometry *>((*it).first);
			const Matrix4x4 &trafo = (*it).second;
			const TriMesh *mesh    = geo->getTriMesh();
			GLuint indexSize       = geo->m_size[GLGeometry::EIndexID];
			GLuint vertexSize      = geo->m_size[GLGeometry::EVertexID];
			GLuint64 indexAddr     = geo->m_addr[GLGeometry::EIndexID];
			GLuint64 vertexAddr    = geo->m_addr[GLGeometry::EVertexID];

			if (trafo != curObjTrafo) {
				loadMatrix(backup * trafo);
				curObjTrafo = trafo;
			}

			int stride = geo->m_stride;
			if (stride != m_stride) {
				glVertexFormatNV(3, GL_FLOAT, stride);
				m_stride = stride;
			}

			glBufferAddressRangeNV(GL_VERTEX_ARRAY_ADDRESS_NV, 0,
				vertexAddr, vertexSize);
			glBufferAddressRangeNV(GL_ELEMENT_ARRAY_ADDRESS_NV, 0,
				indexAddr, indexSize);

			size_t size = mesh->getTriangleCount();

			if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
				/* Draw all triangles */
				glDrawElements(GL_TRIANGLES, (GLsizei) (size * 3),
					GL_UNSIGNED_INT, (GLvoid *) 0);
				m_queuedTriangles += size;
			} else {
				/* Spoon-feed them (keeps the OS responsive) */
				size_t size = mesh->getTriangleCount(), cur = 0;
				while (cur < size) {
					size_t drawAmt = std::min(size - cur,
							MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
					if (drawAmt > 0)
						glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3),
							GL_UNSIGNED_INT, (GLuint *) 0 + cur * 3);
					m_queuedTriangles += drawAmt; cur += drawAmt;
					if (cur < size)
						finish();
				}
			}
		}
	} else {
		for (std::vector<TransformedGPUGeometry>::const_iterator it = allGeometry.begin();
				it != allGeometry.end(); ++it) {
			const GLGeometry *geo  = static_cast<const GLGeometry *>((*it).first);
			const Matrix4x4 &trafo = (*it).second;
			const TriMesh *mesh    = geo->getTriMesh();

			if (trafo != curObjTrafo) {
				loadMatrix(backup * trafo);
				curObjTrafo = trafo;
			}

			glBindBuffer(GL_ARRAY_BUFFER, geo->m_id[GLGeometry::EVertexID]);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geo->m_id[GLGeometry::EIndexID]);

			/* Set up the vertex/normal arrays */
			glVertexPointer(3, GL_FLOAT, geo->m_stride, (GLfloat *) 0);

			size_t size = mesh->getTriangleCount();

			if (EXPECT_TAKEN(m_queuedTriangles + size < MTS_GL_MAX_QUEUED_TRIS)) {
				/* Draw all triangles */
				glDrawElements(GL_TRIANGLES, (GLsizei) (size * 3),
					GL_UNSIGNED_INT, (GLvoid *) 0);
				m_queuedTriangles += size;
			} else {
				/* Spoon-feed them (keeps the OS responsive) */
				size_t size = mesh->getTriangleCount(), cur = 0;
				while (cur < size) {
					size_t drawAmt = std::min(size - cur,
							MTS_GL_MAX_QUEUED_TRIS - m_queuedTriangles);
					if (drawAmt > 0)
						glDrawElements(GL_TRIANGLES, (GLsizei) (drawAmt * 3),
							GL_UNSIGNED_INT, (GLuint *) 0 + cur * 3);
					m_queuedTriangles += drawAmt; cur += drawAmt;
					if (cur < size)
						finish();
				}
			}
		}
	}
	GLRenderer::endDrawingMeshes();
	if (!curObjTrafo.isIdentity())
		loadMatrix(backup);
}

void GLRenderer::blitTexture(const GPUTexture *tex, bool flipVertically,
		bool centerHoriz, bool centerVert, const Vector2i &offset) {
	tex->bind();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	Vector2i scrSize = Vector2i(viewport[2], viewport[3]);
	Vector2i texSize = Vector2i(tex->getSize().x, tex->getSize().y);
	if (scrSize.x == 0 || scrSize.y == 0) {
		tex->unbind();
		return;
	}

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.375f, 0.375f, 0.0f);
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
	glVertex3f((float) upperLeft.x, (float) upperLeft.y, zDepth);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f((float) lowerRight.x, (float) upperLeft.y, zDepth);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f((float) lowerRight.x, (float) lowerRight.y, zDepth);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f((float) upperLeft.x, (float) lowerRight.y, zDepth);
	glEnd();

	tex->unbind();
}

void GLRenderer::blitQuad(bool flipVertically) {
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	Vector2 scrSize((float) viewport[2], (float) viewport[3]);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	const float zDepth = -1.0f;
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, flipVertically ? 1.0f : 0.0f);
	glVertex3f(0.0f, 0.0f, zDepth);
	glTexCoord2f(1.0f, flipVertically ? 1.0f : 0.0f);
	glVertex3f((GLfloat) scrSize.x, 0.0f, zDepth);
	glTexCoord2f(1.0f, flipVertically ? 0.0f : 1.0f);
	glVertex3f((GLfloat) scrSize.x, (GLfloat) scrSize.y, zDepth);
	glTexCoord2f(0.0f, flipVertically ? 0.0f : 1.0f);
	glVertex3f(0.0f, (GLfloat) scrSize.y, zDepth);
	glEnd();
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}

void GLRenderer::drawText(const Point2i &_pos,
		const Font *font, const std::string &text) {
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	Vector2i scrSize = Vector2i(viewport[2], viewport[3]);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, scrSize.x, scrSize.y, 0, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	font->getTexture()->bind();
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
//	glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ONE_MINUS_SRC_COLOR);
	Point2i pos(_pos);
	int initial = pos.x;

	glBegin(GL_QUADS);
	for (size_t i=0; i<text.length(); i++) {
		char character = text[i];
		if (character == '\r')
			continue;
		if (character == '\n') {
			pos.x = initial;
			pos.y += (int) (font->getMaxVerticalBearing()*4.0/3.0);
			continue;
		}

		const Font::Glyph &glyph = font->getGlyph(character);

		Point2 start = Point2(pos + Vector2i(
			glyph.horizontalBearing,
			font->getMaxVerticalBearing() - glyph.verticalBearing
		));
		Point2 end = start + Vector2(glyph.size);
		Point2 txStart = glyph.tx;
		Point2 txEnd = txStart + glyph.ts;

		glTexCoord2f((float) txStart.x, (float) txStart.y);
		glVertex2f(  (float) start.x,   (float)   start.y);
		glTexCoord2f((float) txEnd.x,   (float) txStart.y);
		glVertex2f(  (float) end.x,     (float)   start.y);
		glTexCoord2f((float) txEnd.x,   (float)   txEnd.y);
		glVertex2f(  (float) end.x,     (float)     end.y);
		glTexCoord2f((float) txStart.x, (float)   txEnd.y);
		glVertex2f(  (float) start.x,   (float)     end.y);

		pos.x += glyph.horizontalAdvance;

		if (i+1 < text.length())
			pos.x += font->getKerning(character, text[i+1]);
	}
	glEnd();

	font->getTexture()->unbind();
	glDisable(GL_BLEND);
}

void GLRenderer::setPointSize(Float size) {
	glPointSize((GLfloat) size);
}

void GLRenderer::drawPoint(const Point &p) {
	glBegin(GL_POINTS);
	glVertex3f((float) p.x, (float) p.y, (float) p.z);
	glEnd();
}

void GLRenderer::drawLine(const Point &a, const Point &b) {
	glBegin(GL_LINES);
	glVertex3f((float) a.x, (float) a.y, (float) a.z);
	glVertex3f((float) b.x, (float) b.y, (float) b.z);
	glEnd();
}

void GLRenderer::drawPoint(const Point2 &p) {
	glBegin(GL_POINTS);
	glVertex2f((float) p.x, (float) p.y);
	glEnd();
}

void GLRenderer::drawLine(const Point2 &a, const Point2 &b) {
	glBegin(GL_LINES);
	glVertex2f((float) a.x, (float) a.y);
	glVertex2f((float) b.x, (float) b.y);
	glEnd();
}

void GLRenderer::drawRectangle(const Point2 &a, const Point2 &b) {
	glBegin(GL_LINE_LOOP);
	glVertex2f((float) a.x, (float) a.y);
	glVertex2f((float) b.x, (float) a.y);
	glVertex2f((float) b.x, (float) b.y);
	glVertex2f((float) a.x, (float) b.y);
	glEnd();
}

void GLRenderer::drawFilledRectangle(const Point2 &a, const Point2 &b) {
	glBegin(GL_QUADS);
	glVertex2f((float) a.x, (float) a.y);
	glVertex2f((float) b.x, (float) a.y);
	glVertex2f((float) b.x, (float) b.y);
	glVertex2f((float) a.x, (float) b.y);
	glEnd();
}

void GLRenderer::drawPoint(const Point2i &p) {
	glBegin(GL_POINTS);
	glVertex2i(p.x, p.y);
	glEnd();
}

void GLRenderer::drawLine(const Point2i &a, const Point2i &b) {
	glBegin(GL_LINES);
	glVertex2i(a.x, a.y);
	glVertex2i(b.x, b.y);
	glEnd();
}

void GLRenderer::drawRectangle(const Point2i &a, const Point2i &b) {
	glBegin(GL_LINE_LOOP);
	glVertex2i(a.x, a.y);
	glVertex2i(b.x, a.y);
	glVertex2i(b.x, b.y);
	glVertex2i(a.x, b.y);
	glEnd();
}

void GLRenderer::drawFilledRectangle(const Point2i &a, const Point2i &b) {
	glBegin(GL_QUADS);
	glVertex2i(a.x, a.y);
	glVertex2i(b.x, a.y);
	glVertex2i(b.x, b.y);
	glVertex2i(a.x, b.y);
	glEnd();
}

void GLRenderer::drawEllipse(const Point &center,
		const Vector &axis1, const Vector &axis2) {
	const int nSteps = 100;
	const float stepSize = (float) (2*M_PI/nSteps);
	glBegin(GL_LINE_LOOP);
	for (int i=0; i<100; ++i) {
		Point p = center + axis1 * std::cos(i*stepSize)
			+ axis2 * std::sin(i*stepSize);
		glVertex3f((GLfloat) p.x, (GLfloat) p.y, (GLfloat) p.z);
	}
	glEnd();
}

void GLRenderer::drawAABB(const AABB &aabb) {
	#define V(a,b,c) glVertex3f((GLfloat) aabb.a.x, (GLfloat) aabb.b.y, (GLfloat) aabb.c.z)
	glBegin(GL_LINE_LOOP); V(max,min,max); V(max,min,min); V(max,max,min); V(max,max,max); glEnd();
	glBegin(GL_LINE_LOOP); V(max,max,max); V(max,max,min); V(min,max,min); V(min,max,max); glEnd();
	glBegin(GL_LINE_LOOP); V(max,max,max); V(min,max,max); V(min,min,max); V(max,min,max); glEnd();
	glBegin(GL_LINE_LOOP); V(min,min,max); V(min,max,max); V(min,max,min); V(min,min,min); glEnd();
	glBegin(GL_LINE_LOOP); V(min,min,max); V(min,min,min); V(max,min,min); V(max,min,max); glEnd();
	glBegin(GL_LINE_LOOP); V(min,min,min); V(min,max,min); V(max,max,min); V(max,min,min); glEnd();
	#undef V
}

void GLRenderer::setMatrix(EMatrixType type, const Matrix4x4 &value) {
	glMatrixMode(type == EProjection ? GL_PROJECTION : GL_MODELVIEW);
	loadMatrix(value);
}

Matrix4x4 GLRenderer::getMatrix(EMatrixType type) const {
	return fetchMatrix(type == EProjection ? GL_PROJECTION_MATRIX : GL_MODELVIEW_MATRIX);
}

void GLRenderer::setCamera(const ProjectiveCamera *camera,
		const Point2 &apertureSample, const Point2 &aaSample, Float timeSample) {
	Float time = camera->getShutterOpen() + camera->getShutterOpenTime() * timeSample;

	glMatrixMode(GL_PROJECTION);
	loadMatrix(camera->getProjectionTransform(
		apertureSample, aaSample).getMatrix());
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	/* Apply a rotation to account for the difference in camera
	   conventions. In OpenGL, forward is z=-1, in Mitsuba it is z=+1 */
	glScalef(-1.0f, 1.0f, -1.0f);
	multMatrix(camera->getViewTransform(time).getMatrix());
}

void GLRenderer::setCamera(const Matrix4x4 &proj, const Matrix4x4 &view) {
	glMatrixMode(GL_PROJECTION);
	loadMatrix(proj);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	/* Apply a rotation to account for the difference in camera
	   conventions. In OpenGL, forward is z=-1, in Mitsuba it is z=+1 */
	glScalef(-1.0f, 1.0f, -1.0f);
	multMatrix(view);
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

void GLRenderer::clearTransforms() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GLRenderer::flush() {
	glFlush();
}

void GLRenderer::finish() {
	glFinish();
	m_queuedTriangles = 0;
}

void GLRenderer::setColor(const Color3 &col, Float alpha) {
	glColor4f((GLfloat) col[0], (GLfloat) col[1], (GLfloat) col[2], (GLfloat) alpha);
}

void GLRenderer::setColor(const Spectrum &spec, Float alpha) {
	Float r, g, b;
	spec.toLinearRGB(r, g, b);
	glColor4f((GLfloat) r, (GLfloat) g, (GLfloat) b, (GLfloat) alpha);
}

void GLRenderer::setClearDepth(Float depth) {
	glClearDepth((GLfloat) depth);
}

void GLRenderer::setClearColor(const Color3 &color) {
	glClearColor(
		(GLfloat) color[0],
		(GLfloat) color[1],
		(GLfloat) color[2], 1.0f
	);
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


void GLRenderer::debugString(const std::string &text) {
	if (GLEW_GREMEDY_string_marker)
		glStringMarkerGREMEDY(0, text.c_str());
}

MTS_IMPLEMENT_CLASS(GLRenderer, true, Renderer)
MTS_NAMESPACE_END
