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
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/glgeometry.h>

MTS_NAMESPACE_BEGIN

extern GLEWContext *glewGetContext();

GLGeometry::GLGeometry(const TriMesh *mesh) 
 : GPUGeometry(mesh), m_vertexID(0), m_indexID(0) {
}

void GLGeometry::init() {
	Assert(m_vertexID == 0 && m_indexID == 0);
	glGenBuffers(1, &m_vertexID);
	glGenBuffers(1, &m_indexID);
	refresh();
	
}
void GLGeometry::refresh() {
	Assert(m_vertexID != 0 && m_indexID != 0);
	m_vertexSize = m_mesh->getVertexCount() * sizeof(GLfloat) * 11;
	m_indexSize = m_mesh->getTriangleCount() * sizeof(GLuint) * 3;

	Log(EDebug, "Uploading a GPU geometry object (\"%s\", " SIZE_T_FMT 
		" vertices, " SIZE_T_FMT " triangles, %.1f KiB)",
		getName().c_str(),
		m_mesh->getVertexCount(),
		m_mesh->getTriangleCount(),
		(m_vertexSize + m_indexSize) / 1024.0f);

	GLfloat *vertices = new GLfloat[m_mesh->getVertexCount() * 11];
	GLuint *indices = (GLuint *) m_mesh->getTriangles();
	const Vertex *source = m_mesh->getVertexBuffer();
	int pos = 0;
	for (size_t i=0; i<m_mesh->getVertexCount(); ++i) {
		const Vertex &vtx = source[i];
		vertices[pos++] = (float) vtx.v.x;
		vertices[pos++] = (float) vtx.v.y;
		vertices[pos++] = (float) vtx.v.z;
		vertices[pos++] = (float) vtx.n.x;
		vertices[pos++] = (float) vtx.n.y;
		vertices[pos++] = (float) vtx.n.z;
		vertices[pos++] = (float) vtx.uv.x;
		vertices[pos++] = (float) vtx.uv.y;
		vertices[pos++] = (float) vtx.dpdu.x;
		vertices[pos++] = (float) vtx.dpdu.y;
		vertices[pos++] = (float) vtx.dpdu.y;
	}

	bind();

	bool bindless = glewIsSupported("GL_NV_vertex_buffer_unified_memory");

	glBufferData(GL_ARRAY_BUFFER, m_vertexSize, vertices, GL_STATIC_DRAW);
	if (bindless) {
		glGetBufferParameterui64vNV(GL_ARRAY_BUFFER, GL_BUFFER_GPU_ADDRESS_NV, &m_vertexAddr);
		glMakeBufferResidentNV(GL_ARRAY_BUFFER, GL_READ_ONLY);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indexSize, indices, GL_STATIC_DRAW);
	if (bindless) {
		glGetBufferParameterui64vNV(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_GPU_ADDRESS_NV, &m_indexAddr);
		glMakeBufferResidentNV(GL_ELEMENT_ARRAY_BUFFER, GL_READ_ONLY);
	}
	unbind();

	delete[] vertices;
}

void GLGeometry::bind() {
	glBindBuffer(GL_ARRAY_BUFFER, m_vertexID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexID);
}

void GLGeometry::unbind() {
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void GLGeometry::cleanup() {
	Assert(m_vertexID != 0 && m_indexID != 0);
	Log(EDebug, "Freeing GPU geometry object \"%s\"", getName().c_str());
	glDeleteBuffers(1, &m_vertexID);
	glDeleteBuffers(1, &m_indexID);
	m_vertexID = m_indexID = 0;
}

GLGeometry::~GLGeometry() {
	if (m_vertexID != 0 || m_indexID != 0)
		cleanup();
}

MTS_IMPLEMENT_CLASS(GLGeometry, false, GPUGeometry)
MTS_NAMESPACE_END
