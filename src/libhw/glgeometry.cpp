/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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
#include <mitsuba/hw/glgeometry.h>

MTS_NAMESPACE_BEGIN

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
	m_stride = 3;
	if (m_mesh->getVertexNormals())
		m_stride += 3;
	if (m_mesh->getVertexTexcoords())
		m_stride += 2;
	if (m_mesh->getVertexTangents())
		m_stride += 3;
	if (m_mesh->getVertexColors())
		m_stride += 3;
	m_stride *= sizeof(GLfloat);

	m_vertexSize = (GLuint) (m_mesh->getVertexCount() * m_stride);
	m_indexSize = (GLuint) (m_mesh->getTriangleCount() * sizeof(GLuint) * 3);

	Log(ETrace, "Uploading a GPU geometry object (\"%s\", " SIZE_T_FMT 
		" vertices, " SIZE_T_FMT " triangles, %s)",
		getName().c_str(),
		m_mesh->getVertexCount(),
		m_mesh->getTriangleCount(),
		memString(m_vertexSize + m_indexSize).c_str());

	GLfloat *vertices = new GLfloat[m_mesh->getVertexCount()
		* m_stride/sizeof(GLfloat)];
	GLuint *indices = (GLuint *) m_mesh->getTriangles();
	const Point *sourcePositions = m_mesh->getVertexPositions();
	const Normal *sourceNormals = m_mesh->getVertexNormals();
	const Point2 *sourceTexcoords = m_mesh->getVertexTexcoords();
	const Spectrum *sourceColors = m_mesh->getVertexColors();
	const TangentSpace *sourceTangents = m_mesh->getVertexTangents();

	size_t pos = 0;
	for (size_t i=0; i<m_mesh->getVertexCount(); ++i) {
		vertices[pos++] = (float) sourcePositions[i].x;
		vertices[pos++] = (float) sourcePositions[i].y;
		vertices[pos++] = (float) sourcePositions[i].z;
		if (sourceNormals) {
			vertices[pos++] = (float) sourceNormals[i].x;
			vertices[pos++] = (float) sourceNormals[i].y;
			vertices[pos++] = (float) sourceNormals[i].z;
		}
		if (sourceTexcoords) {
			vertices[pos++] = (float) sourceTexcoords[i].x;
			vertices[pos++] = (float) sourceTexcoords[i].y;
		}
		if (sourceTangents) {
			vertices[pos++] = (float) sourceTangents[i].dpdu.x;
			vertices[pos++] = (float) sourceTangents[i].dpdu.y;
			vertices[pos++] = (float) sourceTangents[i].dpdu.z;
		}
		if (sourceColors) {
			Float r, g, b;
			sourceColors[i].toLinearRGB(r, g, b);
			vertices[pos++] = (float) r;
			vertices[pos++] = (float) g;
			vertices[pos++] = (float) b;
		}
	}
	Assert(pos * sizeof(GLfloat) == m_stride * m_mesh->getVertexCount());

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
	Log(ETrace, "Freeing GPU geometry object \"%s\"", getName().c_str());
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
