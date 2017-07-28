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
#include <mitsuba/hw/glgeometry.h>

MTS_NAMESPACE_BEGIN

GLGeometry::GLGeometry(const TriMesh *mesh)
 : GPUGeometry(mesh) {
    m_id[0] = m_id[1] = 0;
}

void GLGeometry::init() {
    Assert(m_id[0] == 0 && m_id[1] == 0);
    glGenBuffers(2, m_id);
    refresh();
}

void GLGeometry::refresh() {
    Assert(m_id[0] != 0 && m_id[1] != 0);
    m_stride = 3;
    if (m_mesh->hasVertexNormals())
        m_stride += 3;
    if (m_mesh->hasVertexTexcoords())
        m_stride += 2;
    if (m_mesh->hasUVTangents())
        m_stride += 3;
    if (m_mesh->hasVertexColors())
        m_stride += 3;
    m_stride *= sizeof(GLfloat);

    size_t vertexCount = m_mesh->getVertexCount(), triCount = m_mesh->getTriangleCount();
    m_size[EVertexID] = (GLuint) (vertexCount * m_stride);
    m_size[EIndexID] = (GLuint) (triCount * sizeof(GLuint) * 3);

    Log(ETrace, "Uploading a GPU geometry object (\"%s\", " SIZE_T_FMT
        " vertices, " SIZE_T_FMT " triangles, %s)",
        getName().c_str(), vertexCount, triCount,
        memString(m_size[EVertexID] + m_size[EIndexID]).c_str());

    GLfloat *vertices = new GLfloat[vertexCount * m_stride/sizeof(GLfloat)];
    GLuint *indices = (GLuint *) m_mesh->getTriangles();
    const Point *sourcePositions = m_mesh->getVertexPositions();
    const Normal *sourceNormals = m_mesh->getVertexNormals();
    const Point2 *sourceTexcoords = m_mesh->getVertexTexcoords();
    const Color3 *sourceColors = m_mesh->getVertexColors();
    Vector *sourceTangents = NULL;

    if (m_mesh->hasUVTangents()) {
        /* Convert into per-vertex tangents */
        const TangentSpace *triTangents = m_mesh->getUVTangents();
        sourceTangents = new Vector[vertexCount];
        uint32_t *count = new uint32_t[vertexCount];
        memset(sourceTangents, 0, sizeof(Vector)*vertexCount);

        for (size_t i=0; i<triCount; ++i) {
            const Triangle &tri = m_mesh->getTriangles()[i];
            const TangentSpace &tangents = triTangents[i];
            for (int j=0; j<3; ++j) {
                sourceTangents[tri.idx[j]] += tangents.dpdu;
                ++count[tri.idx[j]];
            }
        }

        for (size_t i=0; i<vertexCount; ++i) {
            if (count[i] == 0)
                continue;
            sourceTangents[i] /= (Float) count[i];
        }

        delete[] count;
    }

    size_t pos = 0;
    for (size_t i=0; i<vertexCount; ++i) {
        vertices[pos++] = (GLfloat) sourcePositions[i].x;
        vertices[pos++] = (GLfloat) sourcePositions[i].y;
        vertices[pos++] = (GLfloat) sourcePositions[i].z;
        if (sourceNormals) {
            vertices[pos++] = (GLfloat) sourceNormals[i].x;
            vertices[pos++] = (GLfloat) sourceNormals[i].y;
            vertices[pos++] = (GLfloat) sourceNormals[i].z;
        }
        if (sourceTexcoords) {
            vertices[pos++] = (GLfloat) sourceTexcoords[i].x;
            vertices[pos++] = (GLfloat) sourceTexcoords[i].y;
        }
        if (sourceTangents) {
            vertices[pos++] = (GLfloat) sourceTangents[i].x;
            vertices[pos++] = (GLfloat) sourceTangents[i].y;
            vertices[pos++] = (GLfloat) sourceTangents[i].z;
        }
        if (sourceColors) {
            vertices[pos++] = (GLfloat) sourceColors[i][0];
            vertices[pos++] = (GLfloat) sourceColors[i][1];
            vertices[pos++] = (GLfloat) sourceColors[i][2];
        }
    }
    Assert(pos * sizeof(GLfloat) == m_stride * vertexCount);

    bind();

    glBufferData(GL_ARRAY_BUFFER, m_size[EVertexID], vertices, GL_STATIC_DRAW);
    if (GLEW_NV_vertex_buffer_unified_memory) {
        glGetBufferParameterui64vNV(GL_ARRAY_BUFFER, GL_BUFFER_GPU_ADDRESS_NV, &m_addr[EVertexID]);
        glMakeBufferResidentNV(GL_ARRAY_BUFFER, GL_READ_ONLY);
    }

    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_size[EIndexID], indices, GL_STATIC_DRAW);
    if (GLEW_NV_vertex_buffer_unified_memory) {
        glGetBufferParameterui64vNV(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_GPU_ADDRESS_NV, &m_addr[EIndexID]);
        glMakeBufferResidentNV(GL_ELEMENT_ARRAY_BUFFER, GL_READ_ONLY);
    }
    unbind();

    delete[] vertices;
    if (sourceTangents)
        delete[] sourceTangents;
}

void GLGeometry::bind() {
    glBindBuffer(GL_ARRAY_BUFFER, m_id[EVertexID]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_id[EIndexID]);
}

void GLGeometry::unbind() {
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void GLGeometry::cleanup() {
    Assert(m_id[0] != 0 && m_id[1] != 0);
    Log(ETrace, "Freeing GPU geometry object \"%s\"", getName().c_str());
    glDeleteBuffers(2, m_id);
    m_id[0] = m_id[1] = 0;
}

GLGeometry::~GLGeometry() {
    if (m_id[0] != 0 || m_id[1] != 0)
        cleanup();
}

MTS_IMPLEMENT_CLASS(GLGeometry, false, GPUGeometry)
MTS_NAMESPACE_END
