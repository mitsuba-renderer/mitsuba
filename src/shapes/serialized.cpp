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

#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

/**
 * Serialized model loader
 */
class SerializedMesh : public TriMesh {
public:
	SerializedMesh(const Properties &props) : TriMesh(props) {
		fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		m_name = filePath.stem();

		/* Load the geometry */
		Log(EInfo, "Loading geometry from \"%s\" ..", filePath.leaf().c_str());
		ref<FileStream> stream = new FileStream(filePath, FileStream::EReadOnly);
		stream->setByteOrder(Stream::ENetworkByteOrder);
		ref<TriMesh> mesh = new TriMesh(stream);
		m_triangleCount = mesh->getTriangleCount();
		m_vertexCount = mesh->getVertexCount();
		m_triangles = new Triangle[m_triangleCount];
		if (mesh->hasVertexNormals()) {
			m_normals = new Normal[m_vertexCount];
			memcpy(m_normals, mesh->getVertexNormals(), sizeof(Normal) * m_vertexCount);
		}
		if (mesh->hasVertexTexcoords()) {
			m_texcoords = new Point2[m_vertexCount];
			memcpy(m_texcoords, mesh->getVertexTexcoords(), sizeof(Point2) * m_vertexCount);
		}
		if (mesh->hasVertexColors()) {
			m_colors = new Spectrum[m_vertexCount];
			memcpy(m_colors, mesh->getVertexColors(), sizeof(Spectrum) * m_vertexCount);
		}
		memcpy(m_triangles, mesh->getTriangles(), sizeof(Triangle) * m_triangleCount);

		if (!m_objectToWorld.isIdentity()) {
			for (size_t i=0; i<m_vertexCount; ++i)
				m_positions[i] = m_objectToWorld(m_positions[i]);
			if (m_normals) {
				for (size_t i=0; i<m_vertexCount; ++i)
					m_normals[i] = m_objectToWorld(m_normals[i]);
			}
		}
	}

	SerializedMesh(Stream *stream, InstanceManager *manager) : TriMesh(stream, manager) {
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(SerializedMesh, false, TriMesh)
MTS_EXPORT_PLUGIN(SerializedMesh, "Serialized mesh loader");
MTS_NAMESPACE_END
