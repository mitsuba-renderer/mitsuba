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

MTS_NAMESPACE_BEGIN

/**
 * Serialized model loader
 */
class SerializedMesh : public TriMesh {
public:
	SerializedMesh(const Properties &props) : TriMesh(props) {
		m_name = props.getString("filename");
		std::string filePath = FileResolver::getInstance()->resolve(m_name);

		/* Load the geometry */
		Log(EInfo, "Loading geometry from \"%s\" ..", m_name.c_str());
		ref<FileStream> stream = new FileStream(filePath.c_str(), FileStream::EReadOnly);
		stream->setByteOrder(Stream::ENetworkByteOrder);
		ref<TriMesh> mesh = new TriMesh(stream);
		m_triangleCount = mesh->getTriangleCount();
		m_vertexCount = mesh->getVertexCount();
		m_vertexBuffer = new Vertex[m_vertexCount];
		m_triangles = new Triangle[m_triangleCount];
		memcpy(m_vertexBuffer, mesh->getVertexBuffer(), sizeof(Vertex) * m_vertexCount);
		memcpy(m_triangles, mesh->getTriangles(), sizeof(Triangle) * m_triangleCount);

		if (!m_objectToWorld.isIdentity()) {
			for (size_t i=0; i<m_vertexCount; ++i) {
				Vertex &vertex = m_vertexBuffer[i];
				vertex.v = m_objectToWorld(vertex.v);
				vertex.n = m_objectToWorld(vertex.n);
				vertex.dpdu = m_objectToWorld(vertex.dpdu);
				vertex.dpdv = m_objectToWorld(vertex.dpdv);
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
