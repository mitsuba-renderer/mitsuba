#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

/**
 * Serialized model loader
 */
class SerializedMesh : public TriMesh {
public:
	SerializedMesh(const Properties &props) : TriMesh(props) {
		m_name = FileResolver::getInstance()->resolve(props.getString("filename"));

		/* Load the geometry */
		Log(EInfo, "Loading geometry from \"%s\" ..", m_name.c_str());
		ref<FileStream> stream = new FileStream(m_name, FileStream::EReadOnly);
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
