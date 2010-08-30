#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/zstream.h>

#define MTS_FILEFORMAT_HEADER 0x041C
#define MTS_FILEFORMAT_VERSION_V1 0x01
#define MTS_FILEFORMAT_VERSION_V2 0x02

MTS_NAMESPACE_BEGIN

TriMesh::TriMesh(size_t triangleCount, size_t vertexCount) 
 : Shape(Properties()), m_triangleCount(triangleCount),
 		m_vertexCount(vertexCount), m_flipNormals(false) {
	m_triangles = new Triangle[m_triangleCount];
	m_vertexBuffer = new Vertex[m_vertexCount];
}
	
TriMesh::TriMesh(const std::string &name, Transform worldToObject, Triangle *triangles, 
	size_t triangleCount, Vertex *vertexBuffer, size_t vertexCount) 
	: Shape(Properties()), m_triangles(triangles), m_triangleCount(triangleCount), 
	m_vertexBuffer(vertexBuffer), m_vertexCount(vertexCount), m_flipNormals(false) {
	m_name = name;
	m_worldToObject = worldToObject;
	m_objectToWorld.inverse();
}

TriMesh::TriMesh(const Properties &props) 
 : Shape(props) {
	m_flipNormals = props.getBoolean("flipNormals", false);
	m_triangles = NULL;
	m_vertexBuffer = NULL;
}

TriMesh::TriMesh(Stream *stream, InstanceManager *manager) 
	: Shape(stream, manager) {
	Assert(sizeof(Vertex) == 14*sizeof(Float));
	Assert(sizeof(Triangle) == 3*sizeof(int));
	m_vertexCount = (size_t) stream->readULong();
	m_vertexBuffer = new Vertex[m_vertexCount];
	stream->readFloatArray(reinterpret_cast<Float *>(m_vertexBuffer), 
		m_vertexCount * sizeof(Vertex)/sizeof(Float));
	m_triangleCount = (size_t) stream->readULong();
	m_triangles = new Triangle[m_triangleCount];
	stream->readIntArray(reinterpret_cast<int *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(int));
	m_flipNormals = false;
	configure();
}
	
TriMesh::TriMesh(Stream *_stream) : Shape(Properties()) {
	ref<Stream> stream = _stream;

	Assert(sizeof(Vertex) == 14*sizeof(Float));
	Assert(sizeof(Triangle) == 3*sizeof(int));	
	if (stream->getByteOrder() != Stream::ENetworkByteOrder) 
		Log(EError, "Tried to unserialize a shape from a stream, "
		"which was not previously set to network byte order!");
#if defined(SINGLE_PRECISION)
	bool doublePrecision = false;
#else
	bool doublePrecision = true;
#endif
	
	if (stream->readShort() != MTS_FILEFORMAT_HEADER)
		Log(EError, "Encountered an invalid file format!");

	short version = stream->readShort();

	if (version != MTS_FILEFORMAT_VERSION_V1 &&
		version != MTS_FILEFORMAT_VERSION_V2)
		Log(EError, "Encountered an incompatible file version!");

	if (version == MTS_FILEFORMAT_VERSION_V2)
		stream = new ZStream(stream);

	bool fileDoublePrecision = stream->readBool();
	m_vertexCount = (size_t) stream->readULong();
	m_vertexBuffer = new Vertex[m_vertexCount];
	size_t numEntries = m_vertexCount * sizeof(Vertex)/sizeof(Float);
	Float *target = reinterpret_cast<Float *>(m_vertexBuffer);

	if ((doublePrecision && fileDoublePrecision) ||
		(!doublePrecision && !fileDoublePrecision)) {
		/* Precision matches - load directly into memory */
		stream->readFloatArray(target, numEntries);
	} else if (fileDoublePrecision) {
		/* Double -> Single conversion */
		double *temp = new double[numEntries];
		stream->readDoubleArray(temp, numEntries);
		for (size_t i=0; i<numEntries; ++i)
			target[i] = (Float) temp[i];
		delete[] temp;
	} else {
		/* Single -> Double conversion */
		float *temp = new float[numEntries];
		stream->readSingleArray(temp, numEntries);
		for (size_t i=0; i<numEntries; ++i)
			target[i] = (Float) temp[i];
		delete[] temp;
	}

	m_triangleCount = (size_t) stream->readULong();
	m_triangles = new Triangle[m_triangleCount];
	stream->readIntArray(reinterpret_cast<int *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(int));
	m_flipNormals = false;
	configure();
}

TriMesh::~TriMesh() {
	if (m_triangles)
		delete[] m_triangles;
	if (m_vertexBuffer)
		delete[] m_vertexBuffer;
}

void TriMesh::configure() {
	Shape::configure();

	if (m_areaPDF.isReady())
		return;
	for (size_t i=0; i<m_triangleCount; i++)
		m_areaPDF.put(m_triangles[i].surfaceArea(m_vertexBuffer));
	m_surfaceArea = m_areaPDF.build();
	m_invSurfaceArea = 1.0f / m_surfaceArea;
	/* Generate a bounding sphere */
	m_aabb.reset();
	for (size_t i=0; i<m_vertexCount; i++) 
		m_aabb.expandBy(m_vertexBuffer[i].v);
	m_bsphere.center = m_aabb.getCenter();
	for (size_t i=0; i<m_vertexCount; i++) 
		m_bsphere.expandBy(m_vertexBuffer[i].v);
}

Float TriMesh::sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
	Point2 newSeed = sample;
	int index = m_areaPDF.sampleReuse(newSeed.y);
	sRec.p = m_triangles[index].sample(m_vertexBuffer, sRec.n, newSeed);
	return m_invSurfaceArea;
}

void TriMesh::calculateTangentSpaceBasis(bool hasNormals, bool hasTexCoords, bool complain) {
	/* Calculate smooth normals if there aren't any */
	if (!hasNormals) {
		for (unsigned int i=0; i<m_vertexCount; i++)
			m_vertexBuffer[i].n = Normal(0.0, 0.0f, 0.0f);
		for (unsigned int i=0; i<m_triangleCount; i++) {
			const Point &v0 = m_vertexBuffer[m_triangles[i].idx[0]].v;
			const Point &v1 = m_vertexBuffer[m_triangles[i].idx[1]].v;
			const Point &v2 = m_vertexBuffer[m_triangles[i].idx[2]].v;
			Normal n = Normal(cross(v1 - v0, v2 - v0));
			Float length = n.length();
			if (length != 0)
				n /= length;
			else if (complain)
				Log(EWarn, "Invalid geometry - cannot calculate normals");
			m_vertexBuffer[m_triangles[i].idx[0]].n += n;
			m_vertexBuffer[m_triangles[i].idx[1]].n += n;
			m_vertexBuffer[m_triangles[i].idx[2]].n += n;
		}
		for (unsigned int i=0; i<m_vertexCount; i++) {
			Float length = m_vertexBuffer[i].n.length();
			if (length != 0) {
				m_vertexBuffer[i].n /= length;
			} else {
				/* Choose some bogus value */
				if (complain)
					Log(EWarn, "Could not generate correct mesh normals!");
				m_vertexBuffer[i].n = Normal(1, 0, 0);
			}
		}
	}

	if (m_flipNormals) {
		for (unsigned int i=0; i<m_vertexCount; i++)
			m_vertexBuffer[i].n *= -1;
	}

	if (!hasTexCoords) {
		/* At least create some kind of tangent space basis (fair enough
		   for isotropic BxDFs) */
		for (unsigned int i=0; i<m_vertexCount; i++) {
			if (m_vertexBuffer[i].n.isZero()) {
				if (complain)
					Log(EWarn, "Mesh has zero normals!");
				m_vertexBuffer[i].n = Normal(1, 0, 0);
			}
			coordinateSystem(m_vertexBuffer[i].n, m_vertexBuffer[i].dpdu, m_vertexBuffer[i].dpdv);
		}
	} else {
		/* No. of triangles sharing a vertex */
		int *sharers = new int[m_vertexCount];

		for (unsigned int i=0; i<m_vertexCount; i++) {
			m_vertexBuffer[i].dpdu = Vector(0.0, 0.0f, 0.0f);
			m_vertexBuffer[i].dpdv = Vector(0.0, 0.0f, 0.0f);
			if (m_vertexBuffer[i].n.isZero()) {
				if (complain)
					Log(EWarn, "Mesh has zero normals!");
				m_vertexBuffer[i].n = Normal(1, 0, 0);
			}
			sharers[i] = 0;
		}

		for (unsigned int i=0; i<m_triangleCount; i++) {
			unsigned int idx0 = m_triangles[i].idx[0],
						 idx1 = m_triangles[i].idx[1],
						 idx2 = m_triangles[i].idx[2];
			const Point &v0 = m_vertexBuffer[idx0].v;
			const Point &v1 = m_vertexBuffer[idx1].v;
			const Point &v2 = m_vertexBuffer[idx2].v;
			const Point2 &uv0 = m_vertexBuffer[idx0].uv;
			const Point2 &uv1 = m_vertexBuffer[idx1].uv;
			const Point2 &uv2 = m_vertexBuffer[idx2].uv;

			Vector dP1 = v1 - v0, dP2 = v2 - v0;
			Vector2 dUV1 = uv1 - uv0, dUV2 = uv2 - uv0;

			Float invDet = 1.0f, determinant = dUV1.x * dUV2.y - dUV1.y * dUV2.x;
			if (determinant != 0)
				invDet = 1.0f / determinant;

			Vector dpdu = ( dUV2.y * dP1 - dUV1.y * dP2) * invDet;
			Vector dpdv = (-dUV2.x * dP1 + dUV1.x * dP2) * invDet;

			if (dpdu.length() == 0.0f) {
				/* Recovery - required to recover from invalid geometry */
				Normal n = Normal(cross(v1 - v0, v2 - v0));
				Float length = n.length();
				if (length != 0)
					n /= length;
				else if (complain)
					Log(EWarn, "Mesh contains invalid geometry!");
				dpdu = cross(n, dpdv);

				if (dpdu.length() == 0.0f) {
					/* At least create some kind of tangent space basis 
					   (fair enough for isotropic BxDFs) */
					coordinateSystem(n, dpdu, dpdv);
				}
			}

			if (dpdv.length() == 0.0f) {
				Normal n = Normal(cross(v1 - v0, v2 - v0));
				Float length = n.length();
				if (length != 0)
					n /= length;
				else if (complain)
					Log(EWarn, "Mesh contains invalid geometry!");
				dpdv = cross(dpdu, n);

				if (dpdv.length() == 0.0f) {
					/* At least create some kind of tangent space basis 
					   (fair enough for isotropic BxDFs) */
					coordinateSystem(n, dpdu, dpdv);
				}
			}

			m_vertexBuffer[idx0].dpdu += dpdu;
			m_vertexBuffer[idx1].dpdu += dpdu;
			m_vertexBuffer[idx2].dpdu += dpdu;
			m_vertexBuffer[idx0].dpdv += dpdv;
			m_vertexBuffer[idx1].dpdv += dpdv;
			m_vertexBuffer[idx2].dpdv += dpdv;
			sharers[idx0]++; sharers[idx1]++; sharers[idx2]++;
		}
		/* Orthogonalization + Normalization pass */
		for (unsigned int i=0; i<m_vertexCount; i++) {
			Vector dpdu = m_vertexBuffer[i].dpdu;
			Vector dpdv = m_vertexBuffer[i].dpdv;
 
			if (dpdu.length() == 0.0f || dpdv.length() == 0.0f) {
				/* At least create some kind of tangent space basis 
				   (fair enough for isotropic BxDFs) */
				coordinateSystem(m_vertexBuffer[i].n, dpdu, dpdv);
			} else {
				if (sharers[i] > 0) {
					dpdu /= (Float) sharers[i];
					dpdv /= (Float) sharers[i];
				}
			}

			m_vertexBuffer[i].dpdu = dpdu;
			m_vertexBuffer[i].dpdv = dpdv;
		}
		delete[] sharers;
	}
}

void TriMesh::serialize(Stream *stream, InstanceManager *manager) const {
	Shape::serialize(stream, manager);

	Assert(sizeof(Vertex) == 14*sizeof(Float));
	Assert(sizeof(Triangle) == 3*sizeof(int));
	stream->writeULong(m_vertexCount);
	stream->writeFloatArray(reinterpret_cast<Float *>(m_vertexBuffer), 
		m_vertexCount * sizeof(Vertex)/sizeof(Float));
	stream->writeULong(m_triangleCount);
	stream->writeIntArray(reinterpret_cast<int *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(int));
}

void TriMesh::serialize(Stream *_stream) const {
	ref<Stream> stream = _stream;

	Assert(sizeof(Vertex) == 14*sizeof(Float));
	Assert(sizeof(Triangle) == 3*sizeof(int));

	if (stream->getByteOrder() != Stream::ENetworkByteOrder) 
		Log(EError, "Tried to unserialize a shape from a stream, "
			"which was not previously set to network byte order!");

#if defined(SINGLE_PRECISION)
	bool doublePrecision = false;
#else
	bool doublePrecision = true;
#endif

	stream->writeShort(MTS_FILEFORMAT_HEADER);
	stream->writeShort(MTS_FILEFORMAT_VERSION_V2);

	stream = new ZStream(stream, Z_BEST_COMPRESSION);
	stream->writeBool(doublePrecision);
	stream->writeULong(m_vertexCount);
	stream->writeFloatArray(reinterpret_cast<Float *>(m_vertexBuffer), 
		m_vertexCount * sizeof(Vertex)/sizeof(Float));
	stream->writeULong(m_triangleCount);
	stream->writeIntArray(reinterpret_cast<int *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(int));
}

std::string TriMesh::toString() const {
	std::ostringstream oss;
	oss << getClass()->getName() << "[" << endl
		<< "  name = \"" << m_name<< "\"," << endl
		<< "  triangleCount = " << m_triangleCount << "," << endl
		<< "  vertexCount = " << m_vertexCount << "," << endl
		<< "  flipNormals = " << m_flipNormals << "," << endl
		<< "  surfaceArea = " << m_surfaceArea << "," << endl
		<< "  aabb = " << m_aabb.toString() << "," << endl
		<< "  bsphere = " << m_bsphere.toString() << "," << endl
		<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
		<< "  subsurface = " << indent(m_subsurface.toString()) << "," << endl
		<< "  luminaire = " << indent(m_luminaire.toString()) << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(TriMesh, false, Shape)
MTS_NAMESPACE_END
