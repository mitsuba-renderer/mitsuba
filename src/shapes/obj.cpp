#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/fresolver.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

/**
 * Wavefront OBJ triangle mesh loader
 */
class WavefrontOBJ : public TriMesh {
public:
	struct OBJTriangle {
		unsigned int v[3];
		unsigned int n[3];
		unsigned int uv[3];
	};

	WavefrontOBJ(const Properties &props) : TriMesh(props) {
		m_name = FileResolver::getInstance()->resolve(props.getString("filename"));

		/* Load the geometry */
		Log(EInfo, "Loading geometry from \"%s\" ..", m_name.c_str());
		std::ifstream is(m_name.c_str());
		if (is.bad() || is.fail())
			Log(EError, "Geometry file '%s' not found!", m_name.c_str());

		std::string buf;
		std::vector<Point> vertices;
		std::vector<Normal> normals;
		std::vector<Point2> texcoords;
		std::vector<OBJTriangle> triangles;

		bool hasNormals = false;
		bool hasTexCoords = false;

		while (is >> buf) {
			if (buf == "v") {
				/* Parse + transform vertices */
				Point p;
				is >> p.x >> p.y >> p.z;
				vertices.push_back(m_objectToWorld(p));
			} else if (buf == "vn") {
				Normal n;
				is >> n.x >> n.y >> n.z;
				if (!n.isZero())
					normals.push_back(normalize(m_objectToWorld(n)));
				else
					normals.push_back(n);
				hasNormals = true;
			} else if (buf == "vt") {
				std::string line;
				Float u, v, w;
				std::getline(is, line);
				std::istringstream iss(line);
				iss >> u >> v >> w;
				hasTexCoords = true;
				texcoords.push_back(Point2(u, v));
			} else if (buf == "f") {
				std::string line, tmp;
				std::getline(is, line);
				std::istringstream iss(line);
				OBJTriangle t;
				iss >> tmp; parse(t, 0, tmp);
				iss >> tmp; parse(t, 1, tmp);
				iss >> tmp; parse(t, 2, tmp);
				triangles.push_back(t);
				if (iss >> tmp) {
					parse(t, 1, tmp);
					std::swap(t.v[0], t.v[1]);
					std::swap(t.uv[0], t.uv[1]);
					std::swap(t.n[0], t.n[1]);
					triangles.push_back(t);
				}
			}
		}

		std::vector<Vertex> vertexBuffer(vertices.size());
		std::vector<bool> touched(vertices.size());
		for (unsigned int i=0; i<vertices.size(); i++) {
			vertexBuffer[i].v = vertices[i];
			touched[i] = false;
		}

		/* Collapse the mesh into a more usable form */
		m_triangles = new Triangle[triangles.size()];
		for (unsigned int i=0; i<triangles.size(); i++) {
			Triangle tri;
			for (unsigned int j=0; j<3; j++) {
				unsigned int vertexId = triangles[i].v[j];
				unsigned int normalId = triangles[i].n[j];
				unsigned int uvId = triangles[i].uv[j];

				if (touched[vertexId] == false) {
					if (hasNormals)
						vertexBuffer[vertexId].n = normals[normalId];
					if (hasTexCoords)
						vertexBuffer[vertexId].uv = texcoords[uvId];
					touched[vertexId] = true;
					tri.idx[j] = vertexId;
				} else {
					bool safe = true;
					if (hasNormals && vertexBuffer[vertexId].n 
							!= normals[normalId])
						safe = false;
					if (hasTexCoords && vertexBuffer[vertexId].uv 
							!= texcoords[uvId])
						safe = false;
					if (!safe) {
						/* This vertex was used in conjunction with a different
						   normal / texture coordinate. Now we have to duplicate
						   it. */
						Vertex vertex;
						vertex.v = vertices[vertexId];
						if (hasNormals)
							vertex.n = normals[normalId];
						if (hasTexCoords)
							vertex.uv = texcoords[uvId];
						tri.idx[j] = vertexBuffer.size();
						vertexBuffer.push_back(vertex);
					} else {
						tri.idx[j] = vertexId;
					}
				}
			}
			m_triangles[i] = tri;
		}
		m_vertexBuffer = new Vertex[vertexBuffer.size()];
		for (unsigned int i=0; i<vertexBuffer.size(); i++)
			m_vertexBuffer[i] = vertexBuffer[i];
		m_triangleCount = triangles.size();
		m_vertexCount = vertexBuffer.size();

		calculateTangentSpaceBasis(hasNormals, hasTexCoords);
	}
	
	WavefrontOBJ(Stream *stream, InstanceManager *manager) : TriMesh(stream, manager) {
	}

	void parse(OBJTriangle &t, int i, const std::string &str) {
		std::vector<std::string> tokens = tokenize(str, "/");
		if (tokens.size() == 1) {
			t.v[i] = atoi(tokens[0].c_str())-1;
		} else if (tokens.size() == 2) {
			if (str.find("//") == std::string::npos) {
				t.v[i]  = atoi(tokens[0].c_str())-1;
				t.uv[i] = atoi(tokens[1].c_str())-1;
			} else {
				t.v[i] = atoi(tokens[0].c_str())-1;
				t.n[i] = atoi(tokens[1].c_str())-1;
			}
		} else if (tokens.size() == 3) {
			t.v[i] = atoi(tokens[0].c_str())-1;
			t.uv[i] = atoi(tokens[1].c_str())-1;
			t.n[i] = atoi(tokens[2].c_str())-1;
		} else {
			Log(EError, "Invalid OBJ face format!");
		}
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(WavefrontOBJ, false, TriMesh)
MTS_EXPORT_PLUGIN(WavefrontOBJ, "OBJ triangle mesh loader");
MTS_NAMESPACE_END
