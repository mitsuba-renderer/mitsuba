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
		bool hasNormals = false, hasTexcoords = false;
		bool firstVertex = true;

		std::string name = m_name;

		while (is >> buf) {
			if (buf == "v") {
				/* Parse + transform vertices */
				Point p;
				is >> p.x >> p.y >> p.z;
				vertices.push_back(m_objectToWorld(p));
				if (firstVertex) {
					if (triangles.size() > 0) {
						generateGeometry(name, vertices, normals, texcoords, 
							triangles, hasNormals, hasTexcoords);
						triangles.clear();
					}
					hasNormals = false;
					hasTexcoords = false;
					firstVertex = false;
				}
			} else if (buf == "mtllib") {
				cout << "Got mtllib" << endl;
			} else if (buf == "vn") {
				Normal n;
				is >> n.x >> n.y >> n.z;
				if (!n.isZero())
					normals.push_back(normalize(m_objectToWorld(n)));
				else
					normals.push_back(n);
				hasNormals = true;
			} else if (buf == "g") {
				std::string line;
				std::getline(is, line);
				name = line.substr(1, line.length()-2);
			} else if (buf == "vt") {
				std::string line;
				Float u, v, w;
				std::getline(is, line);
				std::istringstream iss(line);
				iss >> u >> v >> w;
				texcoords.push_back(Point2(u, v));
				hasTexcoords = true;
			} else if (buf == "f") {
				std::string line, tmp;
				std::getline(is, line);
				std::istringstream iss(line);
				firstVertex = true;
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

		generateGeometry(name, vertices, normals, texcoords, 
			triangles, hasNormals, hasTexcoords);
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


	void generateGeometry(const std::string &name,
			const std::vector<Point> &vertices,
			const std::vector<Normal> &normals,
			const std::vector<Point2> &texcoords,
			const std::vector<OBJTriangle> &triangles,
			bool hasNormals, bool hasTexcoords) {
		if (triangles.size() == 0)
			return;
		Log(EInfo, "Creating geometry \"%s\"", name.c_str());
		std::vector<Vertex> vertexBuffer(vertices.size());
		std::vector<bool> touched(vertices.size());
		for (unsigned int i=0; i<vertices.size(); i++) {
			vertexBuffer[i].v = vertices[i];
			vertexBuffer[i].n = Normal(0, 0, 1);
			touched[i] = false;
		}

		/* Collapse the mesh into a more usable form */
		Triangle *triangleArray = new Triangle[triangles.size()];
		for (unsigned int i=0; i<triangles.size(); i++) {
			Triangle tri;
			for (unsigned int j=0; j<3; j++) {
				unsigned int vertexId = triangles[i].v[j];
				unsigned int normalId = triangles[i].n[j];
				unsigned int uvId = triangles[i].uv[j];

				if (touched[vertexId] == false) {
					if (hasNormals)
						vertexBuffer[vertexId].n = normals.at(normalId);
					if (hasTexcoords)
						vertexBuffer[vertexId].uv = texcoords.at(uvId);
					touched[vertexId] = true;
					tri.idx[j] = vertexId;
				} else {
					bool safe = true;
					if (hasNormals && vertexBuffer.at(vertexId).n 
							!= normals.at(normalId))
						safe = false;
					if (hasTexcoords && vertexBuffer.at(vertexId).uv 
							!= texcoords.at(uvId))
						safe = false;
					if (!safe) {
						/* This vertex was used in conjunction with a different
						   normal / texture coordinate. Now we have to duplicate
						   it. */
						Vertex vertex;
						vertex.v = vertices.at(vertexId);
						if (hasNormals)
							vertex.n = normals.at(normalId);
						if (hasTexcoords)
							vertex.uv = texcoords.at(uvId);
						tri.idx[j] = vertexBuffer.size();
						vertexBuffer.push_back(vertex);
					} else {
						tri.idx[j] = vertexId;
					}
				}
			}
			triangleArray[i] = tri;
		}

		Vertex *vertexArray = new Vertex[vertexBuffer.size()];
		for (unsigned int i=0; i<vertexBuffer.size(); i++) 
			vertexArray[i] = vertexBuffer[i];

		ref<TriMesh> mesh = new TriMesh(name, m_worldToObject, triangleArray, 
			triangles.size(), vertexArray, vertexBuffer.size());
		mesh->incRef();
		mesh->calculateTangentSpaceBasis(hasNormals, hasTexcoords);
		m_meshes.push_back(mesh);
	}
	
	WavefrontOBJ(Stream *stream, InstanceManager *manager) : TriMesh(stream, manager) {
	}

	virtual ~WavefrontOBJ() {
		for (size_t i=0; i<m_meshes.size(); ++i)
			m_meshes[i]->decRef();
	}
	
	void configure() {
		Shape::configure();

		m_aabb.reset();
		for (size_t i=0; i<m_meshes.size(); ++i) {
			m_meshes[i]->configure();
			m_aabb.expandBy(m_meshes[i]->getAABB());
		}
	}
	
	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(BSDF::m_theClass)) {
			m_bsdf = static_cast<BSDF *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) 
				m_meshes[i]->addChild(name, child);
		} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
			Assert(m_luminaire == NULL && m_meshes.size() == 1);
			m_luminaire = static_cast<Luminaire *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) {
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else if (cClass->derivesFrom(Subsurface::m_theClass)) {
			Assert(m_subsurface == NULL);
			m_subsurface = static_cast<Subsurface *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) { 
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else {
			Shape::addChild(name, child);
		}
	}

	bool isCompound() const {
		return true;
	}

	Shape *getElement(int index) {
		if (index >= (int) m_meshes.size())
			return NULL;
		return m_meshes[index];
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<TriMesh *> m_meshes;
};

MTS_IMPLEMENT_CLASS_S(WavefrontOBJ, false, TriMesh)
MTS_EXPORT_PLUGIN(WavefrontOBJ, "OBJ triangle mesh loader");
MTS_NAMESPACE_END
