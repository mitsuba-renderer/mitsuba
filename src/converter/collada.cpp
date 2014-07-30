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

#if defined(COLLADA_DOM_2_4)
	#ifndef COLLADA_DOM_SUPPORT141
		#define COLLADA_DOM_SUPPORT141
	#endif
	#ifdef COLLADA_DOM_USING_141
		#undef COLLADA_DOM_USING_141
	#endif
	#define COLLADA_DOM_USING_141
#endif

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/version.h>
#include <mitsuba/core/track.h>
#include <dom/domCOLLADA.h>
#include <dae.h>
#include <dae/daeErrorHandler.h>
#include <dom/domProfile_COMMON.h>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <set>

#if defined(__OSX__)
#include <OpenGL/glu.h>
#else
#if defined(WIN32)
# include <windows.h>
#endif
#include <GL/glu.h>
#endif

#ifndef WIN32
#define __stdcall
#endif

#include "converter.h"

typedef std::map<std::string, std::string> StringMap;
typedef std::map<std::string, int> RefCountMap;
typedef std::multimap<std::string, AbstractAnimationTrack *> AnimationMap;

struct ColladaContext {
	GeometryConverter *cvt;
	RefCountMap refCountMap;
	AnimationMap animations;
	std::set<std::string> serializedGeometry;
	fs::path texturesDirectory;
	fs::path meshesDirectory;
	StringMap idToTexture, fileToId;
	std::ostream &os;
	int trackIndex;

	inline ColladaContext(std::ostream &os) : os(os) { }
};

enum ESourceType {
	EPosition = 0,
	ENormal = 1,
	EUV = 2,
	EVertexColor = 3,
	ELast
};

struct Vec4 {
	Float x, y, z, w;

	inline Vec4(Float x=0, Float y=0, Float z=0, Float w=0)
		: x(x), y(y), z(z), w(w) {
	}

	inline void operator+=(const Vec4 &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
	}

	inline Vec4 operator*(Float f) const {
		return Vec4(x * f, y * f, z * f, w * f);
	}

	inline Float operator[](int i) const {
		return (&x)[i];
	}

	inline Float &operator[](int i) {
		return (&x)[i];
	}

	inline Point toPoint() const {
		return Point(x, y, z);
	}

	inline Vector toVector() const {
		return Vector(x, y, z);
	}

	inline Normal toNormal() const {
		return Normal(x, y, z);
	}

	inline Point2 toPoint2() const {
		return Point2(x, y);
	}
};

struct VertexData {
	size_t nSources;
	bool hasNormals;
	bool hasUVs;
	GLdouble *glPos;
	std::vector<Vec4 *> data;
	std::vector<int> typeToOffset;
	std::vector<int> typeToOffsetInStream;
	std::vector<size_t> typeToCount;

	VertexData() : glPos(NULL) {
	}

	virtual ~VertexData() {
		for (size_t i=0; i<data.size(); ++i) {
			if (data[i])
				delete[] data[i];
		}
		if (glPos)
			delete[] glPos;
	}
};

/* This code is not thread-safe for now */
GLUtesselator *tess = NULL;
std::vector<domUint> tess_data;
std::vector<domUint *> tess_cleanup;
VertexData *tess_vdata = NULL;
size_t tess_nSources;

VertexData *fetchVertexData(Transform transform,
		const domInputLocal_Array &vertInputs,
		const domInputLocalOffset_Array &inputs) {
	VertexData *result = new VertexData();

	result->hasNormals = false;
	result->hasUVs = false;
	result->nSources = inputs.getCount();
	result->data.resize(result->nSources);
	for (size_t i=0; i<result->nSources; ++i)
		result->data[i] = NULL;
	result->typeToOffset.resize(ELast);
	result->typeToCount.resize(ELast);
	result->typeToOffsetInStream.resize(ELast);
	for (int i=0; i<ELast; ++i) {
		result->typeToOffset[i] = result->typeToOffsetInStream[i] = -1;
		result->typeToCount[i] = 0;
	}
	int vertInputIndex = 0;

	for (size_t i=0; i<inputs.getCount(); ++i) {
		int offsetInStream = (int) inputs[i]->getOffset(),
		    offset = offsetInStream;

		daeURI &sourceRef = inputs[i]->getSource();
		sourceRef.resolveElement();
		domSource *source = daeSafeCast<domSource>(sourceRef.getElement());
		std::string semantic = inputs[i]->getSemantic();

		if (semantic == "VERTEX") {
			sourceRef = vertInputs[vertInputIndex]->getSource();
			sourceRef.resolveElement();
			source = daeSafeCast<domSource>(sourceRef.getElement());
			semantic = vertInputs[vertInputIndex]->getSemantic();

			if (vertInputIndex > 0) {
				offset = (int) result->data.size();
				result->data.push_back(NULL);
			}

			if (++vertInputIndex < (int) vertInputs.getCount())
				--i;
		}

		domListOfFloats &floatArray = source->getFloat_array()->getValue();
		domSource::domTechnique_common *techniqueCommon = source->getTechnique_common();
		if (!techniqueCommon)
			SLog(EError, "Data source does not have a <technique_common> tag!");
		domAccessor *accessor = techniqueCommon->getAccessor();
		if (!accessor)
			SLog(EError, "Data source does not have a <accessor> tag!");
		unsigned int nParams = (unsigned int) accessor->getParam_array().getCount(),
			         stride  = (unsigned int) accessor->getStride();
		size_t size = (size_t) accessor->getCount();
		SAssert(nParams <= 4);

		Vec4 *target = new Vec4[size];
		for (size_t j=0; j<size; ++j)
			for (unsigned int k=0; k<nParams; ++k)
				target[j][k] = (Float) floatArray[j*stride+k];

		result->data[offset] = target;

		if (semantic == "POSITION") {
			SAssert(accessor->getStride() == 3);
			SAssert(result->typeToOffset[EPosition] == -1);
			result->typeToOffset[EPosition] = offset;
			result->typeToCount[EPosition] = size;
			result->typeToOffsetInStream[EPosition] = offsetInStream;
			result->glPos = new GLdouble[3*size];
			for (size_t k=0; k<3*size; ++k)
				result->glPos[k] = floatArray[k];
		} else if (semantic == "NORMAL") {
			SAssert(accessor->getStride() == 3);
			SAssert(result->typeToOffset[ENormal] == -1);
			result->hasNormals = true;
			result->typeToOffset[ENormal] = offset;
			result->typeToOffsetInStream[ENormal] = offsetInStream;
			result->typeToCount[ENormal] = size;
		} else if (semantic == "TEXCOORD") {
			SAssert(accessor->getStride() == 2 || accessor->getStride() == 3);
			if (result->typeToOffset[EUV] == -1) {
				result->hasUVs = true;
				result->typeToOffset[EUV] = offset;
				result->typeToOffsetInStream[EUV] = offsetInStream;
				result->typeToCount[EUV] = size;
			} else {
				SLog(EWarn, "Found multiple sets of texture coordinates - ignoring!");
			}
		} else if (semantic == "COLOR") {
			SAssert(accessor->getStride() == 3 || accessor->getStride() == 4);
			SAssert(result->typeToOffset[EVertexColor] == -1);
			result->hasNormals = true;
			result->typeToOffset[EVertexColor] = offset;
			result->typeToOffsetInStream[EVertexColor] = offsetInStream;
			result->typeToCount[EVertexColor] = size;
		} else {
			SLog(EError, "Encountered an unknown source semantic: %s", semantic.c_str());
		}
	}
	SAssert(result->typeToOffset[EPosition] != -1);

	return result;
}

struct Vertex {
	Point p;
	Normal n;
	Point2 uv;
	Vector col;
};

/// For using vertices as keys in an associative structure
struct vertex_key_order : public
	std::binary_function<Vertex, Vertex, bool> {
	static int compare(const Vertex &v1, const Vertex &v2) {
		if (v1.p.x < v2.p.x) return -1;
		else if (v1.p.x > v2.p.x) return 1;
		if (v1.p.y < v2.p.y) return -1;
		else if (v1.p.y > v2.p.y) return 1;
		if (v1.p.z < v2.p.z) return -1;
		else if (v1.p.z > v2.p.z) return 1;
		if (v1.n.x < v2.n.x) return -1;
		else if (v1.n.x > v2.n.x) return 1;
		if (v1.n.y < v2.n.y) return -1;
		else if (v1.n.y > v2.n.y) return 1;
		if (v1.n.z < v2.n.z) return -1;
		else if (v1.n.z > v2.n.z) return 1;
		if (v1.uv.x < v2.uv.x) return -1;
		else if (v1.uv.x > v2.uv.x) return 1;
		if (v1.uv.y < v2.uv.y) return -1;
		else if (v1.uv.y > v2.uv.y) return 1;
		if (v1.col.x < v2.col.x) return -1;
		else if (v1.col.x > v2.col.x) return 1;
		if (v1.col.y < v2.col.y) return -1;
		else if (v1.col.y > v2.col.y) return 1;
		if (v1.col.z < v2.col.z) return -1;
		else if (v1.col.z > v2.col.z) return 1;
		return 0;
	}

	bool operator()(const Vertex &v1, const Vertex &v2) const {
		return compare(v1, v2) < 0;
	}
};

struct SimpleTriangle {
	Point p0, p1, p2;

	inline SimpleTriangle() { }
	inline SimpleTriangle(const Point &p0, const Point &p1, const Point &p2)
		: p0(p0), p1(p1), p2(p2) { }
};

struct triangle_key_order : public std::binary_function<SimpleTriangle, SimpleTriangle, bool> {
	static int compare(const Point &v1, const Point &v2) {
		if (v1.x < v2.x) return -1;
		else if (v1.x > v2.x) return 1;
		if (v1.y < v2.y) return -1;
		else if (v1.y > v2.y) return 1;
		if (v1.z < v2.z) return -1;
		else if (v1.z > v2.z) return 1;
		return 0;
	}

	bool operator()(const SimpleTriangle &t1, const SimpleTriangle &t2) const {
		int result;
		result = compare(t1.p0, t2.p0);
		if (result == -1) return true;
		if (result ==  1) return false;
		result = compare(t1.p1, t2.p1);
		if (result == -1) return true;
		if (result ==  1) return false;
		result = compare(t1.p2, t2.p2);
		if (result == -1) return true;
		if (result ==  1) return false;
		return false;
	}
};

class CustomErrorHandler : public daeErrorHandler {
public:
	void handleError(daeString msg) {
		SLog(EWarn, "Critical COLLADA error: %s", msg);
	}

	void handleWarning(daeString msg) {
		SLog(EWarn, "COLLADA warning: %s", msg);
	}
};

typedef std::map<SimpleTriangle, bool, triangle_key_order> TriangleMap;

static inline Float fromSRGBComponent(Float value) {
	if (value <= (Float) 0.04045)
		return value / (Float) 12.92;
	return std::pow((value + (Float) 0.055)
		/ (Float) (1.0 + 0.055), (Float) 2.4);
}

void writeGeometry(ColladaContext &ctx, const std::string &prefixName, std::string id,
		int geomIndex, std::string matID, Transform transform, VertexData *vData,
		TriangleMap &triMap, bool exportShapeGroup) {
	std::vector<Vertex> vertexBuffer;
	std::vector<Triangle> triangles;
	std::map<Vertex, int, vertex_key_order> vertexMap;
	size_t numMerged = 0, triangleIdx = 0, duplicates = 0;
	Triangle triangle;
	if (tess_data.size() == 0)
		return;

	for (size_t i=0; i<tess_cleanup.size(); ++i)
		delete[] tess_cleanup[i];
	tess_cleanup.clear();

	for (size_t i=0; i<tess_data.size(); i+=tess_nSources) {
		Vertex vertex;
		domUint posRef = tess_data[i+vData->typeToOffsetInStream[EPosition]];
		vertex.p = vData->data[vData->typeToOffset[EPosition]][posRef].toPoint();

		if (vData->typeToOffset[ENormal] != -1) {
			domUint normalRef = tess_data[i+vData->typeToOffsetInStream[ENormal]];
			vertex.n = vData->data[vData->typeToOffset[ENormal]][normalRef].toNormal();
		} else {
			vertex.n = Normal(0.0f);
		}

		if (vData->typeToOffset[EVertexColor] != -1) {
			domUint colorRef = tess_data[i+vData->typeToOffsetInStream[EVertexColor]];
			vertex.col = vData->data[vData->typeToOffset[EVertexColor]][colorRef].toVector();
		} else {
			vertex.col = Normal(0.0f);
		}

		if (vData->typeToOffset[EUV] != -1) {
			domUint uvRef = tess_data[i+vData->typeToOffsetInStream[EUV]];
			vertex.uv = vData->data[vData->typeToOffset[EUV]][uvRef].toPoint2();
			vertex.uv.y = 1-vertex.uv.y; // Invert the V coordinate
		} else {
			vertex.uv = Point2(0.0f);
		}

		int key = -1;
		if (vertexMap.find(vertex) != vertexMap.end()) {
			key = vertexMap[vertex];
			numMerged++;
		} else {
			key = (int) vertexBuffer.size();
			vertexMap[vertex] = (int) key;
			vertexBuffer.push_back(vertex);
		}

		triangle.idx[triangleIdx++] = key;
		if (triangleIdx == 3) {
			Point p0 = vertexBuffer[triangle.idx[0]].p,
				  p1 = vertexBuffer[triangle.idx[1]].p,
				  p2 = vertexBuffer[triangle.idx[2]].p;
			if (triMap.find(SimpleTriangle(p0, p1, p2)) != triMap.end() ||
				triMap.find(SimpleTriangle(p2, p0, p1)) != triMap.end() ||
				triMap.find(SimpleTriangle(p1, p2, p0)) != triMap.end() ||
				triMap.find(SimpleTriangle(p1, p0, p2)) != triMap.end() ||
				triMap.find(SimpleTriangle(p0, p2, p1)) != triMap.end() ||
				triMap.find(SimpleTriangle(p2, p1, p0)) != triMap.end()) {
				/* This triangle is a duplicate from another one which exists
				   in the same geometry group -- we may be dealing with SketchUp,
				   which sometimes exports every face TWICE! */
				duplicates++;
			} else {
				triMap[SimpleTriangle(p0, p1, p2)] = true;
				triangles.push_back(triangle);
			}
			triangleIdx = 0;
		}
	}

	if (duplicates > 0) {
		if (triangles.size() == 0) {
			SLog(EWarn, "\"%s/%s\": Only contains duplicates (%i triangles) of already-existing geometry. Ignoring.",
				prefixName.c_str(), id.c_str(), duplicates);
			ctx.os << "\t<!-- Ignored shape \"" << prefixName << "/"
				<< id << "\" (mat=\"" << matID << "\"), since it only contains duplicate geometry. -->" << endl << endl;
			return;
		} else {
			SLog(EWarn, "Geometry contains %i duplicate triangles!", duplicates);
		}
	}

	SAssert(triangleIdx == 0);
	SLog(EDebug, "\"%s/%s\": Converted " SIZE_T_FMT " triangles, " SIZE_T_FMT
		" vertices (merged " SIZE_T_FMT " vertices).", prefixName.c_str(), id.c_str(),
		triangles.size(), vertexBuffer.size(), numMerged);

	ref<TriMesh> mesh = new TriMesh(prefixName + "/" + id,
		triangles.size(), vertexBuffer.size(),
		vData->typeToOffset[ENormal] != -1,
		vData->typeToOffset[EUV] != -1,
		vData->typeToOffset[EVertexColor] != -1);

	std::copy(triangles.begin(), triangles.end(), mesh->getTriangles());

	Point    *target_positions = mesh->getVertexPositions();
	Normal   *target_normals   = mesh->getVertexNormals();
	Point2   *target_texcoords = mesh->getVertexTexcoords();
	Color3   *target_colors    = mesh->getVertexColors();

	for (size_t i=0; i<vertexBuffer.size(); ++i) {
		*target_positions++ = vertexBuffer[i].p;
		if (target_normals)
			*target_normals ++ = vertexBuffer[i].n;
		if (target_texcoords)
			*target_texcoords++ = vertexBuffer[i].uv;
		if (target_colors) {
			Float r = vertexBuffer[i].col.x;
			Float g = vertexBuffer[i].col.y;
			Float b = vertexBuffer[i].col.z;
			if (!ctx.cvt->m_srgb)
				*target_colors++ = Color3(r, g, b);
			else
				*target_colors++ = Color3(
					fromSRGBComponent(r),
					fromSRGBComponent(g),
					fromSRGBComponent(b)
				);
		}
	}

	id += formatString("_%i", geomIndex);
	std::string filename;

	if (!ctx.cvt->m_geometryFile) {
		filename = id + std::string(".serialized");
		ref<FileStream> stream = new FileStream(ctx.meshesDirectory / filename, FileStream::ETruncReadWrite);
		stream->setByteOrder(Stream::ELittleEndian);
		mesh->serialize(stream);
		stream->close();
		filename = "meshes/" + filename;
	} else {
		ctx.cvt->m_geometryDict.push_back((uint64_t) ctx.cvt->m_geometryFile->getPos());
		mesh->serialize(ctx.cvt->m_geometryFile);
		filename = ctx.cvt->m_geometryFileName.filename().string();
	}

	std::ostringstream matrix;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			matrix << transform.getMatrix().m[i][j] << " ";
	std::string matrixValues = matrix.str();

	if (!exportShapeGroup) {
		ctx.os << "\t<shape id=\"" << id << "\" type=\"serialized\">" << endl;
		ctx.os << "\t\t<string name=\"filename\" value=\"" << filename << "\"/>" << endl;
		if (ctx.cvt->m_geometryFile)
			ctx.os << "\t\t<integer name=\"shapeIndex\" value=\"" << (ctx.cvt->m_geometryDict.size() - 1) << "\"/>" << endl;
		if (!transform.isIdentity()) {
			ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
			ctx.os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
			ctx.os << "\t\t</transform>" << endl;
		}
		if (matID != "")
			ctx.os << "\t\t<ref name=\"bsdf\" id=\"" << matID << "\"/>" << endl;
		ctx.os << "\t</shape>" << endl << endl;
	} else {
		ctx.os << "\t\t<shape type=\"serialized\">" << endl;
		ctx.os << "\t\t\t<string name=\"filename\" value=\"" << filename << "\"/>" << endl;
		if (ctx.cvt->m_geometryFile)
			ctx.os << "\t\t\t<integer name=\"shapeIndex\" value=\"" << (ctx.cvt->m_geometryDict.size() - 1)<< "\"/>" << endl;
		if (matID != "")
			ctx.os << "\t\t\t<ref name=\"bsdf\" id=\"" << matID << "\"/>" << endl;
		ctx.os << "\t\t</shape>" << endl << endl;
	}
}

void exportAnimation(ColladaContext &ctx, const fs::path &path, const std::string &name) {
	AnimationMap::iterator start = ctx.animations.lower_bound(name);
	AnimationMap::iterator end = ctx.animations.upper_bound(name);

	ref<AnimatedTransform> trafo = new AnimatedTransform();

	for (AnimationMap::iterator it = start; it != end; ++it) {
		AbstractAnimationTrack *track = it->second;
		int type = track->getType();
		if (type == FloatTrack::ERotationX
			|| type == FloatTrack::ERotationY
			|| type == FloatTrack::ERotationZ)
			continue;
		trafo->addTrack(track);
	}
	SLog(EDebug, "Writing animation track \"%s\"", path.filename().string().c_str());
	ref<FileStream> fs = new FileStream(path, FileStream::ETruncReadWrite);
	trafo->serialize(fs);
}

void loadGeometry(ColladaContext &ctx, const std::string &instanceName,
		std::string prefixName, Transform transform, domGeometry &geom,
		StringMap &matLookupTable) {
	std::string identifier;
	if (geom.getId() != NULL) {
		identifier = geom.getId();
	} else {
		if (geom.getName() != NULL) {
			identifier = geom.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedGeom_%i", unnamedCtr++);
		}
	}
	TriangleMap triMap;

	SLog(EDebug, "Converting geometry \"%s\" (instantiated by %s)..", identifier.c_str(),
		prefixName == "" ? "/" : prefixName.c_str());
	domMesh *mesh = geom.getMesh().cast();
	if (!mesh)
		SLog(EError, "Invalid geometry type encountered (must be a <mesh>)!");

	bool exportShapeGroup = false;
	if (ctx.serializedGeometry.find(identifier) != ctx.serializedGeometry.end()) {
		if (ctx.animations.find(instanceName) != ctx.animations.end()) {
			ctx.os << "\t<shape id=\"" << prefixName << "/" << identifier << "\" type=\"animatedinstance\">" << endl;
			std::string filename = formatString("animation_%i.serialized", ctx.trackIndex++);
			ctx.os << "\t\t<string name=\"filename\" value=\"" << filename << "\"/>" << endl;
			exportAnimation(ctx, ctx.cvt->m_outputDirectory / filename, instanceName);
		} else {
			ctx.os << "\t<shape id=\"" << prefixName << "/" << identifier << "\" type=\"instance\">" << endl;
		}
		if (!transform.isIdentity()) {
			std::ostringstream matrix;
			for (int i=0; i<4; ++i)
				for (int j=0; j<4; ++j)
					matrix << transform.getMatrix().m[i][j] << " ";
			std::string matrixValues = matrix.str();
			ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
			ctx.os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
			ctx.os << "\t\t</transform>" << endl;
		}
		ctx.os << "\t\t<ref id=\"" << identifier << "\"/>" << endl;
		ctx.os << "\t</shape>" << endl << endl;
		return;
	} else if (ctx.refCountMap[identifier] > 1 || ctx.refCountMap[instanceName]) {
		ctx.os << "\t<shape id=\"" << identifier << "\" type=\"shapegroup\">" << endl;
		exportShapeGroup = true;
	}
	ctx.serializedGeometry.insert(identifier);

	const domInputLocal_Array &vertInputs = mesh->getVertices()->getInput_array();

	int geomIndex = 0;

	domTriangles_Array &trianglesArray = mesh->getTriangles_array();
	for (size_t i=0; i<trianglesArray.getCount(); ++i) {
		domTriangles *triangles = trianglesArray[i];
		domInputLocalOffset_Array &inputs = triangles->getInput_array();
		VertexData *data = fetchVertexData(transform, vertInputs, inputs);
		domListOfUInts &indices = triangles->getP()->getValue();
		tess_data.clear();
		tess_nSources = data->nSources;
		for (size_t j=0; j<indices.getCount(); ++j)
			tess_data.push_back(indices[j]);
		std::string matID;
		if (triangles->getMaterial() == NULL || matLookupTable.find(triangles->getMaterial()) == matLookupTable.end())
			SLog(EWarn, "Referenced material could not be found, substituting a diffuse BRDF.");
		else
			matID = matLookupTable[triangles->getMaterial()];
		writeGeometry(ctx, prefixName, identifier, geomIndex, matID, transform, data, triMap, exportShapeGroup);
		delete data;
		++geomIndex;
	}

	domPolygons_Array &polygonsArray = mesh->getPolygons_array();
	for (size_t i=0; i<polygonsArray.getCount(); ++i) {
		domPolygons *polygons = polygonsArray[i];
		domInputLocalOffset_Array &inputs = polygons->getInput_array();
		VertexData *data = fetchVertexData(transform, vertInputs, inputs);
		domP_Array &indexArray = polygons->getP_array();
		int posOffset = data->typeToOffset[EPosition];

		tess_data.clear();
		tess_nSources = data->nSources;
		for (size_t j=0; j<indexArray.getCount(); ++j) {
			domListOfUInts &indices = indexArray[j]->getValue();
			domUint *temp = new domUint[indices.getCount()];
			bool invalid = false;

			for (int m=0; m<ELast; ++m) {
				int offset = data->typeToOffsetInStream[m];
				if (offset == -1)
					continue;
				size_t count  = data->typeToCount[m];
				for (size_t l = 0; l < indices.getCount(); l+=data->nSources) {
					domUint idx = indices.get(l+offset);
					if (idx >= count) {
						SLog(EWarn, "Encountered the invalid polygon index %i "
							"(must be in [0, " SIZE_T_FMT "]) -- ignoring polygon!", (domInt) idx, count-1);
						invalid = true;
						break;
					}
				}
			}

			if (invalid)
				continue;

			for (size_t l = 0; l<indices.getCount(); ++l)
				temp[l] = indices.get(l);

			tess_vdata = data;
			gluTessBeginPolygon(tess, NULL);
			gluTessBeginContour(tess);

			for (size_t k=0; k<indices.getCount(); k+=data->nSources)
				gluTessVertex(tess, &data->glPos[3*temp[k+posOffset]], (GLvoid *) (temp+k));

			gluTessEndContour(tess);
			gluTessEndPolygon(tess);
			delete[] temp;
		}

		std::string matID;
		if (polygons->getMaterial() == NULL || matLookupTable.find(polygons->getMaterial()) == matLookupTable.end())
			SLog(EWarn, "Referenced material could not be found, substituting a diffuse BRDF.");
		else
			matID = matLookupTable[polygons->getMaterial()];

		writeGeometry(ctx, prefixName, identifier, geomIndex, matID, transform, data, triMap, exportShapeGroup);
		delete data;
		++geomIndex;
	}

	domPolylist_Array &polylistArray = mesh->getPolylist_array();
	for (size_t i=0; i<polylistArray.getCount(); ++i) {
		domPolylist *polylist = polylistArray[i];
		domInputLocalOffset_Array &inputs = polylist->getInput_array();
		VertexData *data = fetchVertexData(transform, vertInputs, inputs);
		domListOfUInts &vcount = polylist->getVcount()->getValue();
		domListOfUInts &indexArray = polylist->getP()->getValue();
		int posOffset = data->typeToOffset[EPosition], indexOffset = 0;
		tess_data.clear();
		tess_nSources = data->nSources;

		for (size_t j=0; j<vcount.getCount(); ++j) {
			size_t vertexCount = (size_t) vcount.get(j);
			domUint *temp = new domUint[vertexCount * data->nSources];

			bool invalid = false;

			for (int m=0; m<ELast; ++m) {
				int offset = data->typeToOffsetInStream[m];
				if (offset == -1)
					continue;
				size_t count  = data->typeToCount[m];
				for (size_t l = 0; l < vertexCount; ++l) {
					domUint idx = indexArray.get(l*data->nSources+offset);
					if (idx >= count) {
						SLog(EWarn, "Encountered the invalid polygon index %i "
							"(must be in [0, " SIZE_T_FMT "]) -- ignoring polygon!", (domInt) idx, count-1);
						invalid = true;
						break;
					}
				}
			}

			if (invalid)
				continue;

			for (size_t l = 0; l<vertexCount * data->nSources; ++l)
				temp[l] = indexArray.get(indexOffset++);

			tess_vdata = data;
			gluTessBeginPolygon(tess, NULL);
			gluTessBeginContour(tess);

			for (size_t k=0; k<vertexCount; k++)
				gluTessVertex(tess, &data->glPos[temp[k*data->nSources+posOffset]*3], (GLvoid *) (temp + k*data->nSources));

			gluTessEndContour(tess);
			gluTessEndPolygon(tess);
			delete[] temp;
		}

		std::string matID;
		if (polylist->getMaterial() == NULL || matLookupTable.find(polylist->getMaterial()) == matLookupTable.end())
			SLog(EWarn, "Referenced material \"%s\" could not be found, substituting a diffuse BRDF.", polylist->getMaterial());
		else
			matID = matLookupTable[polylist->getMaterial()];

		writeGeometry(ctx, prefixName, identifier, geomIndex, matID, transform, data, triMap, exportShapeGroup);
		delete data;
		++geomIndex;
	}
	if (exportShapeGroup) {
		ctx.os << "\t</shape>" << endl << endl;
		if (ctx.animations.find(instanceName) != ctx.animations.end()) {
			ctx.os << "\t<shape id=\"" << prefixName << "/" << identifier << "\" type=\"animatedinstance\">" << endl;
			std::string filename = formatString("animation_%i.serialized", ctx.trackIndex++);
			ctx.os << "\t\t<string name=\"filename\" value=\"" << filename << "\"/>" << endl;
			exportAnimation(ctx, ctx.cvt->m_outputDirectory / filename, instanceName);
		} else {
			ctx.os << "\t<shape id=\"" << prefixName << "/" << identifier << "\" type=\"instance\">" << endl;
			if (!transform.isIdentity()) {
				std::ostringstream matrix;
				for (int i=0; i<4; ++i)
					for (int j=0; j<4; ++j)
						matrix << transform.getMatrix().m[i][j] << " ";
				std::string matrixValues = matrix.str();
				ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
				ctx.os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
				ctx.os << "\t\t</transform>" << endl;
			}
		}
		ctx.os << "\t\t<ref id=\"" << identifier << "\"/>" << endl;
		ctx.os << "\t</shape>" << endl << endl;
	}
}

void loadMaterialParam(ColladaContext &ctx, const std::string &name,
		domCommon_color_or_texture_type *value, bool handleRefs) {
	if (!value)
		return;
	domCommon_color_or_texture_type_complexType::domColor* color =
		value->getColor().cast();
	domCommon_color_or_texture_type_complexType::domTexture* texture =
		value->getTexture().cast();
	if (color && !handleRefs) {
		domFloat4 &colValue = color->getValue();
		if (ctx.cvt->m_srgb)
			ctx.os << "\t\t<srgb name=\"" << name << "\" value=\"";
		else
			ctx.os << "\t\t<rgb name=\"" << name << "\" value=\"";
		ctx.os << colValue.get(0) << " " << colValue.get(1) << " "
		   << colValue.get(2) << "\"/>" << endl;
	} else if (texture && handleRefs) {
		if (ctx.idToTexture.find(texture->getTexture()) == ctx.idToTexture.end()) {
			SLog(EError, "Could not find referenced texture \"%s\"", texture->getTexture());
		} else {
			ctx.os << "\t\t<ref name=\"" << name << "\" id=\""
			   << ctx.idToTexture[texture->getTexture()] << "\"/>" << endl;
		}
	}
}

void loadMaterialParam(ColladaContext &ctx, const std::string &name,
		domCommon_float_or_param_type *value, bool handleRef) {
	if (!value)
		return;
	domCommon_float_or_param_type::domFloat *floatValue = value->getFloat();
	if (!handleRef && floatValue) {
		ctx.os << "\t\t<float name=\"" << name << "\" value=\""
		   << floatValue->getValue() << "\"/>" << endl;
	}
}

void loadMaterial(ColladaContext &ctx, domMaterial &mat) {
	std::string identifier;
	if (mat.getId() != NULL) {
		identifier = mat.getId();
	} else {
		if (mat.getName() != NULL) {
			identifier = mat.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedMat_%i", unnamedCtr++);
		}
	}

	daeURI &effRef = mat.getInstance_effect()->getUrl();
	effRef.resolveElement();
	domEffect *effect = daeSafeCast<domEffect>(effRef.getElement());

	if (!effect)
		SLog(EError, "Referenced effect not found!");

	domProfile_COMMON *commonProfile = daeSafeCast<domProfile_COMMON>
		(effect->getDescendant("profile_COMMON"));

	if (!commonProfile)
		SLog(EError, "Common effect profile not found!");

	/* The following supports a subset of the curious ColladaFX output
       produced by the Blender COLLADA exporter */
	daeTArray<daeSmartRef<domCommon_newparam_type> > &newParamArray = commonProfile->getNewparam_array();
	for (size_t i=0; i<newParamArray.getCount(); ++i) {
		domCommon_newparam_type_complexType *newParam = newParamArray[i];
		domFx_surface_common *surface = newParam->getSurface();
		domFx_sampler2D_common *sampler2D = newParam->getSampler2D();


		if (surface) {
			SAssert(surface->getType() == FX_SURFACE_TYPE_ENUM_2D);
			daeTArray<daeSmartRef<domFx_surface_init_from_common> > &initFromArray
				= surface->getFx_surface_init_common()->getInit_from_array();
			SAssert(initFromArray.getCount() == 1);
			std::string id = initFromArray[0]->getValue().getID();
			if (ctx.idToTexture.find(id) == ctx.idToTexture.end())
				SLog(EError, "Referenced bitmap '%s' not found!", id.c_str());
			ctx.idToTexture[newParam->getSid()] = ctx.idToTexture[id];
		}

		if (sampler2D) {
			std::string id = sampler2D->getSource()->getValue();
			if (ctx.idToTexture.find(id) == ctx.idToTexture.end())
				SLog(EError, "Referenced surface '%s' not found!", id.c_str());
			ctx.idToTexture[newParam->getSid()] = ctx.idToTexture[id];
		}
	}

	domProfile_COMMON::domTechnique *technique = commonProfile->getTechnique();

	if (!technique)
		SLog(EError, "The technique element is missing!");

	domProfile_COMMON::domTechnique::domPhong* phong = technique->getPhong();
	domProfile_COMMON::domTechnique::domBlinn* blinn = technique->getBlinn();
	domProfile_COMMON::domTechnique::domLambert* lambert = technique->getLambert();
	domProfile_COMMON::domTechnique::domConstant* constant = technique->getConstant();

	if (phong) {
		domCommon_color_or_texture_type* diffuse = phong->getDiffuse();
		domCommon_color_or_texture_type* specular = phong->getSpecular();
		domCommon_float_or_param_type* shininess = phong->getShininess();
		bool isDiffuse = false;

		if (specular && specular->getColor().cast()) {
			domFloat4 &colValue = specular->getColor()->getValue();
			if (colValue.get(0) == colValue.get(1) &&
				colValue.get(1) == colValue.get(2) &&
				colValue.get(2) == 0)
				isDiffuse = true;
		}
		if (isDiffuse) {
			ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"diffuse\">" << endl;
			loadMaterialParam(ctx, "reflectance", diffuse, false);
			loadMaterialParam(ctx, "reflectance", diffuse, true);
			ctx.os << "\t</bsdf>" << endl << endl;
		} else {
			ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"phong\">" << endl;
			loadMaterialParam(ctx,  "diffuseReflectance", diffuse, false);
			loadMaterialParam(ctx, "specularReflectance", specular, false);
			loadMaterialParam(ctx, "exponent", shininess, false);
			loadMaterialParam(ctx, "diffuseReflectance", diffuse, true);
			loadMaterialParam(ctx, "specularReflectance", specular, true);
			loadMaterialParam(ctx, "exponent", shininess, true);
			ctx.os << "\t</bsdf>" << endl << endl;
		}
	} else if (lambert) {
		domCommon_color_or_texture_type* diffuse = lambert->getDiffuse();
		ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"diffuse\">" << endl;
		loadMaterialParam(ctx, "reflectance", diffuse, false);
		loadMaterialParam(ctx, "reflectance", diffuse, true);
		ctx.os << "\t</bsdf>" << endl << endl;
	} else if (blinn) {
		domCommon_color_or_texture_type* diffuse = blinn->getDiffuse();
		domCommon_color_or_texture_type* specular = blinn->getSpecular();
		domCommon_float_or_param_type* shininess = blinn->getShininess();
		bool isDiffuse = false;

		if (specular && specular->getColor().cast()) {
			domFloat4 &colValue = specular->getColor()->getValue();
			if (colValue.get(0) == colValue.get(1) &&
				colValue.get(1) == colValue.get(2) &&
				colValue.get(2) == 0)
				isDiffuse = true;
		}
		if (isDiffuse) {
			ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"diffuse\">" << endl;
			loadMaterialParam(ctx, "reflectance", diffuse, false);
			loadMaterialParam(ctx, "reflectance", diffuse, true);
			ctx.os << "\t</bsdf>" << endl << endl;
		} else {
			ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"phong\">" << endl;
			loadMaterialParam(ctx, "diffuseReflectance", diffuse, false);
			loadMaterialParam(ctx, "specularReflectance", specular, false);
			loadMaterialParam(ctx, "exponent", shininess, false);
			loadMaterialParam(ctx, "diffuseReflectance", diffuse, true);
			loadMaterialParam(ctx, "specularReflectance", specular, true);
			loadMaterialParam(ctx, "exponent", shininess, true);
			ctx.os << "\t</bsdf>" << endl << endl;
		}
	} else if (constant) {
		domCommon_float_or_param_type* transparency = constant->getTransparency();
		domCommon_float_or_param_type::domFloat *transparencyValue =
				transparency ? transparency->getFloat() : NULL;
		if (transparencyValue && transparencyValue->getValue() > 0.5) {
			ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"dielectric\"/>" << endl << endl;
		} else {
			SLog(EWarn, "\"%s\": Encountered a \"constant\" COLLADA material, which is currently "
				"unsupported in Mitsuba -- replacing it using a Lambertian material.", identifier.c_str());
			ctx.os << "\t<bsdf id=\"" << identifier << "\" type=\"diffuse\"/>" << endl << endl;
		}
	} else {
		SLog(EError, "Material type not supported! (must be Lambertian/Phong/Blinn/Constant)");
	}
}

void loadLight(ColladaContext &ctx, Transform transform, domLight &light) {
	std::string identifier;
	if (light.getId() != NULL) {
		identifier = light.getId();
	} else {
		if (light.getName() != NULL) {
			identifier = light.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedLight_%i", unnamedCtr++);
		}
	}

	SLog(EDebug, "Converting light \"%s\" ..", identifier.c_str());
	char *end_ptr = NULL;

	// Lights in Mitsuba point along the positive Z axis (COLLADA: neg. Z)
	transform = transform * Transform::scale(Vector(1, 1, -1));

	Point pos = transform(Point(0, 0, 0));
	Point target = transform(Point(0, 0, 1));

	Float intensity = 1;
	const domTechnique_Array &techniques = light.getTechnique_array();
	for (size_t i=0; i<techniques.getCount(); ++i) {
		domTechnique *tech = techniques.get(i);

		daeElement *intensityElement = tech->getChild("intensity");
		if (intensityElement && intensityElement->hasCharData()) {
			std::string charData = intensityElement->getCharData();
			intensity = (Float) strtod(charData.c_str(), &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse the light intensity!");
		}
	}

	domLight::domTechnique_common::domPoint *point = light.getTechnique_common()->getPoint().cast();
	if (point) {
		bool notQuadratic = false;
		if (point->getConstant_attenuation() && point->getConstant_attenuation()->getValue() != 1)
			notQuadratic = true;
		if (point->getLinear_attenuation() && point->getLinear_attenuation()->getValue() != 0)
			notQuadratic = true;
		if (point->getQuadratic_attenuation() && point->getQuadratic_attenuation()->getValue() != 1)
			notQuadratic = true;
		if (notQuadratic)
			SLog(EWarn, "Point light \"%s\" is not a quadratic light! Treating it as one -- expect problems.", identifier.c_str());
		domFloat3 &color = point->getColor()->getValue();
		ctx.os << "\t<emitter id=\"" << identifier << "\" type=\"point\">" << endl;
		ctx.os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl << endl;
		ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
		ctx.os << "\t\t\t<translate x=\"" << pos.x << "\" y=\"" << pos.y << "\" z=\"" << pos.z << "\"/>" << endl;
		ctx.os << "\t\t</transform>" << endl;
		ctx.os << "\t</emitter>" << endl << endl;
	}

	domLight::domTechnique_common::domDirectional *directional = light.getTechnique_common()->getDirectional().cast();
	if (directional) {
		domFloat3 &color = directional->getColor()->getValue();
		ctx.os << "\t<emitter id=\"" << identifier << "\" type=\"directional\">" << endl;
		ctx.os << "\t\t<rgb name=\"irradiance\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl << endl;
		ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
		ctx.os << "\t\t\t<lookat origin=\"" << pos.x << ", " << pos.y << ", " << pos.z << "\" target=\"" << target.x << ", " << target.y << ", " << target.z << "\"/>" << endl;
		ctx.os << "\t\t</transform>" << endl << endl;
		ctx.os << "\t</emitter>" << endl << endl;
	}

	domLight::domTechnique_common::domSpot *spot = light.getTechnique_common()->getSpot().cast();
	if (spot) {
		bool notQuadratic = false;
		if (spot->getConstant_attenuation() && spot->getConstant_attenuation()->getValue() != 1)
			notQuadratic = true;
		if (spot->getLinear_attenuation() && spot->getLinear_attenuation()->getValue() != 0)
			notQuadratic = true;
		if (spot->getQuadratic_attenuation() && spot->getQuadratic_attenuation()->getValue() != 1)
			notQuadratic = true;
		if (notQuadratic)
			SLog(EWarn, "Spot light \"%s\" is not a quadratic light! Treating it as one -- expect problems.", identifier.c_str());
		domFloat3 &color = spot->getColor()->getValue();
		Float falloffAngle = 180.0f;
		if (spot->getFalloff_angle())
			falloffAngle = (Float) spot->getFalloff_angle()->getValue();
		ctx.os << "\t<emitter id=\"" << identifier << "\" type=\"spot\">" << endl;
		ctx.os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl;
		ctx.os << "\t\t<float name=\"cutoffAngle\" value=\"" << falloffAngle/2 << "\"/>" << endl << endl;
		ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
		ctx.os << "\t\t\t<lookat origin=\"" << pos.x << ", " << pos.y << ", " << pos.z << "\" target=\"" << target.x << ", " << target.y << ", " << target.z << "\"/>" << endl;
		ctx.os << "\t\t</transform>" << endl;
		ctx.os << "\t</emitter>" << endl << endl;
	}
	domLight::domTechnique_common::domAmbient *ambient = light.getTechnique_common()->getAmbient().cast();
	if (ambient) {
		domFloat3 &color = ambient->getColor()->getValue();
		ctx.os << "\t<emitter id=\"" << identifier << "\" type=\"constant\">" << endl;
		ctx.os << "\t\t<rgb name=\"radiance\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl;
		ctx.os << "\t</emitter>" << endl << endl;
	}
	if (!point && !spot && !ambient && !directional)
		SLog(EWarn, "Encountered an unknown light type!");
}

void loadImage(ColladaContext &ctx, domImage &image) {
	std::string identifier;
	if (image.getId() != NULL) {
		identifier = image.getId();
	} else {
		if (image.getName() != NULL) {
			identifier = image.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedTexture_%i", unnamedCtr++);
		}
	}

	SLog(EDebug, "Converting texture \"%s\" ..", identifier.c_str());

	std::string uri = image.getInit_from()->getValue().str();
	std::string filename = cdom::uriToFilePath(uri);
	if (filename.empty()) /* When uriToFilePath fails, try to use the path as is */
		filename = uri;

	/* Prevent Linux/OSX fs::path handling issues for DAE files created on Windows */
	for (size_t i=0; i<filename.length(); ++i) {
		if (filename[i] == '\\')
			filename[i] = '/';
	}
	if (ctx.fileToId.find(filename) != ctx.fileToId.end()) {
		ctx.idToTexture[identifier] = ctx.fileToId[filename];
		return;
	}

	ctx.idToTexture[identifier] = identifier;
	ctx.fileToId[filename] = identifier;

	fs::path path = fs::path(filename);
	fs::path targetPath = ctx.texturesDirectory / path.filename();

	std::string extension = boost::to_lower_copy(path.extension().string());
	if (extension == ".rgb")
		SLog(EWarn, "Maya RGB images must be converted to PNG, EXR or JPEG! The 'imgcvt' "
			"utility found in the Maya binary directory can be used to do this.");

	if (!fs::exists(targetPath)) {
		if (!fs::exists(path)) {
			ref<FileResolver> fRes = Thread::getThread()->getFileResolver();
			path = fRes->resolve(path.filename());

			if (!fs::exists(path)) {
				SLog(EWarn, "Found neither \"%s\" nor \"%s\"!", filename.c_str(), path.string().c_str());
				path = ctx.cvt->locateResource(path.filename());
				targetPath = targetPath.parent_path() / path.filename();
				if (path.empty())
					SLog(EError, "Unable to locate a resource -- aborting conversion.");
				else
					fRes->appendPath(path.parent_path());
			}
		}
		if (fs::absolute(path) != fs::absolute(targetPath)) {
			ref<FileStream> input = new FileStream(path, FileStream::EReadOnly);
			ref<FileStream> output = new FileStream(targetPath, FileStream::ETruncReadWrite);
			input->copyTo(output);
			input->close();
			output->close();
		}
	}

	ctx.os << "\t<texture id=\"" << identifier << "\" type=\"bitmap\">" << endl;
	ctx.os << "\t\t<string name=\"filename\" value=\"textures/" << targetPath.filename().string() << "\"/>" << endl;
	ctx.os << "\t</texture>" << endl << endl;
}

void loadCamera(ColladaContext &ctx, Transform transform, domCamera &camera) {
	std::string identifier;
	if (camera.getId() != NULL) {
		identifier = camera.getId();
	} else {
		if (camera.getName() != NULL) {
			identifier = camera.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedCamera_%i", unnamedCtr++);
		}
	}

	SLog(EDebug, "Converting camera \"%s\" ..", identifier.c_str());
	Float aspect = 1.0f;
	int xres=768;

	// Cameras in Mitsuba point along the positive Z axis (COLLADA: neg. Z)
	transform = transform * Transform::scale(Vector(1, 1, -1));

	std::ostringstream matrix;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			matrix << transform.getMatrix().m[i][j] << " ";
	std::string matrixValues = matrix.str();

	domCamera::domOptics::domTechnique_common::domOrthographic* ortho = camera.getOptics()->
		getTechnique_common()->getOrthographic().cast();
	if (ortho) {
		if (ortho->getAspect_ratio().cast() != 0)
			aspect = (Float) ortho->getAspect_ratio()->getValue();
		if (ctx.cvt->m_xres != -1) {
			xres = ctx.cvt->m_xres;
			aspect = (Float) ctx.cvt->m_xres / (Float) ctx.cvt->m_yres;
		}
		ctx.os << "\t<sensor id=\"" << identifier << "\" type=\"orthographic\">" << endl;
	}

	domCamera::domOptics::domTechnique_common::domPerspective* persp = camera.getOptics()->
		getTechnique_common()->getPerspective().cast();
	if (persp) {
		if (ctx.cvt->m_xres != -1) {
			xres = ctx.cvt->m_xres;
			aspect = (Float) ctx.cvt->m_xres / (Float) ctx.cvt->m_yres;
		} else {
			if (persp->getAspect_ratio().cast() != 0)
				aspect = (Float) persp->getAspect_ratio()->getValue();
		}
		ctx.os << "\t<sensor id=\"" << identifier << "\" type=\"perspective\">" << endl;
		if (persp->getXfov().cast()) {
			Float xFov = (Float) persp->getXfov()->getValue();
			Float yFov = radToDeg(2 * std::atan(std::tan(degToRad(xFov)/2) / aspect));
			if (aspect <= 1.0f)
				ctx.os << "\t\t<float name=\"fov\" value=\"" << xFov << "\"/>" << endl;
			else
				ctx.os << "\t\t<float name=\"fov\" value=\"" << yFov << "\"/>" << endl;
		} else if (persp->getYfov().cast()) {
			Float yFov = (Float) persp->getYfov()->getValue();
			Float xFov = radToDeg(2 * std::atan(std::tan(degToRad(yFov)/2) * aspect));
			if (aspect > 1.0f)
				ctx.os << "\t\t<float name=\"fov\" value=\"" << yFov << "\"/>" << endl;
			else
				ctx.os << "\t\t<float name=\"fov\" value=\"" << xFov << "\"/>" << endl;
		}
		ctx.os << "\t\t<float name=\"nearClip\" value=\"" << persp->getZnear()->getValue() << "\"/>" << endl;
		ctx.os << "\t\t<float name=\"farClip\" value=\"" << persp->getZfar()->getValue() << "\"/>" << endl;
		ctx.os << "\t\t<string name=\"fovAxis\" value=\"" << (ctx.cvt->m_mapSmallerSide ? "smaller" : "larger") << "\"/>" << endl;
	}

	ctx.os << endl;
	ctx.os << "\t\t<transform name=\"toWorld\">" << endl;
	ctx.os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
	ctx.os << "\t\t</transform>" << endl << endl;
	ctx.os << "\t\t<sampler id=\"" << identifier << "_sampler\" type=\"ldsampler\">" << endl;
	ctx.os << "\t\t\t<integer name=\"sampleCount\" value=\"4\"/>" << endl;
	ctx.os << "\t\t</sampler>" << endl << endl;
	ctx.os << "\t\t<film id=\"" << identifier << "_film\" type=\"" << ctx.cvt->m_filmType << "\">" << endl;
	ctx.os << "\t\t\t<integer name=\"width\" value=\"" << xres << "\"/>" << endl;
	ctx.os << "\t\t\t<integer name=\"height\" value=\"" << (int) (xres/aspect) << "\"/>" << endl;
	ctx.os << "\t\t\t<rfilter type=\"gaussian\"/>" << endl;
	ctx.os << "\t\t</film>" << endl;
	ctx.os << "\t</sensor>" << endl << endl;
}

void loadNode(ColladaContext &ctx, Transform transform, domNode &node, std::string prefixName) {
	std::string identifier;
	if (node.getId() != NULL) {
		identifier = node.getId();
	} else {
		if (node.getName() != NULL) {
			identifier = node.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedNode_%i", unnamedCtr);
		}
	}
	prefixName = prefixName + std::string("/") + identifier;
	SLog(EDebug, "Converting node \"%s\" ..", identifier.c_str());

	daeTArray<daeSmartRef<daeElement> > children = node.getChildren();
	/* Parse transformations */
	for (size_t i=0; i<children.getCount(); ++i) {
		daeElement *element = children.get(i);
		if (element->typeID() == domRotate::ID()) {
			/* Skip rotations labeled as "post-rotationY". Maya and 3ds max export these with some
			   cameras, which introduces an incorrect 90 degree rotation unless ignored */
			if (element->hasAttribute("sid") && element->getAttribute("sid") == "post-rotationY")
				continue;
			daeTArray<double> value = daeSafeCast<domRotate>(element)->getValue();
			Vector axis((Float) value.get(0), (Float) value.get(1), (Float) value.get(2));
			Float angle = (Float) value.get(3);
			if (angle != 0) {
				if (axis.isZero()) {
					SLog(EWarn, "Encountered a rotation around a zero vector -- ignoring!");
				} else {
					transform = transform *
						Transform::rotate(axis, (Float) value.get(3));
				}
			}
		} else if (element->typeID() == domTranslate::ID()) {
			daeTArray<double> value = daeSafeCast<domTranslate>(element)->getValue();
			transform = transform *
				Transform::translate(Vector((Float) value.get(0), (Float) value.get(1), (Float) value.get(2)));
		} else if (element->typeID() == domScale::ID()) {
			daeTArray<double> value = daeSafeCast<domScale>(element)->getValue();
			transform = transform *
				Transform::scale(Vector((Float) value.get(0), (Float) value.get(1), (Float) value.get(2)));
		} else if (element->typeID() == domLookat::ID()) {
			daeTArray<double> value = daeSafeCast<domLookat>(element)->getValue();
			transform = transform *
				Transform::lookAt(
					Point((Float) value.get(0), (Float) value.get(1), (Float) value.get(2)),
					Point((Float) value.get(3), (Float) value.get(4), (Float) value.get(5)),
					Vector((Float) value.get(6), (Float) value.get(7), (Float) value.get(8))
			);
		} else if (element->typeID() == domMatrix::ID()) {
			daeTArray<double> value = daeSafeCast<domMatrix>(element)->getValue();
			Matrix4x4 matrix(
				(Float) value.get(0), (Float) value.get(1), (Float) value.get(2), (Float) value.get(3),
				(Float) value.get(4), (Float) value.get(5), (Float) value.get(6), (Float) value.get(7),
				(Float) value.get(8), (Float) value.get(9), (Float) value.get(10), (Float) value.get(11),
				(Float) value.get(12), (Float) value.get(13), (Float) value.get(14), (Float) value.get(15)
			);
			transform = transform * Transform(matrix);
		}
	}


	/* Iterate over all geometry references */
	domInstance_geometry_Array &instanceGeometries = node.getInstance_geometry_array();
	for (size_t i=0; i<instanceGeometries.getCount(); ++i) {
		domInstance_geometry *inst = instanceGeometries[i];
		domGeometry *geom = daeSafeCast<domGeometry>(inst->getUrl().getElement());
		domBind_material *bmat = inst->getBind_material();
		StringMap matLookupTable;
		if (bmat) {
			domBind_material::domTechnique_common *technique = bmat->getTechnique_common();
			if (!technique)
				SLog(EError, "bind_material does not contain a <technique_common> element!");
			domInstance_material_Array &instMaterials = technique->getInstance_material_array();

			for (size_t i=0; i<instMaterials.getCount(); ++i) {
				domInstance_material *instMat = instMaterials[i];
				daeURI &matRef = instMat->getTarget();
				matRef.resolveElement();
				domMaterial *material = daeSafeCast<domMaterial>(matRef.getElement());
				matLookupTable[instMat->getSymbol()] = material->getId();
			}
		} else {
			SLog(EWarn, "instance_geometry does not contain a <bind_material> element!");
		}

		if (!geom)
			SLog(EError, "Could not find a referenced geometry object!");
		loadGeometry(ctx, identifier, prefixName, transform, *geom, matLookupTable);
	}

	/* Iterate over all light references */
	domInstance_light_Array &instanceLights = node.getInstance_light_array();
	for (size_t i=0; i<instanceLights.getCount(); ++i) {
		domInstance_light *inst = instanceLights[i];
		domLight *light = daeSafeCast<domLight>(inst->getUrl().getElement());
		if (!light)
			SLog(EError, "Could not find a referenced light!");
		loadLight(ctx, transform, *light);
	}

	/* Iterate over all camera references */
	domInstance_camera_Array &instanceCameras = node.getInstance_camera_array();
	for (size_t i=0; i<instanceCameras.getCount(); ++i) {
		domInstance_camera *inst = instanceCameras[i];
		domCamera *camera = daeSafeCast<domCamera>(inst->getUrl().getElement());
		if (camera == NULL)
			SLog(EError, "Could not find a referenced camera!");
		loadCamera(ctx, transform, *camera);
	}

	/* Recursively iterate through sub-nodes */
	domNode_Array &nodes = node.getNode_array();
	for (size_t i=0; i<nodes.getCount(); ++i)
		loadNode(ctx, transform, *nodes[i], prefixName);

	/* Recursively iterate through <instance_node> elements */
	domInstance_node_Array &instanceNodes = node.getInstance_node_array();
	for (size_t i=0; i<instanceNodes.getCount(); ++i) {
		domNode *node = daeSafeCast<domNode>(instanceNodes[i]->getUrl().getElement());
		if (!node)
			SLog(EError, "Could not find a referenced node!");
		loadNode(ctx, transform, *node, prefixName + std::string("/") + identifier);
	}
}

void computeRefCounts(ColladaContext &ctx, domNode &node) {
	domInstance_geometry_Array &instanceGeometries = node.getInstance_geometry_array();
	for (size_t i=0; i<instanceGeometries.getCount(); ++i) {
		domInstance_geometry *inst = instanceGeometries[i];
		domGeometry *geom = daeSafeCast<domGeometry>(inst->getUrl().getElement());
		if (geom->getId() != NULL) {
			if (ctx.refCountMap.find(geom->getId()) == ctx.refCountMap.end())
				ctx.refCountMap[geom->getId()] = 1;
			else
				ctx.refCountMap[geom->getId()]++;
		}
	}

	/* Recursively iterate through sub-nodes */
	domNode_Array &nodes = node.getNode_array();
	for (size_t i=0; i<nodes.getCount(); ++i)
		computeRefCounts(ctx, *nodes[i]);

	/* Recursively iterate through <instance_node> elements */
	domInstance_node_Array &instanceNodes = node.getInstance_node_array();
	for (size_t i=0; i<instanceNodes.getCount(); ++i) {
		domNode *node = daeSafeCast<domNode>(instanceNodes[i]->getUrl().getElement());
		computeRefCounts(ctx, *node);
	}
}

void loadAnimation(ColladaContext &ctx, domAnimation &anim) {
	std::string identifier;
	if (anim.getId() != NULL) {
		identifier = anim.getId();
	} else {
		if (anim.getName() != NULL) {
			identifier = anim.getName();
		} else {
			static int unnamedCtr = 0;
			identifier = formatString("unnamedAnimation_%i", unnamedCtr);
		}
	}
	SLog(EDebug, "Loading animation \"%s\" ..", identifier.c_str());

	domChannel_Array &channels = anim.getChannel_array();
	for (size_t i=0; i<channels.getCount(); ++i) {
		domChannel *channel = channels[i];
		std::vector<std::string> target = tokenize(channel->getTarget(), "./");
		SAssert(target.size() >= 2);
		if (ctx.refCountMap.find(target[0]) == ctx.refCountMap.end())
			ctx.refCountMap[target[0]] = 1;
		else
			ctx.refCountMap[target[0]]++;

		daeURI &sourceRef = channel->getSource();
		sourceRef.resolveElement();
		domSampler *sampler = daeSafeCast<domSampler>(sourceRef.getElement());
		if (!sampler)
			SLog(EError, "Referenced animation sampler not found!");
		const domInputLocal_Array &inputs = sampler->getInput_array();
		AbstractAnimationTrack *track = NULL;
		AbstractAnimationTrack::EType trackType = AbstractAnimationTrack::EInvalid;
		boost::to_lower(target[1]);
		if (target.size() > 2)
			boost::to_lower(target[2]);
		if (target[1] == "location" || target[1] == "translate") {
			if (target.size() == 2) {
				trackType = VectorTrack::ETranslationXYZ;
			} else if (target[2] == "x") {
				trackType = FloatTrack::ETranslationX;
			} else if (target[2] == "y") {
				trackType = FloatTrack::ETranslationY;
			} else if (target[2] == "z") {
				trackType = FloatTrack::ETranslationZ;
			}
		} else if (target[1] == "scale") {
			if (target.size() == 2) {
				trackType = VectorTrack::EScaleXYZ;
			} else if (target[2] == "x") {
				trackType = FloatTrack::EScaleX;
			} else if (target[2] == "y") {
				trackType = FloatTrack::EScaleY;
			} else if (target[2] == "z") {
				trackType = FloatTrack::EScaleZ;
			}
		} else if ((target[1] == "rotationx" || target[1] == "rotatex") && target.size() == 3 && target[2] == "angle") {
				trackType = FloatTrack::ERotationX;
		} else if ((target[1] == "rotationy" || target[1] == "rotatey") && target.size() == 3 && target[2] == "angle") {
				trackType = FloatTrack::ERotationY;
		} else if ((target[1] == "rotationz" || target[1] == "rotatez") && target.size() == 3 && target[2] == "angle") {
				trackType = FloatTrack::ERotationZ;
		}

		if (trackType == FloatTrack::EInvalid) {
			SLog(EWarn, "Skipping unsupported animation track of type %s.%s",
				target[1].c_str(), target.size() > 2 ? target[2].c_str() : "");
			continue;
		}

		for (size_t j=0; j<inputs.getCount(); ++j) {
			sourceRef = inputs[j]->getSource();
			sourceRef.resolveElement();
			domSource *source = daeSafeCast<domSource>(sourceRef.getElement());
			if (!source)
				SLog(EError, "Referenced animation source not found!");
			std::string semantic = inputs[j]->getSemantic();
			domSource::domTechnique_common *techniqueCommon = source->getTechnique_common();
			if (!techniqueCommon)
				SLog(EError, "Data source does not have a <technique_common> tag!");
			domAccessor *accessor = techniqueCommon->getAccessor();
			if (!accessor)
				SLog(EError, "Data source does not have a <accessor> tag!");
			unsigned int stride  = (unsigned int) accessor->getStride();
			size_t size = (size_t) accessor->getCount();

			if (!track) {
				if (trackType == VectorTrack::EScaleXYZ || trackType == VectorTrack::ETranslationXYZ)
					track = new VectorTrack(trackType, size);
				else
					track = new FloatTrack(trackType, size);
				track->incRef();
			} else {
				SAssert(track->getSize() == size);
			}

			if (semantic == "INPUT") {
				domListOfFloats &floatArray = source->getFloat_array()->getValue();
				SAssert(stride == 1);
				for (size_t i=0; i<size; ++i)
					track->setTime(i, (Float) floatArray[i]);
			} else if (semantic == "OUTPUT") {
				domListOfFloats &floatArray = source->getFloat_array()->getValue();
				if (trackType == VectorTrack::ETranslationXYZ || trackType == VectorTrack::EScaleXYZ) {
					SAssert(stride == 3);
					for (size_t i=0; i<size; ++i)
						((VectorTrack *) track)->setValue(i,
							Vector((Float) floatArray[i*3+0], (Float) floatArray[i*3+1],
								(Float) floatArray[i*3+2]));
				} else {
					SAssert(stride == 1);
					for (size_t i=0; i<size; ++i)
						((FloatTrack *) track)->setValue(i, (Float) floatArray[i]);
				}
			} else if (semantic == "INTERPOLATION") {
				/// Ignored for now
			} else {
				SLog(EWarn, "Encountered an unsupported semantic: \"%s\"", semantic.c_str());
			}
		}

		if (track)
			ctx.animations.insert(AnimationMap::value_type(target[0], track));
	}
}

void mergeRotations(ColladaContext &ctx) {
	for (AnimationMap::iterator it = ctx.animations.begin();
		it != ctx.animations.end();) {
		std::string key = it->first;
		AnimationMap::iterator start = ctx.animations.lower_bound(key);
		AnimationMap::iterator end = ctx.animations.upper_bound(key);
		FloatTrack *tracks[] = { NULL, NULL, NULL };

		for (AnimationMap::iterator it2 = start; it2 != end; ++it2) {
			if (it2->second->getType() == FloatTrack::ERotationX)
				tracks[0] = (FloatTrack *) it2->second;
			else if (it2->second->getType() == FloatTrack::ERotationY)
				tracks[1] = (FloatTrack *) it2->second;
			else if (it2->second->getType() == FloatTrack::ERotationZ)
				tracks[2] = (FloatTrack *) it2->second;
		}

		if (!tracks[0] && !tracks[1] && !tracks[2]) {
			it = end;
			continue;
		}

		SLog(EDebug, "Converting rotation track of \"%s\" to quaternions ..",
				key.c_str());

		std::set<Float> times;
		for (size_t i=0; i<3; ++i)
			for (size_t j=0; j<(tracks[i] ? tracks[i]->getSize() : (size_t) 0); ++j)
				times.insert(tracks[i]->getTime(j));

		QuatTrack *newTrack = new QuatTrack(QuatTrack::ERotationQuat, times.size());
		size_t idx = 0;

		for (std::set<Float>::iterator it2 = times.begin();
			it2 != times.end(); ++it2) {
			Float time = *it2, rot[3];
			for (int i=0; i<3; ++i)
				rot[i] = tracks[i] ? (tracks[i]->eval(time) * (M_PI/180)) : (Float) 0;

			newTrack->setTime(idx, time);
			newTrack->setValue(idx, Quaternion::fromEulerAngles(
				Quaternion::EEulerXYZ, rot[0], rot[1], rot[2]));
			idx++;
		}

		newTrack->incRef();
		ctx.animations.insert(AnimationMap::value_type(key, newTrack));

		it = ctx.animations.upper_bound(key);
	}
}


GLvoid __stdcall tessBegin(GLenum type) {
	SAssert(type == GL_TRIANGLES);
}

GLvoid __stdcall tessVertex(void *data) {
	const domUint *raw = (domUint *) data;
	for (size_t i=0; i<tess_nSources; ++i)
		tess_data.push_back(raw[i]);
}

GLvoid __stdcall tessCombine(GLdouble coords[3], void *vertex_data[4],
	GLfloat weight[4], void **outData) {
	domUint *result = new domUint[tess_nSources], size = 0;
	tess_cleanup.push_back(result);

	for (size_t i=0; i<tess_nSources; ++i) {
		int offset = 0;
		bool found = false;
		for (int j=0; j<ELast; ++j) {
			if (tess_vdata->typeToOffsetInStream[j] == (int) i) {
				offset = tess_vdata->typeToOffset[j];
				size = tess_vdata->typeToCount[j];
				found = true;
				break;
			}
		}
		SAssert(found);

		// this will be very slow -- let's hope that it happens rarely
		Vec4 *oldVec = tess_vdata->data[offset];
		Vec4 *newVec = new Vec4[(size_t) size + 1];
		memcpy(newVec, oldVec, (size_t) size * sizeof(Vec4));

		newVec[size] = Vec4(0.0f);
		for (int j=0; j<4; ++j) {
			void *ptr = vertex_data[j];
			if (!ptr)
				continue;
			domUint idx = ((domUint *) ptr)[i];
			newVec[size] += newVec[idx] * weight[j];
		}
		tess_vdata->data[offset] = newVec;
		result[i] = size;

		delete oldVec;
	}
	*outData = result;
}

GLvoid __stdcall tessEnd() { }
GLvoid __stdcall tessEdgeFlag(GLboolean) { }
GLvoid __stdcall tessError(GLenum error) {
	SLog(EError, "The GLU tesselator generated an error: %s!", gluErrorString(error));
}

void GeometryConverter::convertCollada(const fs::path &inputFile,
	std::ostream &os,
	const fs::path &texturesDirectory,
	const fs::path &meshesDirectory) {
	CustomErrorHandler errorHandler;
	daeErrorHandler::setErrorHandler(&errorHandler);
	SLog(EInfo, "Loading \"%s\" ..", inputFile.filename().string().c_str());
	#if defined(COLLADA_DOM_2_4)
		DAE *dae = new DAE(NULL, NULL, "1.4.1");
		if (dae->load(inputFile.string().c_str()) != DAE_OK)
			SLog(EError, "Could not load \"%s\"!", inputFile.string().c_str());
		ColladaDOM141::domCOLLADA *document = dae->getDom141(inputFile.string().c_str());
	#else
		DAE *dae = new DAE();
		if (dae->load(inputFile.string().c_str()) != DAE_OK)
			SLog(EError, "Could not load \"%s\"!", inputFile.string().c_str());
		domCOLLADA *document = dae->getDom(inputFile.string().c_str());
	#endif
	domVisual_scene *visualScene = daeSafeCast<domVisual_scene>
		(document->getDescendant("visual_scene"));
	if (!visualScene)
		SLog(EError, "Could not find a visual_scene!");

	/* Configure the GLU tesselator */
	tess = gluNewTess();
	if (!tess)
		SLog(EError, "Could not allocate a GLU tesselator!");

	gluTessCallback(tess, GLU_TESS_VERTEX, reinterpret_cast<GLvoid (__stdcall *)()>(&tessVertex));
	gluTessCallback(tess, GLU_TESS_BEGIN, reinterpret_cast<GLvoid (__stdcall *)()>(&tessBegin));
	gluTessCallback(tess, GLU_TESS_END, reinterpret_cast<GLvoid (__stdcall *)()>(&tessEnd));
	gluTessCallback(tess, GLU_TESS_ERROR, reinterpret_cast<GLvoid (__stdcall *)()>(&tessError));
	gluTessCallback(tess, GLU_TESS_COMBINE, reinterpret_cast<GLvoid (__stdcall *)()>(&tessCombine));
	gluTessCallback(tess, GLU_TESS_EDGE_FLAG, reinterpret_cast<GLvoid (__stdcall *)()>(&tessEdgeFlag));

	domNode_Array &nodes = visualScene->getNode_array();
	os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl << endl;
	os << "<!--" << endl << endl;
	os << "\tAutomatically converted from COLLADA" << endl << endl;
	os << "-->" << endl << endl;
	os << "<scene version=\"" << MTS_VERSION << "\">" << endl;
	os << "\t<integrator id=\"integrator\" type=\"direct\"/>" << endl << endl;

	ColladaContext ctx(os);
	ctx.meshesDirectory = meshesDirectory;
	ctx.texturesDirectory = texturesDirectory;
	ctx.cvt = this;
	ctx.trackIndex = 0;

	ref<Timer> timer = new Timer();

	if (m_importMaterials) {
		SLog(EInfo, "Importing materials ..");
		domLibrary_images_Array &libraryImages = document->getLibrary_images_array();
		for (size_t i=0; i<libraryImages.getCount(); ++i) {
			domImage_Array &images = libraryImages[i]->getImage_array();
			for (size_t j=0; j<images.getCount(); ++j)
				loadImage(ctx, *images.get(j));
		}

		domLibrary_materials_Array &libraryMaterials = document->getLibrary_materials_array();
		for (size_t i=0; i<libraryMaterials.getCount(); ++i) {
			domMaterial_Array &materials = libraryMaterials[i]->getMaterial_array();
			for (size_t j=0; j<materials.getCount(); ++j)
				loadMaterial(ctx, *materials.get(j));
		}
	}

	if (m_importAnimations) {
		SLog(EInfo, "Importing animations ..");
		domLibrary_animations_Array &libraryAnimations = document->getLibrary_animations_array();
		for (size_t i=0; i<libraryAnimations.getCount(); ++i) {
			domAnimation_Array &animations = libraryAnimations[i]->getAnimation_array();
			for (size_t j=0; j<animations.getCount(); ++j)
				loadAnimation(ctx, *animations[j]);
		}
		mergeRotations(ctx);
	}

	SLog(EInfo, "Importing scene ..");
	for (size_t i=0; i<nodes.getCount(); ++i)
		computeRefCounts(ctx, *nodes[i]);

	for (size_t i=0; i<nodes.getCount(); ++i)
		loadNode(ctx, Transform(), *nodes[i], "");

	for (AnimationMap::iterator it = ctx.animations.begin();
		it != ctx.animations.end(); ++it)
		it->second->decRef();

	SLog(EInfo, "Done, took %s", timeString(timer->getMilliseconds()/1000.0f).c_str());

	os << "</scene>" << endl;

	gluDeleteTess(tess);
	daeErrorHandler::setErrorHandler(NULL);
	delete dae;
}

