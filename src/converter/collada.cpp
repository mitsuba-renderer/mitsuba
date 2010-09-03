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

#define BOOST_FILESYSTEM_NO_LIB 
#define BOOST_SYSTEM_NO_LIB 

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/fresolver.h>
#include <dae.h>
#include <dom/domCOLLADA.h>
#include <dom/domProfile_COMMON.h>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#if defined(__OSX__)
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include "converter.h"

typedef std::map<std::string, std::string> StringMap;

enum ESourceType {
	EPosition = 0,
	ENormal = 1,
	EUV = 2,
	EVertexColors = 3,
	ELast
};
		
struct Vec4 {
	Float x, y, z, w;

	inline Vec4(Float x=0, Float y=0, Float z=0, Float w=0) 
		: x(x), y(y), z(z), w(w) {
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
size_t tess_nSources;

VertexData *fetchVertexData(Transform transform, std::ostream &os, 
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
	result->typeToOffsetInStream.resize(ELast);
	for (int i=0; i<ELast; ++i)
		result->typeToOffset[i] = result->typeToOffsetInStream[i] = -1;
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
				offset = result->data.size();
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
		int nParams = (int) accessor->getParam_array().getCount(),
			stride  = (int) accessor->getStride(),
			size    = (int) accessor->getCount();
		SAssert(nParams <= 4);

		Vec4 *target = new Vec4[size];
		for (int j=0; j<size; ++j)
			for (int k=0; k<nParams; ++k)
				target[j][k] = (Float) floatArray[j*stride+k];

		result->data[offset] = target;

		if (semantic == "POSITION") {
			SAssert(accessor->getStride() == 3);
			SAssert(result->typeToOffset[EPosition] == -1);
			result->typeToOffset[EPosition] = offset;
			result->typeToOffsetInStream[EPosition] = offsetInStream;
			result->glPos = new GLdouble[3*size];
			for (int k=0; k<3*size; ++k)
				result->glPos[k] = floatArray[k];
		} else if (semantic == "NORMAL") {
			SAssert(accessor->getStride() == 3);
			SAssert(result->typeToOffset[ENormal] == -1);
			result->hasNormals = true;
			result->typeToOffset[ENormal] = offset;
			result->typeToOffsetInStream[ENormal] = offsetInStream;
		} else if (semantic == "TEXCOORD") {
			SAssert(accessor->getStride() == 2 || accessor->getStride() == 3);
			if (result->typeToOffset[EUV] == -1) {
				result->hasUVs = true;
				result->typeToOffset[EUV] = offset;
				result->typeToOffsetInStream[EUV] = offsetInStream;
			} else {
				SLog(EWarn, "Found multiple sets of texture coordinates - ignoring!");
			}
		} else if (semantic == "COLOR") {
			SLog(EWarn, "Found per-vertex colors - ignoring. Please bake into a texture "
				"(Lighting/shading -> Batch Bake in Maya)");
			result->typeToOffset[EVertexColors] = offset;
			result->typeToOffsetInStream[EVertexColors] = offsetInStream;
		} else {
			SLog(EError, "Encountered an unknown source semantic: %s", semantic.c_str());
		}
	}
	SAssert(result->typeToOffset[EPosition] != -1);

	return result;
}

/// For using vertices as keys in an associative structure
struct vertex_key_order : public 
	std::binary_function<Vertex, Vertex, bool> {
	static int compare(const Vertex &v1, const Vertex &v2) {
		if (v1.v.x < v2.v.x) return -1;
		else if (v1.v.x > v2.v.x) return 1;
		if (v1.v.y < v2.v.y) return -1;
		else if (v1.v.y > v2.v.y) return 1;
		if (v1.v.z < v2.v.z) return -1;
		else if (v1.v.z > v2.v.z) return 1;
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

typedef std::map<SimpleTriangle, bool, triangle_key_order> TriangleMap;

void writeGeometry(std::string prefixName, std::string id, int geomIndex, std::string matID, Transform transform, 
		std::ostream &os, VertexData *vData, TriangleMap &triMap, const std::string &meshesDirectory) {
	std::vector<Vertex> vertexBuffer;
	std::vector<Triangle> triangles;
	std::map<Vertex, int, vertex_key_order> vertexMap;
	size_t numMerged = 0, triangleIdx = 0, duplicates = 0;
	Triangle triangle;
	if (tess_data.size() == 0)
		return;

	id += formatString("_%i", geomIndex);

	for (size_t i=0; i<tess_data.size(); i+=tess_nSources) {
		Vertex vertex;
		domUint posRef = tess_data[i+vData->typeToOffsetInStream[EPosition]];
		vertex.v = vData->data[vData->typeToOffset[EPosition]][posRef].toPoint();

		if (vData->typeToOffset[ENormal] != -1) {
			domUint normalRef = tess_data[i+vData->typeToOffsetInStream[ENormal]];
			vertex.n = vData->data[vData->typeToOffset[ENormal]][normalRef].toNormal();
		}

		if (vData->typeToOffset[EUV] != -1) {
			domUint uvRef = tess_data[i+vData->typeToOffsetInStream[EUV]];
			vertex.uv = vData->data[vData->typeToOffset[EUV]][uvRef].toPoint2();
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
			Point p0 = vertexBuffer[triangle.idx[0]].v,
				p1 = vertexBuffer[triangle.idx[1]].v,
				p2 = vertexBuffer[triangle.idx[2]].v;
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
			SLog(EWarn, "%s: Only contains duplicates of already-existing geometry. Ignoring.");
			os << "\t<!-- Ignored shape \"" << prefixName << "/" 
				<< id << "\" (mat=\"" << matID << "\"), since it only contains duplicate geometry. -->" << endl << endl;
			return;
		} else {
			SLog(EWarn, "Geometry contains %i duplicate triangles!", duplicates);
		}
	}

	SAssert(triangleIdx == 0);
	SLog(EInfo, "%s: Converted " SIZE_T_FMT " triangles, " SIZE_T_FMT 
		" vertices (merged " SIZE_T_FMT " vertices).", id.c_str(),
		triangles.size(), vertexBuffer.size(), numMerged);
	
	ref<TriMesh> mesh = new TriMesh(triangles.size(), vertexBuffer.size());
	std::copy(triangles.begin(), triangles.end(), mesh->getTriangles());
	std::copy(vertexBuffer.begin(), vertexBuffer.end(), mesh->getVertexBuffer());
	mesh->calculateTangentSpaceBasis(vData->typeToOffset[ENormal]!=-1, vData->typeToOffset[EUV]!=-1);

	std::string filename = meshesDirectory + id + std::string(".serialized");
	ref<FileStream> stream = new FileStream(filename, FileStream::ETruncReadWrite);
	stream->setByteOrder(Stream::ENetworkByteOrder);
	mesh->serialize(stream);
	stream->close();

	std::ostringstream matrix;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			matrix << transform.getMatrix()->m[i][j] << " ";
	std::string matrixValues = matrix.str();

	os << "\t<shape id=\"" << prefixName << "/" << id << "\" type=\"serialized\">" << endl;
	os << "\t\t<string name=\"filename\" value=\"" << filename << "\"/>" << endl;
	if (!transform.isIdentity()) {
		os << "\t\t<transform name=\"toWorld\">" << endl;
		os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
		os << "\t\t</transform>" << endl;
	}
	if (matID != "") 
		os << "\t\t<ref name=\"bsdf\" id=\"" << matID << "\"/>" << endl;
	os << "\t</shape>" << endl << endl;
}

void loadGeometry(std::string prefixName, Transform transform, std::ostream &os, domGeometry &geom, 
		StringMap &matLookupTable, const std::string &meshesDir) {
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

	SLog(EInfo, "Converting geometry \"%s\" (instantiated by %s)..", identifier.c_str(), 
		prefixName == "" ? "/" : prefixName.c_str());
	domMesh *mesh = geom.getMesh().cast();
	if (!mesh)
		SLog(EError, "Invalid geometry type encountered (must be a <mesh>)!");

	const domInputLocal_Array &vertInputs = mesh->getVertices()->getInput_array();

	int geomIndex = 0;

	domTriangles_Array &trianglesArray = mesh->getTriangles_array();
	for (size_t i=0; i<trianglesArray.getCount(); ++i) {
		domTriangles *triangles = trianglesArray[i];
		domInputLocalOffset_Array &inputs = triangles->getInput_array();
		VertexData *data = fetchVertexData(transform, os, vertInputs, inputs);
		domListOfUInts &indices = triangles->getP()->getValue();
		tess_data.clear();
		tess_nSources = data->nSources;
		for (size_t j=0; j<indices.getCount(); ++j)
			tess_data.push_back(indices[j]);
		std::string matID;
		if (triangles->getMaterial() == NULL || matLookupTable.find(triangles->getMaterial()) == matLookupTable.end())
			SLog(EWarn, "Referenced material could not be found, substituting a lambertian BRDF.");
		else
			matID = matLookupTable[triangles->getMaterial()];
		writeGeometry(prefixName, identifier, geomIndex, matID, transform, os, data, triMap, meshesDir);
		delete data;
		++geomIndex;
	}

	domPolygons_Array &polygonsArray = mesh->getPolygons_array();
	for (size_t i=0; i<polygonsArray.getCount(); ++i) {
		domPolygons *polygons = polygonsArray[i];
		domInputLocalOffset_Array &inputs = polygons->getInput_array();
		VertexData *data = fetchVertexData(transform, os, vertInputs, inputs);
		domP_Array &indexArray = polygons->getP_array();
		int posOffset = data->typeToOffset[EPosition];

		tess_data.clear();
		tess_nSources = data->nSources;
		for (size_t j=0; j<indexArray.getCount(); ++j) {
			domListOfUInts &indices = indexArray[j]->getValue();
			domUint *temp = new domUint[indices.getCount()];
			for (size_t l = 0; l<indices.getCount(); ++l)
				temp[l] = indices.get(l);

			gluTessBeginPolygon(tess, NULL);
			gluTessBeginContour(tess);

			for (size_t k=0; k<indices.getCount(); k+=data->nSources) 
				gluTessVertex(tess, &data->glPos[temp[k+posOffset]*3], (GLvoid *) (k+temp));

			gluTessEndContour(tess);
			gluTessEndPolygon(tess);
			delete[] temp;
		}

		std::string matID;
		if (polygons->getMaterial() == NULL || matLookupTable.find(polygons->getMaterial()) == matLookupTable.end())
			SLog(EWarn, "Referenced material could not be found, substituting a lambertian BRDF.");
		else
			matID = matLookupTable[polygons->getMaterial()];

		writeGeometry(prefixName, identifier, geomIndex, matID, transform, os, data, triMap, meshesDir);
		delete data;
		++geomIndex;
	}

	domPolylist_Array &polylistArray = mesh->getPolylist_array();
	for (size_t i=0; i<polylistArray.getCount(); ++i) {
		domPolylist *polylist = polylistArray[i];
		domInputLocalOffset_Array &inputs = polylist->getInput_array();
		VertexData *data = fetchVertexData(transform, os, vertInputs, inputs);
		domListOfUInts &vcount = polylist->getVcount()->getValue();
		domListOfUInts &indexArray = polylist->getP()->getValue();
		int posOffset = data->typeToOffset[EPosition], indexOffset = 0;
		tess_data.clear();
		tess_nSources = data->nSources;

		for (size_t j=0; j<vcount.getCount(); ++j) {
			size_t vertexCount = (size_t) vcount.get(j);

			domUint *temp = new domUint[vertexCount * data->nSources];
			for (size_t l = 0; l<vertexCount * data->nSources; ++l)
				temp[l] = indexArray.get(indexOffset++);

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
			SLog(EWarn, "Referenced material could not be found, substituting a lambertian BRDF.");
		else
			matID = matLookupTable[polylist->getMaterial()];

		writeGeometry(prefixName, identifier, geomIndex, matID, transform, os, data, triMap, meshesDir);
		delete data;
		++geomIndex;
	}
}

void loadMaterialParam(GeometryConverter *cvt, std::ostream &os, const std::string &name, StringMap &idToTexture, 
		domCommon_color_or_texture_type *value, bool handleRefs) {
	if (!value)
		return;
	domCommon_color_or_texture_type_complexType::domColor* color = 
		value->getColor().cast();
	domCommon_color_or_texture_type_complexType::domTexture* texture = 
		value->getTexture().cast();
	if (color && !handleRefs) {
		domFloat4 &colValue = color->getValue();
		if (cvt->m_srgb)
			os << "\t\t<srgb name=\"" << name << "\" value=\"";
		else
			os << "\t\t<rgb name=\"" << name << "\" value=\"";
		os << colValue.get(0) << " " << colValue.get(1) << " " 
		   << colValue.get(2) << "\"/>" << endl;
	} else if (texture && handleRefs) {
		if (idToTexture.find(texture->getTexture()) == idToTexture.end()) {
			SLog(EError, "Could not find referenced texture \"%s\"", texture->getTexture());
		} else {
			os << "\t\t<ref name=\"" << name << "\" id=\""
			   << idToTexture[texture->getTexture()] << "\"/>" << endl;
		}
	}
}

void loadMaterialParam(GeometryConverter *cvt, std::ostream &os, const std::string &name, StringMap &,
		domCommon_float_or_param_type *value, bool handleRef) {
	if (!value)
		return;
	domCommon_float_or_param_type::domFloat *floatValue = value->getFloat();
	if (!handleRef && floatValue) {
		os << "\t\t<float name=\"" << name << "\" value=\""
		   << floatValue->getValue() << "\"/>" << endl;
	}
}

void loadMaterial(GeometryConverter *cvt, std::ostream &os, domMaterial &mat, StringMap &_idToTexture) {
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
	StringMap idToTexture = _idToTexture;

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
			if (idToTexture.find(id) == idToTexture.end())
				SLog(EError, "Referenced bitmap '%s' not found!", id.c_str());
			idToTexture[newParam->getSid()] = idToTexture[id];
		}

		if (sampler2D) {
			std::string id = sampler2D->getSource()->getValue();
			if (idToTexture.find(id) == idToTexture.end())
				SLog(EError, "Referenced surface '%s' not found!", id.c_str());
			idToTexture[newParam->getSid()] = idToTexture[id];
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

		if (specular->getColor().cast()) {
			domFloat4 &colValue = specular->getColor()->getValue();
			if (colValue.get(0) == colValue.get(1) &&
				colValue.get(1) == colValue.get(2) &&
				colValue.get(2) == 0)
				isDiffuse = true;
		}
		if (isDiffuse) {
			os << "\t<bsdf id=\"" << identifier << "\" type=\"lambertian\">" << endl;
			loadMaterialParam(cvt, os, "reflectance", idToTexture, diffuse, false);
			loadMaterialParam(cvt, os, "reflectance", idToTexture, diffuse, true);
			os << "\t</bsdf>" << endl << endl;
		} else {
			os << "\t<bsdf id=\"" << identifier << "\" type=\"phong\">" << endl;
			os << "\t\t<float name=\"specularReflectance\" value=\"1\"/>" << endl;
			os << "\t\t<float name=\"diffuseReflectance\" value=\"1\"/>" << endl;
			loadMaterialParam(cvt, os, "diffuseColor", idToTexture, diffuse, false);
			loadMaterialParam(cvt, os, "specularColor", idToTexture, specular, false);
			loadMaterialParam(cvt, os, "exponent", idToTexture, shininess, false);
			loadMaterialParam(cvt, os, "diffuseColor", idToTexture, diffuse, true);
			loadMaterialParam(cvt, os, "specularColor", idToTexture, specular, true);
			loadMaterialParam(cvt, os, "exponent", idToTexture, shininess, true);
			os << "\t</bsdf>" << endl << endl;
		}
	} else if (lambert) {
		domCommon_float_or_param_type* transparency = lambert->getTransparency();
		domCommon_float_or_param_type::domFloat *transparencyValue = 
				transparency ? transparency->getFloat() : NULL;
		if (transparencyValue && transparencyValue->getValue() > 0.5) {
			os << "\t<bsdf id=\"" << identifier << "\" type=\"dielectric\"/>" << endl << endl;
		} else {
			domCommon_color_or_texture_type* diffuse = lambert->getDiffuse();
			os << "\t<bsdf id=\"" << identifier << "\" type=\"lambertian\">" << endl;
			loadMaterialParam(cvt, os, "reflectance", idToTexture, diffuse, false);
			loadMaterialParam(cvt, os, "reflectance", idToTexture, diffuse, true);
			os << "\t</bsdf>" << endl << endl;
		}
	} else if (blinn) {
		SLog(EWarn, "\"%s\": Encountered a \"blinn\" COLLADA material, which is currently "
			"unsupported in Mitsuba -- replacing it using a Phong material.", identifier.c_str());
		domCommon_color_or_texture_type* diffuse = blinn->getDiffuse();
		domCommon_color_or_texture_type* specular = blinn->getSpecular();
		domCommon_float_or_param_type* shininess = blinn->getShininess();
		bool isDiffuse = false;

		if (specular->getColor().cast()) {
			domFloat4 &colValue = specular->getColor()->getValue();
			if (colValue.get(0) == colValue.get(1) &&
				colValue.get(1) == colValue.get(2) &&
				colValue.get(2) == 0)
				isDiffuse = true;
		}
		if (isDiffuse) {
			os << "\t<bsdf id=\"" << identifier << "\" type=\"lambertian\">" << endl;
			loadMaterialParam(cvt, os, "reflectance", idToTexture, diffuse, false);
			loadMaterialParam(cvt, os, "reflectance", idToTexture, diffuse, true);
			os << "\t</bsdf>" << endl << endl;
		} else {
			os << "\t<bsdf id=\"" << identifier << "\" type=\"blinn\">" << endl;
			os << "\t\t<float name=\"specularReflectance\" value=\"1\"/>" << endl;
			os << "\t\t<float name=\"diffuseReflectance\" value=\"1\"/>" << endl;
			loadMaterialParam(cvt, os, "diffuseColor", idToTexture, diffuse, false);
			loadMaterialParam(cvt, os, "specularColor", idToTexture, specular, false);
			loadMaterialParam(cvt, os, "exponent", idToTexture, shininess, false);
			loadMaterialParam(cvt, os, "diffuseColor", idToTexture, diffuse, true);
			loadMaterialParam(cvt, os, "specularColor", idToTexture, specular, true);
			loadMaterialParam(cvt, os, "exponent", idToTexture, shininess, true);
			os << "\t</bsdf>" << endl << endl;
		}
	} else if (constant) {
		domCommon_float_or_param_type* transparency = constant->getTransparency();
		domCommon_float_or_param_type::domFloat *transparencyValue = 
				transparency ? transparency->getFloat() : NULL;
		if (transparencyValue && transparencyValue->getValue() > 0.5) {
			os << "\t<bsdf id=\"" << identifier << "\" type=\"dielectric\"/>" << endl << endl;
		} else {
			SLog(EWarn, "\"%s\": Encountered a \"constant\" COLLADA material, which is currently "
				"unsupported in Mitsuba -- replacing it using a Lambertian material.", identifier.c_str());
			os << "\t<bsdf id=\"" << identifier << "\" type=\"lambertian\"/>" << endl << endl;
		}
	} else {
		SLog(EError, "Material type not supported! (must be Lambertian/Phong/Blinn/Constant)");
	}
}

void loadLight(Transform transform, std::ostream &os, domLight &light) {
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

	SLog(EInfo, "Converting light \"%s\" ..", identifier.c_str());
	char *end_ptr = NULL;

	// Lights in Mitsuba point along the positive Z axis (COLLADA: neg. Z)
	transform = transform * Transform::scale(Vector(1, 1, -1));

	Point pos = transform(Point(0, 0, 0));
	Vector target = transform(Point(0, 0, 1));

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
		os << "\t<luminaire id=\"" << identifier << "\" type=\"point\">" << endl;
		os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl << endl;
		os << "\t\t<transform name=\"toWorld\">" << endl;
		os << "\t\t\t<translate x=\"" << pos.x << "\" y=\"" << pos.y << "\" z=\"" << pos.z << "\"/>" << endl;
		os << "\t\t</transform>" << endl;
		os << "\t</luminaire>" << endl << endl;
	}

	domLight::domTechnique_common::domDirectional *directional = light.getTechnique_common()->getDirectional().cast();
	if (directional) {
		domFloat3 &color = directional->getColor()->getValue();
		os << "\t<luminaire id=\"" << identifier << "\" type=\"directional\">" << endl;
		os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl << endl;
		os << "\t\t<transform name=\"toWorld\">" << endl;
		os << "\t\t\t<lookAt ox=\"" << pos.x << "\" oy=\"" << pos.y << "\" oz=\"" << pos.z << "\" tx=\"" << target.x << "\" ty=\"" << target.y << "\" tz=\"" << target.z << "\"/>" << endl;
		os << "\t\t</transform>" << endl << endl;
		os << "\t</luminaire>" << endl << endl;
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
		os << "\t<luminaire id=\"" << identifier << "\" type=\"spot\">" << endl;
		os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl;
		os << "\t\t<float name=\"cutoffAngle\" value=\"" << falloffAngle/2 << "\"/>" << endl << endl;
		os << "\t\t<transform name=\"toWorld\">" << endl;
		os << "\t\t\t<lookAt ox=\"" << pos.x << "\" oy=\"" << pos.y << "\" oz=\"" << pos.z << "\" tx=\"" << target.x << "\" ty=\"" << target.y << "\" tz=\"" << target.z << "\"/>" << endl;
		os << "\t\t</transform>" << endl;
		os << "\t</luminaire>" << endl << endl;
	}
	domLight::domTechnique_common::domAmbient *ambient = light.getTechnique_common()->getAmbient().cast();
	if (ambient) {
		domFloat3 &color = ambient->getColor()->getValue();
		os << "\t<luminaire id=\"" << identifier << "\" type=\"constant\">" << endl;
		os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl;
		os << "\t</luminaire>" << endl << endl;
	}
	if (!point && !spot && !ambient && !directional)
		SLog(EWarn, "Encountered an unknown light type!");
}

void loadImage(GeometryConverter *cvt, std::ostream &os, const std::string &textureDir, 
		domImage &image, StringMap &idToTexture, StringMap &fileToId) {
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

	SLog(EInfo, "Converting texture \"%s\" ..", identifier.c_str());

	std::string filename = cdom::uriToFilePath(image.getInit_from()->getValue().str());
	if (fileToId.find(filename) != fileToId.end()) {
		idToTexture[identifier] = fileToId[filename];
		return;
	}

	idToTexture[identifier] = identifier;
	fileToId[filename] = identifier;

	boost::filesystem::path path = boost::filesystem::path(filename, boost::filesystem::native);
	std::string targetPath = textureDir + path.leaf();

	if (endsWith(filename, ".rgb")) 
		SLog(EWarn, "Maya RGB images must be converted to PNG, EXR or JPEG! The 'imgcvt' "
		"utility found in the Maya binary directory can be used to do this.");

	if (!FileStream::exists(targetPath)) {
		ref<FileResolver> fRes = FileResolver::getInstance();
		std::string resolved = fRes->resolve(path.leaf());
		if (!FileStream::exists(filename)) {
			if (!FileStream::exists(resolved)) {
				SLog(EWarn, "Found neither \"%s\" nor \"%s\"!", filename.c_str(), resolved.c_str());
				filename = cvt->locateResource(filename);
				if (filename == "")
					SLog(EError, "Unable to locate a resource -- aborting conversion.");
			} else {
				filename = resolved;
			}
		}
		ref<FileStream> input = new FileStream(filename, FileStream::EReadOnly);
		ref<FileStream> output = new FileStream(targetPath, FileStream::ETruncReadWrite);
		input->copyTo(output);
		input->close();
		output->close();
	}

	os << "\t<texture id=\"" << identifier << "\" type=\"ldrtexture\">" << endl;
	os << "\t\t<string name=\"filename\" value=\"" << textureDir + path.leaf() << "\"/>" << endl;
	os << "\t</texture>" << endl << endl;
}

void loadCamera(GeometryConverter *cvt, Transform transform, std::ostream &os, domCamera &camera) {
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

	SLog(EInfo, "Converting camera \"%s\" ..", identifier.c_str());
	Float aspect = 1.0f;
	int xres=768;

	// Cameras in Mitsuba point along the positive Z axis (COLLADA: neg. Z)
	transform = transform * Transform::scale(Vector(1,1,-1));

	std::ostringstream matrix;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			matrix << transform.getMatrix()->m[i][j] << " ";
	std::string matrixValues = matrix.str();

	domCamera::domOptics::domTechnique_common::domOrthographic* ortho = camera.getOptics()->
		getTechnique_common()->getOrthographic().cast();
	if (ortho) {
		if (ortho->getAspect_ratio().cast() != 0)
			aspect = (Float) ortho->getAspect_ratio()->getValue();
		if (cvt->m_xres != -1) {
			xres = cvt->m_xres;
			aspect = (Float) cvt->m_xres / (Float) cvt->m_yres;
		}
		os << "\t<camera id=\"" << identifier << "\" type=\"orthographic\">" << endl;
	}

	domCamera::domOptics::domTechnique_common::domPerspective* persp = camera.getOptics()->
		getTechnique_common()->getPerspective().cast();
	if (persp) {
		if (cvt->m_xres != -1) {
			xres = cvt->m_xres;
			aspect = (Float) cvt->m_xres / (Float) cvt->m_yres;
		} else {
			if (persp->getAspect_ratio().cast() != 0) {
				aspect = (Float) persp->getAspect_ratio()->getValue();
				if (std::abs(aspect-0.1) < Epsilon) {
					SLog(EWarn, "Found the suspicious aspect ratio \"0.1\", which is likely due to a bug in Blender 2.5"
						" - setting to 1.0. Please use the \"-r\" parameter to override the resolution.");
					aspect = 1.0f;
				}
			}
		}
		os << "\t<camera id=\"" << identifier << "\" type=\"perspective\">" << endl;
		if (persp->getXfov().cast()) {
			Float xFov = (Float) persp->getXfov()->getValue();
			if (std::abs(xFov-1.0f) < Epsilon && cvt->m_fov == -1) {
				SLog(EWarn, "Found the suspicious field of view value \"1.0\", which is likely due to a bug in Blender 2.5"
					" - setting to 45deg. Please use the \"-f\" parameter to override this.");
				xFov = 45.0f;
			}
			Float yFov = radToDeg(2 * std::atan(std::tan(degToRad(xFov)/2) / aspect));
			if (cvt->m_fov != -1)
				xFov = yFov = cvt->m_fov;
			if (aspect <= 1.0f)
				os << "\t\t<float name=\"fov\" value=\"" << xFov << "\"/>" << endl;
			else
				os << "\t\t<float name=\"fov\" value=\"" << yFov << "\"/>" << endl;
		} else if (persp->getYfov().cast()) {
			Float yFov = (Float) persp->getYfov()->getValue();
			if (std::abs(yFov-1.0) < Epsilon && cvt->m_fov == -1) {
				SLog(EWarn, "Found the suspicious field of view value \"1.0\", which is likely due to a bug in Blender 2.5"
					" - setting to 45deg. Please use the \"-f\" parameter to override this.");
				yFov = 45.0f;
			}
			Float xFov = radToDeg(2 * std::atan(std::tan(degToRad(yFov)/2) * aspect));
			if (cvt->m_fov != -1)
				xFov = yFov = cvt->m_fov;
			if (aspect > 1.0f)
				os << "\t\t<float name=\"fov\" value=\"" << yFov << "\"/>" << endl;
			else
				os << "\t\t<float name=\"fov\" value=\"" << xFov << "\"/>" << endl;
		}
		os << "\t\t<float name=\"nearClip\" value=\"" << persp->getZnear()->getValue() << "\"/>" << endl;
		os << "\t\t<float name=\"farClip\" value=\"" << persp->getZfar()->getValue() << "\"/>" << endl;
		os << "\t\t<boolean name=\"mapSmallerSide\" value=\"" << (cvt->m_mapSmallerSide ? "true" : "false") << "\"/>" << endl;
	}

	os << endl;
	os << "\t\t<transform name=\"toWorld\">" << endl;
	os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
	os << "\t\t</transform>" << endl << endl;
	os << "\t\t<sampler type=\"ldsampler\">" << endl;
	os << "\t\t\t<integer name=\"sampleCount\" value=\"" << cvt->m_samplesPerPixel << "\"/>" << endl;
	os << "\t\t</sampler>" << endl << endl;
	os << "\t\t<film type=\"exrfilm\">" << endl;
	os << "\t\t\t<integer name=\"width\" value=\"" << xres << "\"/>" << endl;
	os << "\t\t\t<integer name=\"height\" value=\"" << (int) (xres/aspect) << "\"/>" << endl;
	os << "\t\t\t<rfilter type=\"gaussian\"/>" << endl;
	os << "\t\t</film>" << endl;
	os << "\t</camera>" << endl << endl;
}

void loadNode(GeometryConverter *cvt, Transform transform, std::ostream &os, 
		domNode &node, std::string prefixName, const std::string &meshesDir) {
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
	SLog(EInfo, "Converting node \"%s\" ..", identifier.c_str());

	daeTArray<daeSmartRef<daeElement> > children = node.getChildren();
	/* Parse transformations */
	for (size_t i=0; i<children.getCount(); ++i) {
		daeElement *element = children.get(i);
		if (element->typeID() == domRotate::ID()) {
			/* Skip rotations labeled as "post-rotationY". Maya exports these with some cameras,
			   which introduces an incorrect 90 degree rotation unless ignored */
			if (element->hasAttribute("sid") && element->getAttribute("sid") == "post-rotationY") 
				continue;
			daeTArray<double> value = daeSafeCast<domRotate>(element)->getValue();
			transform = transform *
				Transform::rotate(Vector((Float) value.get(0), (Float) value.get(1), (Float) value.get(2)), (Float) value.get(3));
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
			ref<Matrix4x4> matrix = new Matrix4x4(
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
		loadGeometry(prefixName, transform, os, *geom, matLookupTable, meshesDir);
	}

	/* Iterate over all light references */
	domInstance_light_Array &instanceLights = node.getInstance_light_array();
	for (size_t i=0; i<instanceLights.getCount(); ++i) {
		domInstance_light *inst = instanceLights[i];
		domLight *light = daeSafeCast<domLight>(inst->getUrl().getElement());
		if (!light)
			SLog(EError, "Could not find a referenced light!");
		loadLight(transform, os, *light);
	}

	/* Iterate over all camera references */
	domInstance_camera_Array &instanceCameras = node.getInstance_camera_array();
	for (size_t i=0; i<instanceCameras.getCount(); ++i) {
		domInstance_camera *inst = instanceCameras[i];
		domCamera *camera = daeSafeCast<domCamera>(inst->getUrl().getElement());
		if (camera == NULL)
			SLog(EError, "Could not find a referenced camera!");
		loadCamera(cvt, transform, os, *camera);
	}

	/* Recursively iterate through sub-nodes */
	domNode_Array &nodes = node.getNode_array();
	for (size_t i=0; i<nodes.getCount(); ++i) 
		loadNode(cvt, transform, os, *nodes[i], prefixName + std::string("/") + identifier, meshesDir);

	/* Recursively iterate through <instance_node> elements */
	domInstance_node_Array &instanceNodes = node.getInstance_node_array();
	for (size_t i=0; i<instanceNodes.getCount(); ++i) {
		domNode *node = daeSafeCast<domNode>(instanceNodes[i]->getUrl().getElement());
		if (!node)
			SLog(EError, "Could not find a referenced node!");
		loadNode(cvt, transform, os, *node, prefixName + std::string("/") + identifier, meshesDir);
	}
}

#ifndef WIN32
#define __stdcall
#endif

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
	SLog(EWarn, "Detected a self-intersecting polygon!");
}

GLvoid __stdcall tessEnd() {
}

GLvoid __stdcall tessError(GLenum error) {
	SLog(EError, "The GLU tesselator generated an error: %s!", gluErrorString(error));
}

GLvoid __stdcall tessEdgeFlag(GLboolean) {
}

void GeometryConverter::convertCollada(const std::string &inputFile, 
	std::ostream &os,
	const std::string &textureDirectory,
	const std::string &meshesDirectory) {
#if defined(__LINUX__)
	std::string path = std::string("file://") +
		FileResolver::getInstance()->resolveAbsolute(inputFile);
#else
	std::string path = inputFile;
#endif

	DAE *dae = new DAE();
	SLog(EInfo, "Loading \"%s\" ..", path.c_str());
	if (dae->load(path.c_str()) != DAE_OK) 
		SLog(EError, "Could not load \"%s\"!", path.c_str());

	domCOLLADA *document = dae->getDom(path.c_str());
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
	os << "<scene>" << endl;
	os << "\t<integrator id=\"integrator\" type=\"direct\"/>" << endl << endl;

	StringMap idToTexture, fileToId;
	domLibrary_images_Array &libraryImages = document->getLibrary_images_array();
	for (size_t i=0; i<libraryImages.getCount(); ++i) {
		domImage_Array &images = libraryImages[i]->getImage_array();
		for (size_t j=0; j<images.getCount(); ++j) 
			loadImage(this, os, textureDirectory, *images.get(j), idToTexture, fileToId);
	}

	domLibrary_materials_Array &libraryMaterials = document->getLibrary_materials_array();
	for (size_t i=0; i<libraryMaterials.getCount(); ++i) {
		domMaterial_Array &materials = libraryMaterials[i]->getMaterial_array();
		for (size_t j=0; j<materials.getCount(); ++j) 
			loadMaterial(this, os, *materials.get(j), idToTexture);
	}

	for (size_t i=0; i<nodes.getCount(); ++i) 
		loadNode(this, Transform(), os, *nodes[i], "", meshesDirectory);

	os << "</scene>" << endl;

	gluDeleteTess(tess);
	delete dae;
}

