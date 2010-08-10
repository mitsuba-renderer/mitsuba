/**
 * Mitsuba COLLADA -> XML converter
 * 
 * Takes a DAE file and turns it into a scene description and separate mesh files
 * using a compact binary format. All associated files are copied into newly created 
 * 'textures' and 'meshes' directories
 *
 * Currently supports the following subset of the COLLADA specification:
 * - Arbitrary polygonal meshes
 * - Lambert and Phong materials (allowed to be textured)
 * - Cameras
 * - Spot and Point and Ambient lights
 *
 * When exporting using Maya/FBX, be sure to have it convert all NURBS surfaces into
 * "Software Render Meshes". Triangulation is not required (the code below does this
 * automatically for arbitrary polygonal meshes). The Light and camera export options
 * should be activated, since they are off by default. While modeling the scene, it is
 * advisable to use light sources with an inverse square falloff. Otherwise, the 
 * illumination will be competely different when rendering in Mitsuba (the image might 
 * be pitch black). Note that most BRDFs in Mitsuba treat surfaces as one-sided, thus they
 * will appear black when seen from the back.
 *
 * The conversion barfs when it gets more than 10MB in one single XML string 
 * (error: xmlSAX2Characters: huge text node: out of memory). In this case, split the
 * mesh into smaller pieces or recompile libxml with a higher limit.
 *
 * Since Mitsuba does not support per-vertex colors and prefers textures, any vertex colors 
 * part of the input file are not converted and should instead be baked to textures beforehand 
 * (e.g. using Lighting/shading -> Batch bake in Maya).
 */

#define BOOST_FILESYSTEM_NO_LIB 
#define BOOST_SYSTEM_NO_LIB 

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/hw/glrenderer.h>
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

using namespace mitsuba;

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
	std::vector<ESourceType> offsetToType;
	std::vector<int> typeToOffset;

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
	
/* Globals */
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
	result->data.resize(inputs.getCount());
	for (size_t i=0; i<inputs.getCount(); ++i)
		result->data[i] = NULL;
	result->offsetToType.resize(inputs.getCount());
	result->typeToOffset.resize(ELast);
	for (int i=0; i<ELast; ++i)
		result->typeToOffset[i] = -1;

	for (size_t i=0; i<inputs.getCount(); ++i) {
		int offset = (int) inputs[i]->getOffset();
		daeURI &sourceRef = inputs[i]->getSource();
		sourceRef.resolveElement();
		domSource *source = daeSafeCast<domSource>(sourceRef.getElement());

		if (!strcmp(inputs[i]->getSemantic(), "VERTEX")) {
			SAssert(vertInputs.getCount() == 1);
			sourceRef = vertInputs[0]->getSource();
			sourceRef.resolveElement();
			source = daeSafeCast<domSource>(sourceRef.getElement());
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
		for (int j=0; j<size; ++j) {
			for (int k=0; k<nParams; ++k)
				target[j][k] = (Float) floatArray[j*stride+k];
		}

		result->data[offset] = target;

		if (!strcmp(inputs[i]->getSemantic(), "VERTEX")) {
			SAssert(accessor->getStride() == 3);
			SAssert(result->typeToOffset[EPosition] == -1);
			result->offsetToType[offset] = EPosition;
			result->typeToOffset[EPosition] = offset;
			result->glPos = new GLdouble[3*size];
			for (int k=0; k<3*size; ++k)
				result->glPos[k] = floatArray[k];
		} else if (!strcmp(inputs[i]->getSemantic(), "NORMAL")) {
			SAssert(accessor->getStride() == 3);
			SAssert(result->typeToOffset[ENormal] == -1);
			result->hasNormals = true;
			result->offsetToType[offset] = ENormal;
			result->typeToOffset[ENormal] = offset;
		} else if (!strcmp(inputs[i]->getSemantic(), "TEXCOORD")) {
			SAssert(accessor->getStride() == 2 || accessor->getStride() == 3);
			if (result->typeToOffset[EUV] == -1) {
				result->hasUVs = true;
				result->offsetToType[offset] = EUV;
				result->typeToOffset[EUV] = offset;
			} else {
				SLog(EWarn, "Found multiple texture coordinate records - ignoring!");
			}
		} else if (!strcmp(inputs[i]->getSemantic(), "COLOR")) {
			SLog(EWarn, "Found per-vertex colors - ignoring. Please bake into a texture (Lighting/shading -> Batch Bake in Maya)");
			result->offsetToType[offset] = EVertexColors;
			result->typeToOffset[EVertexColors] = offset;
		} else {
			SLog(EError, "Encountered an unknown source semantic: %s", 
				inputs[i]->getSemantic());
		}
	}
	SAssert(result->typeToOffset[EPosition] != -1);

	return result;
}

/// For using vertices as keys in an associative structure
struct vertex_key_order : public 
	std::binary_function<Vertex, Vertex, bool> {
public:
	bool operator()(const Vertex &v1, const Vertex &v2) const {
		if (v1.v.x < v2.v.x) return true;
		else if (v1.v.x > v2.v.x) return false;
		if (v1.v.y < v2.v.y) return true;
		else if (v1.v.y > v2.v.y) return false;
		if (v1.v.z < v2.v.z) return true;
		else if (v1.v.z > v2.v.z) return false;
		if (v1.n.x < v2.n.x) return true;
		else if (v1.n.x > v2.n.x) return false;
		if (v1.n.y < v2.n.y) return true;
		else if (v1.n.y > v2.n.y) return false;
		if (v1.n.z < v2.n.z) return true;
		else if (v1.n.z > v2.n.z) return false;
		if (v1.uv.x < v2.uv.x) return true;
		else if (v1.uv.x > v2.uv.x) return false;
		if (v1.uv.y < v2.uv.y) return true;
		else if (v1.uv.y > v2.uv.y) return false;
		return false;
	}
};

void writeGeometry(std::string id, int geomIndex, std::string matID, Transform transform, std::ostream &os, VertexData *vData) {
	std::vector<Vertex> vertexBuffer;
	std::vector<Triangle> triangles;
	std::map<Vertex, int, vertex_key_order> vertexMap;
	size_t numMerged = 0, triangleIdx = 0;
	Triangle triangle;

	id += formatString("_%i", geomIndex);

	for (size_t i=0; i<tess_data.size(); i+=tess_nSources) {
		Vertex vertex;
		domUint posRef = tess_data[i+vData->typeToOffset[EPosition]];
		vertex.v = vData->data[vData->typeToOffset[EPosition]][posRef].toPoint();

		if (vData->typeToOffset[ENormal] != -1) {
			domUint normalRef = tess_data[i+vData->typeToOffset[ENormal]];
			vertex.n = vData->data[vData->typeToOffset[ENormal]][normalRef].toNormal();
		}
		
		if (vData->typeToOffset[EUV] != -1) {
			domUint uvRef = tess_data[i+vData->typeToOffset[EUV]];
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
			triangles.push_back(triangle);
			triangleIdx = 0;
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

	std::string filename = formatString("meshes/%s.serialized", id.c_str());
	ref<FileStream> stream = new FileStream(filename, FileStream::ETruncReadWrite);
	stream->setByteOrder(Stream::ENetworkByteOrder);
	mesh->serialize(stream);
	stream->close();

	std::ostringstream matrix;
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			matrix << transform.getMatrix()->m[i][j] << " ";
	std::string matrixValues = matrix.str();

	os << "\t<shape id=\"" << id << "\" type=\"serialized\">" << endl;
	os << "\t\t<string name=\"filename\" value=\"" << filename.c_str() << "\"/>" << endl;
	if (!transform.isIdentity()) {
		os << "\t\t<transform name=\"toWorld\">" << endl;
		os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
		os << "\t\t</transform>" << endl;
	}
	os << "\t\t<ref name=\"bsdf\" id=\"" << matID << "\"/>" << endl;
	os << "\t</shape>" << endl << endl;
}

void loadGeometry(std::string nodeName, Transform transform, std::ostream &os, domGeometry &geom, 
		StringMap &matLookupTable) {
	SLog(EInfo, "Converting geometry \"%s\" (instantiated by %s)..", geom.getName(), nodeName.c_str());

	domMesh *mesh = geom.getMesh().cast();
	if (!mesh)
		SLog(EError, "Invalid geometry type encountered (must be a <mesh>)!");

	const domInputLocal_Array &vertInputs = mesh->getVertices()->getInput_array();

	std::string xmlName = nodeName + std::string("_") + geom.getId();
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
		if (matLookupTable.find(triangles->getMaterial()) == matLookupTable.end())
			SLog(EError, "Referenced material could not be found!");
		writeGeometry(xmlName, geomIndex, matLookupTable[triangles->getMaterial()], transform, os, data);
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
	
		if (polygons->getMaterial() == NULL) 
			SLog(EError, "No material reference specified!");

		if (matLookupTable.find(polygons->getMaterial()) == matLookupTable.end())
			SLog(EError, "Referenced material could not be found!");

		writeGeometry(xmlName, geomIndex, matLookupTable[polygons->getMaterial()], transform, os, data);
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
			size_t vertexCount = vcount.get(i);

			domUint *temp = new domUint[vertexCount * data->nSources];
			for (size_t l = 0; l<vertexCount * data->nSources; ++l)
				temp[l] = indexArray.get(indexOffset++);

			gluTessBeginPolygon(tess, NULL);
			gluTessBeginContour(tess);

			for (size_t k=0; k<vertexCount*data->nSources; k+=data->nSources) 
				gluTessVertex(tess, &data->glPos[temp[k+posOffset]*3], (GLvoid *) (k+temp));

			gluTessEndContour(tess);
			gluTessEndPolygon(tess);
			delete[] temp;
		}

		if (polylist->getMaterial() == NULL) 
			SLog(EError, "No material reference specified!");

		if (matLookupTable.find(polylist->getMaterial()) == matLookupTable.end())
			SLog(EError, "Referenced material could not be found!");

		writeGeometry(xmlName, geomIndex, matLookupTable[polylist->getMaterial()], transform, os, data);
		delete data;
		++geomIndex;
	}
}

void loadMaterialParam(std::ostream &os, const std::string &name, StringMap &idToTexture, 
		domCommon_color_or_texture_type *value, bool handleRefs) {
	if (!value)
		return;
	domCommon_color_or_texture_type_complexType::domColor* color = 
		value->getColor().cast();
	domCommon_color_or_texture_type_complexType::domTexture* texture = 
		value->getTexture().cast();
	if (color && !handleRefs) {
		domFloat4 &colValue = color->getValue();
		os << "\t\t<rgb name=\"" << name << "\" value=\""
		   << colValue.get(0) << " " << colValue.get(1) << " " 
		   << colValue.get(2) << "\"/>" << endl;
	} else if (texture && handleRefs) {
		os << "\t\t<ref name=\"" << name << "\" id=\""
		   << idToTexture[texture->getTexture()] << "\"/>" << endl;
	}
}

void loadMaterialParam(std::ostream &os, const std::string &name, StringMap &,
		domCommon_float_or_param_type *value, bool handleRef) {
	if (!value)
		return;
	domCommon_float_or_param_type::domFloat *floatValue = value->getFloat();
	if (!handleRef && floatValue) {
		os << "\t\t<float name=\"" << name << "\" value=\""
		   << floatValue->getValue() << "\"/>" << endl;
	}
}

void loadMaterial(std::ostream &os, domMaterial &mat, StringMap &idToTexture) {
	SLog(EInfo, "Converting material \"%s\" ..", mat.getName());

	daeURI &effRef = mat.getInstance_effect()->getUrl();
	effRef.resolveElement();
	domEffect *effect = daeSafeCast<domEffect>(effRef.getElement());
	
	if (!effect)
		SLog(EError, "Referenced effect not found!");

	domProfile_COMMON *commonProfile = daeSafeCast<domProfile_COMMON>
		(effect->getDescendant("profile_COMMON"));
	
	if (!commonProfile)
		SLog(EError, "Common effect profile not found!");

	domProfile_COMMON::domTechnique *technique = commonProfile->getTechnique();
	
	if (!technique)
		SLog(EError, "The technique element is missing!");

	domProfile_COMMON::domTechnique::domPhong* phong = technique->getPhong();
	domProfile_COMMON::domTechnique::domLambert* lambert = technique->getLambert();

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
			os << "\t<bsdf id=\"" << mat.getId() << "\" type=\"lambertian\">" << endl;
			loadMaterialParam(os, "reflectance", idToTexture, diffuse, false);
			loadMaterialParam(os, "reflectance", idToTexture, diffuse, true);
			os << "\t</bsdf>" << endl << endl;
		} else {
			os << "\t<bsdf id=\"" << mat.getId() << "\" type=\"phong\">" << endl;
			loadMaterialParam(os, "diffuseReflectance", idToTexture, diffuse, false);
			loadMaterialParam(os, "specularReflectance", idToTexture, specular, false);
			loadMaterialParam(os, "exponent", idToTexture, shininess, false);
			loadMaterialParam(os, "diffuseReflectance", idToTexture, diffuse, true);
			loadMaterialParam(os, "specularReflectance", idToTexture, specular, true);
			loadMaterialParam(os, "exponent", idToTexture, shininess, true);
			os << "\t</bsdf>" << endl << endl;
		}
	} else if (lambert) {
		domCommon_color_or_texture_type* diffuse = lambert->getDiffuse();
		os << "\t<bsdf id=\"" << mat.getId() << "\" type=\"lambertian\">" << endl;
		loadMaterialParam(os, "reflectance", idToTexture, diffuse, false);
		loadMaterialParam(os, "reflectance", idToTexture, diffuse, true);
		os << "\t</bsdf>" << endl << endl;
	} else {
		SLog(EError, "Material type not supported! (must be Lambertian/Phong)");
	}
}

void loadLight(Transform transform, std::ostream &os, domLight &light) {
	SLog(EInfo, "Converting light \"%s\" ..", light.getName());
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
		domFloat3 &color = point->getColor()->getValue();
		os << "\t<luminaire id=\"" << light.getId() << "\" type=\"point\">" << endl;
		os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl << endl;
		os << "\t\t<transform name=\"toWorld\">" << endl;
		os << "\t\t\t<translate x=\"" << pos.x << "\" y=\"" << pos.y << "\" z=\"" << pos.z << "\"/>" << endl;
		os << "\t\t</transform>" << endl;
		os << "\t</luminaire>" << endl << endl;
	}
	domLight::domTechnique_common::domSpot *spot = light.getTechnique_common()->getSpot().cast();
	if (spot) {
		domFloat3 &color = spot->getColor()->getValue();
		Float falloffAngle = 180.0f;
		if (spot->getFalloff_angle())
			falloffAngle = (Float) spot->getFalloff_angle()->getValue();
		os << "\t<luminaire id=\"" << light.getId() << "\" type=\"spot\">" << endl;
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
		os << "\t<luminaire id=\"" << light.getId() << "\" type=\"constant\">" << endl;
		os << "\t\t<rgb name=\"intensity\" value=\"" << color[0]*intensity << " " << color[1]*intensity << " " << color[2]*intensity << "\"/>" << endl;
		os << "\t</luminaire>" << endl << endl;
	}
	if (!point && !spot && !ambient)
		SLog(EWarn, "Encountered an unknown light type!");
}

void loadImage(std::ostream &os, domImage &image, StringMap &idToTexture, StringMap &fileToId) {
	SLog(EInfo, "Converting texture \"%s\" ..", image.getName());

	std::string filename = cdom::uriToFilePath(image.getInit_from()->getValue().str());
	if (fileToId.find(filename) != fileToId.end()) {
		idToTexture[image.getId()] = fileToId[filename];
		return;
	}

	idToTexture[image.getId()] = image.getId();
	fileToId[filename] = image.getId();

	boost::filesystem::path path = boost::filesystem::path(filename, boost::filesystem::native);
	std::string leaf = path.leaf();

	if (endsWith(leaf, ".rgb")) 
		SLog(EWarn, "Maya RGB images must be converted to PNG, EXR or JPEG! The 'imgcvt' "
		"utility found in the Maya binary directory can be used to do this.");

	if (!FileStream::exists(filename)) {
		if (!FileStream::exists(leaf))
			SLog(EError, "Found neither \"%s\" nor \"%s\"!", filename.c_str(), leaf.c_str());
		filename = leaf;
	}

	ref<FileStream> input = new FileStream(filename, FileStream::EReadOnly);
	ref<FileStream> output = new FileStream(formatString("textures/%s", leaf.c_str()), FileStream::ETruncReadWrite);
	input->copyTo(output);

	os << "\t<texture id=\"" << image.getId() << "\" type=\"ldrtexture\">" << endl;
	os << "\t\t<string name=\"filename\" value=\"textures/" << leaf << "\"/>" << endl;
	os << "\t</texture>" << endl << endl;
}

void loadCamera(Transform transform, std::ostream &os, domCamera &camera) {
	SLog(EInfo, "Converting camera \"%s\" ..", camera.getName());
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
		os << "\t<camera id=\"" << camera.getId() << "\" type=\"orthographic\">" << endl;
	}

	domCamera::domOptics::domTechnique_common::domPerspective* persp = camera.getOptics()->
		getTechnique_common()->getPerspective().cast();
	if (persp) {
		if (persp->getAspect_ratio().cast() != 0)
			aspect = (Float) persp->getAspect_ratio()->getValue();
		os << "\t<camera id=\"" << camera.getId() << "\" type=\"perspective\">" << endl;
		if (persp->getXfov().cast()) {
			Float yFov = radToDeg(2 * std::atan(std::tan(degToRad((Float) persp->getXfov()->getValue())/2) / aspect));
			if (aspect <= 1.0f)
				os << "\t\t<float name=\"fov\" value=\"" << persp->getXfov()->getValue() << "\"/>" << endl;
			else
				os << "\t\t<float name=\"fov\" value=\"" << yFov << "\"/>" << endl;
			os << "\t\t<float name=\"nearClip\" value=\"" << persp->getZnear()->getValue() << "\"/>" << endl;
			os << "\t\t<float name=\"farClip\" value=\"" << persp->getZfar()->getValue() << "\"/>" << endl;
		} else if (persp->getYfov().cast()) {
			Float xFov = radToDeg(2 * std::atan(std::tan(degToRad((Float) persp->getYfov()->getValue())/2) * aspect));
			if (aspect > 1.0f)
				os << "\t\t<float name=\"fov\" value=\"" << persp->getYfov()->getValue() << "\"/>" << endl;
			else
				os << "\t\t<float name=\"fov\" value=\"" << xFov << "\"/>" << endl;
			os << "\t\t<float name=\"nearClip\" value=\"" << persp->getZnear()->getValue() << "\"/>" << endl;
			os << "\t\t<float name=\"farClip\" value=\"" << persp->getZfar()->getValue() << "\"/>" << endl;
		}
	}

	os << endl;
	os << "\t\t<transform name=\"toWorld\">" << endl;
	os << "\t\t\t<matrix value=\"" << matrixValues.substr(0, matrixValues.length()-1) << "\"/>" << endl;
	os << "\t\t</transform>" << endl << endl;

	os << "\t\t<sampler type=\"stratified\"/>" << endl << endl;
	os << "\t\t<film type=\"exrfilm\">" << endl;
	os << "\t\t\t<integer name=\"width\" value=\"" << xres << "\"/>" << endl;
	os << "\t\t\t<integer name=\"height\" value=\"" << (int) (xres/aspect) << "\"/>" << endl;
	os << "\t\t\t<rfilter type=\"gaussian\"/>" << endl;
	os << "\t\t</film>" << endl;
	os << "\t</camera>" << endl << endl;
}

void loadNode(Transform transform, std::ostream &os, domNode &node) {
	SLog(EInfo, "Converting node \"%s\" ..", node.getName());

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
		}
	}
	
	/* Iterate over all geometry references */
	domInstance_geometry_Array &instanceGeometries = node.getInstance_geometry_array();
	for (size_t i=0; i<instanceGeometries.getCount(); ++i) {
		domInstance_geometry *inst = instanceGeometries[i];
		domGeometry *geom = daeSafeCast<domGeometry>(inst->getUrl().getElement());
		domBind_material *bmat = inst->getBind_material();
		if (!bmat)
			SLog(EError, "instance_geometry does not contain a <bind_material> element!");
		domBind_material::domTechnique_common *technique = bmat->getTechnique_common();
		if (!technique)
			SLog(EError, "bind_material does not contain a <technique_common> element!");
		domInstance_material_Array &instMaterials = technique->getInstance_material_array();

		StringMap matLookupTable;
		for (size_t i=0; i<instMaterials.getCount(); ++i) {
			domInstance_material *instMat = instMaterials[i];
			daeURI &matRef = instMat->getTarget();
			matRef.resolveElement();
			domMaterial *material = daeSafeCast<domMaterial>(matRef.getElement());
			matLookupTable[instMat->getSymbol()] = material->getId();
		}

		if (!geom)
			SLog(EError, "Could not find a referenced geometry object!");
		loadGeometry(node.getName(), transform, os, *geom, matLookupTable);
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
		loadCamera(transform, os, *camera);
	}

	/* Recursively iterate through sub-nodes */
	domNode_Array &nodes = node.getNode_array();
	for (size_t i=0; i<nodes.getCount(); ++i) 
		loadNode(transform, os, *nodes[i]);

	/* Recursively iterate through <instance_node> elements */
	domInstance_node_Array &instanceNodes = node.getInstance_node_array();
	for (size_t i=0; i<instanceNodes.getCount(); ++i) {
		domNode *node = daeSafeCast<domNode>(instanceNodes[i]->getUrl().getElement());
		if (!node)
			SLog(EError, "Could not find a referenced node!");
		loadNode(transform, os, *node);
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

int colladaMain(int argc, char **argv) {
	if (argc < 3) {
		cout << "Syntax: mtsimport <DAE source file URL> <XML destination file>" << endl;
		return -1;
	}

#if defined(__LINUX__)
	std::string path = std::string("file://") +
		FileResolver::getInstance()->resolveAbsolute(argv[1]);
#else
	std::string path = argv[1];
#endif

	DAE *dae = new DAE();
	SLog(EInfo, "Loading \"%s\" ..", path.c_str());
	if (dae->load(path.c_str()) != DAE_OK) 
		SLog(EError, "Could not load \"%s\"!", argv[1]);

	std::ofstream os(argv[2]);
	if (os.fail())
		SLog(EError, "Could not write to \"%s\"!", argv[2]);

	domCOLLADA *document = dae->getDom(path.c_str());
	domVisual_scene *visualScene = daeSafeCast<domVisual_scene>
		(document->getDescendant("visual_scene"));
	if (!visualScene)
		SLog(EError, "Could not find a visual_scene!");

	/* Configure the GLU tesselator */
	tess = gluNewTess();
	if (!tess) 
		SLog(EError, "Could not allocate a GLU tesselator!");

#ifndef WIN32
	int status = mkdir("textures", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status != 0 && errno != EEXIST)
		SLog(EError, "Could not create the directory \"textures\"");
	status = mkdir("meshes", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status != 0 && errno != EEXIST)
		SLog(EError, "Could not create the directory \"meshes\"");
#else
	int status = CreateDirectory("textures", NULL);
	if (status == 0 && GetLastError() != ERROR_ALREADY_EXISTS)
		SLog(EError, "Could not create the directory \"textures\"");
	status = CreateDirectory("meshes", NULL);
	if (status == 0 && GetLastError() != ERROR_ALREADY_EXISTS)
		SLog(EError, "Could not create the directory \"meshes\"");
#endif

	gluTessCallback(tess, GLU_TESS_VERTEX, reinterpret_cast<GLvoid (__stdcall *)()>(&tessVertex));
	gluTessCallback(tess, GLU_TESS_BEGIN, reinterpret_cast<GLvoid (__stdcall *)()>(&tessBegin));
	gluTessCallback(tess, GLU_TESS_END, reinterpret_cast<GLvoid (__stdcall *)()>(&tessEnd));
	gluTessCallback(tess, GLU_TESS_ERROR, reinterpret_cast<GLvoid (__stdcall *)()>(&tessError));
	gluTessCallback(tess, GLU_TESS_COMBINE, reinterpret_cast<GLvoid (__stdcall *)()>(&tessCombine));
	gluTessCallback(tess, GLU_TESS_EDGE_FLAG, reinterpret_cast<GLvoid (__stdcall *)()>(&tessEdgeFlag));

	domNode_Array &nodes = visualScene->getNode_array();
	os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl << endl;
	os << "<!--" << endl;
	os << "\tAutomatically converted from COLLADA" << endl;
	os << "-->" << endl << endl;
	os << "<scene>" << endl;
	os << "\t<integrator type=\"direct\"/>" << endl << endl;

	SLog(EInfo, "Converting to \"%s\" ..", argv[2]);

	StringMap idToTexture, fileToId;
	domLibrary_images_Array &libraryImages = document->getLibrary_images_array();
	for (size_t i=0; i<libraryImages.getCount(); ++i) {
		domImage_Array &images = libraryImages[i]->getImage_array();
		for (size_t j=0; j<images.getCount(); ++j) 
			loadImage(os, *images.get(j), idToTexture, fileToId);
	}

	domLibrary_materials_Array &libraryMaterials = document->getLibrary_materials_array();
	for (size_t i=0; i<libraryMaterials.getCount(); ++i) {
		domMaterial_Array &materials = libraryMaterials[i]->getMaterial_array();
		for (size_t j=0; j<materials.getCount(); ++j) 
			loadMaterial(os, *materials.get(j), idToTexture);
	}

	for (size_t i=0; i<nodes.getCount(); ++i) 
		loadNode(Transform(), os, *nodes[i]);

	os << "</scene>" << endl;

	gluDeleteTess(tess);
	delete dae;
	return 0;
}

int ubi_main(int argc, char **argv) {
	int retval;
	/* Initialize the core framework */
	Class::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();

	Thread::getThread()->getLogger()->setLogLevel(EInfo);

#ifdef WIN32
	char lpFilename[1024];
	if (GetModuleFileNameA(NULL,
		lpFilename, sizeof(lpFilename))) {
		FileResolver *resolver = FileResolver::getInstance();
		resolver->addPathFromFile(lpFilename);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}
#else
	char exePath[PATH_MAX];
	if (getcwd(exePath, PATH_MAX)) {
		FileResolver *resolver = FileResolver::getInstance();
		resolver->addPathFromFile(exePath);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}
#endif

	try {
		/* An OpenGL context may be required for the GLU tesselator */
		ref<Session> session = Session::create();
		ref<Device> device = Device::create(session);
		ref<Renderer> renderer = Renderer::create(session);

		session->init();
		device->init();
		renderer->init(device);

		device->makeCurrent(renderer);

		retval = colladaMain(argc, argv);

		if (retval != -1)
			cout << "Finished conversion" << endl;

		renderer->shutdown();
		device->shutdown();
		session->shutdown();
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
		retval = -1;
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type!" << endl;
		retval = -1;
	}

	/* Shutdown the core framework */
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	Class::staticShutdown();
	
	return retval;
}

#if !defined(__OSX__)
int main(int argc, char **argv) {
	return ubi_main(argc, argv);
}
#endif


