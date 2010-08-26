#define BOOST_FILESYSTEM_NO_LIB 
#define BOOST_SYSTEM_NO_LIB 

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/trimesh.h>
#include <boost/filesystem.hpp>
#include <fstream>
#include "converter.h"

std::string copyTexture(GeometryConverter *cvt, const std::string &textureDir, std::string filename) {
	SLog(EInfo, "Copying texture \"%s\" ..", filename.c_str());

#if defined(WIN32)
	for (size_t i=0; i<filename.length(); ++i)
		if (filename[i] == '/')
			filename[i] = '\\';
#else
	for (size_t i=0; i<filename.length(); ++i)
		if (filename[i] == '\\')
			filename[i] = '/';
#endif
	
	boost::filesystem::path path = boost::filesystem::path(filename, boost::filesystem::native);
	std::string targetPath = textureDir + path.leaf();

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
		output->close();
		input->close();
	}

	return targetPath;
}

void addMaterial(GeometryConverter *cvt, std::ostream &os, const std::string &mtlName,
		const std::string &texturesDir, const Spectrum &diffuseValue, 
		const std::string &diffuseMap, const std::string maskMap) {
	if (mtlName == "") 
		return;
	SLog(EInfo, "Copying material \"%s\" ..", mtlName.c_str());
	std::string indent = "";

	if (maskMap != "") {
		indent = "\t";
		os << "\t<bsdf id=\"" << mtlName << "\" type=\"mask\">" << endl;
		os << "\t\t<texture name=\"opacity\" type=\"ldrtexture\">" << endl;
		os << "\t\t\t<string name=\"filename\" value=\"" << copyTexture(cvt, texturesDir, maskMap) << "\"/>" << endl;
		os << "\t\t</texture>" << endl;
		os << "\t\t<bsdf type=\"lambertian\">" << endl;
	} else {
		os << "\t<bsdf id=\"" << mtlName << "\" type=\"lambertian\">" << endl;
	}

	if (diffuseMap == "") {
		Float r, g, b;
		diffuseValue.toLinearRGB(r, g, b);
		os << indent << "\t\t<rgb name=\"reflectance\" value=\"" 
			<< r << " " << g << " " << b << "\"/>" << endl;
	} else {
		os << indent << "\t\t<texture name=\"reflectance\" type=\"ldrtexture\">" << endl
		   << indent << "\t\t\t<string name=\"filename\" value=\"" << copyTexture(cvt, texturesDir, diffuseMap) << "\"/>" << endl
		   << indent << "\t\t</texture>" << endl;
	}

	os << indent << "\t</bsdf>" << endl << endl;

	if (maskMap != "") 
		os << "\t</bsdf>" << endl;
}

void parseMaterials(GeometryConverter *cvt, std::ostream &os, const std::string &texturesDir, 
		const std::string &mtlFileName) {
	SLog(EInfo, "Loading OBJ materials from \"%s\" ..", mtlFileName.c_str());
	std::ifstream is(mtlFileName.c_str());
	if (is.bad() || is.fail())
		SLog(EError, "Unexpected I/O error while accessing material file '%s'!", mtlFileName.c_str());
	std::string buf, line;
	std::string mtlName;
	Spectrum diffuse(0.0f);
	std::string diffuseMap, maskMap;

	while (is >> buf) {
		if (buf == "newmtl") {
			addMaterial(cvt, os, mtlName, texturesDir, diffuse, diffuseMap, maskMap);
			std::getline(is, line);
			mtlName = trim(line.substr(1, line.length()-1));
			diffuse = Spectrum(0.0f);
			diffuseMap = "";
			maskMap = "";
		} else if (buf == "Kd") {
			Float r, g, b;
			is >> r >> g >> b;
			if (cvt->m_srgb)
				diffuse.fromSRGB(r, g, b);
			else
				diffuse.fromLinearRGB(r, g, b);
		} else if (buf == "map_Kd") {
			std::getline(is, line);
			diffuseMap = trim(line.substr(1, line.length()-1));
		} else if (buf == "map_d") {
			std::getline(is, line);
			maskMap = trim(line.substr(1, line.length()-1));
		} else {
			/* Ignore */
			std::getline(is, line);
		}
	}
	addMaterial(cvt, os, mtlName, texturesDir, diffuse, diffuseMap, maskMap);
}

void GeometryConverter::convertOBJ(const std::string &inputFile, 
	std::ostream &os,
	const std::string &textureDirectory,
	const std::string &meshesDirectory) {

	std::ifstream is(inputFile.c_str());
	if (is.bad() || is.fail())
		SLog(EError, "Could not open OBJ file '%s'!", inputFile.c_str());
	
	os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl << endl;
	os << "<!--" << endl << endl;
	os << "\tAutomatically converted from Wavefront OBJ" << endl << endl;
	os << "-->" << endl << endl;
	os << "<scene>" << endl;
	os << "\t<integrator id=\"integrator\" type=\"direct\"/>" << endl << endl;
	
	std::string buf, line;
	while (is >> buf) {
		if (buf == "mtllib") {
			std::getline(is, line);
			std::string mtlName = trim(line.substr(1, line.length()-1));
			ref<FileResolver> fRes = FileResolver::getInstance()->clone();
			fRes->addPathFromFile(fRes->resolveAbsolute(inputFile));
			std::string fullMtlName = fRes->resolve(mtlName);
			if (FileStream::exists(fullMtlName))
				parseMaterials(this, os, textureDirectory, fullMtlName);
			else
				SLog(EWarn, "Could not find referenced material library '%s'", mtlName.c_str());
		} else {
			/* Ignore */
			std::getline(is, line);
		}
	}

	Properties objProps("obj");
	objProps.setString("filename", inputFile);

	ref<Shape> rootShape = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(Shape::m_theClass, objProps));
	SAssert(rootShape->isCompound());

	int ctr = 0;
	while (true) {
		TriMesh *mesh = static_cast<TriMesh *>(rootShape->getElement(ctr++));
		if (!mesh)
			break;
		std::string filename = meshesDirectory + mesh->getName() + std::string(".serialized");
		SLog(EInfo, "Saving \"%s\"", filename.c_str());
		ref<FileStream> stream = new FileStream(filename, FileStream::ETruncReadWrite);
		stream->setByteOrder(Stream::ENetworkByteOrder);
		mesh->serialize(stream);
		stream->close();
		os << "\t<shape id=\"" << mesh->getName() << "\" type=\"serialized\">" << endl;
		os << "\t\t<string name=\"filename\" value=\"" << filename.c_str() << "\"/>" << endl;
		if (mesh->getBSDF() != NULL) { 
			const std::string &matID = mesh->getBSDF()->getName();
			os << "\t\t<ref name=\"bsdf\" id=\"" << matID << "\"/>" << endl;
		}
		os << "\t</shape>" << endl << endl;
	}
	os << "</scene>" << endl;
}
