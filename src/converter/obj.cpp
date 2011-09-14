/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#define BOOST_FILESYSTEM_NO_LIB 
#define BOOST_SYSTEM_NO_LIB 

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/version.h>
#include <mitsuba/render/scene.h>
#include <boost/filesystem/fstream.hpp>
#include "converter.h"

std::string copyTexture(GeometryConverter *cvt, const fs::path &textureDir, std::string filename) {
	SLog(EInfo, "Copying texture \"%s\" ..", filename.c_str());

	boost::filesystem::path path = boost::filesystem::path(filename);
	fs::path targetPath = textureDir / path.leaf();
	fs::path resolved = filename;

	if (!fs::exists(targetPath)) {
		ref<FileResolver> fRes = Thread::getThread()->getFileResolver();
		if (!fs::exists(resolved)) {
			resolved = fRes->resolve(path.leaf());
			if (!fs::exists(resolved)) {
				SLog(EWarn, "Found neither \"%s\" nor \"%s\"!", filename.c_str(), resolved.file_string().c_str());
				resolved = cvt->locateResource(path.leaf());
				targetPath = targetPath.parent_path() / resolved.leaf();
				if (resolved.empty())
					SLog(EError, "Unable to locate a resource -- aborting conversion.");
			}
		}	

		if (fs::complete(resolved) != fs::complete(targetPath)) {
			ref<FileStream> input = new FileStream(resolved, FileStream::EReadOnly);
			ref<FileStream> output = new FileStream(targetPath, FileStream::ETruncReadWrite);
			input->copyTo(output);
			output->close();
			input->close();
		}
	}

#if BOOST_FILESYSTEM_VERSION == 3
	return targetPath.leaf().string();
#else
	return targetPath.leaf();
#endif
}

void addMaterial(GeometryConverter *cvt, std::ostream &os, const std::string &mtlName,
		const fs::path &texturesDir, const Spectrum &diffuseValue, 
		const std::string &diffuseMap, const std::string maskMap) {
	if (mtlName == "") 
		return;
	SLog(EInfo, "Copying material \"%s\" ..", mtlName.c_str());
	std::string indent = "";

	if (maskMap != "") {
		indent = "\t";
		os << "\t<bsdf id=\"" << mtlName << "\" type=\"mask\">" << endl;
		os << "\t\t<texture name=\"opacity\" type=\"bitmap\">" << endl;
		os << "\t\t\t<string name=\"filename\" value=\"textures/" << copyTexture(cvt, texturesDir, maskMap) << "\"/>" << endl;
		os << "\t\t</texture>" << endl;
		os << "\t\t<bsdf type=\"diffuse\">" << endl;
	} else {
		os << "\t<bsdf id=\"" << mtlName << "\" type=\"diffuse\">" << endl;
	}

	if (diffuseMap == "") {
		Float r, g, b;
		diffuseValue.toLinearRGB(r, g, b);
		os << indent << "\t\t<rgb name=\"reflectance\" value=\"" 
			<< r << " " << g << " " << b << "\"/>" << endl;
	} else {
		os << indent << "\t\t<texture name=\"reflectance\" type=\"bitmap\">" << endl
		   << indent << "\t\t\t<string name=\"filename\" value=\"textures/" << copyTexture(cvt, texturesDir, diffuseMap) << "\"/>" << endl
		   << indent << "\t\t</texture>" << endl;
	}

	os << indent << "\t</bsdf>" << endl << endl;

	if (maskMap != "") 
		os << "\t</bsdf>" << endl;
}

void parseMaterials(GeometryConverter *cvt, std::ostream &os, const fs::path &texturesDir, 
		const fs::path &mtlFileName, std::set<std::string> &mtlList) {
	SLog(EInfo, "Loading OBJ materials from \"%s\" ..", mtlFileName.file_string().c_str());
	fs::ifstream is(mtlFileName);
	if (is.bad() || is.fail())
		SLog(EError, "Unexpected I/O error while accessing material file '%s'!", 
			mtlFileName.file_string().c_str());
	std::string buf, line;
	std::string mtlName;
	Spectrum diffuse(0.0f);
	std::string diffuseMap, maskMap;

	while (is >> buf) {
		if (buf == "newmtl") {
			mtlList.insert(mtlName);
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

void GeometryConverter::convertOBJ(const fs::path &inputFile, 
	std::ostream &os,
	const fs::path &textureDirectory,
	const fs::path &meshesDirectory) {

	fs::ifstream is(inputFile);
	if (is.bad() || is.fail())
		SLog(EError, "Could not open OBJ file '%s'!", inputFile.file_string().c_str());

	os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl << endl;
	os << "<!--" << endl << endl;
	os << "\tAutomatically converted from Wavefront OBJ" << endl << endl;
	os << "-->" << endl << endl;
	os << "<scene version=\"" << MTS_VERSION << "\">" << endl;
	os << "\t<integrator id=\"integrator\" type=\"direct\"/>" << endl << endl;

	std::string buf, line;
	std::set<std::string> mtlList;
	while (is >> buf) {
		if (buf == "mtllib" && m_importMaterials) {
			std::getline(is, line);
			std::string mtlName = trim(line.substr(1, line.length()-1));
			ref<FileResolver> fRes = Thread::getThread()->getFileResolver()->clone();
			fRes->addPath(fs::complete(fRes->resolve(inputFile)).parent_path());
			fs::path fullMtlName = fRes->resolve(mtlName);
			if (fs::exists(fullMtlName))
				parseMaterials(this, os, textureDirectory, fullMtlName, mtlList);
			else
				SLog(EWarn, "Could not find referenced material library '%s'", mtlName.c_str());
		} else {
			/* Ignore */
			std::getline(is, line);
		}
	}

	Properties objProps("obj");
	objProps.setString("filename", inputFile.file_string());

	ref<Shape> rootShape = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Shape), objProps));
	SAssert(rootShape->isCompound());

	int ctr = 0;
	while (true) {
		TriMesh *mesh = static_cast<TriMesh *>(rootShape->getElement(ctr++));
		if (!mesh)
			break;
		os << "\t<shape id=\"" << mesh->getName() << "\" type=\"serialized\">" << endl;

		if (!m_geometryFile) {
			std::string filename = mesh->getName() + std::string(".serialized");
			SLog(EInfo, "Saving \"%s\"", filename.c_str());
			ref<FileStream> stream = new FileStream(meshesDirectory / filename, FileStream::ETruncReadWrite);
			stream->setByteOrder(Stream::ELittleEndian);
			mesh->serialize(stream);
			stream->close();
			os << "\t\t<string name=\"filename\" value=\"meshes/" << filename.c_str() << "\"/>" << endl;
		} else {
			m_geometryDict.push_back((uint32_t) m_geometryFile->getPos());
			SLog(EInfo, "Saving mesh \"%s\"", mesh->getName().c_str());
			mesh->serialize(m_geometryFile);
			os << "\t\t<string name=\"filename\" value=\"" << m_geometryFileName.filename() << "\"/>" << endl;
			os << "\t\t<integer name=\"shapeIndex\" value=\"" << (m_geometryDict.size()-1) << "\"/>" << endl;
		}

		if (mesh->getBSDF() != NULL && 
				mtlList.find(mesh->getBSDF()->getName()) != mtlList.end()) { 
			const std::string &matID = mesh->getBSDF()->getName();
			os << "\t\t<ref name=\"bsdf\" id=\"" << matID << "\"/>" << endl;
		} else {
			os << "\t\t<bsdf type=\"diffuse\"/>" << endl;
		}
		os << "\t</shape>" << endl << endl;
	}
	os << "</scene>" << endl;
}
