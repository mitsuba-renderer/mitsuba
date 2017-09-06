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

#define BOOST_FILESYSTEM_NO_LIB
#define BOOST_SYSTEM_NO_LIB

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/version.h>
#include <mitsuba/render/scene.h>
#include <boost/filesystem/fstream.hpp>
#include "converter.h"

std::set<std::string> availableTextures;

void referenceTexture(GeometryConverter *cvt, std::ostream &os, const std::string &parameter,
        const std::string &indent, const fs::path &textureDir, std::string filename) {
    /* Prevent Linux/OSX fs::path handling issues for OBJ files created on Windows */
    for (size_t i=0; i<filename.length(); ++i) {
        if (filename[i] == '\\')
            filename[i] = '/';
    }

    fs::path path = fs::path(filename);
    fs::path targetPath = textureDir / path.filename();
    std::string textureName = targetPath.filename().string();

    if (availableTextures.find(textureName) == availableTextures.end()) {
        SLog(EInfo, "Copying texture \"%s\" ..", textureName.c_str());

        if (!fs::exists(targetPath)) {
            if (!fs::exists(path)) {
                ref<FileResolver> fRes = Thread::getThread()->getFileResolver();
                path = fRes->resolve(path.filename());
                if (!fs::exists(path)) {
                    SLog(EWarn, "Found neither \"%s\" nor \"%s\"!", filename.c_str(), path.string().c_str());
                    path = cvt->locateResource(path.filename());
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
                output->close();
                input->close();
            }
        }

        os << indent << "<texture name=\"" << parameter << "\" type=\"bitmap\" id=\"" << textureName << "\">" << endl
            << indent << "\t<string name=\"filename\" value=\"textures/" << textureName << "\"/>" << endl
            << indent << "</texture>" << endl;
        availableTextures.insert(textureName);
    } else {
        os << indent << "<ref name=\"" << parameter << "\" id=\"" << textureName << "\"/>" << endl;
    }
}

void addMaterial(GeometryConverter *cvt, std::ostream &os, const std::string &mtlName,
        const fs::path &texturesDir, const Spectrum &diffuseValue,
        const std::string &diffuseMap, const std::string maskMap) {
    if (mtlName == "")
        return;
    SLog(EInfo, "Importing material \"%s\" ..", mtlName.c_str());
    std::string indent = "";

    if (maskMap != "") {
        os << "\t<bsdf id=\"" << mtlName << "_material\" type=\"mask\">" << endl;
        referenceTexture(cvt, os, "opacity", "\t\t", texturesDir, maskMap);
        os << "\t\t<bsdf type=\"diffuse\">" << endl;
    } else {
        os << "\t<bsdf id=\"" << mtlName << "_material\" type=\"diffuse\">" << endl;
    }

    if (diffuseMap == "") {
        Float r, g, b;
        diffuseValue.toLinearRGB(r, g, b);
        os << indent << "\t\t<rgb name=\"reflectance\" value=\""
            << r << " " << g << " " << b << "\"/>" << endl;
    } else {
        referenceTexture(cvt, os, "reflectance", "\t\t", texturesDir, diffuseMap);
    }

    os << indent << "\t</bsdf>" << endl << endl;

    if (maskMap != "")
        os << "\t</bsdf>" << endl;
}

void parseMaterials(GeometryConverter *cvt, std::ostream &os, const fs::path &texturesDir,
        const fs::path &mtlFileName, std::set<std::string> &mtlList) {
    SLog(EInfo, "Loading OBJ materials from \"%s\" ..", mtlFileName.string().c_str());
    fs::ifstream is(mtlFileName);
    if (is.bad() || is.fail())
        SLog(EError, "Unexpected I/O error while accessing material file '%s'!",
            mtlFileName.string().c_str());
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
        SLog(EError, "Could not open OBJ file '%s'!", inputFile.string().c_str());

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
            fRes->prependPath(fs::absolute(fRes->resolve(inputFile)).parent_path());
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
    objProps.setString("filename", inputFile.string());

    ref<Shape> rootShape = static_cast<Shape *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Shape), objProps));
    SAssert(rootShape->isCompound());

    int ctr = 0;
    while (true) {
        TriMesh *mesh = static_cast<TriMesh *>(rootShape->getElement(ctr++));
        if (!mesh)
            break;
        os << "\t<shape id=\"" << mesh->getName() << "_mesh\" type=\"serialized\">" << endl;

        if (!m_geometryFile) {
            std::string filename = mesh->getName() + std::string(".serialized");
            SLog(EInfo, "Saving \"%s\"", filename.c_str());
            ref<FileStream> stream = new FileStream(meshesDirectory / filename, FileStream::ETruncReadWrite);
            stream->setByteOrder(Stream::ELittleEndian);
            mesh->serialize(stream);
            stream->close();
            os << "\t\t<string name=\"filename\" value=\"meshes/" << filename.c_str() << "\"/>" << endl;
        } else {
            m_geometryDict.push_back((uint64_t) m_geometryFile->getPos());
            SLog(EInfo, "Saving mesh \"%s\" ..", mesh->getName().c_str());
            mesh->serialize(m_geometryFile);
            os << "\t\t<string name=\"filename\" value=\"" << m_geometryFileName.filename().string() << "\"/>" << endl;
            os << "\t\t<integer name=\"shapeIndex\" value=\"" << (m_geometryDict.size()-1) << "\"/>" << endl;
        }

        if (mesh->getBSDF() != NULL &&
                mtlList.find(mesh->getBSDF()->getID()) != mtlList.end()) {
            const std::string &matID = mesh->getBSDF()->getID();
            os << "\t\t<ref name=\"bsdf\" id=\"" << matID << "_material\"/>" << endl;
        } else {
            os << "\t\t<bsdf type=\"diffuse\"/>" << endl;
        }
        os << "\t</shape>" << endl << endl;
    }
    os << "</scene>" << endl;
}
