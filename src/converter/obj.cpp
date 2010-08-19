#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/trimesh.h>
#include "converter.h"

void GeometryConverter::convertOBJ(const std::string &inputFile, 
	std::ostream &os,
	const std::string &textureDirectory,
	const std::string &meshesDirectory) {

	Properties objProps("obj");
	objProps.setString("filename", inputFile);

	ref<Shape> rootShape = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(Shape::m_theClass, objProps));
	SAssert(rootShape->isCompound());
	os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl << endl;
	os << "<!--" << endl << endl;
	os << "\tAutomatically converted from Wavefront OBJ" << endl << endl;
	os << "-->" << endl << endl;
	os << "<scene>" << endl;
	os << "\t<integrator id=\"integrator\" type=\"direct\"/>" << endl << endl;

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
