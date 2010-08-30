#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/framework/Wrapper4InputSource.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fresolver.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <set>

XERCES_CPP_NAMESPACE_USE

#include "converter.h"

class ImporterDOMErrorHandler : public DOMErrorHandler {
public:
	inline ImporterDOMErrorHandler() { }

	bool handleError(const DOMError& domError) {
		ELogLevel logLevel;

		if (domError.getSeverity() == DOMError::DOM_SEVERITY_WARNING)
			logLevel = EWarn;
		else
			logLevel = EError;

		SLog(logLevel, "%s (line %i, char %i): %s",
			XMLString::transcode(domError.getLocation()->getURI()),
			domError.getLocation()->getLineNumber(),
			domError.getLocation()->getColumnNumber(),
			XMLString::transcode(domError.getMessage()));
		return true;
	}
};

void findRemovals(DOMNode *node, std::set<std::string> &removals) {
	if (node) {
		char *nodeName = XMLString::transcode(node->getNodeName());
		if (strcmp(nodeName, "ref") == 0) {
			XMLString::release(&nodeName);
			return;
		}
		XMLString::release(&nodeName);
		if (node->getNodeType() == DOMNode::ELEMENT_NODE && node->hasAttributes()) {
			DOMNamedNodeMap *attributes = node->getAttributes();
			for (size_t i=0; i<attributes->getLength(); ++i) {
				DOMAttr *attribute = (DOMAttr*) attributes->item(i);
				char *name = XMLString::transcode(attribute->getName());
				char *value = XMLString::transcode(attribute->getValue());

				if (strcmp(name, "id") == 0)
					removals.insert(value);

				XMLString::release(&name);
				XMLString::release(&value);
			}
		}
		for (DOMNode *child = node->getFirstChild(); child != 0; child=child->getNextSibling())
			findRemovals(child, removals);
	}
}

bool cleanupPass(DOMNode *node, const std::set<std::string> &removals) {
	if (node) {
		char *nodeName = XMLString::transcode(node->getNodeName());
		if (strcmp(nodeName, "ref") == 0) {
			XMLString::release(&nodeName);
			return false;
		}

		if (node->getNodeType() == DOMNode::ELEMENT_NODE && node->hasAttributes()) {
			DOMNamedNodeMap *attributes = node->getAttributes();
			for (size_t i=0; i<attributes->getLength(); ++i) {
				DOMAttr *attribute = (DOMAttr*) attributes->item(i);
				char *name = XMLString::transcode(attribute->getName());
				char *value = XMLString::transcode(attribute->getValue());

				if (strcmp(name, "id") == 0 && removals.find(value) != removals.end()) {
					XMLString::release(&name);
					XMLString::release(&value);
					return true; /* Remove this node */
				}

				XMLString::release(&name);
				XMLString::release(&value);
			}
			XMLString::release(&nodeName);
		} else if (node->getNodeType() == DOMNode::TEXT_NODE) {
			return true;
		}
		DOMNode *child = node->getFirstChild();
		while (child) {
			DOMNode *next = child->getNextSibling();
			bool doRemove = cleanupPass(child, removals);
			if (doRemove)
				node->removeChild(child);
			child = next;
		}
	}
	return false;
}

void GeometryConverter::convert(const std::string &inputFile, 
	const std::string &outputDirectory, 
	const std::string &sceneName,
	const std::string &adjustmentFile) {

	std::string textureDirectory = "textures";
	std::string meshesDirectory = "meshes";
	std::string outputFile = sceneName;

	SLog(EInfo, "Creating directories ..");
#if !defined(WIN32)
	if (outputDirectory != "") {
		textureDirectory = outputDirectory + "/textures";
		meshesDirectory = outputDirectory + "/meshes";
		outputFile = outputDirectory + std::string("/") + sceneName;
	}

	int status = mkdir(textureDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status != 0 && errno != EEXIST)
		SLog(EError, "Could not create the directory \"%s\"", textureDirectory.c_str());
	status = mkdir(meshesDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status != 0 && errno != EEXIST)
		SLog(EError, "Could not create the directory \"%s\"", meshesDirectory.c_str());
	
	textureDirectory += "/";
	meshesDirectory += "/";
#else
	if (outputDirectory != "") {
		textureDirectory = outputDirectory + "\\textures";
		meshesDirectory = outputDirectory + "\\meshes";
		outputFile = outputDirectory + std::string("\\") + sceneName;
	}

	int status = CreateDirectory(textureDirectory.c_str(), NULL);
	if (status == 0 && GetLastError() != ERROR_ALREADY_EXISTS)
		SLog(EError, "Could not create the directory \"%s\"", textureDirectory.c_str());
	status = CreateDirectory(meshesDirectory.c_str(), NULL);
	if (status == 0 && GetLastError() != ERROR_ALREADY_EXISTS)
		SLog(EError, "Could not create the directory \"%s\"", meshesDirectory.c_str());

	textureDirectory += "\\";
	meshesDirectory += "\\";
#endif

	std::ostringstream os;

	SLog(EInfo, "Beginning conversion ..");
	if (endsWith(toLowerCase(inputFile), ".dae") || endsWith(toLowerCase(inputFile), ".zae")) {
		convertCollada(inputFile, os, textureDirectory, meshesDirectory);
	} else if (endsWith(toLowerCase(inputFile), ".obj")) {
		convertOBJ(inputFile, os, textureDirectory, meshesDirectory);
	} else {
		SLog(EError, "Unknown input format (must end in either .DAE, .ZAE or .OBJ)");
	}

	if (adjustmentFile != "") {
		SLog(EInfo, "Applying adjustments ..");
		static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
		DOMImplementationLS *impl = DOMImplementationRegistry::getDOMImplementation(gLS);

		DOMLSParser *parser = impl->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
		DOMConfiguration *conf(parser->getDomConfig());
		ImporterDOMErrorHandler errorHandler;
		conf->setParameter(XMLUni::fgDOMErrorHandler, &errorHandler);

		std::string xmlString = os.str();
		MemBufInputSource* memBufIS = new MemBufInputSource((const XMLByte*) xmlString.c_str(), 
			xmlString.length(), "bufID", false);
		Wrapper4InputSource *wrapper = new Wrapper4InputSource(memBufIS, false);
		XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *doc = parser->parse(wrapper);
		XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *adj = parser->parseURI(adjustmentFile.c_str());
	
		std::set<std::string> removals, emptyList;
		cleanupPass(adj, emptyList);
		findRemovals(adj, removals);
		cleanupPass(doc, removals);

		DOMElement *docRoot = doc->getDocumentElement();
		DOMElement *adjRoot = adj->getDocumentElement();

		DOMNode *insertBeforeNode = NULL;
		for (DOMNode *child = docRoot->getFirstChild(); child != 0; child=child->getNextSibling()) {
			char *name = XMLString::transcode(child->getNodeName());
			if (strcmp(name, "shape") == 0) {
				insertBeforeNode = child;
				break;
			}
			XMLString::release(&name);
		}

		if (insertBeforeNode == NULL) {
			/* No shape node found, use the camera node instead */
			for (DOMNode *child = docRoot->getFirstChild(); child != 0; child=child->getNextSibling()) {
				char *name = XMLString::transcode(child->getNodeName());
				if (strcmp(name, "camera") == 0) {
					insertBeforeNode = child;
					break;
				}
				XMLString::release(&name);
			}
			SAssertEx(insertBeforeNode != NULL, "Internal error while applying adjustments: cannot find shape/camera node");
		}

		for (DOMNode *child = adjRoot->getFirstChild(); child != 0; child=child->getNextSibling())
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				docRoot->insertBefore(doc->importNode(child, true), insertBeforeNode);

		DOMLSSerializer *serializer = impl->createLSSerializer();
		DOMConfiguration *serConf(serializer->getDomConfig());
		DOMLSOutput *output = impl->createLSOutput();
		serConf->setParameter(XMLUni::fgDOMErrorHandler, &errorHandler);
		serConf->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
		XMLFormatTarget *target = new LocalFileFormatTarget(outputFile.c_str());
		output->setByteStream(target);
		serializer->write(doc, output);
		delete output;
		delete target;
		delete wrapper;
		delete memBufIS;
		delete serializer;
		parser->release();
	} else {
		std::ofstream ofile(outputFile.c_str());
		if (ofile.fail())
			SLog(EError, "Could not write to \"%s\"!", outputFile.c_str());
		ofile << os.str();
		ofile.close();
	}
	m_filename = outputFile;
}

