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

// Mitsuba's "Assert" macro conflicts with Xerces' XSerializeEngine::Assert(...).
// This becomes a problem when using a PCH which contains mitsuba/core/logger.h
#if defined(Assert)
# undef Assert
#endif
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
#include <xercesc/framework/MemBufFormatTarget.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <set>

XERCES_CPP_NAMESPACE_USE

#include "converter.h"

class ImporterDOMErrorHandler : public DOMErrorHandler {
public:
    inline ImporterDOMErrorHandler() { }

    bool handleError(const DOMError& domError) {
        SLog(EWarn, "%s (line %i, char %i): %s",
            XMLString::transcode(domError.getLocation()->getURI()),
            domError.getLocation()->getLineNumber(),
            domError.getLocation()->getColumnNumber(),
            XMLString::transcode(domError.getMessage()));

        if (domError.getSeverity() != DOMError::DOM_SEVERITY_WARNING)
            SLog(EError, "Encountered a critical DOM error -- giving up!");

        return true;
    }
};

void createNodeMap(DOMNode *node, std::map<std::string, DOMNode *> &nodes) {
    if (node) {
        char *nodeName = XMLString::transcode(node->getNodeName());
        if (strcmp(nodeName, "ref") == 0) {
            XMLString::release(&nodeName);
            return;
        }
        XMLString::release(&nodeName);
        if (node->getNodeType() == DOMNode::ELEMENT_NODE && node->hasAttributes()) {
            DOMNamedNodeMap *attributes = node->getAttributes();
            XMLCh *idString = XMLString::transcode("id");
            DOMAttr *attribute = (DOMAttr *) attributes->getNamedItem(idString);
            XMLString::release(&idString);
            if (attribute != NULL) {
                char *value = XMLString::transcode(attribute->getValue());
                nodes[value] = node;
                XMLString::release(&value);
            }
        }
        for (DOMNode *child = node->getFirstChild(); child != 0; child=child->getNextSibling())
            createNodeMap(child, nodes);
    }
}

void cleanup(DOMNode *node) {
    DOMNode *child = node->getFirstChild();
    while (child) {
        DOMNode *next = child->getNextSibling();
        if (child->getNodeType() == DOMNode::TEXT_NODE)
            node->removeChild(child);
        else
            cleanup(child);
        child = next;
    }
}

void GeometryConverter::convert(const fs::path &inputFile,
    const fs::path &outputDirectory,
    const fs::path &sceneName,
    const fs::path &adjustmentFile) {

    fs::path textureDirectory = "textures";
    fs::path meshesDirectory = "meshes";
    fs::path outputFile = sceneName;

    this->m_outputDirectory = outputDirectory;

    if (!outputDirectory.empty()) {
        textureDirectory = outputDirectory / "textures";
        meshesDirectory = outputDirectory / "meshes";
        outputFile = outputDirectory / sceneName;
    }

    if (m_packGeometry) {
        m_geometryFileName = outputDirectory / sceneName;
        m_geometryFileName.replace_extension(".serialized");
        m_geometryFile = new FileStream(m_geometryFileName, FileStream::ETruncReadWrite);
        m_geometryFile->setByteOrder(Stream::ELittleEndian);
    }

    if (!fs::exists(textureDirectory)) {
        SLog(EInfo, "Creating directory \"%s\" ..", textureDirectory.string().c_str());
        fs::create_directory(textureDirectory);
    }

    if (!fs::exists(meshesDirectory) && !m_packGeometry) {
        SLog(EInfo, "Creating directory \"%s\" ..", meshesDirectory.string().c_str());
        fs::create_directory(meshesDirectory);
    }

    std::ostringstream os;
    SLog(EInfo, "Beginning conversion ..");

    std::string extension = boost::to_lower_copy(inputFile.extension().string());

    if (extension == ".dae" || extension == ".zae") {
        convertCollada(inputFile, os, textureDirectory, meshesDirectory);
    } else if (extension == ".obj") {
        convertOBJ(inputFile, os, textureDirectory, meshesDirectory);
    } else {
        SLog(EError, "Unknown input format (must end in either .DAE, .ZAE or .OBJ)");
    }

    if (!adjustmentFile.empty()) {
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
        if (adj == NULL)
            SLog(EError, "Could not parse adjustments file!");

        std::map<std::string, DOMNode *> nodeMap;
        createNodeMap(doc, nodeMap);

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

        for (DOMNode *child = adjRoot->getFirstChild(); child != 0; child=child->getNextSibling()) {
            if (child->getNodeType() == DOMNode::ELEMENT_NODE) {
                char *temp = XMLString::transcode(child->getNodeName());
                std::string nodeName(temp);
                XMLString::release(&temp);
                std::string id;
                if (child->getNodeType() == DOMNode::ELEMENT_NODE && child->hasAttributes()) {
                    DOMNamedNodeMap *attributes = child->getAttributes();
                    XMLCh *idString = XMLString::transcode("id");
                    DOMAttr *attribute = (DOMAttr *) attributes->getNamedItem(idString);
                    XMLString::release(&idString);
                    if (attribute) {
                        char *value = XMLString::transcode(attribute->getValue());
                        id = value;
                        XMLString::release(&value);
                    }
                }
                if (id != "" && nodeMap.find(id) != nodeMap.end()) {
                    DOMNode *node = nodeMap[id], *parent = node->getParentNode();
                    if (nodeName == "append") {
                        for (DOMNode *child2 = child->getFirstChild(); child2 != 0; child2=child2->getNextSibling())
                            node->insertBefore(doc->importNode(child2, true), NULL);
                    } else if (nodeName == "prepend") {
                        for (DOMNode *child2 = child->getFirstChild(); child2 != 0; child2=child2->getNextSibling())
                            node->insertBefore(doc->importNode(child2, true), node->getFirstChild());
                    } else if (nodeName == "remove") {
                        parent->removeChild(node);
                    } else if (parent == insertBeforeNode->getParentNode()) {
                        parent->removeChild(node);
                        docRoot->insertBefore(doc->importNode(child, true), insertBeforeNode);
                    } else {
                        parent->replaceChild(doc->importNode(child, true), node);
                    }
                } else {
                    if (nodeName == "append" || nodeName == "prepend" || nodeName == "remove")
                        SLog(EError, "Adjustments file: Found an append/prepend/remove element, "
                            " which refers to a nonexistant node");
                    docRoot->insertBefore(doc->importNode(child, true), insertBeforeNode);
                }
            }
        }
        cleanup(doc);

        DOMLSSerializer *serializer = impl->createLSSerializer();
        DOMConfiguration *serConf = serializer->getDomConfig();
        serConf->setParameter(XMLUni::fgDOMErrorHandler, &errorHandler);
        if (serConf->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serConf->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
        DOMLSOutput *output = impl->createLSOutput();
        MemBufFormatTarget *target = new MemBufFormatTarget();
        output->setByteStream(target);
        serializer->write(doc, output);
        const XMLByte *content = target->getRawBuffer();
        std::ostringstream oss;

        /* Turn leading spaces into tabs */
        bool newline = true;
        int numSpaces = 0;
        for (size_t i=0; i<target->getLen(); ++i) {
            char data = content[i];
            switch (data) {
                case ' ':
                    if (newline)
                        numSpaces++;
                    else
                        oss << data;
                    break;
                case '\r':
                case '\n':
                    oss << data;
                    newline = true;
                    numSpaces = 0;
                    break;
                default:
                    if (newline) {
                        for (int i=0; i<numSpaces/2; ++i)
                            oss << '\t';
                    }
                    oss << data;
                    newline = false;
                    numSpaces = 0;
                    break;
            }
        }

        fs::ofstream os(outputFile);
        os << oss.str();
        os.close();
        delete output;
        delete target;
        delete wrapper;
        delete memBufIS;
        delete serializer;
        parser->release();
    } else {
        fs::ofstream ofile(outputFile);
        if (ofile.fail())
            SLog(EError, "Could not write to \"%s\"!", outputFile.string().c_str());
        ofile << os.str();
        ofile.close();
    }
    if (m_geometryFile) {
        for (size_t i=0; i<m_geometryDict.size(); ++i)
            m_geometryFile->writeULong((uint64_t) m_geometryDict[i]);
        m_geometryFile->writeUInt((uint32_t) m_geometryDict.size());
        m_geometryFile->close();
    }

    m_filename = outputFile;
}

