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

#include <mitsuba/core/platform.h>
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scenehandler.h>

MTS_NAMESPACE_BEGIN

ref<Scene> Utility::loadScene(const std::string &filename,
		const ParameterMap &params) {
	/* Prepare for parsing scene descriptions */
	FileResolver *resolver = Thread::getThread()->getFileResolver();
	SAXParser* parser = new SAXParser();
	fs::path schemaPath = resolver->resolveAbsolute("schema/scene.xsd");

	/* Check against the 'scene.xsd' XML Schema */
	parser->setDoSchema(true);
	parser->setValidationSchemaFullChecking(true);
	parser->setValidationScheme(SAXParser::Val_Always);
	parser->setExternalNoNamespaceSchemaLocation(schemaPath.file_string().c_str());
	parser->setCalculateSrcOfs(true);

	SceneHandler *handler = new SceneHandler(parser, params);
	parser->setDoNamespaces(true);
	parser->setDocumentHandler(handler);
	parser->setErrorHandler(handler);
		
	parser->parse(filename.c_str());
	ref<Scene> scene = handler->getScene();

	delete parser;
	delete handler;

	return scene;
}

MTS_IMPLEMENT_CLASS(Utility, true, Object)
MTS_NAMESPACE_END
