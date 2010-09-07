#include <mitsuba/core/platform.h>
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/render/util.h>
#include <mitsuba/render/shandler.h>

MTS_NAMESPACE_BEGIN

ref<Scene> Utility::loadScene(const std::string &filename) {
	/* Prepare for parsing scene descriptions */
	FileResolver *resolver = FileResolver::getInstance();
	SAXParser* parser = new SAXParser();
	std::string schemaPath = resolver->resolveAbsolute("schema/scene.xsd");

	/* Check against the 'scene.xsd' XML Schema */
	parser->setDoSchema(true);
	parser->setValidationSchemaFullChecking(true);
	parser->setValidationScheme(SAXParser::Val_Always);
	parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

	std::map<std::string, std::string> parameters;
	SceneHandler *handler = new SceneHandler(parameters);
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
