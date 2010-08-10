#include <xercesc/parsers/SAXParser.hpp>
#include "glwidget.h"
#include "sceneloader.h"
#include "../mitsuba/shandler.h"

SceneLoader::SceneLoader(FileResolver *resolver, const std::string &filename) 
	: Thread("load"), m_resolver(resolver), m_filename(filename) {
	m_wait = new WaitFlag();
}

SceneLoader::~SceneLoader() {
}

void SceneLoader::run() {
	FileResolver::setInstance(m_resolver);
	SAXParser* parser = new SAXParser();
	std::string lowerCase = m_filename;
	for(size_t i=0; i<m_filename.size();++i)
		lowerCase[i] = std::tolower(m_filename[i]);

	SceneHandler *handler = new SceneHandler();
	m_result = new SceneContext();
	try {
		QSettings settings("mitsuba-renderer.org", "qtgui");
		m_result->srgb = settings.value("preview_sRGB", true).toBool();
		m_result->gamma = (Float) settings.value("preview_gamma", 2.2).toDouble();
		m_result->reinhardKey = (Float) settings.value("preview_reinhardKey", 0.18).toDouble();
		m_result->reinhardBurn = (Float) settings.value("preview_reinhardBurn", 0.0).toDouble();
		m_result->exposure = (Float) settings.value("preview_exposure", 0).toDouble();
		m_result->shadowMapResolution = settings.value("preview_shadowMapResolution", 256).toInt();
		m_result->clamping = (Float) settings.value("preview_clamping", 0.1f).toDouble();
		m_result->previewMethod = (EPreviewMethod) settings.value("preview_previewMethod", EOpenGL).toInt();
		m_result->toneMappingMethod = (EToneMappingMethod) settings.value("preview_toneMappingMethod", EGamma).toInt();

		if (endsWith(lowerCase, ".exr")) {
			/* This is an image, not a scene */
			ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
			ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, fs);

			m_result->mode = ERender;
			m_result->framebuffer = bitmap;
			m_result->fileName = QString(m_filename.c_str());
			m_result->shortName = QFileInfo(m_filename.c_str()).fileName();
			m_result->pathLength = 2;
		} else {
			std::string schemaPath = m_resolver->resolveAbsolute("schema/scene.xsd");

			/* Check against the 'scene.xsd' XML Schema */
			parser->setDoSchema(true);
			parser->setValidationSchemaFullChecking(true);
			parser->setValidationScheme(SAXParser::Val_Always);
			parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

			/* Set the SAX handler */
			parser->setDoNamespaces(true);
			parser->setDocumentHandler(handler);
			parser->setErrorHandler(handler);

			std::string filenameWithoutExtension = m_resolver->resolveDest( 
				FileResolver::getFilenameWithoutExtension(m_filename));

			SLog(EInfo, "Parsing scene description from \"%s\" ..", m_filename.c_str());
			parser->parse(m_filename.c_str());
			ref<Scene> scene = handler->getScene();

			scene->setSourceFile(m_filename.c_str());
			scene->setDestinationFile(filenameWithoutExtension);
			scene->initialize();

			if (scene->getIntegrator() == NULL)
				SLog(EError, "The scene contains no integrator! Aborting..");
			if (scene->getCamera() == NULL)
				SLog(EError, "The scene contains no camera! Aborting..");
			if (scene->getCamera()->getFilm() == NULL)
				SLog(EError, "The scene contains no film! Aborting..");
			if (scene->getLuminaires().size() == 0)
				SLog(EError, "The scene contains no light sources! Aborting..");
		
			Vector2i size = scene->getFilm()->getSize();
			Camera *camera = scene->getCamera();

			m_result->scene = scene;
			m_result->sceneResID = Scheduler::getInstance()->registerResource(scene);
			m_result->renderJob = NULL;
			m_result->movementScale = scene->getBSphere().radius / 2000.0f;
			m_result->mode = EPreview;
			m_result->framebuffer = new Bitmap(size.x, size.y, 128);
			m_result->framebuffer->clear();
			m_result->fileName = QString(m_filename.c_str());
			m_result->shortName = QFileInfo(m_filename.c_str()).fileName();
			m_result->up = camera->getInverseViewTransform()(Vector(0, 1, 0));
			m_result->scrollOffset = Vector2i(0, 0);
			m_result->pathLength = m_result->detectPathLength();
		}
	} catch (const std::exception &e) {
		m_error = e.what();
		delete m_result;
		m_result = NULL;
	} catch (...) {
		m_error = "An unknown type of error occurred!";
		delete m_result;
		m_result = NULL;
	}
	m_wait->set(true);
	delete parser;
	delete handler;
}

