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
#include <xercesc/parsers/SAXParser.hpp>
#include "glwidget.h"
#include "sceneloader.h"
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <boost/algorithm/string.hpp>

XERCES_CPP_NAMESPACE_USE

SceneLoader::SceneLoader(FileResolver *resolver, const fs::path &filename,
        const fs::path &destFile, const std::map<std::string, std::string, SimpleStringOrdering> &parameters)
    : Thread("load"), m_resolver(resolver), m_filename(fromFsPath(filename)), m_destFile(destFile),
      m_parameters(parameters) {
    m_wait = new WaitFlag();
    m_versionError = false;
}

SceneLoader::~SceneLoader() {
}

void SceneLoader::run() {
    Thread::getThread()->setFileResolver(m_resolver);
    SAXParser* parser = new SAXParser();
    QFileInfo fileInfo(m_filename);
    QString suffix = fileInfo.suffix().toLower();

    SceneHandler *handler = new SceneHandler(m_parameters);
    m_result = new SceneContext();
    try {
        QSettings settings;
        m_result->srgb = settings.value("preview_sRGB", true).toBool();
        m_result->gamma = (Float) settings.value("preview_gamma", 2.2).toDouble();
        m_result->reinhardKey = (Float) settings.value("preview_reinhardKey", 0.18).toDouble();
        m_result->reinhardBurn = (Float) settings.value("preview_reinhardBurn", -10.0).toDouble();
        m_result->exposure = (Float) settings.value("preview_exposure", 0).toDouble();
        m_result->shadowMapResolution = settings.value("preview_shadowMapResolution", 256).toInt();
        m_result->clamping = (Float) settings.value("preview_clamping", 0.1f).toDouble();
        m_result->previewMethod = (EPreviewMethod) settings.value("preview_method", EOpenGL).toInt();
        if (m_result->previewMethod != EOpenGL && m_result->previewMethod != EDisabled)
            m_result->previewMethod = EOpenGL;
        m_result->toneMappingMethod = (EToneMappingMethod) settings.value("preview_toneMappingMethod", EGamma).toInt();
        m_result->diffuseSources = settings.value("preview_diffuseSources", true).toBool();
        m_result->diffuseReceivers = settings.value("preview_diffuseReceivers", false).toBool();
        m_result->currentLayer = 0;

        if (suffix == "exr" || suffix == "png"  || suffix == "jpg" || suffix == "jpeg" ||
            suffix == "hdr" || suffix == "rgbe" || suffix == "pfm" || suffix == "ppm") {
            /* This is an image, not a scene */
            ref<FileStream> fs = new FileStream(toFsPath(m_filename), FileStream::EReadOnly);
            ref<Bitmap> bitmap = new Bitmap(Bitmap::EAuto, fs, (suffix == "exr") ? "*" : "");

            if (suffix == "exr" && bitmap->getLayers().size() > 1) {
                std::map<std::string, Bitmap*> bitmaps = bitmap->split();
                for (std::map<std::string, Bitmap *>::iterator it = bitmaps.begin(); it != bitmaps.end(); ++it) {
                    std::string name = it->first;
                    Bitmap *layer = it->second;
                    layer->incRef();
                    m_result->layers.push_back(std::make_pair(name, layer));
                }
                bitmap = m_result->layers[m_result->currentLayer].second;
            }

            bitmap = bitmap->convert(Bitmap::ERGBA, Bitmap::EFloat32);

            m_result->mode = ERender;
            m_result->framebuffer = bitmap;
            m_result->fileName = m_filename;
            m_result->shortName = fileInfo.fileName();
            m_result->scrollOffset = Vector2i(0, 0);
            m_result->pathLength = 2;
        } else {
            fs::path schemaPath = m_resolver->resolveAbsolute("data/schema/scene.xsd");

            /* Check against the 'scene.xsd' XML Schema */
            parser->setDoSchema(true);
            parser->setValidationSchemaFullChecking(true);
            parser->setValidationScheme(SAXParser::Val_Always);
            parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

            /* Set the SAX handler */
            parser->setDoNamespaces(true);
            parser->setDocumentHandler(handler);
            parser->setErrorHandler(handler);

            fs::path
                filename = toFsPath(m_filename),
                filePath = fs::absolute(filename).parent_path(),
                baseName = filename.stem();

            SLog(EInfo, "Parsing scene description from \"%s\" ..", filename.string().c_str());

            if (!fs::exists(filename))
                SLog(EError, "Unable to load scene \"%s\": file not found!",
                    filename.string().c_str());

            try {
                parser->parse(filename.c_str());
            } catch (const VersionException &ex) {
                m_versionError = true;
                m_version = ex.getVersion();
                throw;
            }

            ref<Scene> scene = handler->getScene();

            scene->setSourceFile(filename);
            scene->setDestinationFile(m_destFile.empty() ? (filePath / baseName) : m_destFile);
            scene->initialize();

            if (scene->getIntegrator() == NULL)
                SLog(EError, "Unable to load scene: no integrator found!");
            if (scene->getSensor() == NULL)
                SLog(EError, "Unable to load scene: no sensor found!");
            if (scene->getSensor()->getFilm() == NULL)
                SLog(EError, "Unable to load scene: no film found!");
            if (scene->getEmitters().size() == 0)
                SLog(EError, "Unable to load scene: no light sources found!");
            Vector2i size = scene->getFilm()->getCropSize();
            Sensor *sensor = scene->getSensor();

            /* Also generate a DOM representation for the Qt-based GUI */
            QFile file(m_filename);
            if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
                Log(EError, "Unable to open the file \"%s\"", filename.string().c_str());
            QString errorMsg;
            int line, column;
            if (!m_result->doc.setContent(&file, &errorMsg, &line, &column))
                SLog(EError, "Unable to parse file: error %s at line %i, colum %i",
                    qPrintable(errorMsg), line, column);

            m_result->scene = scene;
            m_result->sceneResID = Scheduler::getInstance()->registerResource(scene);
            m_result->renderJob = NULL;
            m_result->movementScale = scene->getBSphere().radius / 2000.0f;
            m_result->mode = EPreview;
            m_result->framebuffer = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, size);
            m_result->framebuffer->clear();
            m_result->fileName = m_filename;
            m_result->shortName = fileInfo.fileName();
            if (sensor->getClass()->derivesFrom(MTS_CLASS(PerspectiveCamera))) {
                m_result->up = static_cast<PerspectiveCamera *>(sensor)->getWorldTransform(
                    sensor->getShutterOpen() + 0.5f * sensor->getShutterOpenTime())(Vector(0, 1, 0));
            } else {
                m_result->up = Vector(0.0f);
            }
            m_result->scrollOffset = Vector2i(0, 0);
            m_result->originalSize = m_result->scene->getFilm()->getCropSize();
            m_result->pathLength = m_result->detectPathLength();
            m_result->showKDTree = false;
            m_result->shownKDTreeLevel = 0;
        }
    } catch (const std::exception &e) {
        m_error = e.what();
        delete m_result;
        m_result = NULL;
    } catch (...) {
        cout << "Don't know what" << endl;
        m_error = "An unknown type of error occurred!";
        delete m_result;
        m_result = NULL;
    }
    m_wait->set(true);
    delete parser;
    delete handler;
}

