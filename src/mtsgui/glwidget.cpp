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

#include <mitsuba/mitsuba.h>
#include <QtGui/QtGui>
#include "glwidget.h"
#include "preview.h"
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/qmc.h>
#include <mitsuba/hw/font.h>
#include <boost/tuple/tuple.hpp>

GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent), m_context(NULL) {
    setAutoBufferSwap(false);
    setFocusPolicy(Qt::StrongFocus);
    m_clock = new Timer();
    m_wheelTimer = new Timer();
    m_animationTimer = new Timer();
    m_statusTimer = new Timer();
    m_movementTimer = new QTimer(parent);
    m_movementTimer->setInterval(20);
    m_movementTimer->setSingleShot(false);
    m_redrawTimer = new QTimer(parent);
    m_redrawTimer->setInterval(200);
    connect(m_movementTimer, SIGNAL(timeout()), this, SLOT(timerImpulse()));
    connect(m_redrawTimer, SIGNAL(timeout()), this, SLOT(updateGL()));
    m_renderer = Renderer::create(NULL);
    m_device = new QtDevice(this);
    m_preview = new PreviewThread(m_device, m_renderer);
    connect(m_preview, SIGNAL(caughtException(const QString &)),
        this, SLOT(onException(const QString &)), Qt::QueuedConnection);
    connect(m_preview, SIGNAL(statusMessage(const QString &)),
        this, SIGNAL(statusMessage(const QString &)), Qt::QueuedConnection);
    m_invertMouse = false;
    m_navigationMode = EFlythrough;
    m_ignoreMouseEvent = QPoint(0, 0);
    m_didSetCursor = false;
#if defined(MTS_GUI_SOFTWARE_FALLBACK)
    m_softwareFallback = true;
#else
    m_softwareFallback = false;
#endif
    m_ignoreResizeEvents = false;
    m_ignoreScrollEvents = false;
    m_animation = false;
    m_cropping = false;
    setAcceptDrops(true);
}

GLWidget::~GLWidget() {
    shutdown();
}

void GLWidget::shutdown() {
    if (m_preview)
        m_preview->quit();
}

void GLWidget::onException(const QString &what) {
    QString errorString("A critical exception occurred in the preview rendering thread. "
            "Please make sure that you are using the most recent graphics drivers. "
            "Mitsuba will now switch to a slow software fallback mode, which only "
            "supports the rendering preview but no tonemapping and no real-time "
            "preview/navigation. The encountered error was:\n\n%1");
    QMessageBox::critical(this, tr("Critical exception"),
        errorString.arg(what), QMessageBox::Ok);
    m_softwareFallback = true;
}

void GLWidget::setSourceFromResource(GPUProgram *program,
    GPUProgram::EType type, const QString &resourceName) {
    QFile file(resourceName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        SLog(EError, "Could not open the shader resource: %s",
            resourceName.toLocal8Bit().constData());
        return;
    }

    QTextStream in(&file);
    const QString source = in.readAll();
    program->setSource(type, source.toStdString());
}

GPUProgram* GLWidget::createGPUProgram(const std::string &name,
        const QString &vertexResource, const QString &fragmentResource) {
    GPUProgram *prog = m_renderer->createGPUProgram(name);
    if (prog == NULL) {
        SLog(EError, "Could not create the GPU program");
        return NULL;
    }
    setSourceFromResource(prog, GPUProgram::EVertexProgram,   vertexResource);
    setSourceFromResource(prog, GPUProgram::EFragmentProgram, fragmentResource);
    return prog;
}

void GLWidget::initializeGL() {
    /* Load the Mitsuba logo into a texture */
    QResource res("/resources/logo.png");
    SAssert(res.isValid());
    ref<Stream> mStream = new MemoryStream(res.size());
    mStream->write(res.data(), res.size());
    mStream->seek(0);
    ref<Bitmap> bitmap = new Bitmap(Bitmap::EPNG, mStream);
    m_logoSize = Vector2(bitmap->getWidth(), bitmap->getHeight());
    m_device->init();
    m_renderer->init(m_device);
    m_logoTexture = m_renderer->createGPUTexture("Logo", bitmap);
    m_logoTexture->setFilterType(GPUTexture::ENearest);
    m_logoTexture->setMipMapped(false);
    m_font = new Font(Font::EBitstreamVeraMono14);

    std::vector<std::string> missingExtensions;

    if (!m_renderer->getCapabilities()->isSupported(
        RendererCapabilities::EShadingLanguage))
        missingExtensions.push_back("GLSL");
    if (!m_renderer->getCapabilities()->isSupported(
        RendererCapabilities::ERenderToTexture))
        missingExtensions.push_back("Render-To-Texture");
    if (!m_renderer->getCapabilities()->isSupported(
        RendererCapabilities::EFloatingPointTextures))
        missingExtensions.push_back("Floating point textures");
    if (!m_renderer->getCapabilities()->isSupported(
        RendererCapabilities::EFloatingPointBuffer))
        missingExtensions.push_back("Floating point render buffers");
    if (!m_renderer->getCapabilities()->isSupported(
        RendererCapabilities::EVertexBufferObjects))
        missingExtensions.push_back("Vertex buffer objects");
    if (!m_renderer->getCapabilities()->isSupported(
        RendererCapabilities::EBufferBlit))
        missingExtensions.push_back("Fast blitting");

    if (!missingExtensions.empty() || m_softwareFallback) {
#if !defined(MTS_GUI_SOFTWARE_FALLBACK)
        /* Show a warning message unless the fallback mode
           was explicitly requested */
        std::ostringstream oss;
        oss << "You machine is missing the following required "
            "OpenGL capabilities: ";
        for (size_t i=0; i<missingExtensions.size(); ++i) {
            oss << missingExtensions[i];
            if (i+1 < missingExtensions.size())
                oss << ", ";
        }
#if MTS_SSE
        oss << ". Please make sure that you are using the most "
            << "recent graphics drivers.\n\nMitsuba will now switch "
            << "to a software fallback mode, which supports "
            << "the rendering preview but no real-time preview/navigation.";
#else
        oss << ". Please make sure that you are using the most "
            << "recent graphics drivers.\n\nMitsuba will now switch "
            << "to a slow software fallback mode, which only supports "
            << "the rendering preview but no tonemapping and no "
            << "real-time preview/navigation.";
#endif
        m_errorString = QString(oss.str().c_str());
        m_softwareFallback = true;
#endif
        // Don't redraw as often, since this is now quite costly
        m_redrawTimer->setInterval(1000);
    } else {
        m_gammaTonemap        = createGPUProgram("Tonemapper [Gamma]",
            ":/shaders/gamma.vert",        ":/shaders/gamma.frag");
        m_reinhardTonemap     = createGPUProgram("Tonemapper [Reinhard et al. 2002]",
            ":/shaders/reinhard.vert",     ":/shaders/reinhard.frag");
        m_luminanceProgram    = createGPUProgram("Log-luminance program",
            ":/shaders/logluminance.vert", ":/shaders/logluminance.frag");
        m_downsamplingProgram = createGPUProgram("Downsampling program",
            ":/shaders/downsampling.vert", ":/shaders/downsampling.frag");

        if (!m_preview->isRunning()) {
#if defined(__WINDOWS__)
            wglMakeCurrent(NULL, NULL);
#endif
            m_preview->start();
            m_preview->waitUntilStarted();
#if defined(__WINDOWS__)
            makeCurrent();
#endif
        }

        m_gammaTonemap->init();
        m_reinhardTonemap->init();
        m_downsamplingProgram->init();
        m_luminanceProgram->init();
    }
    m_logoTexture->init();
    m_font->init(m_renderer);
    m_redrawTimer->start();
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}

void GLWidget::setScene(SceneContext *context) {
    if (m_movementTimer->isActive())
        m_movementTimer->stop();
    m_context = context;

    if (context && context->layers.size() > 0) {
        m_statusMessage =
            formatString("Showing layer \"%s\" (%i/%i); use '[' and ']' to switch.",
                m_context->layers[m_context->currentLayer].first.c_str(),
                m_context->currentLayer+1, m_context->layers.size());
        m_statusTimer->reset();
    }

    if (context && context->scene == NULL)
        context = NULL;

    m_preview->setSceneContext(context, true, false);
    m_framebufferChanged = true;
    m_mouseDrag = m_animation = m_cropping = false;
    m_leftKeyDown = m_rightKeyDown = m_upKeyDown = m_downKeyDown = false;
    m_aabb.reset();
    setCursor(Qt::ArrowCursor);
    updateGeometry();
    updateScrollBars();
    updateGL();
}

void GLWidget::refreshScene() {
    m_framebufferChanged = true;
    resetPreview();
    updateGeometry();
    updateScrollBars();
    if (m_context->mode == EPreview)
        resumePreview();
}

void GLWidget::setPathLength(int length) {
    if (m_context->pathLength != length) {
        m_context->pathLength = length;
        resetPreview();
    }
}

void GLWidget::setShadowMapResolution(int res) {
    if (res != m_context->shadowMapResolution) {
        m_context->shadowMapResolution = res;
        resetPreview();
    }
}

void GLWidget::setDiffuseSources(bool value) {
    if (value != m_context->diffuseSources) {
        m_context->diffuseSources = value;
        resetPreview();
    }
}

void GLWidget::setDiffuseReceivers(bool value) {
    if (value != m_context->diffuseReceivers) {
        m_context->diffuseReceivers = value;
        resetPreview();
    }
}

void GLWidget::setPreviewMethod(EPreviewMethod method) {
    if (method != m_context->previewMethod) {
        m_context->previewMethod = method;
        resetPreview();
    }
}

void GLWidget::setClamping(Float clamping) {
    if (clamping != m_context->clamping) {
        m_context->clamping = clamping;
        resetPreview();
    }
}

void GLWidget::setGamma(bool srgb, Float gamma) {
    if (srgb != m_context->srgb || gamma != m_context->gamma) {
        m_context->srgb = srgb;
        m_context->gamma = gamma;
        updateGL();
    }
}

void GLWidget::setToneMappingMethod(EToneMappingMethod method) {
    if (method != m_context->toneMappingMethod) {
        m_context->toneMappingMethod = method;
        updateGL();
    }
}

void GLWidget::setExposure(Float exposure) {
    if (exposure != m_context->exposure) {
        m_context->exposure = exposure;
        updateGL();
    }
}

void GLWidget::setReinhardKey(Float value) {
    if (value != m_context->reinhardKey) {
        m_context->reinhardKey = value;
        updateGL();
    }
}

void GLWidget::setReinhardBurn(Float value) {
    if (value != m_context->reinhardBurn) {
        m_context->reinhardBurn = value;
        updateGL();
    }
}

void GLWidget::downloadFramebuffer() {
    bool createdFramebuffer = false;

    if (!m_preview->isRunning()
        || m_context->previewMethod == EDisabled) {
        m_context->framebuffer->clear();
        return;
    }

    makeCurrent();
    if (m_framebuffer == NULL ||
        m_framebuffer->getBitmap() != m_context->framebuffer) {
        if (m_framebuffer)
            m_framebuffer->cleanup();
        m_framebuffer = m_renderer->createGPUTexture("Framebuffer",
            m_context->framebuffer);
        m_framebuffer->setMipMapped(false);
        m_framebuffer->setFilterType(GPUTexture::ENearest);
        createdFramebuffer = true;
    }

    PreviewQueueEntry entry = m_preview->acquireBuffer(1000);
    if (entry.buffer == NULL) {
        if (createdFramebuffer)
            m_framebuffer->init();
        return;
    }

    Point3i size = entry.buffer->getSize();
    m_renderer->finish();
    entry.buffer->download(m_context->framebuffer);

    // Scale by the number of developed VPLs
    float *targetData = m_context->framebuffer->getFloat32Data();
    float factor = 1.0f / entry.vplSampleOffset;

    for (size_t pos=0, total = size.x*size.y; pos<total; ++pos) {
        for (int i=0; i<3; ++i)
            targetData[4*pos+i] *= factor;
        targetData[4*pos+3] = 1.0f;
    }

    if (createdFramebuffer)
        m_framebuffer->init();
    else
        m_framebuffer->refresh();
    m_preview->releaseBuffer(entry);
}

QSize GLWidget::sizeHint() const {
    float ratio = windowHandle()->devicePixelRatio();
    QSize minimumSize(440, 170);
    if (m_context) {
        return QSize(std::max(minimumSize.width(), int(m_context->framebuffer->getWidth()/ratio)),
                     std::max(minimumSize.height(), int(m_context->framebuffer->getHeight()/ratio)));
    } else {
        return minimumSize;
    }
}

void GLWidget::focusOutEvent(QFocusEvent *event) {
    m_leftKeyDown = m_rightKeyDown = m_upKeyDown = m_downKeyDown = false;
    if (m_movementTimer->isActive())
        m_movementTimer->stop();
}

void GLWidget::timerImpulse() {
    if (!m_context || !m_context->scene || !m_preview->isRunning()) {
        m_movementTimer->stop();
        return;
    }
    if (m_animation) {
        Float x = std::min(m_animationTimer->getMilliseconds() / 500.0f, 1.0f);
        Float t = x*x*x*(x*(x*6-15)+10); // smootherstep by Ken Perlin

        Point origin = (1-t) * m_animationOrigin0 + t * m_animationOrigin1;
        Point target = (1-t) * m_animationTarget0 + t * m_animationTarget1;

        setWorldTransform(
            Transform::lookAt(origin, target, m_context->up));

        if (x == 1.0f)
            m_animation = false;
    }

    if (!(m_mouseDrag && m_navigationMode == EStandard) && !m_animation) {
        ProjectiveCamera *camera = getProjectiveCamera();

        Float delta = m_context->movementScale
            * m_clock->getMilliseconds();

        if (m_leftKeyDown)
            setWorldTransform(getWorldTransform() *
                Transform::translate(Vector(delta,0,0)));
        if (m_rightKeyDown)
            setWorldTransform(getWorldTransform() *
                Transform::translate(Vector(-delta,0,0)));
        if (m_downKeyDown) {
            setWorldTransform(getWorldTransform() *
                Transform::translate(Vector(0,0,-delta)));
            camera->setFocusDistance(camera->getFocusDistance() + delta);
        }
        if (m_upKeyDown) {
            setWorldTransform(getWorldTransform() *
                Transform::translate(Vector(0,0,delta)));
            camera->setFocusDistance(camera->getFocusDistance() - delta);
        }
    }

    m_clock->reset();

    if (m_context->renderJob) {
        m_context->renderJob->cancel();
        m_context->cancelled = true;
        m_context->cancelMode = EPreview;
        return;
    } else if (m_context->mode != EPreview) {
        m_context->mode = EPreview;
        // causes updateUI to be called in the main window
        emit stopRendering();
    }

    if (!(m_leftKeyDown || m_rightKeyDown || m_upKeyDown || m_downKeyDown
        || m_wheelTimer->getMilliseconds() < 200 || m_animation)) {
        if (m_movementTimer->isActive())
            m_movementTimer->stop();
    }

    resetPreview();
}

bool GLWidget::askReallyCancelRendering() {
    try {
        Float renderTime = m_context->renderJob->getRenderTime();

        if (renderTime < 10) /* Only ask for jobs that have been rendering for a bit */
            return true;
    } catch (...) {
        return true;
    }

    bool cancel = QMessageBox::question(this,
        "Really cancel rendering?", "Camera motion detected. Do you really want to cancel the rendering?",
        QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes;

    return cancel;
}

void GLWidget::resetPreview() {
    if (!m_context || !m_context->scene || !m_preview->isRunning())
        return;
    bool motion = m_leftKeyDown || m_rightKeyDown ||
        m_upKeyDown || m_downKeyDown || m_mouseDrag ||
        m_wheelTimer->getMilliseconds() < 200 || m_animation;
    m_preview->setSceneContext(m_context, false, motion);
    updateGL();
}

void GLWidget::startCrop(ECropType type) {
    if (!m_context || !m_context->scene || m_context->renderJob)
        return;
    m_cropType = type;
    m_cropping = true;
    m_cropStart = m_cropEnd = Point2i(0);
    setCursor(Qt::CrossCursor);
}

void GLWidget::keyPressEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_Escape) {
        if (m_cropping) {
            m_cropping = false;
            setCursor(Qt::ArrowCursor);
            updateGL();
        } else {
            emit quit();
        }
    } else if (event->key() == Qt::Key_R) {
        emit beginRendering();
    }
    else if (event->key() == Qt::Key_S) {
        emit stopRendering();
    }

    if (!m_context ||
        (event->isAutoRepeat() && event->key() != Qt::Key_BracketLeft &&
         event->key() != Qt::Key_BracketRight))
        return;

    switch (event->key()) {
        case Qt::Key_PageUp:   m_context->movementScale *= 2; break;
        case Qt::Key_PageDown: m_context->movementScale /= 2; break;
        case Qt::Key_Left:
            if (event->modifiers() & Qt::AltModifier) {
                emit switchTab(-1);
                return;
            }
            m_leftKeyDown = true; break;
        case Qt::Key_Right:
            if (event->modifiers() & Qt::AltModifier) {
                emit switchTab(1);
                return;
            }
            m_rightKeyDown = true; break;
        case Qt::Key_Up:
            m_upKeyDown = true; break;
        case Qt::Key_Down:
            m_downKeyDown = true; break;
        case Qt::Key_BracketLeft:
            if (m_context->layers.size() > 0) {
                m_context->currentLayer = std::max(0, m_context->currentLayer-1);
                m_context->framebuffer = m_context->layers[m_context->currentLayer].second->convert(Bitmap::ERGBA, Bitmap::EFloat32);
                m_statusMessage =
                    formatString("Showing layer \"%s\" (%i/%i); use '[' and ']' to switch.",
                        m_context->layers[m_context->currentLayer].first.c_str(),
                        m_context->currentLayer+1, m_context->layers.size());
                m_statusTimer->reset();
            } else if (m_context->showKDTree) {
                m_context->shownKDTreeLevel
                    = std::max(0, m_context->shownKDTreeLevel - 1);
                updateGL();
                return;
            }
            break;
        case Qt::Key_BracketRight:
            if (m_context->layers.size() > 0) {
                m_context->currentLayer = std::min((int) m_context->layers.size() - 1, m_context->currentLayer+1);
                m_context->framebuffer = m_context->layers[m_context->currentLayer].second->convert(Bitmap::ERGBA, Bitmap::EFloat32);
                m_statusMessage =
                    formatString("Showing layer \"%s\" (%i/%i); use '[' and ']' to switch.",
                        m_context->layers[m_context->currentLayer].first.c_str(),
                        m_context->currentLayer+1, m_context->layers.size());
                m_statusTimer->reset();
            } else if (m_context->showKDTree) {
                m_context->shownKDTreeLevel++;
                updateGL();
                return;
            }
            break;
        case Qt::Key_M: startCrop(ECropAndMagnify); break;
        case Qt::Key_C: startCrop(ECrop); break;
        case Qt::Key_A: {
            if (m_context->selectionMode == EScene) {
                m_context->selectionMode = ENothing;
                m_aabb.reset();
            } else if (m_context->scene) {
                m_context->selectionMode = EScene;
                m_aabb = m_context->scene->getKDTree()->getAABB();
            }
            m_context->selectedShape = NULL;
            emit selectionChanged();
        }
        // break stmt. intentionally missing
        case Qt::Key_F: {
            reveal(m_aabb);
        };
        break;
    }
    if (!m_movementTimer->isActive() && (m_leftKeyDown || m_rightKeyDown
            || m_upKeyDown || m_downKeyDown)) {

        if (m_context->renderJob && !askReallyCancelRendering()) {
            m_leftKeyDown = m_rightKeyDown = m_upKeyDown = m_downKeyDown = false;
            return;
        }

        m_clock->reset();
        m_movementTimer->start();
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event) {
    if (event->isAutoRepeat() || !m_context || !m_preview->isRunning())
        return;
    switch (event->key()) {
        case Qt::Key_Left:
            m_leftKeyDown = false; break;
        case Qt::Key_Right:
            m_rightKeyDown = false; break;
        case Qt::Key_Up:
            m_upKeyDown = false; break;
        case Qt::Key_Down:
            m_downKeyDown = false; break;
    }
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
    float ratio = windowHandle()->devicePixelRatio();
    QPoint rel = event->pos() - m_mousePos / ratio;
    if (!m_context || !m_context->scene
     || !m_mouseDrag || rel == QPoint(0,0) || m_animation)
        return;

    m_mousePos = event->pos();

    //  if (m_ignoreMouseEvent == rel) {
    if (m_ignoreMouseEvent != QPoint(0, 0)) {
        m_ignoreMouseEvent = QPoint(0, 0);
        return;
    }

    if (m_cropping) {
        Point2i offset = upperLeft();
        Vector2i maxCrop = m_context->scene->getFilm()->getCropSize();
        m_cropEnd.x = std::min(std::max(0, m_mousePos.x()-offset.x), maxCrop.x-1);
        m_cropEnd.y = std::min(std::max(0, m_mousePos.y()-offset.y), maxCrop.y-1);

        if (event->modifiers().testFlag(Qt::ShiftModifier)) {
            int dx = m_cropEnd.x-m_cropStart.x, dy = m_cropEnd.y-m_cropStart.y;

            if (std::abs(dx) > std::abs(dy))
                m_cropEnd.y = std::min(std::max(0, m_cropStart.y + std::abs(dx)*(dy<0 ? -1 : 1)), maxCrop.y-1);
            else
                m_cropEnd.x = std::min(std::max(0, m_cropStart.x + std::abs(dy)*(dx<0 ? -1 : 1)), maxCrop.x-1);
        }

        m_statusMessage =
            formatString("%s: %i x %i pixels",
                m_cropType == ECrop ? "Crop" : "Crop & Magnify",
                std::abs(m_cropEnd.x-m_cropStart.x), std::abs(m_cropEnd.y-m_cropStart.y));
        m_statusTimer->reset();

        updateGL();
        return;
    }

    PerspectiveCamera *camera = getPerspectiveCamera();
    if (!camera || !m_preview->isRunning())
        return;

    Transform invView = getWorldTransform();
    Point p = invView(Point(0,0,0));
    Vector d = invView(Vector(0,0,1));
    bool didMove = false;

    Float focusDistance = camera->getFocusDistance(),
        nearClip = camera->getNearClip(),
        farClip = camera->getFarClip();

    bool setCursor = true;
    if (m_context->renderJob) {
        if (!askReallyCancelRendering())
            return;
        setCursor = false;
    }
    if (focusDistance <= nearClip || focusDistance >= farClip) {
        focusDistance = autoFocus();
        camera->setFocusDistance(focusDistance);
    }

    Point target = p + d * focusDistance;
    Vector up = m_context->up;

    if (!m_didSetCursor && setCursor) {
        QApplication::setOverrideCursor(Qt::BlankCursor);
        m_didSetCursor = true;
    }

    if (event->buttons() & Qt::LeftButton) {
        if (m_navigationMode == EStandard) {
            Frame frame(up);
            Point2 coords = toSphericalCoordinates(frame.toLocal(normalize(p-target)));
            coords.y -=  0.001f * rel.x() * m_mouseSensitivity * (isRightHanded() ? 1 : -1);
            coords.x -=  0.001f * rel.y() * m_mouseSensitivity * (m_invertMouse ? -1.0f : 1.0f);
            p = target + focusDistance * frame.toWorld(sphericalDirection(coords.x, coords.y));

            if (coords.x < 0 || coords.x > M_PI)
                m_context->up *= -1;

            setWorldTransform(Transform::lookAt(p, target, m_context->up));
        } else {
            Float yaw = -.03f * rel.x() * m_mouseSensitivity;
            Float pitch = .03f * rel.y() * m_mouseSensitivity;
            if (m_invertMouse)
                pitch *= -1;

            Transform trafo = invView
                    * Transform::rotate(Vector(1,0,0), pitch)
                    * Transform::rotate(Vector(0,1,0), yaw);
            d = trafo(Vector(0,0,1));

            setWorldTransform(Transform::lookAt(p, p+d, up));
        }
        didMove = true;
    } else if (event->buttons() & Qt::MidButton) {
        setWorldTransform(invView *
            Transform::translate(Vector((Float) rel.x(), (Float) rel.y(), 0)
                * m_mouseSensitivity * .6f * m_context->movementScale));
        didMove = true;
    } else if (event->buttons() & Qt::RightButton) {
        if (event->modifiers() & Qt::ShiftModifier) {
            Float roll = rel.x() * m_mouseSensitivity * .02f;
            Float fovChange = rel.y() * m_mouseSensitivity * .03f;

            m_context->up = Transform::rotate(d, isRightHanded() ? -roll : roll)(up);
            setWorldTransform(Transform::lookAt(p, p+d, m_context->up));

            camera->setXFov(std::min(std::max((Float) 1.0f, camera->getXFov()
                + fovChange), (Float) 160.0f));
            m_statusMessage =
                formatString("Field of view: %.2f degrees", camera->getXFov());
            m_statusTimer->reset();
        } else {
            Float focusDistance = camera->getFocusDistance(),
                nearClip = camera->getNearClip(),
                farClip = camera->getFarClip();

            if (focusDistance <= nearClip || focusDistance >= farClip)
                focusDistance = autoFocus();

            Float oldFocusDistance = focusDistance;
            focusDistance = std::min(std::max(focusDistance * std::pow((Float) (1 - 2e-3f),
                    (Float) -rel.y() * m_mouseSensitivity * m_context->movementScale),
                    1.2f*nearClip), farClip/1.2f);

            camera->setFocusDistance(focusDistance);
            p = p + (oldFocusDistance - focusDistance) * d;

            setWorldTransform(Transform::lookAt(p, p+d, up));
        }
        didMove = true;
    }
    camera->configure();

    /* Re-center cursor as needed */
    QPoint global = mapToGlobal(m_mousePos);
    QDesktopWidget *desktop = QApplication::desktop();
    QRect geo = desktop->screenGeometry(global);

    if (global.x() < 50 || global.y() < 50 ||
        global.x() > desktop->width()-50 ||
        global.y() > desktop->height()-50) {
        QPoint target = geo.center();
        m_ignoreMouseEvent = target - global;
        QCursor::setPos(target);
    }

    if (didMove)
        timerImpulse();
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
    if (m_context == NULL || m_context->scene == NULL)
        return;

    float ratio = windowHandle()->devicePixelRatio();
    m_mousePos = event->pos()*ratio;

    m_initialMousePos = mapToGlobal(event->pos());
    m_mouseDrag = true;

    if (m_cropping) {
        Point2i offset = upperLeft();
        m_cropStart = m_cropEnd = Point2i(m_mousePos.x() - offset.x,
            m_mousePos.y() - offset.y);
        Vector2i maxCrop = m_context->scene->getFilm()->getCropSize();

        if (event->buttons() != Qt::LeftButton
            || m_cropStart.x < 0 || m_cropStart.y < 0 ||
               m_cropStart.x >= maxCrop.x ||
               m_cropStart.y >= maxCrop.y) {
            m_cropping = false;
            setCursor(Qt::ArrowCursor);
        } else {
            setCursor(Qt::SizeFDiagCursor);
        }
        return;
    }

    const PerspectiveCamera *camera = getPerspectiveCamera();
    if (!camera || !m_preview->isRunning())
        return;

    if (event->buttons() == Qt::LeftButton && m_navigationMode == EStandard) {
        Point2i offset = upperLeft();
        Point2 sample = Point2(m_mousePos.x() - offset.x,
                m_mousePos.y() - offset.y);
        Intersection its;
        Ray ray;

        camera->sampleRay(ray, sample, Point2(0.5f), 0.5f);

        if (m_context->scene->rayIntersect(ray, its)) {
            m_statusMessage =
                formatString("Selected shape \"%s\"", its.shape->getName().c_str());
            m_statusTimer->reset();
            m_context->selectedShape = its.instance ? its.instance : its.shape;
            AABB aabb(m_context->selectedShape->getAABB());
            bool newSelection = (m_aabb != aabb);
            m_context->selectionMode = EShape;
            if (newSelection) {
                m_aabb = aabb;
                emit selectionChanged();
            }
        } else if (m_aabb.isValid()) {
            m_aabb.reset();
            m_context->selectionMode = ENothing;
            m_context->selectedShape = NULL;
            emit selectionChanged();
        }
    }
}

void GLWidget::mouseDoubleClickEvent(QMouseEvent *event) {
    if (!m_preview->isRunning())
        return;
    if (m_navigationMode == EStandard && m_aabb.isValid()
        && event->buttons() & Qt::LeftButton)
        reveal(m_aabb);
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
    if (event->buttons() == 0) {
        if (m_cropping) {
            setCursor(Qt::ArrowCursor);
            if (m_cropEnd.x < m_cropStart.x)
                std::swap(m_cropEnd.x, m_cropStart.x);
            if (m_cropEnd.y < m_cropStart.y)
                std::swap(m_cropEnd.y, m_cropStart.y);

            int width  = m_cropEnd.x - m_cropStart.x,
                height = m_cropEnd.y - m_cropStart.y;

            if (width > 1 && height > 1)
                emit crop(m_cropType, m_cropStart.x,
                    m_cropStart.y, width, height);

            m_cropping = false;
        }
        m_mouseDrag = false;

        if (m_didSetCursor) {
            resetPreview();
            QApplication::restoreOverrideCursor();
            QCursor::setPos(m_initialMousePos);
            m_didSetCursor = false;
        }
    }
}

void GLWidget::wheelEvent(QWheelEvent *event) {
    QScrollBar *bar = event->orientation() == Qt::Vertical
        ? m_vScroll : m_hScroll;

    if (bar->isVisible()) {
        int oldStep = bar->singleStep();
        bar->setSingleStep(event->delta()/4);
#if defined(__OSX__)
        /* Support two-finger swipes */
        if (std::abs(event->delta()) < 8*15)
            bar->setSingleStep(event->delta());
#endif

        bar->triggerAction(event->delta() > 0 ?
            QAbstractSlider::SliderSingleStepSub
            : QAbstractSlider::SliderSingleStepAdd);
        bar->setSingleStep(oldStep);
    } else {
        if (!m_preview->isRunning() || m_context == NULL || m_context->scene == NULL || m_animation)
            return;
        PerspectiveCamera *camera = getPerspectiveCamera();
        if (!camera)
            return;

        if (m_context->renderJob && !askReallyCancelRendering())
            return;

        Float focusDistance = camera->getFocusDistance(),
            nearClip = camera->getNearClip(),
            farClip = camera->getFarClip();

        if (focusDistance <= nearClip || focusDistance >= farClip)
            focusDistance = autoFocus();

        Float oldFocusDistance = focusDistance;
        focusDistance = std::min(std::max(focusDistance * std::pow((Float) (1 - 1e-3f),
                (Float) event->delta()), 1.2f*nearClip), farClip/1.2f);

        camera->setFocusDistance(focusDistance);

        Transform invView = getWorldTransform();
        Point p  = invView(Point(0,0,0));
        Vector d = invView(Vector(0,0,1));
        Vector up = invView(Vector(0,1,0));

        p = p + (oldFocusDistance - focusDistance) * d;

        setWorldTransform(Transform::lookAt(p, p+d, up));

        m_wheelTimer->reset();
        if (!m_movementTimer->isActive())
            m_movementTimer->start();
        resetPreview();
    }
    event->accept();
}

Float GLWidget::autoFocus() const {
    if (m_context == NULL || m_context->scene == NULL)
        return std::numeric_limits<Float>::infinity();
    const Scene *scene = m_context->scene;
    const ProjectiveCamera *camera = getProjectiveCamera();
    if (!camera)
        return 0.0f;
    Float variance = 0.0625f; // (0.25f ^ 2)
    Float t, avgDistance = 0, weightSum = 0;
    Vector2i size = camera->getFilm()->getCropSize();
    ConstShapePtr ptr;
    Normal n;
    Ray ray;
    Point2 uv;

    for (size_t sampleIndex=0; sampleIndex<200; ++sampleIndex) {
        Point2 sample(
            radicalInverse(2, sampleIndex),
            radicalInverse(3, sampleIndex));
        camera->sampleRay(ray, Point2(sample.x * size.x, sample.y*size.y), Point2(0.5f), 0.5f);
        if (scene->rayIntersect(ray, t, ptr, n, uv)) {
            Float weight = math::fastexp(-0.5 / variance * (
                std::pow(sample.x - 0.5f, (Float) 2) +
                std::pow(sample.y - 0.5f, (Float) 2)));
            avgDistance += t * weight;
            weightSum += weight;
        }
    }

    if (weightSum == 0)
        return 0.5f * (camera->getNearClip() + camera->getFarClip());
    else
        return avgDistance / weightSum;
}

void GLWidget::paintGL() {
    m_renderer->setDepthTest(true);
    m_renderer->setDepthMask(true);
    if (m_context == NULL) {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        m_renderer->setBlendMode(Renderer::EBlendAlpha);
        m_renderer->blitTexture(m_logoTexture);
        m_renderer->setBlendMode(Renderer::EBlendNone);
    } else if (m_context != NULL) {
        Vector2i size = m_context->framebuffer->getSize();
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, m_device->getSize().x, m_device->getSize().y, 0, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.375f, 0.375f, 0.0f);

        PreviewQueueEntry entry;
        GPUTexture *buffer = NULL;
        Point2i upperLeft = this->upperLeft();
        Point2i lowerRight = upperLeft +
                m_context->framebuffer->getSize();

        float ratio = windowHandle()->devicePixelRatio();
        if (width() * ratio > size.x || height() * ratio > size.y) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            /* Draw a border to highlight the region occupied by the image */
            glColor4f(0.4f, 0.4f, 0.4f, 1.0f);
            glRecti(upperLeft.x-1, upperLeft.y-1,
                    lowerRight.x, lowerRight.y);
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if (m_context->mode == EPreview) {
            if (!m_preview->isRunning() || m_context->previewMethod == EDisabled) {
                /* No preview thread running - just show a grey screen */
                if (m_cropping && m_cropStart != m_cropEnd) {
                    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glEnable(GL_COLOR_LOGIC_OP);
                    glLogicOp(GL_XOR);
                    Point2i cropStart = m_cropStart + upperLeft;
                    Point2i cropEnd = m_cropEnd + upperLeft;
                    m_renderer->setDepthTest(false);
                    glRecti(cropStart.x, cropStart.y,
                            cropEnd.x, cropEnd.y);
                    glDisable(GL_COLOR_LOGIC_OP);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                }
                swapBuffers();
                return;
            }
            bool motion = m_leftKeyDown || m_rightKeyDown ||
                m_upKeyDown || m_downKeyDown || m_mouseDrag ||
                m_wheelTimer->getMilliseconds() < 200 || m_animation;
            entry = m_preview->acquireBuffer(motion ? -1 : 50);
            if (entry.buffer == NULL) {
                /* Unsuccessful at acquiring a buffer in the
                   alloted time - just leave and keep the current display */
                return;
            }
            size = Vector2i(entry.buffer->getSize().x, entry.buffer->getSize().y);
            buffer = entry.buffer;
        } else if (m_context->mode == ERender) {
            if (m_framebuffer == NULL ||
                m_framebuffer->getBitmap()->getWidth() != m_context->framebuffer->getWidth() ||
                m_framebuffer->getBitmap()->getHeight() != m_context->framebuffer->getHeight() ||
                (m_softwareFallback && m_framebuffer->getBitmap() != m_fallbackBitmap) ||
                (!m_softwareFallback && m_framebuffer->getBitmap() != m_context->framebuffer)) {
                if (m_framebuffer)
                    m_framebuffer->cleanup();
                if (!m_softwareFallback) {
                    m_framebuffer = m_renderer->createGPUTexture("Framebuffer",
                        m_context->framebuffer);
                } else {
                    m_fallbackBitmap = new Bitmap(Bitmap::ERGBA, Bitmap::EUInt8,
                        m_context->framebuffer->getSize());
                    m_fallbackBitmap->clear();
                    m_framebuffer = m_renderer->createGPUTexture("Framebuffer",
                        m_fallbackBitmap);
#if MTS_SSE
                    m_cpuTonemap = new TonemapCPU;
#endif
                }
                m_framebuffer->setMipMapped(false);
                m_framebuffer->setFilterType(GPUTexture::ENearest);
                m_framebuffer->init();
            }

            if (m_framebufferChanged) {
                if (m_softwareFallback) {
#if MTS_SSE
                    Float mult = 1.0;
                    if (m_context->mode == EPreview) {
                        mult /= entry.vplSampleOffset;
                    }
                    m_cpuTonemap->setLuminanceInfo(m_context->framebuffer,mult);
#else
                    /* Manually generate a gamma-corrected image
                       on the CPU (with gamma=2.2) - this will be slow! */
                    Bitmap *source = m_context->framebuffer;
                    float *sourceData = source->getFloat32Data();
                    uint8_t *targetData = (uint8_t *) m_fallbackBitmap->getData();
                    for (int y=0; y<source->getHeight(); ++y) {
                        for (int x=0; x<source->getWidth(); ++x) {
                            const float invGammaValue = 0.45455f;
                            *targetData++ = (uint8_t) std::max(std::min(
                                std::pow(*sourceData++, invGammaValue) * 255.0f, 255.0f), 0.0f);
                            *targetData++ = (uint8_t) std::max(std::min(
                                std::pow(*sourceData++, invGammaValue) * 255.0f, 255.0f), 0.0f);
                            *targetData++ = (uint8_t) std::max(std::min(
                                std::pow(*sourceData++, invGammaValue) * 255.0f, 255.0f), 0.0f);
                            ++sourceData;
                            *targetData++ = 255;
                        }
                    }
#endif
                }
                m_framebuffer->refresh();
                m_framebufferChanged = false;
            }

            size = Vector2i(m_framebuffer->getSize().x, m_framebuffer->getSize().y);
            buffer = m_framebuffer;
        } else {
            return;
        }
        bool hasDepth = m_context->mode == EPreview
            && m_context->previewMethod == EOpenGL;

        if (m_softwareFallback) {
#if MTS_SSE
            m_cpuTonemap->setSRGB(m_context->srgb);
            m_cpuTonemap->setInvGamma(static_cast<Float>(1) / m_context->gamma);

            if (m_context->toneMappingMethod == EGamma) {
                Float invWhitePoint = std::pow((Float) 2.0f,
                    m_context->exposure);
                if (m_context->mode == EPreview)
                    invWhitePoint /= entry.vplSampleOffset;

                m_cpuTonemap->setInvWhitePoint(invWhitePoint);
                m_cpuTonemap->gammaTonemap(m_context->framebuffer,
                    m_fallbackBitmap);
            } else if (m_context->toneMappingMethod == EReinhard) {
                Float mult = 1.0;
                if (m_context->mode == EPreview)
                    mult /= entry.vplSampleOffset;

                /* Getting the luminance info is rather expensive, avoid if the
                   multiplier has not changed */
                if (mult != m_cpuTonemap->multiplier()) {
                    m_cpuTonemap->setLuminanceInfo(m_context->framebuffer,mult);
                }
                Float burn = std::min((Float) 1, std::max((Float)1e-8f,
                    1-(m_context->reinhardBurn + 10) / 20.0f));
                Float scale = m_context->reinhardKey /
                              m_cpuTonemap->logAvgLuminance();
                Float Lwhite = m_cpuTonemap->maxLuminance() * scale;

                m_cpuTonemap->setScale(scale);
                m_cpuTonemap->setMultiplier(mult);
                m_cpuTonemap->setInvWhitePoint(1 / (Lwhite * (burn*burn)));
                m_cpuTonemap->reinhardTonemap(m_context->framebuffer,
                    m_fallbackBitmap);
            }
            m_framebuffer->refresh();
#endif
            buffer->bind();
            m_renderer->setColor(Spectrum(1.0f));
            m_renderer->blitTexture(buffer, false,
                !m_hScroll->isVisible(), !m_vScroll->isVisible(),
                -m_context->scrollOffset);
            buffer->unbind();
        } else if (m_context->toneMappingMethod == EGamma) {
            Float invWhitePoint = std::pow((Float) 2.0f, m_context->exposure);
            if (m_context->mode == EPreview)
                invWhitePoint /= entry.vplSampleOffset;

            if (hasDepth)
                buffer->bind(1, 1);
            m_gammaTonemap->bind();
            m_gammaTonemap->setParameter("colorSource", 0);
            m_gammaTonemap->setParameter("depthSource", 1);
            m_gammaTonemap->setParameter("invWhitePoint", invWhitePoint);
            m_gammaTonemap->setParameter("invGamma", 1/m_context->gamma);
            m_gammaTonemap->setParameter("sRGB", m_context->srgb);
            m_gammaTonemap->setParameter("hasDepth", hasDepth);
            m_renderer->blitTexture(buffer, m_context->mode == EPreview,
                !m_hScroll->isVisible(), !m_vScroll->isVisible(),
                -m_context->scrollOffset);
            m_gammaTonemap->unbind();
        } else if (m_context->toneMappingMethod == EReinhard) {
            if (m_luminanceBuffer[0] == NULL || m_luminanceBuffer[0]->getSize() != Point3i(size.x, size.y, 1)) {
                for (int i=0; i<2; ++i) {
                    m_luminanceBuffer[i] = m_renderer->createGPUTexture(formatString("Luminance buffer %i", i),
                            new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, size));
                    m_luminanceBuffer[i]->setComponentFormat(GPUTexture::EFloat32);
                    m_luminanceBuffer[i]->setPixelFormat(GPUTexture::ERGB);
                    m_luminanceBuffer[i]->setSize(Point3i(size.x, size.y, 1));
                    m_luminanceBuffer[i]->setFilterType(GPUTexture::ENearest);
                    m_luminanceBuffer[i]->setFrameBufferType(GPUTexture::EColorBuffer);
                    m_luminanceBuffer[i]->setMipMapped(false);
                    m_luminanceBuffer[i]->init();
                }
            }

            Float multiplier = 1.0f;
            if (m_context->mode == EPreview)
                multiplier /= entry.vplSampleOffset;

            /* Compute the luminace & log luminance */
            m_luminanceBuffer[0]->activateTarget();
            m_luminanceBuffer[0]->clear();
            buffer->bind();
            m_luminanceProgram->bind();
            m_luminanceProgram->setParameter("source", buffer);
            m_luminanceProgram->setParameter("multiplier", multiplier);
            m_renderer->blitQuad(m_context->mode == EPreview);
            m_luminanceProgram->unbind();
            buffer->unbind();
            m_luminanceBuffer[0]->releaseTarget();

            /* Iteratively downsample the image until left with a 1x1 pixel sum over
               the whole image */
            int source = 0, target = 1;
            Vector2i sourceSize(size.x, size.y);
            m_downsamplingProgram->bind();
            m_downsamplingProgram->setParameter("invSourceSize", Vector2(1.0f/size.x, 1.0f/size.y));
            while (sourceSize != Vector2i(1,1)) {
                target = 1-source;
                Vector2i targetSize((int) std::ceil(sourceSize.x/2.0f), (int) std::ceil(sourceSize.y/2.0f));
                m_luminanceBuffer[target]->activateTarget();
                m_luminanceBuffer[target]->clear();
                m_luminanceBuffer[source]->bind();
                m_downsamplingProgram->setParameter("source", m_luminanceBuffer[source]);
                m_downsamplingProgram->setParameter("sourceSize", Vector2(sourceSize.x, sourceSize.y));
                m_downsamplingProgram->setParameter("targetSize", Vector2(targetSize.x, targetSize.y));
                m_luminanceBuffer[target]->setTargetRegion(Point2i(0, 0), targetSize);
                m_renderer->blitQuad(true);
                m_luminanceBuffer[source]->unbind();
                m_luminanceBuffer[target]->releaseTarget();
                sourceSize = targetSize;
                source = target;
            }
            m_downsamplingProgram->unbind();
            Float logAvgLuminance, maxLuminance;
            Color3 result = m_luminanceBuffer[target]->getPixel(0, 0);
            logAvgLuminance = result[0];
            maxLuminance = result[1];
            logAvgLuminance = math::fastexp(logAvgLuminance / (size.x*size.y));

            if (!std::isfinite(logAvgLuminance)) {
                SLog(EWarn, "Could not determine the average log-luminance, since the image contains NaNs/infs/negative values");
                logAvgLuminance = 1;
            }

            Float burn = std::min((Float) 1, std::max((Float) 1e-8f, 1-(m_context->reinhardBurn + 10) / 20.0f)),
                  scale = m_context->reinhardKey / logAvgLuminance,
                  Lwhite = maxLuminance * scale;

            if (hasDepth)
                buffer->bind(1, 1);

            m_reinhardTonemap->bind();
            m_reinhardTonemap->setParameter("colorSource", 0);
            m_reinhardTonemap->setParameter("depthSource", 1);
            m_reinhardTonemap->setParameter("scale", scale);
            m_reinhardTonemap->setParameter("multiplier", multiplier);
            m_reinhardTonemap->setParameter("invWp2", 1 / (Lwhite * Lwhite * std::pow(burn, (Float) 4)));
            m_reinhardTonemap->setParameter("invGamma", 1/m_context->gamma);
            m_reinhardTonemap->setParameter("sRGB", m_context->srgb);
            m_reinhardTonemap->setParameter("hasDepth", hasDepth);
            m_renderer->blitTexture(buffer, m_context->mode == EPreview,
                !m_hScroll->isVisible(), !m_vScroll->isVisible(),
                -m_context->scrollOffset);
            m_reinhardTonemap->unbind();
        }

        if (m_context->mode == EPreview) {
            m_preview->releaseBuffer(entry);
            if (m_context->showKDTree) {
                std::string message = "kd-tree visualization mode\nPress '[' "
                    "and ']' to change the shown level";
                Vector2i size = m_font->getSize(message);
                int pos = m_statusMessage.empty() ? 10 : 30;
                m_renderer->setBlendMode(Renderer::EBlendAlpha);
                m_renderer->setColor(Spectrum(0.0f), 0.5f);
                m_renderer->drawFilledRectangle(Point2i(5, pos-2), Point2i(15, pos+5)+size);
                m_renderer->setBlendMode(Renderer::EBlendNone);
                m_renderer->setColor(Spectrum(1.0f));
                m_renderer->setColor(Spectrum(1.0f));
                m_renderer->drawText(Point2i(10, pos), m_font, message);
            }
        }

        if (m_cropping && m_cropStart != m_cropEnd) {
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            glEnable(GL_COLOR_LOGIC_OP);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glLogicOp(GL_XOR);
            Point2i cropStart = m_cropStart + upperLeft;
            Point2i cropEnd = m_cropEnd + upperLeft;
            m_renderer->setDepthTest(false);
            glRecti(cropStart.x, cropStart.y,
                    cropEnd.x, cropEnd.y);
            glDisable(GL_COLOR_LOGIC_OP);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        const ProjectiveCamera *camera = NULL;
        if (m_context->scene) {
            camera = getProjectiveCamera();
            if (!camera)
                m_statusMessage = "Camera type is incompatible with the OpenGL preview!";
        }

        if (!m_statusMessage.empty() && (m_context->mode == EPreview
            || (m_context->mode == ERender && m_context->scene == NULL))) {
            m_renderer->setDepthTest(false);
            Vector2i size = m_font->getSize(m_statusMessage);
            m_renderer->setBlendMode(Renderer::EBlendAlpha);
            m_renderer->setColor(Spectrum(0.0f), 0.5f);
            m_renderer->drawFilledRectangle(Point2i(5, 7), Point2i(15, 15)+size);
            m_renderer->setBlendMode(Renderer::EBlendNone);
            m_renderer->setColor(Spectrum(1.0f));
            m_renderer->drawText(Point2i(10, 10), m_font, m_statusMessage);
            if (m_statusTimer->getMilliseconds() > 2000)
                m_statusMessage = "";
        }

        if (m_context->mode == EPreview && camera) {
            m_renderer->setDepthTest(true);
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            m_renderer->setCamera(camera);
            glPushAttrib(GL_VIEWPORT_BIT);
            upperLeft = this->upperLeft(true);
            glViewport(upperLeft.x, upperLeft.y, size.x, size.y);
            m_renderer->setDepthMask(false);
            m_renderer->setDepthTest(true);
            m_renderer->setBlendMode(Renderer::EBlendAdditive);

            if (m_context->showKDTree) {
                oglRenderKDTree(m_context->scene->getKDTree());
                const ref_vector<Shape> &shapes = m_context->scene->getShapes();
                for (size_t j=0; j<shapes.size(); ++j)
                    if (shapes[j]->getKDTree())
                        oglRenderKDTree(shapes[j]->getKDTree());
            }

            if (m_navigationMode == EStandard && m_aabb.isValid())
                m_renderer->drawAABB(m_aabb);

            m_renderer->setBlendMode(Renderer::EBlendNone);
            glPopAttrib();
        }
    }
    swapBuffers();
}

void GLWidget::oglRenderKDTree(const KDTreeBase<AABB> *kdtree) {
    std::stack<boost::tuple<const KDTreeBase<AABB>::KDNode *, AABB, uint32_t> > stack;

    stack.push(boost::make_tuple(kdtree->getRoot(), kdtree->getTightAABB(), 0));
    Float brightness = 0.1f;

    glEnable(GL_LINE_SMOOTH);
    glLineWidth(0.6f);
    while (!stack.empty()) {
        const KDTreeBase<AABB>::KDNode *node = boost::get<0>(stack.top());
        AABB aabb = boost::get<1>(stack.top());
        int level = boost::get<2>(stack.top());
        stack.pop();
        m_renderer->setColor(Spectrum(brightness));
        m_renderer->drawAABB(aabb);

        if (!node->isLeaf()) {
            int axis = node->getAxis();
            float split = node->getSplit();
            if (level + 1 <= m_context->shownKDTreeLevel) {
                Float tmp = aabb.max[axis];
                aabb.max[axis] = split;
                stack.push(boost::make_tuple(node->getLeft(), aabb, level+1));
                aabb.max[axis] = tmp;
                aabb.min[axis] = split;
                stack.push(boost::make_tuple(node->getRight(), aabb, level+1));
            } else {
                aabb.min[axis] = split;
                aabb.max[axis] = split;
                Spectrum color;
                color.fromLinearRGB(0, 0, 4*brightness);
                m_renderer->setColor(color);
                m_renderer->drawAABB(aabb);
            }
        }
    }
    glDisable(GL_LINE_SMOOTH);
}


void GLWidget::resizeEvent(QResizeEvent *event) {
    if (m_context && m_ignoreResizeEvents) {
        event->accept();
        return;
    }
    QGLWidget::resizeEvent(event);
}

void GLWidget::updateScrollBars() {
    if (!m_context) {
        m_hScroll->hide();
        m_vScroll->hide();
        return;
    }

    float ratio = windowHandle()->devicePixelRatio();
    int reqWidth = m_context->framebuffer->getWidth(),
        reqHeight = m_context->framebuffer->getHeight(),
        width = size().width()*ratio,
        height = size().height()*ratio;

    if (m_hScroll->isVisible())
        height += m_hScroll->size().height();
    if (m_vScroll->isVisible())
        width += m_vScroll->size().width();

    bool needsHScroll = false, needsVScroll = false;

    m_ignoreScrollEvents = true;
    for (int i=0; i<2; ++i) {
        if (width < reqWidth) {
            if (!m_hScroll->isVisible())
                m_hScroll->show();
            if (!needsHScroll) {
                height -= m_hScroll->sizeHint().height();
                needsHScroll = true;
            }
            m_hScroll->setMaximum(reqWidth-width);
            m_hScroll->setPageStep(width);
            if (m_context->scrollOffset.x + width > reqWidth)
                m_context->scrollOffset.x = reqWidth - width;
        } else if (m_hScroll->isVisible()) {
            m_hScroll->setMaximum(0);
            m_hScroll->hide();
            m_context->scrollOffset.x = 0;
        }

        if (height < reqHeight) {
            if (!m_vScroll->isVisible()) {
                m_vScroll->show();
            }
            if (!needsVScroll) {
                width -= m_vScroll->sizeHint().width();
                needsVScroll = true;
            }
            m_vScroll->setMaximum(reqHeight-height);
            m_vScroll->setPageStep(height);
            if (m_context->scrollOffset.y + height > reqHeight)
                m_context->scrollOffset.y = reqHeight - height;
        } else if (m_vScroll->isVisible()) {
            m_vScroll->setMaximum(0);
            m_vScroll->hide();
            m_context->scrollOffset.y = 0;
        }
    }

    m_hScroll->setValue(m_context->scrollOffset.x);
    m_vScroll->setValue(m_context->scrollOffset.y);
    m_ignoreScrollEvents = false;
    updateGeometry();
}

void GLWidget::setScrollBars(QScrollBar *hScroll, QScrollBar *vScroll) {
    m_hScroll = hScroll;
    m_vScroll = vScroll;
    connect(m_hScroll, SIGNAL(valueChanged(int)), this, SLOT(onScroll()));
    connect(m_vScroll, SIGNAL(valueChanged(int)), this, SLOT(onScroll()));
}

void GLWidget::onScroll() {
    if (!m_context || m_ignoreScrollEvents)
        return;
    m_context->scrollOffset.x = m_hScroll->value();
    m_context->scrollOffset.y = m_vScroll->value();
    updateGL();
}

void GLWidget::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasFormat("text/uri-list")) {
        event->setDropAction(Qt::CopyAction);
        event->acceptProposedAction();
    }
}

void GLWidget::dropEvent(QDropEvent *event) {
    QList<QUrl> urls = event->mimeData()->urls();
    for (int i=0; i<urls.size(); ++i)
        emit loadFileRequest(urls[i].toLocalFile());
    event->acceptProposedAction();
}

void GLWidget::resumePreview() {
    if (m_preview->isRunning())
        m_preview->resume();
}

void GLWidget::onUpdateView() {
    m_framebufferChanged = true;
}

Point2i GLWidget::upperLeft(bool flipY) const {
    if (!m_context)
        return Point2i(0, 0);
    Vector2i outputSize(m_context->framebuffer->getWidth(),
            m_context->framebuffer->getHeight());
    Vector2i deviceSize = m_device->getSize();

    return Point2i(
        deviceSize.x < outputSize.x ? (-m_context->scrollOffset.x) : (deviceSize.x - outputSize.x)/2,
        deviceSize.y < outputSize.y ? (flipY ? (deviceSize.y-outputSize.y
            + m_context->scrollOffset.y) : -m_context->scrollOffset.y) : (deviceSize.y - outputSize.y)/2);
}

void GLWidget::reveal(const AABB &aabb) {
    if (!m_context || m_animation || !aabb.isValid())
        return;
    BSphere bsphere = aabb.getBSphere();
    PerspectiveCamera *camera = getPerspectiveCamera();
    if (!camera)
        return;

    Transform invView = getWorldTransform();
    Point p  = invView(Point(0,0,0));
    Vector d = invView(Vector(0,0,1));

    Float fov = std::min(camera->getXFov(), camera->getYFov())*0.9f/2;
    Float distance = bsphere.radius/std::tan(fov * M_PI/180.0f);
    camera->setFocusDistance(distance);

    m_animationTimer->reset();
    m_animationOrigin0 = p;
    m_animationOrigin1 = bsphere.center - distance*d;
    m_animationTarget0 = m_animationOrigin0 + d;
    m_animationTarget1 = m_animationOrigin1 + d;
    m_animation = true;

    if (!m_movementTimer->isActive())
        m_movementTimer->start();
}

void GLWidget::resizeGL(int width, int height) {
    glViewport(0, 0, (GLint) width, (GLint) height);
    m_device->setSize(Vector2i(width, height));
}

