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

#include <QtGui>
#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
#include "glwidget.h"
#include "preview.h"
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/hw/font.h>

GLWidget::GLWidget(QWidget *parent) :
	QGLWidget(parent), m_context(NULL) {
	setAutoBufferSwap(false);
	setFocusPolicy(Qt::StrongFocus);
	m_clock = new Timer();
	m_wheelTimer = new Timer();
	m_animationTimer = new Timer();
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

void GLWidget::initializeGL() {
	/* Load the Mitsuba logo into a texture */
	QResource res("/resources/logo.png");
	SAssert(res.isValid());
	ref<Stream> mStream = new MemoryStream(res.size());
	mStream->write(res.data(), res.size());
	mStream->setPos(0);
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

	if (missingExtensions.size() > 0 || m_softwareFallback) {
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
		oss << ". Please make sure that you are using the most "
			<< "recent graphics drivers.\n\nMitsuba will now switch "
			<< "to a slow software fallback mode, which only supports "
			<< "the rendering preview but no tonemapping and no "
			<< "real-time preview/navigation.";
		m_errorString = QString(oss.str().c_str());
		m_softwareFallback = true;
#endif
		// Don't redraw as often, since this is now quite costly
		m_redrawTimer->setInterval(1000);
	} else {
		m_gammaTonemap = m_renderer->createGPUProgram("Tonemapper [Gamma]");
		m_reinhardTonemap = m_renderer->createGPUProgram("Tonemapper [Reinhard et al. 2002]");

		m_gammaTonemap->setSource(GPUProgram::EVertexProgram,
			"void main() {\n"
			"	gl_Position = ftransform();\n"
			"   gl_TexCoord[0]  = gl_MultiTexCoord0;\n"
			"}\n"
		);

		m_gammaTonemap->setSource(GPUProgram::EFragmentProgram,
			"#version 120\n"
			"uniform sampler2D colorSource, depthSource;\n"
			"uniform float invWhitePoint, invGamma;\n"
			"uniform bool sRGB, hasDepth;\n"
			"\n"
			"float toSRGB(float value) {\n"
			"	if (value < 0.0031308)\n"
			"		return 12.92 * value;\n"
			"	return 1.055 * pow(value, 0.41666) - 0.055;\n"
			"}\n"
			"\n"
			"void main() {\n"
			"	vec4 color = texture2D(colorSource, gl_TexCoord[0].xy) * invWhitePoint;\n"
			"	if (sRGB)\n"
			"		gl_FragColor = vec4(toSRGB(color.r), toSRGB(color.g), toSRGB(color.b), 1);\n"
			"	else\n"
			"		gl_FragColor = vec4(pow(color.rgb, vec3(invGamma)), 1);\n"
			"	gl_FragDepth = hasDepth ? texture2D(depthSource, gl_TexCoord[0].xy).r : 0.5;\n"
			"}\n"
		);

		m_reinhardTonemap->setSource(GPUProgram::EVertexProgram,
			"void main() {\n"
			"	gl_Position = ftransform();\n"
			"   gl_TexCoord[0]  = gl_MultiTexCoord0;\n"
			"}\n"
		);

		m_reinhardTonemap->setSource(GPUProgram::EFragmentProgram,
			"#version 120\n"
			"uniform sampler2D colorSource, depthSource;\n"
			"uniform float key, invWpSqr, invGamma, multiplier;\n"
			"uniform bool sRGB, hasDepth;\n"
			"\n"
			"float toSRGB(float value) {\n"
			"	if (value < 0.0031308)\n"
			"		return 12.92 * value;\n"
			"	return 1.055 * pow(value, 0.41666) - 0.055;\n"
			"}\n"
			"\n"
			"void main() {\n"
			"	const mat3 rgb2xyz = mat3(0.412453,  0.357580,  0.180423,\n"
			"                             0.212671,  0.715160,  0.072169,\n"
			"					          0.019334,  0.119193,  0.950227);\n"
			"\n"
			"	const mat3 xyz2rgb = mat3(3.240479, -1.537150, -0.498535,\n"
			"                            -0.969256,  1.875991,  0.041556,\n"
			"					          0.055648, -0.204043,  1.057311);\n"
			"\n"
			"	vec4 color = texture2D(colorSource, gl_TexCoord[0].xy)*multiplier;\n"
			"   vec3 xyz = rgb2xyz * color.rgb;\n"
			"   float normalization = 1.0/(xyz.x + xyz.y + xyz.z);\n"
			"   vec3 Yxy = vec3(xyz.x*normalization, xyz.y*normalization, xyz.y);\n"
			"   float Lp = Yxy.z*key;\n"
			"   Yxy.z = Lp * (1.0 + Lp*invWpSqr) / (1.0+Lp);\n"
			"   xyz = vec3(Yxy.x * (Yxy.z/Yxy.y), Yxy.z, (Yxy.z/Yxy.y) * (1.0 - Yxy.x - Yxy.y));\n"
			"	color.rgb = xyz2rgb * xyz;\n"
			"	if (sRGB)\n"
			"		gl_FragColor = vec4(toSRGB(color.r), toSRGB(color.g), toSRGB(color.b), 1);\n"
			"	else\n"
			"		gl_FragColor = vec4(pow(color.rgb, vec3(invGamma)), 1);\n"
			"	gl_FragDepth = hasDepth ? texture2D(depthSource, gl_TexCoord[0].xy).r : 0.5;\n"
			"}\n"
		);

		m_luminanceProgram = m_renderer->createGPUProgram("Log-luminance program");
		m_luminanceProgram->setSource(GPUProgram::EVertexProgram,
			"void main() {\n"
			"	gl_Position = ftransform();\n"
			"   gl_TexCoord[0]  = gl_MultiTexCoord0;\n"
			"}\n"
		);

		m_luminanceProgram->setSource(GPUProgram::EFragmentProgram,
			"uniform sampler2D source;\n"
			"uniform float multiplier;\n"
			"\n"
			"void main() {\n"
			"	vec4 color = texture2D(source, gl_TexCoord[0].xy);\n"
			"	float luminance = multiplier * (color.r * 0.212671 + color.g * 0.715160 + color.b * 0.072169);\n"
			"	if (luminance < 0.0 || luminance != luminance) luminance = 0.0; // catch NaNs and negative numbers\n"
			"	float logLuminance = log(0.001+luminance);\n"
			"	gl_FragColor = vec4(logLuminance, luminance, 0.0, 1.0);"
			"}\n"
		);

		m_downsamplingProgram = m_renderer->createGPUProgram("Downsampling program");
		m_downsamplingProgram->setSource(GPUProgram::EVertexProgram,
			"uniform vec2 targetSize;\n"
			"void main() {\n"
			"	gl_Position = ftransform();\n"
			"   gl_TexCoord[0].xy = vec2(gl_MultiTexCoord0.x * targetSize.x, \n"
			"                            gl_MultiTexCoord0.y * targetSize.y);\n"
			"}\n"
		);

		m_downsamplingProgram->setSource(GPUProgram::EFragmentProgram,
			"uniform sampler2D source;\n"
			"uniform vec2 activeRegionSize;\n"
			"uniform vec2 invSourceSize;\n"
			"\n"
			"/* Perform a texture lookup by pixel coordinates */\n"
			"vec4 lookupPixel(vec2 coords) {\n"
			"   coords = coords + vec2(0.5, 0.5);\n"
			"	if (coords.x < 0.0 || coords.y < 0.0 ||\n"
			"		coords.x > activeRegionSize.x || coords.y > activeRegionSize.y)\n"
			"		return vec4(0);\n"
			"	else\n"
			"		return texture2D(source, coords*invSourceSize);\n"
			"}\n"
			"\n"
			"/* Find the max. luminance and the sum of all log-luminance values */\n"
			"void main() {\n"
			"   vec2 pos = (gl_TexCoord[0].xy-vec2(.5, .5))*2.0;\n"
			"	vec2 pixel0 = lookupPixel(pos).rg,\n"
			"        pixel1 = lookupPixel(pos + vec2(1, 0)).rg,\n"
			"        pixel2 = lookupPixel(pos + vec2(0, 1)).rg,\n"
			"        pixel3 = lookupPixel(pos + vec2(1, 1)).rg;\n"
			"	gl_FragColor.r = pixel0.r + pixel1.r + pixel2.r + pixel3.r;\n"
			"	gl_FragColor.g = max(pixel0.g, max(pixel1.g, max(pixel2.g, pixel3.g)));\n"
			"	gl_FragColor.ba = vec2(0, 1);\n"
			"}\n"
		);

		if (!m_preview->isRunning()) {
#if defined(WIN32)
			wglMakeCurrent(NULL, NULL);
#endif
			m_preview->start();
			m_preview->waitUntilStarted();
#if defined(WIN32)
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
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(0.6f);
}

void GLWidget::setScene(SceneContext *context) {
	if (m_movementTimer->isActive())
		m_movementTimer->stop();
	m_context = context;

	if (context && context->scene == NULL)
		context = NULL;
	m_preview->setSceneContext(context, true, false);
	m_framebufferChanged = true;
	m_mouseDrag = m_animation = false;
	m_leftKeyDown = m_rightKeyDown = m_upKeyDown = m_downKeyDown = false;
	m_aabb.reset();
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
	ref<Bitmap> sourceBitmap = new Bitmap(size.x, size.y, 128);
	m_renderer->finish();
	entry.buffer->download(sourceBitmap);

	// Scale by the number of developed VPLs 
	const float *sourceData = sourceBitmap->getFloatData();
	float *targetData = m_context->framebuffer->getFloatData();
	float factor = 1.0f / entry.vplSampleOffset;

	for (size_t pos=0, total = size.x*size.y*4; pos<total; ++pos)
		targetData[pos] = sourceData[pos] * factor;

	if (m_context->mode == EPreview && (m_context->previewMethod == ERayTrace
			|| m_context->previewMethod == ERayTraceCoherent)) {
		/* Set alpha channel to 1 */
		for (size_t pos=0, total = size.x*size.y*4; pos<total; pos+=4)
			targetData[pos+3] = 1.0f;
	}

	if (createdFramebuffer)
		m_framebuffer->init();
	else
		m_framebuffer->refresh();
	m_preview->releaseBuffer(entry);
}

QSize GLWidget::sizeHint() const {
	QSize minimumSize(440, 170);
	if (m_context) {
		return QSize(std::max(minimumSize.width(), m_context->framebuffer->getWidth()), 
					 std::max(minimumSize.height(), m_context->framebuffer->getHeight()));
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
	ProjectiveCamera *camera 
		= static_cast<ProjectiveCamera *>(m_context->scene->getCamera());

	if (m_animation) {
		Float x = std::min(m_animationTimer->getMilliseconds() / 500.0f, 1.0f);
		Float t = x*x*x*(x*(x*6-15)+10); // smootherstep by Ken Perlin

		Point origin = (1-t) * m_animationOrigin0 + t * m_animationOrigin1;
		Point target = (1-t) * m_animationTarget0 + t * m_animationTarget1;

		if (camera->getViewTransform().det3x3() > 0) 
			camera->setInverseViewTransform(
			Transform::lookAt(origin, target, m_context->up));
		else 
			camera->setInverseViewTransform(
				Transform::lookAt(origin, target, m_context->up) *
				Transform::scale(Vector(-1,1,1)));

		if (x == 1.0f)
			m_animation = false;
	}

	if (!(m_mouseDrag && m_navigationMode == EStandard) && !m_animation) {
		Float moveSpeed = m_context->movementScale
			* m_clock->getMilliseconds();

		if (m_leftKeyDown)
			camera->setViewTransform(
			Transform::translate(Vector(moveSpeed,0,0))
			* camera->getViewTransform());
		if (m_rightKeyDown)
			camera->setViewTransform(
			Transform::translate(Vector(-moveSpeed,0,0))
			* camera->getViewTransform());
		if (m_downKeyDown)
			camera->setViewTransform(
			Transform::translate(Vector(0,0,moveSpeed))
			* camera->getViewTransform());
		if (m_upKeyDown)
			camera->setViewTransform(
			Transform::translate(Vector(0,0,-moveSpeed))
			* camera->getViewTransform());
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

void GLWidget::resetPreview() {
	if (!m_context || !m_context->scene || !m_preview->isRunning())
		return;
	bool motion = m_leftKeyDown || m_rightKeyDown || 
		m_upKeyDown || m_downKeyDown || m_mouseDrag ||
		m_wheelTimer->getMilliseconds() < 200 || m_animation;
	m_preview->setSceneContext(m_context, false, motion);
	updateGL();
}

void GLWidget::keyPressEvent(QKeyEvent *event) {
	if (event->key() == Qt::Key_Escape)
		emit quit();
	else if (event->key() == Qt::Key_R)
		emit beginRendering();

	if (event->isAutoRepeat() || !m_context)
		return;
	
	switch (event->key()) {
		case Qt::Key_PageUp:   m_context->movementScale *= 2; break;
		case Qt::Key_PageDown: m_context->movementScale /= 2; break;
		case Qt::Key_Left:
			m_leftKeyDown = true; break;
		case Qt::Key_Right:
			m_rightKeyDown = true; break;
		case Qt::Key_Up:
			m_upKeyDown = true; break;
		case Qt::Key_Down:
			m_downKeyDown = true; break;
		case Qt::Key_A: {
			if (m_context->selectionMode == EScene) {
				m_context->selectionMode = ENothing;
				m_aabb.reset();
			} else {
				m_context->selectionMode = EScene;
				m_aabb = m_context->scene->getAABB();
			}
			m_context->selectedShape = NULL;
		}
		// break intentionally missing
		case Qt::Key_F: {
			reveal(m_aabb);
		};
		break;
	}
	if (!m_movementTimer->isActive() && (m_leftKeyDown || m_rightKeyDown || m_upKeyDown || m_downKeyDown)) {
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
		case Qt::Key_BracketLeft:
			if (m_context->showKDTree) {
				m_context->shownKDTreeLevel
					= std::max(0, m_context->shownKDTreeLevel - 1);
				resetPreview();
			}
			break;
		case Qt::Key_BracketRight:
			if (m_context->showKDTree) {
				m_context->shownKDTreeLevel++;
				resetPreview();
			}
			break;
	}
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
	if (!m_preview->isRunning() || m_animation)
		return;

	QPoint rel = event->pos() - m_mousePos;
	m_mousePos = event->pos();

	if (!m_context || !m_context->scene || !m_mouseDrag || rel == QPoint(0,0))
		return;

	//	if (m_ignoreMouseEvent == rel) {
	if (m_ignoreMouseEvent != QPoint(0, 0)) {
		m_ignoreMouseEvent = QPoint(0, 0);
		return;
	}

	if (!m_didSetCursor) {
		QApplication::setOverrideCursor(Qt::BlankCursor);
		m_didSetCursor = true;
	}

	PerspectiveCamera *camera = static_cast<PerspectiveCamera *>(m_context->scene->getCamera());
	Point p = camera->getInverseViewTransform()(Point(0,0,0));
	Vector d = camera->getInverseViewTransform()(Vector(0,0,1));
	bool didMove = false;
	
	Float focusDepth = camera->getFocusDepth(),
		nearClip = camera->getNearClip(),
		farClip = camera->getFarClip();

	if (focusDepth <= nearClip || focusDepth >= farClip) {
		focusDepth = autoFocus();
		camera->setFocusDepth(focusDepth);
	}

	Point target = p + d * focusDepth;
	Vector up = m_context->up;

	if (event->buttons() & Qt::LeftButton) {
		if (m_navigationMode == EStandard) {
			Frame frame(up);
			Point2 coords = toSphericalCoordinates(frame.toLocal(normalize(p-target)));
			coords.y -=  0.001f * rel.x() * m_mouseSensitivity * (camera->getViewTransform().det3x3() < 0 ? 1 : -1);
			coords.x -=  0.001f * rel.y() * m_mouseSensitivity * (m_invertMouse ? -1.0f : 1.0f);
			p = target + focusDepth * frame.toWorld(sphericalDirection(coords.x, coords.y));

			if (coords.x < 0 || coords.x > M_PI) 
				m_context->up *= -1;

			if (camera->getViewTransform().det3x3() > 0) 
				camera->setInverseViewTransform(Transform::lookAt(p, target, m_context->up));
			else
				camera->setInverseViewTransform(
					Transform::lookAt(p, target, m_context->up) *
					Transform::scale(Vector(-1,1,1))
				);
		} else {
			Float yaw = -.03f * rel.x() * m_mouseSensitivity;
			Float pitch = -.03f * rel.y() * m_mouseSensitivity;
			if (m_invertMouse) 
				pitch *= -1;

			Transform trafo = Transform::rotate(Vector(0,1,0), yaw)
					* Transform::rotate(Vector(1,0,0), pitch)
					* camera->getViewTransform();
			d = trafo.inverse()(Vector(0,0,1));

			if (camera->getViewTransform().det3x3() > 0) 
				camera->setInverseViewTransform(Transform::lookAt(p, p+d, up));
			else
				camera->setInverseViewTransform(
					Transform::lookAt(p, p+d, up) *
					Transform::scale(Vector(-1,1,1))
				);
		}
		didMove = true;
	} else if (event->buttons() & Qt::MidButton) {
		camera->setViewTransform(
			Transform::translate(Vector((Float) rel.x(), (Float) -rel.y(), 0) 
				* m_mouseSensitivity * .6f * m_context->movementScale)
			* camera->getViewTransform());
		didMove = true;
	} else if (event->buttons() & Qt::RightButton) {
		if (event->modifiers() & Qt::ShiftModifier) {
			Float roll = rel.x() * m_mouseSensitivity * .02f;
			Float fovChange = rel.y() * m_mouseSensitivity * .03f;

			if (camera->getViewTransform().det3x3() > 0) {
				m_context->up = Transform::rotate(d, roll)(up);
				camera->setInverseViewTransform(Transform::lookAt(p, p+d, m_context->up));
			} else {
				m_context->up = Transform::rotate(d, -roll)(up);
				camera->setInverseViewTransform(
					Transform::lookAt(p, p+d, m_context->up) *
					Transform::scale(Vector(-1,1,1))
				);
			}

			camera->setFov(std::min(std::max((Float) 1.0f, camera->getFov() 
				+ fovChange), (Float) 150.0f));
		} else {
			Float focusDepth = camera->getFocusDepth(),
				nearClip = camera->getNearClip(),
				farClip = camera->getFarClip();

			if (focusDepth <= nearClip || focusDepth >= farClip) 
				focusDepth = autoFocus();

			Float oldFocusDepth = focusDepth;
			focusDepth = std::min(std::max(focusDepth * std::pow((Float) (1 - 2e-3f),
					(Float) -rel.y() * m_mouseSensitivity * m_context->movementScale),
					1.2f*nearClip), farClip/1.2f);

			camera->setFocusDepth(focusDepth);
			Vector d = Vector(camera->getImagePlaneNormal());
			p = p + (oldFocusDepth - focusDepth) * d;

			if (camera->getViewTransform().det3x3() > 0) 
				camera->setInverseViewTransform(Transform::lookAt(p, p+d, up));
			else
				camera->setInverseViewTransform(
					Transform::lookAt(p, p+d, up) *
					Transform::scale(Vector(-1,1,1))
				);
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
	
void GLWidget::wheelEvent(QWheelEvent *event) {
	if (!m_preview->isRunning() || m_context == NULL || m_context->scene == NULL || m_animation)
		return;

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
		PerspectiveCamera *camera 
			= static_cast<PerspectiveCamera *>(m_context->scene->getCamera());
		Float focusDepth = camera->getFocusDepth(),
			nearClip = camera->getNearClip(),
			farClip = camera->getFarClip();

		if (focusDepth <= nearClip || focusDepth >= farClip) 
			focusDepth = autoFocus();

		Float oldFocusDepth = focusDepth;
		focusDepth = std::min(std::max(focusDepth * std::pow(1 - 1e-3f, 
				event->delta()), 1.2f*nearClip), farClip/1.2f);

		camera->setFocusDepth(focusDepth);

		Vector up = camera->getInverseViewTransform()(Vector(0,1,0));
		Vector d = Vector(camera->getImagePlaneNormal());
		Point o = camera->getPosition() + (oldFocusDepth - focusDepth) * d;

		if (camera->getViewTransform().det3x3() > 0) 
			camera->setInverseViewTransform(Transform::lookAt(o, o+d, up));
		else
			camera->setInverseViewTransform(
				Transform::lookAt(o, o+d, up) *
				Transform::scale(Vector(-1,1,1))
			);

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
	Vector2i filmSize(scene->getFilm()->getSize());
	const PerspectiveCamera *camera = static_cast<const PerspectiveCamera *>(
		scene->getCamera());
	Float t, avgDistance = 0, importance = 0;
	ConstShapePtr ptr;
	Normal n;
	Ray ray;
	Float variance = std::pow(filmSize.x / (Float) 4, (Float) 2);

	for (size_t sampleIndex=0; sampleIndex<200; ++sampleIndex) {
		Point2 sample(
			radicalInverse(2, sampleIndex) * filmSize.x,
			radicalInverse(3, sampleIndex) * filmSize.y);
		scene->getCamera()->generateRay(sample, Point2(0, 0), 0, ray);
		if (scene->rayIntersect(ray, t, ptr, n)) {
			Float weight = std::fastexp(-0.5 / variance * (
				std::pow(sample.x - filmSize.x / (Float) 2, (Float) 2) +
				std::pow(sample.y - filmSize.y / (Float) 2, (Float) 2)));
			avgDistance += t * weight;
			importance += weight;
		}
	}

	if (importance == 0)
		return 0.5f * (camera->getNearClip() + camera->getFarClip());
	else
		return avgDistance / importance;
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
	if (!m_preview->isRunning() || m_context == NULL || m_context->scene == NULL)
		return;
	m_mousePos = event->pos();
	m_initialMousePos = mapToGlobal(m_mousePos);
	m_mouseDrag = true;
	ProjectiveCamera *camera = static_cast<ProjectiveCamera *>(m_context->scene->getCamera());

	if (m_navigationMode == EStandard && event->buttons() & Qt::LeftButton) {
		Point2i offset = upperLeft();
		Point2 sample = Point2(m_mousePos.x() - offset.x, m_mousePos.y() - offset.y);
		Intersection its;
		Float nearT, farT;
		Ray ray;

		camera->generateRay(sample, Point2(0, 0), 0, ray);
		if (!m_aabb.isValid() || !m_aabb.rayIntersect(ray, nearT, farT))
			nearT = std::numeric_limits<Float>::infinity();

		if (m_context->scene->rayIntersect(ray, its)) {
			if (its.t < nearT) {
				SLog(EInfo, "Selected shape \"%s\"", its.shape->getName().c_str());
				m_context->selectedShape = its.shape;
				bool newSelection = (m_aabb != its.shape->getAABB());
				m_aabb = its.shape->getAABB();
				m_context->selectionMode = EShape;
				if (newSelection) 
					return;
			}
		} else {
			if (nearT == std::numeric_limits<Float>::infinity()) {
				m_aabb.reset();
				m_context->selectionMode = ENothing;
				m_context->selectedShape = NULL;
				return;
			}
		}

		m_storedViewTransform = camera->getViewTransform();
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
	if (!m_preview->isRunning())
		return;
	if (event->buttons() == 0) {
		m_mouseDrag = false;
		if (m_didSetCursor) {
			resetPreview();
			QApplication::restoreOverrideCursor();
			QCursor::setPos(m_initialMousePos);
			m_didSetCursor = false;
		}
	}
}

void GLWidget::paintGL() {
	m_renderer->setDepthTest(true);
	m_renderer->setDepthMask(true);
	if (m_context == NULL) {
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		m_renderer->clear();
		m_renderer->setBlendMode(Renderer::EBlendAlpha);
		m_renderer->blitTexture(m_logoTexture);
		m_renderer->setBlendMode(Renderer::EBlendNone);
	} else if (m_context != NULL) {
		Vector2i size;
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		PreviewQueueEntry entry;
		GPUTexture *buffer = NULL;

		if (m_context->mode == EPreview) {
			if (!m_preview->isRunning() || m_context->previewMethod == EDisabled) {
				/* No preview thread running - just show a grey screen */
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
					m_fallbackBitmap = new Bitmap(m_context->framebuffer->getWidth(),
						m_context->framebuffer->getHeight(), 24);
					m_fallbackBitmap->clear();
					m_framebuffer = m_renderer->createGPUTexture("Framebuffer", 
						m_fallbackBitmap);
				}
				m_framebuffer->setMipMapped(false);
				m_framebuffer->setFilterType(GPUTexture::ENearest);
				m_framebuffer->init();
			}

			if (m_framebufferChanged) {
				if (m_softwareFallback) {
					/* Manually generate a gamma-corrected image 
					   on the CPU (with gamma=2.2) - this will be slow! */
					Bitmap *source = m_context->framebuffer;
					float *sourceData = source->getFloatData();
					uint8_t *targetData = m_fallbackBitmap->getData();
					for (int y=0; y<source->getHeight(); ++y) {
						for (int x=0; x<source->getWidth(); ++x) {
							const float invGammaValue = 0.45455f;
							*targetData++ = (uint8_t) std::max(std::min(std::pow(*sourceData++, invGammaValue) * 255.0f, 255.0f), 0.0f);
							*targetData++ = (uint8_t) std::max(std::min(std::pow(*sourceData++, invGammaValue) * 255.0f, 255.0f), 0.0f);
							*targetData++ = (uint8_t) std::max(std::min(std::pow(*sourceData++, invGammaValue) * 255.0f, 255.0f), 0.0f);
							sourceData++;
						}
					}
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
			&& (m_context->previewMethod == EOpenGL ||
				m_context->previewMethod == EOpenGLSinglePass);

		if (m_softwareFallback) {
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
					m_luminanceBuffer[i] = m_renderer->createGPUTexture(formatString("Luminance buffer %i", i)); 
					m_luminanceBuffer[i]->setFormat(GPUTexture::EFloat32RGB);
					m_luminanceBuffer[i]->setSize(Point3i(size.x, size.y, 1));
					m_luminanceBuffer[i]->setFilterType(GPUTexture::ENearest);
					m_luminanceBuffer[i]->setFrameBufferType(GPUTexture::EColorBuffer);
					m_luminanceBuffer[i]->setMipMapped(false);
					m_luminanceBuffer[i]->init();
				}
			}

			Float multiplier = 1.0;
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
			#if 0
				/* For debugging .. */
				ref<Bitmap> bitmap = new Bitmap(size.x, size.y, 128);
				m_luminanceBuffer[0]->download(bitmap);
				ref<FileStream> fs = new FileStream("test.exr", FileStream::ETruncReadWrite);
				bitmap->save(Bitmap::EEXR, fs);
				fs->close();
			#endif

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
				m_luminanceBuffer[source]->bind();
				m_downsamplingProgram->setParameter("source", m_luminanceBuffer[source]);
				m_downsamplingProgram->setParameter("activeRegionSize", Vector2(sourceSize.x, sourceSize.y));
				m_downsamplingProgram->setParameter("targetSize", Vector2(targetSize.x, targetSize.y));
				m_luminanceBuffer[target]->setTargetRegion(Point2i(0, 0), targetSize);
				m_renderer->blitQuad(true);
				m_luminanceBuffer[source]->unbind();
				m_luminanceBuffer[target]->releaseTarget();
				sourceSize = targetSize;
				source = target;
			}
			m_downsamplingProgram->unbind();
			Float logLuminance, maxLuminance, unused;
			m_luminanceBuffer[target]->getPixel(0, 0).toLinearRGB(logLuminance, maxLuminance, unused);
			logLuminance = std::fastexp(logLuminance / (size.x*size.y));
			if (mts_isnan(logLuminance) || std::isinf(logLuminance)) {
				SLog(EWarn, "Could not determine the average log-luminance, since the image contains NaNs/infs/negative values");
				logLuminance = 1;
			}

			if (hasDepth)
				buffer->bind(1, 1);
			m_reinhardTonemap->bind();
			m_reinhardTonemap->setParameter("colorSource", 0);
			m_reinhardTonemap->setParameter("depthSource", 1);
			m_reinhardTonemap->setParameter("key", m_context->reinhardKey/logLuminance);
			m_reinhardTonemap->setParameter("multiplier", multiplier);
			m_reinhardTonemap->setParameter("invWpSqr", std::pow((Float) 2, m_context->reinhardBurn));
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
				m_renderer->setColor(Spectrum(1.0f));
				m_renderer->drawText(Point2i(10, 10), m_font, 
					formatString("kd-tree visualization mode\nPress '[' "
					"and ']' to change the shown level"));
			}
		}

		GLint viewport[4];	
		glGetIntegerv(GL_VIEWPORT, viewport);
		if (viewport[2] > size.x || viewport[3] > size.y) {
			/* Draw a border to highlight the region occupied by the image */
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			Vector2i upperLeft((viewport[2] - size.x)/2, (viewport[3] - size.y)/2);
			Vector2i lowerRight = upperLeft + Vector2i(size.x, size.y);
			glColor4f(0.4f, 0.4f, 0.4f, 1.0f);
			glBegin(GL_LINE_LOOP);
			glVertex2f(upperLeft.x, upperLeft.y);
			glVertex2f(lowerRight.x, upperLeft.y);
			glVertex2f(lowerRight.x, lowerRight.y);
			glVertex2f(upperLeft.x, lowerRight.y);
			glEnd();
			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		}

		if (m_context->mode == EPreview) {
			const ProjectiveCamera *camera = static_cast<const ProjectiveCamera *>
				(m_context->scene->getCamera());
			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
			m_renderer->setCamera(camera);
			glPushAttrib(GL_VIEWPORT_BIT);
			Point2i offset = upperLeft(true);
			glViewport(offset.x, offset.y, size.x, size.y);
			m_renderer->setDepthMask(false);
			m_renderer->setDepthTest(true);
			m_renderer->setBlendMode(Renderer::EBlendAdditive);

			if (m_context->showKDTree) {
				oglRenderKDTree(m_context->scene->getKDTree());
				const std::vector<Shape *> shapes = m_context->scene->getShapes();
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

	int reqWidth = m_context->framebuffer->getWidth(),
		reqHeight = m_context->framebuffer->getHeight(),
		width = size().width(),
		height = size().height();

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
	Vector2i filmSize = m_context->scene->getCamera()->getFilm()->getSize();
	Vector2i deviceSize = m_device->getSize();

	return Point2i(
		deviceSize.x < filmSize.x ? (-m_context->scrollOffset.x) : (deviceSize.x - filmSize.x)/2,
		deviceSize.y < filmSize.y ? (flipY ? (deviceSize.y-filmSize.y 
			+ m_context->scrollOffset.y) : -m_context->scrollOffset.y) : (deviceSize.y - filmSize.y)/2);
}

void GLWidget::reveal(const AABB &aabb) {
	if (!m_context || m_animation || !aabb.isValid())
		return;
	BSphere bsphere = aabb.getBSphere();
	PerspectiveCamera *camera = static_cast<PerspectiveCamera *>(m_context->scene->getCamera());

	Vector d = Vector(camera->getImagePlaneNormal());
	Float fov = std::min(camera->getXFov(), camera->getYFov())*0.9f/2;
	Float distance = bsphere.radius/std::tan(fov * M_PI/180.0f);
	camera->setFocusDepth(distance);

	m_animationTimer->reset();
	m_animationOrigin0 = camera->getPosition();
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

