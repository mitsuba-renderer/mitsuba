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

#include <mitsuba/core/statistics.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/hw/session.h>
#include <mitsuba/hw/device.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

/**
 * Rasterization-based global illuminated technique using hardware
 * accelerated renderings of the scene under point source illumination. Based on 
 * "Instant Radiosity" by Alexander Keller in Computer Graphics Proceedings, 
 * Annual Conference Series, SIGGRAPH 97, pp. 49-56. 
 */
class VPLIntegrator : public Integrator {
public:
	VPLIntegrator(const Properties &props) : Integrator(props) {
		/* Number of virtual point lights */
		m_vplCount = props.getInteger("vplCount", 1000);
		/* Shadow map resolution (e.g. 512x512) */
		m_shadowMapResolution = props.getInteger("shadowMapResolution", 512);
		/* Max. depth (expressed as path length) */
		m_maxDepth = props.getInteger("maxDepth", 5);
		/* Relative clamping factor (0=no clamping, 1=full clamping) */
		m_clamping = props.getFloat("clamping", 0.1f);

		m_session = Session::create();
		m_device = Device::create(m_session);
		m_renderer = Renderer::create(m_session);
	
		m_random = new Random();
	}

	/// Draw the full scene using additive blending and shadow maps
	void drawShadowedScene(const Scene *scene, const VPL &vpl) {
	const std::vector<std::pair<const TriMesh *, Transform> > meshes = m_shaderManager->getMeshes();
	const std::vector<std::pair<const TriMesh *, Transform> > transpMeshes = m_shaderManager->getTransparentMeshes();
		const ProjectiveCamera *camera = static_cast<const ProjectiveCamera *>(scene->getCamera());
		Point2 jitter(0.5f - m_random->nextFloat(), 0.5f - m_random->nextFloat());
		Transform projectionTransform = camera->getGLProjectionTransform(jitter);
		m_renderer->setCamera(projectionTransform.getMatrix(), camera->getViewTransform().getMatrix());
		Transform clipToWorld = camera->getViewTransform().inverse() * projectionTransform.inverse();
		Point camPos = scene->getCamera()->getPosition();

		m_renderer->beginDrawingMeshes();
		for (unsigned int j=0; j<meshes.size(); j++) {
			const TriMesh *mesh = meshes[j].first;
			bool hasTransform = !meshes[j].second.isIdentity();
			m_shaderManager->configure(vpl, mesh->getBSDF(), mesh->getLuminaire(),
				camPos, !mesh->hasVertexNormals());
			if (hasTransform)
				m_renderer->pushTransform(meshes[j].second);
			m_renderer->drawTriMesh(mesh);
			if (hasTransform)
				m_renderer->popTransform();
			m_shaderManager->unbind();
		}
		m_renderer->endDrawingMeshes();

		m_shaderManager->drawBackground(clipToWorld, camPos, 1.0f);

		if (transpMeshes.size() > 0) {
			m_renderer->setDepthMask(false);
			m_renderer->setBlendMode(Renderer::EBlendAlpha);
			m_renderer->beginDrawingMeshes();
			for (unsigned int j=0; j<transpMeshes.size(); j++) {
				const TriMesh *mesh = transpMeshes[j].first;
				bool hasTransform = !transpMeshes[j].second.isIdentity();
				m_shaderManager->configure(vpl, mesh->getBSDF(), mesh->getLuminaire(),
					camPos, !mesh->hasVertexNormals());
				if (hasTransform)
					m_renderer->pushTransform(transpMeshes[j].second);
				m_renderer->drawTriMesh(mesh);
				if (hasTransform)
					m_renderer->popTransform();
				m_shaderManager->unbind();
			}
			m_renderer->endDrawingMeshes();
			m_renderer->setDepthMask(true);
			m_renderer->setBlendMode(Renderer::EBlendNone);
		}
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

		if (m_vpls.size() == 0) {
			Float normalization = (Float) 1 / generateVPLs(scene, m_random,
					0, m_vplCount, m_maxDepth, true, m_vpls);
			for (size_t i=0; i<m_vpls.size(); ++i)
				m_vpls[i].P *= normalization;
			Log(EInfo, "Generated %i virtual point lights", m_vpls.size());
		}
		return true;
	}

	void cancel() {
		m_cancel = true;
	}

	bool render(Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID, int samplerResID) {
		ref<Camera> camera = scene->getCamera();
		ref<Film> film = camera->getFilm();
		m_cancel = false;

		if (!camera->getClass()->derivesFrom(MTS_CLASS(ProjectiveCamera)))
			Log(EError, "The VPL renderer requires a projective camera!");

		/* Initialize hardware rendering */
		m_framebuffer = m_renderer->createGPUTexture("Framebuffer", NULL);
		m_framebuffer->setFrameBufferType(GPUTexture::EColorBuffer);
		m_framebuffer->setFormat(GPUTexture::EFloat32RGB);
		m_framebuffer->setSize(Point3i(film->getSize().x, film->getSize().y, 1));
		m_framebuffer->setFilterType(GPUTexture::ENearest);
		m_framebuffer->setMipMapped(false);

		m_accumBuffer = m_renderer->createGPUTexture("Accumulation buffer",
			new Bitmap(film->getSize().x, film->getSize().y, 128));
		m_accumBuffer->setFrameBufferType(GPUTexture::EColorBuffer);
		m_framebuffer->setFormat(GPUTexture::EFloat32RGB);
		m_accumBuffer->setMipMapped(false);

		m_session->init();
		m_device->setSize(film->getSize());
		m_device->init();
		m_device->setVisible(false);
		m_renderer->init(m_device);

		if (!m_renderer->getCapabilities()->isSupported(
			RendererCapabilities::EShadingLanguage))
			Log(EError, "Support for GLSL is required!");
		if (!m_renderer->getCapabilities()->isSupported(
			RendererCapabilities::ERenderToTexture))
			Log(EError, "Render-to-texture support is required!");
		if (!m_renderer->getCapabilities()->isSupported(
			RendererCapabilities::EFloatingPointTextures))
			Log(EError, "Floating point texture support is required!");
		if (!m_renderer->getCapabilities()->isSupported(
			RendererCapabilities::EFloatingPointBuffer))
			Log(EError, "Floating point render buffer support is required!");
		if (!m_renderer->getCapabilities()->isSupported(
			RendererCapabilities::EVertexBufferObjects))
			Log(EError, "Vertex buffer object support is required!");
		if (!m_renderer->getCapabilities()->isSupported(
			RendererCapabilities::EGeometryShaders))
			Log(EError, "Geometry shader support is required!");

		/* Initialize and clear the framebuffer */
		m_framebuffer->init();
		m_accumBuffer->init();
		m_accumBuffer->activateTarget();
		m_accumBuffer->clear();
		m_accumBuffer->releaseTarget();

		m_shaderManager = new VPLShaderManager(scene, m_renderer);
		m_shaderManager->setShadowMapResolution(m_shadowMapResolution);
		m_shaderManager->setClamping(m_clamping);
		m_shaderManager->init();

		ProgressReporter progress("Rendering", m_vpls.size(), job);
		for (size_t i=0; i<m_vpls.size() && !m_cancel; ++i) {
			const VPL &vpl = m_vpls[i];

			m_renderer->setDepthMask(true);
			m_renderer->setDepthTest(true);
			m_renderer->setBlendMode(Renderer::EBlendNone);
			m_shaderManager->setVPL(vpl);

			m_framebuffer->activateTarget();
			m_framebuffer->clear();
			drawShadowedScene(scene, vpl);
			m_framebuffer->releaseTarget();

			m_renderer->setDepthMask(false);
			m_renderer->setDepthTest(false);
			m_renderer->setBlendMode(Renderer::EBlendAdditive);
			m_accumBuffer->activateTarget();
			m_renderer->blitTexture(m_framebuffer, true);
			m_accumBuffer->releaseTarget();

			if ((i%20) == 0) { 
				m_renderer->flush();
				m_renderer->checkError();
			}
			progress.update(i);
		}

		m_accumBuffer->download();
		film->fromBitmap(m_accumBuffer->getBitmap());

		m_shaderManager->cleanup();
		m_shaderManager = NULL;

		m_framebuffer->cleanup();
		m_accumBuffer->cleanup();
		m_renderer->shutdown();
		m_device->shutdown();
		m_session->shutdown();
		return !m_cancel;
	}

	MTS_DECLARE_CLASS()
private:
	ref<Session> m_session;
	ref<Device> m_device;
	ref<Renderer> m_renderer;
	ref<GPUTexture> m_framebuffer, m_accumBuffer;
	ref<VPLShaderManager> m_shaderManager;
	std::deque<VPL> m_vpls;
	ref<Random> m_random;
	int m_vplCount, m_maxDepth;
	int m_shadowMapResolution;
	Float m_clamping;
	bool m_cancel;
};

MTS_IMPLEMENT_CLASS(VPLIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(VPLIntegrator, "VPL-based integrator");
MTS_NAMESPACE_END
