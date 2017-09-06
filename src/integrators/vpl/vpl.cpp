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

#include <mitsuba/core/statistics.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/hw/session.h>
#include <mitsuba/hw/device.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{vpl}{Virtual Point Light integrator}
 * \order{14}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{2} will lead to direct-only illumination.
 *         \default{\code{5}}
 *     }
 *     \parameter{shadowMap\showbreak Resolution}{\Integer}{
 *       Resolution of the shadow maps that are used
 *       to compute the point-to-point visibility \default{512}
 *     }
 *     \parameter{clamping}{\Float}{
 *       A relative clamping factor between $[0,1]$ that is
 *       used to control the rendering artifact discussed below.
 *       \default{0.1}
 *     }
 * }
 *
 * This integrator implements a hardware-accelerated global illumination
 * rendering technique based on the Instant Radiosity method by Keller
 * \cite{Keller1997Instant}. This is the same approach that is also used in
 * Mitsuba's real-time preview; the reason for providing it as a separate
 * integrator plugin is to enable automated (e.g. scripted) usage.
 *
 * The method roughly works as follows: during a pre-process pass, any present direct
 * and indirect illumination is converted into a set of \emph{virtual point light}
 * sources (VPLs). The scene is then separately rendered many times, each time using
 * a different VPL as a source of illumination. All of the renderings created in this
 * manner are accumulated to create the final output image.
 *
 * Because the individual rendering steps can be exectuted on a
 * graphics card, it is possible to render many (i.e. 100-1000) VPLs
 * per second. The method is not without problems, however. In particular,
 * it performs poorly when rendering glossy materials, and it produces
 * artifacts in corners and creases . Mitsuba automatically limits
 * the ``glossyness'' of materials to reduce the effects of the former
 * problem. A \code{clamping} parameter is provided to control the latter
 * (see the figure below). The number of samples per pixel specified to
 * the sampler is interpreted as the number of VPLs that should be rendered.
 *
 * \renderings{
 *     \rendering{\code{clamping=0}: With clamping fully disabled, bright
 *     blotches appear in corners and creases.}{integrator_vpl_clamping0}
 *     \rendering{\code{clamping=0.3}: Higher clamping factors remove these
 *     artifacts, but they lead to visible energy loss (the rendering
 *     is too dark in certain areas). The default of \code{0.1} is
 *     usually reasonable.}{integrator_vpl_clamping03}
 * }
 */
class VPLIntegrator : public Integrator {
public:
    VPLIntegrator(const Properties &props) : Integrator(props) {
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
        const ProjectiveCamera *sensor = static_cast<const ProjectiveCamera *>(scene->getSensor());

        Point2 aaSample = Point2(m_random->nextFloat(), m_random->nextFloat());
        Point2 apertureSample(0.5f);
        if (sensor->needsApertureSample())
            apertureSample = Point2(m_random->nextFloat(), m_random->nextFloat());

        Transform projTransform = sensor->getProjectionTransform(apertureSample, aaSample);
        Transform worldTransform = sensor->getWorldTransform()->eval(
            sensor->getShutterOpen() +
            m_random->nextFloat() * sensor->getShutterOpenTime());
        m_shaderManager->setVPL(vpl);
        m_framebuffer->activateTarget();
        m_framebuffer->clear();
        m_renderer->setCamera(projTransform.getMatrix(), worldTransform.getInverseMatrix());
        m_shaderManager->drawAllGeometryForVPL(vpl, sensor);
        m_shaderManager->drawBackground(sensor, projTransform, vpl.emitterScale);
        m_framebuffer->releaseTarget();
    }

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

        if (!(scene->getSensor()->getType() & Sensor::EProjectiveCamera))
            Log(EError, "The VPL integrator requires a projective camera "
                "(e.g. perspective/thinlens/orthographic/telecentric)!");

        m_vpls.clear();
        size_t sampleCount = scene->getSampler()->getSampleCount();
        Float normalization = (Float) 1 / generateVPLs(scene, m_random,
                0, sampleCount, m_maxDepth, true, m_vpls);
        for (size_t i=0; i<m_vpls.size(); ++i) {
            m_vpls[i].P *= normalization;
            m_vpls[i].emitterScale *= normalization;
        }
        Log(EInfo, "Generated %i virtual point lights", m_vpls.size());

        return true;
    }

    void cancel() {
        m_cancel = true;
    }

    bool render(Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID, int samplerResID) {
        ref<Sensor> sensor = scene->getSensor();
        ref<Film> film = sensor->getFilm();
        m_cancel = false;

        if (!sensor->getClass()->derivesFrom(MTS_CLASS(ProjectiveCamera)))
            Log(EError, "The VPL renderer requires a projective camera!");

        /* Initialize hardware rendering */
        m_framebuffer = m_renderer->createGPUTexture("Framebuffer", NULL);
        m_framebuffer->setFrameBufferType(GPUTexture::EColorBuffer);
        m_framebuffer->setComponentFormat(GPUTexture::EFloat32);
        m_framebuffer->setPixelFormat(GPUTexture::ERGB);
        m_framebuffer->setSize(Point3i(film->getSize().x, film->getSize().y, 1));
        m_framebuffer->setFilterType(GPUTexture::ENearest);
        m_framebuffer->setMipMapped(false);

        m_accumBuffer = m_renderer->createGPUTexture("Accumulation buffer",
            new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, film->getSize()));
        m_accumBuffer->setFrameBufferType(GPUTexture::EColorBuffer);
        m_framebuffer->setComponentFormat(GPUTexture::EFloat32);
        m_framebuffer->setPixelFormat(GPUTexture::ERGB);
        m_accumBuffer->setMipMapped(false);

        MTS_AUTORELEASE_BEGIN()
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

        m_shaderManager = new VPLShaderManager(m_renderer);
        m_shaderManager->setShadowMapResolution(m_shadowMapResolution);
        m_shaderManager->setClamping(m_clamping);
        m_shaderManager->init();
        m_shaderManager->setScene(scene);

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
        progress.finish();

        m_accumBuffer->download();
        film->setBitmap(m_accumBuffer->getBitmap());

        m_shaderManager->cleanup();
        m_shaderManager = NULL;

        m_framebuffer->cleanup();
        m_accumBuffer->cleanup();
        m_renderer->shutdown();
        m_device->shutdown();
        m_session->shutdown();
        MTS_AUTORELEASE_END()
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
    int m_maxDepth;
    int m_shadowMapResolution;
    Float m_clamping;
    bool m_cancel;
};

MTS_IMPLEMENT_CLASS(VPLIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(VPLIntegrator, "VPL-based integrator");
MTS_NAMESPACE_END
