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

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{ao}{Ambient occlusion integrator}
 * \order{0}
 * \parameters{
 *     \parameter{shadingSamples}{\Integer}{Specifies the number of
 *         shading samples that should be computed per primary ray
 *         \default{1}}
 *
 *     \parameter{rayLength}{\Float}{Specifies the world-space length of the
 *         ambient occlusion rays that will be cast. \default{\code{-1}, i.e. automatic}}.
 * }
 * \renderings{
 *    \rendering{A view of the scene on page \pageref{fig:rungholt}, rendered using
 *               the Ambient Occlusion integrator}{integrator_ao}
 *    \rendering{A corresponding rendering created using the standard \pluginref{path} tracer}{integrator_ao_path}
 * }
 * Ambient Occlusion is a simple non-photorealistic rendering technique that simulates the exposure of an object
 * to uniform illumination incident from all direction. It produces approximate shadowing between closeby
 * objects, as well as darkening in corners, creases, and cracks. The scattering models associated with objects
 * in the scene are ignored.
 */

class AmbientOcclusionIntegrator : public SamplingIntegrator {
public:
    AmbientOcclusionIntegrator(const Properties &props) : SamplingIntegrator(props) {
        m_shadingSamples = props.getSize("shadingSamples", 1);
        m_rayLength = props.getFloat("rayLength", -1);
    }

    /// Unserialize from a binary data stream
    AmbientOcclusionIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_shadingSamples = stream->readSize();
        m_rayLength = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeSize(m_shadingSamples);
        stream->writeFloat(m_rayLength);
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);

        if (m_shadingSamples > 1)
            sampler->request2DArray(m_shadingSamples);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);
        if (m_rayLength < 0) {
            m_rayLength = scene->getAABB().getBSphere().radius * 0.5f;
            Log(EInfo, "Setting occlusion ray length to %f", m_rayLength);
        }
        return true;
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        Spectrum Li(0.0f);
        Point2 sample;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        if (!rRec.rayIntersect(ray)) {
            /* If no intersection could be found, possibly return
               radiance from a background emitter */
            return Spectrum(1.0f);
        }

        /* Figure out how many shading samples to take, and where the
           required random numbers should come from */
        Point2 *sampleArray = &sample;
        size_t numShadingSamples = m_shadingSamples;

        bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

        if (numShadingSamples > 1 && rRec.depth == 1 && !adaptiveQuery) {
            sampleArray = rRec.sampler->next2DArray(numShadingSamples);
        } else {
            /* This integrator is used recursively by another integrator.
               Be less accurate as this sample will not directly be observed. */
            numShadingSamples = 1;
            sample = rRec.nextSample2D();
        }

        const Intersection &its = rRec.its;
        for (size_t i=0; i<numShadingSamples; ++i) {
            Vector d = its.toWorld(warp::squareToCosineHemisphere(sampleArray[i]));

            Ray shadowRay(its.p, d, Epsilon, m_rayLength, ray.time);
            if (!rRec.scene->rayIntersect(shadowRay))
                Li += Spectrum(1.0f);
        }

        Li /= static_cast<Float>(numShadingSamples);

        return Li;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "AmbientOcclusionIntegrator[" << endl
            << "  shadingSamples = " << m_shadingSamples << "," << endl
            << "  rayLength = " << m_rayLength << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    size_t m_shadingSamples;
    Float m_rayLength;
};

MTS_IMPLEMENT_CLASS_S(AmbientOcclusionIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(AmbientOcclusionIntegrator, "Ambient occlusion integrator");
MTS_NAMESPACE_END
