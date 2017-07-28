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

#include <mitsuba/core/plugin.h>
#include "irrcache_proc.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{irrcache}{Irradiance caching integrator}
 * \order{15}
 * \parameters{
 *     \parameter{resolution}{\Integer}{Elevational resolution of the stratified
 *      final gather hemisphere. The azimuthal resolution is two times this value. \default{14, i.e. $2\cdot14^2$=392 samples in total}}
 *     \parameter{quality}{\Float}{Quality factor (the $\kappa$ parameter of
 *     Tabellion et al. \cite{Tabellion2004Approximate})\default{1.0, which is adequate for most cases}}
 *     \parameter{gradients}{\Boolean}{Use irradiance gradients \cite{Ward1992Irradiance}?\default{\code{true}}}
 *     \parameter{clampNeighbor}{\Boolean}{Use neighbor clamping \cite{Krivanek2006Making}?\default{\code{true}}}
 *     \parameter{clampScreen}{\Boolean}{Use a screen-space clamping criterion \cite{Tabellion2004Approximate}? \default{\code{true}}}
 *     \parameter{overture}{\Boolean}{Do an overture pass before starting the main rendering process?
 *      Usually a good idea.\default{\code{true}}}
 *     \parameter{quality\showbreak Adjustment}{\Float}{When an overture pass is used, Mitsuba subsequently reduces
 *      the quality parameter by this amount to interpolate amongst more samples, creating a visually
 *      smoother result. \default{0.5}}
 *     \parameter{indirectOnly}{\Boolean}{Only show the indirect illumination? This can be useful to check
 *      the interpolation quality. \default{\code{false}}}
 *     \parameter{debug}{\Boolean}{Visualize the sample placement? \default{\code{false}}}
 * }
 * \renderings{
 *  \unframedbigrendering{Illustration of the effect of the different optimizatations
 *   that are provided by this plugin}{integrator_irrcache}
 * }
 * This ``meta-integrator'' implements \emph{irradiance caching} by Ward and Heckbert
 * \cite{Ward1988Ray}. This method computes and caches irradiance information at a sparse
 * set of scene locations and efficiently determines approximate values at other
 * locations using interpolation.
 *
 * This plugin only provides the caching and interpolation part---another plugin
 * is still needed to do the actual computation of irradiance values at cache points.
 * This is done using nesting, e.g. as follows:
 *
 * \begin{xml}[caption={Instantiation of a photon mapper with irradiance caching}]
 * <integrator type="irrcache">
 *     <integrator type="photonmapper"/>
 * </integrator>
 * \end{xml}
 *
 * When a radiance query involves a non-diffuse material, all computation
 * is forwarded to the sub-integrator, i.e. \pluginref{irrcache} is passive.
 * Otherwise, the existing cache points are interpolated to approximate the
 * emitted radiance, or a new cache point is created if the resulting accuracy
 * would be too low.
 * By default, this integrator also performs a distributed overture pass before
 * rendering, which is recommended to avoid artifacts resulting from the
 * addition of samples as rendering proceeds.
 *
 * Note that wrapping an integrator into \pluginref{irrcache} adds one extra
 * light bounce. For instance, the method resulting from using \pluginref{direct}
 * in an irradiance cache renders two-bounce direct illumination.
 *
 * The generality of this implementation allows it to be used in conjunction
 * with photon mapping (the most likely application) as well as all other
 * sampling-based integrators in Mitsuba. Several optimizations are used to
 * improve the achieved interpolation quality, namely irradiance gradients
 * \cite{Ward1992Irradiance}, neighbor clamping \cite{Krivanek2006Making}, a screen-space
 * clamping metric and an improved error function \cite{Tabellion2004Approximate}.
 */

class IrradianceCacheIntegrator : public SamplingIntegrator {
public:
    IrradianceCacheIntegrator(const Properties &props) : SamplingIntegrator(props) {
        /* Elevational resolution of the stratified final gather hemisphere.
           The azimuthal resolution is three times this value. Default:
           2*14^2=392 samples */
        m_resolution = props.getInteger("resolution", 14);
        /* If set to true, the irradiance cache will be filled by a
           parallel overture pass before the main rendering process starts.
           This is strongly recommended. */
        m_overture = props.getBoolean("overture", true);
        /* Quality setting (\kappa in the [Tabellion et al.] paper).
           A value of 1 should be adequate in most cases. */
        m_quality = props.getFloat("quality", 1.0f);
        /* Multiplicative factor for the quality parameter following an
           overture pass. This can be used to interpolate amongst more
           samples, creating a visually smoother result. Must be
           1 or less. */
        m_qualityAdjustment = props.getFloat("qualityAdjustment", .5f);
        /* If set to true, sample locations will be visually highlighted */
        m_debug = props.getBoolean("debug", false);
        /* Should irradiance gradients be used? Generally, this will
           significantly improve the interpolation quality.*/
        m_gradients = props.getBoolean("gradients", true);
        /* Should neighbor clamping [Krivanek et al.] be used? This
           propagates geometry information amongst close-by samples
           and generally leads to better sample placement. */
        m_clampNeighbor = props.getBoolean("clampNeighbor", true);
        /* If set to true, the influence region of samples will be clamped
           using the screen-space metric by [Tabellion et al.]?
           Turning this off may lead to excessive sample placement. */
        m_clampScreen = props.getBoolean("clampScreen", true);
        /* If set to true, direct illumination will be suppressed -
           useful for checking the interpolation quality */
        m_indirectOnly = props.getBoolean("indirectOnly", false);

        if (m_debug)
            m_overture = false;

        Assert(m_qualityAdjustment > 0 && m_qualityAdjustment <= 1);
    }

    IrradianceCacheIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_irrCache = static_cast<IrradianceCache *>(manager->getInstance(stream));
        m_subIntegrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
        m_resolution = stream->readInt();
        m_quality = stream->readFloat();
        m_qualityAdjustment = stream->readFloat();
        m_diffScaleFactor = stream->readFloat();
        m_clampScreen = stream->readBool();
        m_clampNeighbor = stream->readBool();
        m_overture = stream->readBool();
        m_gradients = stream->readBool();
        m_debug = stream->readBool();
        m_indirectOnly = stream->readBool();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        manager->serialize(stream, m_irrCache.get());
        manager->serialize(stream, m_subIntegrator.get());
        stream->writeInt(m_resolution);
        stream->writeFloat(m_quality);
        stream->writeFloat(m_qualityAdjustment);
        stream->writeFloat(m_diffScaleFactor);
        stream->writeBool(m_clampScreen);
        stream->writeBool(m_clampNeighbor);
        stream->writeBool(m_overture);
        stream->writeBool(m_gradients);
        stream->writeBool(m_debug);
        stream->writeBool(m_indirectOnly);
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        m_subIntegrator->configureSampler(scene, sampler);
        m_diffScaleFactor = std::sqrt((Float) sampler->getSampleCount());
    }

    void bindUsedResources(ParallelProcess *proc) const {
        m_subIntegrator->bindUsedResources(proc);
    }

    void wakeup(ConfigurableObject *parent,
            std::map<std::string, SerializableObject *> &params) {
        m_subIntegrator->wakeup(this, params);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();

        if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
            if (!cClass->derivesFrom(MTS_CLASS(SamplingIntegrator)))
                Log(EError, "The sub-integrator must be derived from the class SamplingIntegrator");
            m_subIntegrator = static_cast<SamplingIntegrator *>(child);
            m_subIntegrator->setParent(this);
        } else {
            Integrator::addChild(name, child);
        }
    }

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;

        if (m_subIntegrator == NULL)
            Log(EError, "No sub-integrator was specified!");

        if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;

        ref<Scheduler> sched = Scheduler::getInstance();
        m_irrCache = new IrradianceCache(scene->getAABB());
        m_irrCache->clampNeighbor(m_clampNeighbor);
        m_irrCache->clampScreen(m_clampScreen);
        m_irrCache->useGradients(m_gradients);
        m_irrCache->setQuality(m_quality);

        std::string irrCacheStatus;
        if (m_overture)
            irrCacheStatus += "overture, ";
        if (m_debug)
            irrCacheStatus += "debug, ";
        if (m_gradients)
            irrCacheStatus += "gradients, ";
        if (m_clampNeighbor)
            irrCacheStatus += "clampNeighbor, ";
        if (m_clampScreen)
            irrCacheStatus += "clampScreen, ";

        Log(EDebug, "Irradiance cache status : %s", irrCacheStatus.c_str());
        Log(EDebug, "  - Gather resolution   : %ix%i = %i samples", m_resolution, 2*m_resolution, 2*m_resolution*m_resolution);
        Log(EDebug, "  - Quality setting     : %.2f (adjustment: %.2f)", m_quality, m_qualityAdjustment);

        if (m_overture) {
            int subIntegratorResID = sched->registerResource(m_subIntegrator);
            ref<OvertureProcess> proc = new OvertureProcess(job, m_resolution, m_gradients,
                m_clampNeighbor, m_clampScreen, m_quality);
            m_proc = proc;
            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("subIntegrator", subIntegratorResID);
            bindUsedResources(proc);
            sched->schedule(proc);
            sched->unregisterResource(subIntegratorResID);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
                Log(EWarn, "The overture pass did not complete sucessfully!");
                return false;
            }

            ref<const IrradianceRecordVector> vec = proc->getSamples();
            Log(EDebug, "Overture pass generated %i irradiance samples", vec->size());
            for (size_t i=0; i<vec->size(); ++i)
                m_irrCache->insert(new IrradianceCache::Record((*vec)[i]));

            m_irrCache->setQuality(m_quality * m_qualityAdjustment);
        }
        return true;
    }

    void cancel() {
        if (m_proc) {
            Scheduler::getInstance()->cancel(m_proc);
        } else {
            SamplingIntegrator::cancel();
            m_subIntegrator->cancel();
        }
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        Intersection &its = rRec.its;

        if (m_indirectOnly && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance))
            rRec.type ^= RadianceQueryRecord::EDirectSurfaceRadiance;

        if (rRec.rayIntersect(ray)) {
            const BSDF *bsdf = its.getBSDF(ray);

            if (bsdf && (bsdf->getType() & BSDF::EAll) == BSDF::EDiffuseReflection &&
                    (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
                Spectrum E;

                if (Frame::cosTheta(its.wi) <= 0) {
                    E = Spectrum(0.0f);
                } else if (!m_irrCache->get(its, E)) {
                    handleMiss(ray, rRec, E);

                    if (m_debug)
                        E.fromLinearRGB(1e3, 0, 0);
                }

                rRec.type ^= RadianceQueryRecord::EIndirectSurfaceRadiance;

                return E * bsdf->getDiffuseReflectance(its) * INV_PI +
                    m_subIntegrator->Li(ray, rRec);
            }
        }
        return m_subIntegrator->Li(ray, rRec);
    }

    void handleMiss(RayDifferential ray, const RadianceQueryRecord &rRec,
            Spectrum &E) const {
        /* Handle an irradiance cache miss */
        HemisphereSampler *hs = m_hemisphereSampler.get();
        Sampler *sampler = m_sampleGenerator.get();
        RadianceQueryRecord rRec2;
        if (hs == NULL) {
            Properties props("independent");
            sampler = static_cast<Sampler *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(Sampler), props));
            hs = new HemisphereSampler(m_resolution, 2 * m_resolution);
            m_hemisphereSampler.set(hs);
            m_sampleGenerator.set(sampler);
        }

        /* Generate stratified cosine-weighted samples and compute
           rotational + translational gradients */
        hs->generateDirections(rRec.its);
        sampler->generate(Point2i(0));

        for (unsigned int j=0; j<hs->getM(); j++) {
            for (unsigned int k=0; k<hs->getN(); k++) {
                HemisphereSampler::SampleEntry &entry = (*hs)(j, k);
                    entry.dist = std::numeric_limits<Float>::infinity();
                rRec2.recursiveQuery(rRec,
                    RadianceQueryRecord::ERadianceNoEmission | RadianceQueryRecord::EDistance);
                rRec2.extra = 1;
                rRec2.sampler = sampler;
                entry.L = m_subIntegrator->Li(RayDifferential(rRec.its.p, entry.d, ray.time), rRec2);
                entry.dist = rRec2.dist;
                sampler->advance();
            }
        }

        hs->process(rRec.its);
        /* Undo ray differential scaling done by the integrator */
        if (ray.hasDifferentials)
            ray.scaleDifferential(m_diffScaleFactor);
        m_irrCache->put(ray, rRec.its, *hs);
        E = hs->getIrradiance();
    }

    Spectrum E(const Scene *scene, const Intersection &its, const Medium *medium,
            Sampler *sampler, int nSamples, bool handleIndirect) const {
        Spectrum EDir(0.0f), EIndir(0.0f);
        DirectSamplingRecord dRec(its);

        /* Sample the direct illumination component */
        for (int i=0; i<nSamples; i++) {
            int maxIntermediateInteractions = -1;
            Spectrum directRadiance = scene->sampleAttenuatedEmitterDirect(
                dRec, its, medium, maxIntermediateInteractions, sampler->next2D());

            if (!directRadiance.isZero()) {
                Float dp = dot(dRec.d, its.shFrame.n);
                if (dp > 0)
                    EDir += directRadiance * dp;
            }
        }

        if (handleIndirect) {
            RadianceQueryRecord rRec(scene, sampler);
            rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission, medium);
            rRec.its = its;
            if (!m_irrCache->get(rRec.its, EIndir))
                handleMiss(RayDifferential(), rRec, EIndir);
        }

        return (EDir / (Float) nSamples) + EIndir;
    }

    const Integrator *getSubIntegrator(int idx) const {
        if (idx != 0)
            return NULL;
        return m_subIntegrator.get();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "IrradianceCacheIntegrator[" << endl
            << "  subIntegrator = " << indent(m_subIntegrator->toString()) << "," << endl
            << "  resolution = " << m_resolution << "," << endl
            << "  irrCache = " << indent(m_irrCache->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    mutable ThreadLocal<HemisphereSampler> m_hemisphereSampler;
    mutable ThreadLocal<Sampler> m_sampleGenerator;
    mutable ref<IrradianceCache> m_irrCache;
    ref<SamplingIntegrator> m_subIntegrator;
    ref<ParallelProcess> m_proc;
    Float m_quality, m_qualityAdjustment, m_diffScaleFactor;
    bool m_clampScreen, m_clampNeighbor;
    bool m_overture, m_gradients, m_debug, m_indirectOnly;
    int m_resolution;
};

MTS_IMPLEMENT_CLASS_S(IrradianceCacheIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(IrradianceCacheIntegrator, "Irradiance cache");
MTS_NAMESPACE_END
