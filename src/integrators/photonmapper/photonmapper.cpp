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
#include <mitsuba/render/common.h>
#include <mitsuba/render/gatherproc.h>
#include "bre.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{photonmapper}{Photon map integrator}
 * \order{6}
 * \parameters{
 *     \parameter{directSamples}{\Integer}{Number of samples used for the
 *        direct illumination component \default{16}}
 *     \parameter{glossySamples}{\Integer}{Number of samples used for the indirect
 *        illumination component of glossy materials \default{32}}
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{1} will only render directly visible light sources.
 *         \code{2} will lead to single-bounce (direct-only) illumination,
 *         and so on. \default{\code{-1}}
 *     }
 *     \parameter{globalPhotons}{\Integer}{Number of photons that will be collected for the global photon map\default{250000}}
 *     \parameter{causticPhotons}{\Integer}{Number of photons that will be collected for the caustic photon map\default{250000}}
 *     \parameter{volumePhotons}{\Integer}{Number of photons that will be collected for the volumetric photon map\default{250000}}
 *     \parameter{globalLookup\showbreak Radius}{\Float}{Maximum radius of photon lookups in the global photon map (relative to the scene size)\default{0.05}}
 *     \parameter{causticLookup\showbreak Radius}{\Float}{Maximum radius of photon lookups in the caustic photon map (relative to the scene size)\default{0.0125}}
 *     \parameter{lookupSize}{\Integer}{Number of photons that should be fetched in photon map queries\default{120}}
 *     \parameter{granularity}{\Integer}{
 *      Granularity of photon tracing work units for the purpose
 *      of parallelization (in \# of shot particles) \default{0, i.e. decide automatically}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 *     \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *        which the implementation will start to use the ``russian roulette''
 *        path termination criterion. \default{\code{5}}
 *     }
 * }
 * This plugin implements the two-pass photon mapping algorithm as proposed by Jensen \cite{Jensen1996Global}.
 * The implementation partitions the illumination into three different classes (diffuse, caustic, and volumetric),
 * and builds a separate photon map for each class.
 *
 * Following this, a standard recursive ray tracing pass is started which performs kernel density estimation
 * using these photon maps. Since the photon maps are visualized directly, the result will appear ``blotchy'' (\figref{pmap-blotchy})
 * unless an extremely large number of photons is used. A simple remedy is to combine the photon mapper with
 * an irradiance cache, which performs \emph{final gathering} to remove these artifacts. Due to its caching nature,
 * the rendering process will be faster as well.
 * \begin{xml}[caption={Instantiation of a photon mapper with irradiance caching}]
 * <integrator type="irrcache">
 *     <integrator type="photonmapper"/>
 * </integrator>
 * \end{xml}
 * \renderings{
 *   \rendering{Rendered using plain photon mapping}{integrator_photonmapper_1}
 *   \rendering{Rendered using photon mapping together with irradiance caching}{integrator_photonmapper_2}
 *   \vspace{-2mm}
 *   \caption{\label{fig:pmap-blotchy}Sponza atrium illuminated by a point light and rendered using 5 million photons.
 *   Irradiance caching significantly accelerates the rendering time and eliminates the ``blotchy''
 *   kernel density estimation artifacts. Model courtesy of Marko Dabrovic.}
 * }
 *
 * When the scene contains participating media, the Beam Radiance Estimate \cite{Jarosz2008Beam}
 * by Jarosz et al. is used to estimate the illumination due to volumetric scattering.
 *
 * \remarks{
 *     \item Currently, only homogeneous participating media are supported by this implementation
 * }
 */
class PhotonMapIntegrator : public SamplingIntegrator {
public:
    PhotonMapIntegrator(const Properties &props) : SamplingIntegrator(props),
          m_parentIntegrator(NULL) {
        /* Number of lsamples for direct illumination */
        m_directSamples = props.getInteger("directSamples", 16);
        /* Number of BSDF samples when intersecting a glossy material */
        m_glossySamples = props.getInteger("glossySamples", 32);
        /* Depth to start using russian roulette when tracing photons */
        m_rrDepth = props.getInteger("rrDepth", 5);
        /* Longest visualized path length (\c -1 = infinite).
           A value of \c 1 will visualize only directly visible light sources.
           \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
        m_maxDepth = props.getInteger("maxDepth", -1);
        /**
         * When encountering an ideally specular material, the photon mapper places
         * a sample on each lobe (e.g. reflection *and* transmission). This leads
         * to an exponential growth in running time but greatly reduces variance and
         * is therefore usually worth it. This parameter specifies after how many
         * bounces this behavior should be stopped.
         */
        m_maxSpecularDepth = props.getInteger("maxSpecularDepth", 4);
        /* Granularity of photon tracing work units (in shot particles, 0 => decide automatically) */
        m_granularity = props.getInteger("granularity", 0);
        /* Number of photons to collect for the global photon map */
        m_globalPhotons = props.getSize("globalPhotons", 250000);
        /* Number of photons to collect for the caustic photon map */
        m_causticPhotons = props.getSize("causticPhotons", 250000);
        /* Number of photons to collect for the volumetric photon map */
        m_volumePhotons = props.getSize("volumePhotons", 250000);
        /* Max. radius of lookups in the global photon map (relative to the scene size) */
        m_globalLookupRadiusRel = props.getFloat("globalLookupRadius", 0.05f);
        /* Max. radius of lookups in the caustic photon map (relative to the scene size) */
        m_causticLookupRadiusRel = props.getFloat("causticLookupRadius", 0.0125f);
        /* Minimum amount of photons to consider a photon map lookup valid */
        int lookupSize = props.getInteger("lookupSize", 120);
        /* Minimum amount of photons to consider a volumetric photon map lookup valid */
        m_globalLookupSize = props.getInteger("globalLookupSize", lookupSize);
        /* Maximum number of results for caustic photon map lookups */
        m_causticLookupSize = props.getInteger("causticLookupSize", lookupSize);
        /* Approximate number of volume photons to be used in a lookup */
        m_volumeLookupSize = props.getInteger("volumeLookupSize", lookupSize);
        /* Should photon gathering steps exclusively run on the local machine? */
        m_gatherLocally = props.getBoolean("gatherLocally", true);
        /* Indicates if the gathering steps should be canceled if not enough photons are generated. */
        m_autoCancelGathering = props.getBoolean("autoCancelGathering", true);
        /* When this flag is set to true, contributions from directly
         * visible emitters will not be included in the rendered image */
        m_hideEmitters = props.getBoolean("hideEmitters", false);

        if (m_maxDepth == 0) {
            Log(EError, "maxDepth must be greater than zero!");
        } else if (m_maxDepth == -1) {
            /**
             * An infinite depth is currently not supported, since
             * the photon tracing step uses a Halton sequence
             * that is based on a finite-sized prime number table
             */
            m_maxDepth = 128;
        }

        m_causticPhotonMapID = m_globalPhotonMapID = m_breID = 0;
    }

    /// Unserialize from a binary data stream
    PhotonMapIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager), m_parentIntegrator(NULL) {
        m_directSamples = stream->readInt();
        m_glossySamples = stream->readInt();
        m_maxDepth = stream->readInt();
        m_maxSpecularDepth = stream->readInt();
        m_rrDepth = stream->readInt();
        m_globalPhotons = stream->readSize();
        m_causticPhotons = stream->readSize();
        m_volumePhotons = stream->readSize();
        m_globalLookupRadius = stream->readFloat();
        m_causticLookupRadius = stream->readFloat();
        m_globalLookupSize = stream->readInt();
        m_causticLookupSize = stream->readInt();
        m_volumeLookupSize = stream->readInt();
        m_gatherLocally = stream->readBool();
        m_autoCancelGathering = stream->readBool();
        m_hideEmitters = stream->readBool();
        m_causticPhotonMapID = m_globalPhotonMapID = m_breID = 0;
        configure();
    }

    virtual ~PhotonMapIntegrator() {
        ref<Scheduler> sched = Scheduler::getInstance();
        if (m_globalPhotonMapID)
            sched->unregisterResource(m_globalPhotonMapID);
        if (m_causticPhotonMapID)
            sched->unregisterResource(m_causticPhotonMapID);
        if (m_breID)
            sched->unregisterResource(m_breID);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeInt(m_directSamples);
        stream->writeInt(m_glossySamples);
        stream->writeInt(m_maxDepth);
        stream->writeInt(m_maxSpecularDepth);
        stream->writeInt(m_rrDepth);
        stream->writeSize(m_globalPhotons);
        stream->writeSize(m_causticPhotons);
        stream->writeSize(m_volumePhotons);
        stream->writeFloat(m_globalLookupRadius);
        stream->writeFloat(m_causticLookupRadius);
        stream->writeInt(m_globalLookupSize);
        stream->writeInt(m_causticLookupSize);
        stream->writeInt(m_volumeLookupSize);
        stream->writeBool(m_gatherLocally);
        stream->writeBool(m_autoCancelGathering);
        stream->writeBool(m_hideEmitters);
    }

    /// Configure the sampler for a specified amount of direct illumination samples
    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        if (m_directSamples > 1)
            sampler->request2DArray(m_directSamples);
        int bsdfSamples = std::max(m_directSamples, m_glossySamples);
        if (bsdfSamples > 1)
            sampler->request2DArray(bsdfSamples);

        bool hasDelta = false;
        const ref_vector<Shape> &shapes = scene->getShapes();
        for (size_t i=0; i<shapes.size(); ++i) {
            const BSDF *bsdf = shapes[i]->getBSDF();
            if (bsdf && bsdf->getType() & BSDF::EDelta)
                hasDelta = true;
        }

        if (!hasDelta)
            m_causticPhotons = 0;
    }

    void configure() {
        m_invEmitterSamples = 1.0f / m_directSamples;
        m_invGlossySamples = 1.0f / m_glossySamples;
    }

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
        /* Create a deterministic sampler for the photon gathering step */
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("halton")));
        /* Create a sampler instance for every core */
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }
        int qmcSamplerID = sched->registerMultiResource(samplers);
        for (size_t i=0; i<samplers.size(); ++i)
            samplers[i]->decRef();

        const ref_vector<Medium> &media = scene->getMedia();
        for (ref_vector<Medium>::const_iterator it = media.begin(); it != media.end(); ++it) {
            if (!(*it)->isHomogeneous())
                Log(EError, "Inhomogeneous media are currently not supported by the photon mapper!");
        }

        if (m_globalPhotonMap.get() == NULL && m_globalPhotons > 0) {
            /* Generate the global photon map */
            ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
                GatherPhotonProcess::ESurfacePhotons, m_globalPhotons,
                m_granularity, m_maxDepth-1, m_rrDepth, m_gatherLocally,
                m_autoCancelGathering, job);

            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", qmcSamplerID);

            m_proc = proc;
            sched->schedule(proc);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess)
                return false;

            ref<PhotonMap> globalPhotonMap = proc->getPhotonMap();
            if (globalPhotonMap->isFull()) {
                Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
                    SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

                m_globalPhotonMap = globalPhotonMap;
                m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
                m_globalPhotonMap->build();
                m_globalPhotonMapID = sched->registerResource(m_globalPhotonMap);
            }
        }

        if (m_causticPhotonMap.get() == NULL && m_causticPhotons > 0) {
            /* Generate the caustic photon map */
            ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
                GatherPhotonProcess::ECausticPhotons, m_causticPhotons,
                m_granularity, m_maxDepth-1, m_rrDepth, m_gatherLocally,
                m_autoCancelGathering, job);

            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", qmcSamplerID);

            m_proc = proc;
            sched->schedule(proc);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess)
                return false;

            ref<PhotonMap> causticPhotonMap = proc->getPhotonMap();
            if (causticPhotonMap->isFull()) {
                Log(EDebug, "Caustic photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
                    SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

                m_causticPhotonMap = causticPhotonMap;
                m_causticPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
                m_causticPhotonMap->build();
                m_causticPhotonMapID = sched->registerResource(m_causticPhotonMap);
            }
        }

        size_t volumePhotons = scene->getMedia().size() == 0 ? 0 : m_volumePhotons;
        if (m_volumePhotonMap.get() == NULL && volumePhotons > 0) {
            /* Generate the volume photon map */
            ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
                GatherPhotonProcess::EVolumePhotons, volumePhotons,
                m_granularity, m_maxDepth-1, m_rrDepth, m_gatherLocally,
                m_autoCancelGathering, job);

            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", qmcSamplerID);

            m_proc = proc;
            sched->schedule(proc);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess)
                return false;

            ref<PhotonMap> volumePhotonMap = proc->getPhotonMap();
            if (volumePhotonMap->isFull()) {
                Log(EDebug, "Volume photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
                    SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

                volumePhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
                volumePhotonMap->build();
                m_bre = new BeamRadianceEstimator(volumePhotonMap, m_volumeLookupSize);
                m_breID = sched->registerResource(m_bre);
            }
        }

        /* Adapt to scene extents */
        m_globalLookupRadius = m_globalLookupRadiusRel * scene->getBSphere().radius;
        m_causticLookupRadius = m_causticLookupRadiusRel * scene->getBSphere().radius;

        sched->unregisterResource(qmcSamplerID);

        return true;
    }

    void setParent(ConfigurableObject *parent) {
        if (parent->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)))
            m_parentIntegrator = static_cast<SamplingIntegrator *>(parent);
        else
            m_parentIntegrator = this;
    }

    /// Specify globally shared resources
    void bindUsedResources(ParallelProcess *proc) const {
        if (m_globalPhotonMap.get())
            proc->bindResource("globalPhotonMap", m_globalPhotonMapID);
        if (m_causticPhotonMap.get())
            proc->bindResource("causticPhotonMap", m_causticPhotonMapID);
        if (m_bre.get())
            proc->bindResource("bre", m_breID);
    }

    /// Connect to globally shared resources
    void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
        if (!m_globalPhotonMap.get() && params.find("globalPhotonMap") != params.end())
            m_globalPhotonMap = static_cast<PhotonMap *>(params["globalPhotonMap"]);
        if (!m_causticPhotonMap.get() && params.find("causticPhotonMap") != params.end())
            m_causticPhotonMap = static_cast<PhotonMap *>(params["causticPhotonMap"]);
        if (!m_bre.get() && params.find("bre") != params.end())
            m_bre = static_cast<BeamRadianceEstimator *>(params["bre"]);

        if (parent && parent->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)))
            m_parentIntegrator = static_cast<SamplingIntegrator *>(parent);
        else
            m_parentIntegrator = this;
    }

    void cancel() {
        SamplingIntegrator::cancel();
        if (m_proc)
            Scheduler::getInstance()->cancel(m_proc);
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        Spectrum LiSurf(0.0f), LiMedium(0.0f), transmittance(1.0f);
        Intersection &its = rRec.its;
        const Scene *scene = rRec.scene;

        bool cacheQuery = (rRec.extra & RadianceQueryRecord::ECacheQuery);
        bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);

        if (rRec.medium) {
            Ray mediumRaySegment(ray, 0, its.t);
            transmittance = rRec.medium->evalTransmittance(mediumRaySegment);
            mediumRaySegment.mint = ray.mint;
            if (rRec.type & RadianceQueryRecord::EVolumeRadiance &&
                    (rRec.depth < m_maxDepth || m_maxDepth < 0) && m_bre.get() != NULL)
                LiMedium = m_bre->query(mediumRaySegment, rRec.medium);
        }

        if (!its.isValid()) {
            /* If no intersection could be found, possibly return
               attenuated radiance from a background luminaire */
            if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
                LiSurf = scene->evalEnvironment(ray);
            return LiSurf * transmittance + LiMedium;
        }

        /* Possibly include emitted radiance if requested */
        if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
            LiSurf += its.Le(-ray.d);

        /* Include radiance from a subsurface scattering model if requested */
        if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
            LiSurf += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

        const BSDF *bsdf = its.getBSDF(ray);

        if (rRec.depth >= m_maxDepth && m_maxDepth > 0)
            return LiSurf * transmittance + LiMedium;

        unsigned int bsdfType = bsdf->getType() & BSDF::EAll;

        /* Irradiance cache query -> treat as diffuse */
        bool isDiffuse = (bsdfType == BSDF::EDiffuseReflection) || cacheQuery;

        bool hasSpecular = bsdfType & BSDF::EDelta;

        /* Exhaustively recurse into all specular lobes? */
        bool exhaustiveSpecular = rRec.depth < m_maxSpecularDepth && !cacheQuery;

        if (isDiffuse && (dot(its.shFrame.n, ray.d) < 0 || (bsdf->getType() & BSDF::EBackSide))) {
            /* 1. Diffuse indirect */
            int maxDepth = m_maxDepth == -1 ? INT_MAX : (m_maxDepth-rRec.depth);
            if (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance && m_globalPhotonMap.get())
                LiSurf += m_globalPhotonMap->estimateIrradiance(its.p,
                    its.shFrame.n, m_globalLookupRadius, maxDepth,
                    m_globalLookupSize) * bsdf->getDiffuseReflectance(its) * INV_PI;
            if (rRec.type & RadianceQueryRecord::ECausticRadiance && m_causticPhotonMap.get())
                LiSurf += m_causticPhotonMap->estimateIrradiance(its.p,
                    its.shFrame.n, m_causticLookupRadius, maxDepth,
                    m_causticLookupSize) * bsdf->getDiffuseReflectance(its) * INV_PI;
        }

        if (hasSpecular && exhaustiveSpecular
            && (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
            /* 1. Specular indirect */
            int compCount = bsdf->getComponentCount();
            RadianceQueryRecord rRec2;
            for (int i=0; i<compCount; i++) {
                unsigned int type = bsdf->getType(i);
                if (!(type & BSDF::EDelta))
                    continue;
                /* Sample the BSDF and recurse */
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                bRec.component = i;
                Spectrum bsdfVal = bsdf->sample(bRec, Point2(0.5f));
                if (bsdfVal.isZero())
                    continue;

                rRec2.recursiveQuery(rRec, RadianceQueryRecord::ERadiance);
                RayDifferential bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);
                if (its.isMediumTransition())
                    rRec2.medium = its.getTargetMedium(bsdfRay.d);

                LiSurf += bsdfVal * m_parentIntegrator->Li(bsdfRay, rRec2);
            }
        }

        /* Estimate the direct illumination if this is requested */
        int numEmitterSamples = m_directSamples, numBSDFSamples;
        Float weightLum, weightBSDF;
        Point2 *sampleArray;
        Point2 sample;

        if (rRec.depth > 1 || cacheQuery || adaptiveQuery) {
            /* This integrator is used recursively by another integrator.
               Be less accurate as this sample will not directly be observed. */
            numBSDFSamples = numEmitterSamples = 1;
            weightLum = weightBSDF = 1.0f;
        } else {
            if (isDiffuse) {
                numBSDFSamples = m_directSamples;
                weightBSDF = weightLum = m_invEmitterSamples;
            } else {
                numBSDFSamples = m_glossySamples;
                weightLum = m_invEmitterSamples;
                weightBSDF = m_invGlossySamples;
            }
        }

        if ((bsdfType & BSDF::ESmooth) && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
            DirectSamplingRecord dRec(its);

            if (numEmitterSamples > 1) {
                sampleArray = rRec.sampler->next2DArray(m_directSamples);
            } else {
                sample = rRec.nextSample2D(); sampleArray = &sample;
            }

            for (int i=0; i<numEmitterSamples; ++i) {
                int interactions = m_maxDepth - rRec.depth - 1;
                Spectrum value = scene->sampleAttenuatedEmitterDirect(
                        dRec, its, rRec.medium, interactions,
                        sampleArray[i], rRec.sampler);

                /* Estimate the direct illumination if this is requested */
                if (!value.isZero()) {
                    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    if (!bsdfVal.isZero()) {
                        /* Calculate prob. of having sampled that direction
                           using BSDF sampling */

                        if (!hasSpecular || exhaustiveSpecular)
                            bRec.typeMask = BSDF::ESmooth;

                        Float bsdfPdf = (emitter->isOnSurface()
                                && dRec.measure == ESolidAngle
                                && interactions == 0)
                                ? bsdf->pdf(bRec) : (Float) 0.0f;

                        /* Weight using the power heuristic */
                        const Float weight = miWeight(dRec.pdf * numEmitterSamples,
                                bsdfPdf * numBSDFSamples) * weightLum;

                        LiSurf += value * bsdfVal * weight;
                    }
                }
            }
        }

        /* ==================================================================== */
        /*                            BSDF sampling                             */
        /* ==================================================================== */

        /* Sample direct compontent via BSDF sampling if this is generally requested AND
             the BSDF is smooth, or there is a delta component that was not handled by the
             exhaustive sampling loop above */
        bool bsdfSampleDirect = (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
            ((bsdfType & BSDF::ESmooth) || (hasSpecular && !exhaustiveSpecular));

        /* Sample indirect component via BSDF sampling if this is generally requested AND
            the BSDF is non-diffuse (diffuse is handled by the global photon map)
            or there is a delta component that was not handled by the exhaustive sampling loop
            above. */
        bool bsdfSampleIndirect = (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance) &&
            !isDiffuse && ((bsdfType & BSDF::ESmooth) || (hasSpecular && !exhaustiveSpecular));

        if (bsdfSampleDirect || bsdfSampleIndirect) {
            if (numBSDFSamples > 1) {
                sampleArray = rRec.sampler->next2DArray(
                    std::max(m_directSamples, m_glossySamples));
            } else {
                sample = rRec.nextSample2D(); sampleArray = &sample;
            }

            RadianceQueryRecord rRec2;
            Intersection &bsdfIts = rRec2.its;

            DirectSamplingRecord dRec(its);
            for (int i=0; i<numBSDFSamples; ++i) {
                /* Sample BSDF * cos(theta) */
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                if (!hasSpecular || exhaustiveSpecular)
                    bRec.typeMask = BSDF::ESmooth;

                Float bsdfPdf;
                Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
                if (bsdfVal.isZero())
                    continue;

                /* Trace a ray in this direction */
                RayDifferential bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);

                Spectrum value;
                bool hitEmitter = false;
                if (scene->rayIntersect(bsdfRay, bsdfIts)) {
                    /* Intersected something - check if it was a luminaire */
                    if (bsdfIts.isEmitter() && bsdfSampleDirect) {
                        value = bsdfIts.Le(-bsdfRay.d);
                        dRec.setQuery(bsdfRay, bsdfIts);
                        hitEmitter = true;
                    }
                } else if (bsdfSampleDirect) {
                    /* Intersected nothing -- perhaps there is an environment map? */
                    const Emitter *env = scene->getEnvironmentEmitter();

                    if (env) {
                        value = env->evalEnvironment(bsdfRay);
                        if (env->fillDirectSamplingRecord(dRec, bsdfRay))
                            hitEmitter = true;
                    }
                }

                if (hitEmitter) {
                    const Float emitterPdf = scene->pdfEmitterDirect(dRec);

                    Spectrum transmittance = rRec2.medium ?
                        rRec2.medium->evalTransmittance(Ray(bsdfRay, 0, bsdfIts.t)) : Spectrum(1.0f);

                    const Float weight = miWeight(bsdfPdf * numBSDFSamples,
                        emitterPdf * numEmitterSamples) * weightBSDF;

                    LiSurf += value * bsdfVal * weight * transmittance;
                }

                /* Recurse */
                if (bsdfSampleIndirect) {
                    rRec2.recursiveQuery(rRec,
                        RadianceQueryRecord::ERadianceNoEmission);
                    rRec2.type ^= RadianceQueryRecord::EIntersection;

                    if (its.isMediumTransition())
                        rRec2.medium = its.getTargetMedium(bsdfRay.d);

                    LiSurf += bsdfVal * m_parentIntegrator->Li(bsdfRay, rRec2) * weightBSDF;
                }
            }
        }

        return LiSurf * transmittance + LiMedium;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "PhotonMapIntegrator[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  maxSpecularDepth = " << m_maxSpecularDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  directSamples = " << m_directSamples << "," << endl
            << "  glossySamples = " << m_glossySamples << "," << endl
            << "  globalPhotons = " << m_globalPhotons << "," << endl
            << "  causticPhotons = " << m_causticPhotons << "," << endl
            << "  volumePhotons = " << m_volumePhotons << "," << endl
            << "  gatherLocally = " << m_gatherLocally << "," << endl
            << "  globalLookupRadius = " << m_globalLookupRadius << "," << endl
            << "  causticLookupRadius = " << m_causticLookupRadius << "," << endl
            << "  globalLookupSize = " << m_globalLookupSize << "," << endl
            << "  causticLookupSize = " << m_causticLookupSize << "," << endl
            << "  volumeLookupSize = " << m_volumeLookupSize << endl
            << "]";
        return oss.str();
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA; pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    MTS_DECLARE_CLASS()
private:
    ref<PhotonMap> m_globalPhotonMap;
    ref<PhotonMap> m_causticPhotonMap;
    ref<PhotonMap> m_volumePhotonMap;
    ref<BeamRadianceEstimator> m_bre;
    ref<ParallelProcess> m_proc;
    SamplingIntegrator *m_parentIntegrator;
    int m_globalPhotonMapID, m_causticPhotonMapID, m_breID;
    size_t m_globalPhotons, m_causticPhotons, m_volumePhotons;
    int m_globalLookupSize, m_causticLookupSize, m_volumeLookupSize;
    Float m_globalLookupRadiusRel, m_globalLookupRadius;
    Float m_causticLookupRadiusRel, m_causticLookupRadius;
    Float m_invEmitterSamples, m_invGlossySamples;
    int m_granularity, m_directSamples, m_glossySamples;
    int m_rrDepth, m_maxDepth, m_maxSpecularDepth;
    bool m_gatherLocally, m_autoCancelGathering;
    bool m_hideEmitters;
};

MTS_IMPLEMENT_CLASS_S(PhotonMapIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(PhotonMapIntegrator, "Photon map integrator");
MTS_NAMESPACE_END
