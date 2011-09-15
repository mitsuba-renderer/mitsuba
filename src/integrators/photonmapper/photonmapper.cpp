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

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/gatherproc.h>
#include "bre.h"

MTS_NAMESPACE_BEGIN

/**
 * Basic two-pass integrator, which visualizes a computed photon map 
 */
class PhotonMapIntegrator : public SampleIntegrator {
public:
	PhotonMapIntegrator(const Properties &props) : SampleIntegrator(props) {
		/* Number of luminaire samples for direct illumination */
		m_directSamples = props.getInteger("directSamples", 16);
		/* Number of BSDF samples when intersecting a glossy material */
		m_glossySamples = props.getInteger("glossySamples", 32);
		/* Depth to start using russian roulette when tracing photons */
		m_rrDepth = props.getInteger("rrDepth", 10);
		/* Longest visualized path length (\c -1 = infinite). 
		   A value of \c 1 will visualize only directly visible light sources.
		   \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
		m_maxDepth = props.getInteger("maxDepth", 5);
		/**
		 * When encountering an ideally specular material within the first few
		 * boundes, the photon mapper places a sample on each lobe
		 * (e.g. reflection *and* transmission). This greatly reduces
		 * variance and is therefore usually worth it
		 */
		m_maxSpecularDepth = props.getInteger("maxSpecularDepth", 4);
		/* Granularity of photon tracing work units (in shot particles, 0 => decide automatically) */
		m_granularity = props.getInteger("granularity", 0);
		/* Number of photons to collect for the global photon map */
		m_globalPhotons = props.getSize("globalPhotons", 200000);
		/* Number of photons to collect for the caustic photon map */
		m_causticPhotons = props.getSize("causticPhotons", 0);
		/* Number of photons to collect for the volumetric photon map */
		m_volumePhotons = props.getSize("volumePhotons", 0);
		/* Max. radius of lookups in the global photon map (relative to the scene size) */
		m_globalLookupRadiusRel = props.getFloat("globalLookupRadius", 0.05f);
		/* Max. radius of lookups in the caustic photon map (relative to the scene size) */
		m_causticLookupRadiusRel = props.getFloat("causticLookupRadius", 0.0125f);
		/* Minimum amount of photons to consider a volumetric photon map lookup valid */
		m_globalLookupSize = props.getInteger("globalLookupSize", 120);
		/* Maximum number of results for caustic photon map lookups */
		m_causticLookupSize = props.getInteger("causticLookupSize", 120);
		/* Approximate number of volume photons to be used in a lookup */
		m_volumeLookupSize = props.getInteger("volumeLookupSize", 120);
		/* Should photon gathering steps exclusively run on the local machine? */
		m_gatherLocally = props.getBoolean("gatherLocally", true);
		/* Indicates if the gathering steps should be canceled if not enough photons are generated. */
		m_autoCancelGathering = props.getBoolean("autoCancelGathering", true);

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
	 : SampleIntegrator(stream, manager) {
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
		SampleIntegrator::serialize(stream, manager);
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
	}

	/// Configure the sampler for a specified amount of direct illumination samples
	void configureSampler(Sampler *sampler) {
		if (m_directSamples > 1)
			sampler->request2DArray(m_directSamples);
		int glossySamples = std::max(m_directSamples, m_glossySamples);
		if (glossySamples > 1)
			sampler->request2DArray(glossySamples);
	}

	void configure() {
		m_invLuminaireSamples = 1.0f / m_directSamples;
		m_invGlossySamples = 1.0f / m_glossySamples;
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, 
			int sceneResID, int cameraResID, int samplerResID) {
		SampleIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);
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
		int qmcSamplerID = sched->registerManifoldResource(samplers); 

		const std::set<Medium *> &media = scene->getMedia();
		for (std::set<Medium *>::const_iterator it = media.begin(); it != media.end(); ++it) {
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
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

			ref<PhotonMap> globalPhotonMap = proc->getPhotonMap();
			if (globalPhotonMap->isFull()) {
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
				m_granularity, 3, m_rrDepth, m_gatherLocally,
				m_autoCancelGathering, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;
	
			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			Log(EDebug, "Caustic photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

			ref<PhotonMap> causticPhotonMap = proc->getPhotonMap();
			if (causticPhotonMap->isFull()) {
				m_causticPhotonMap = causticPhotonMap;
				m_causticPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
				m_causticPhotonMap->build();
				m_causticPhotonMapID = sched->registerResource(m_causticPhotonMap);
			}
		}

		if (m_volumePhotonMap.get() == NULL && m_volumePhotons > 0) {
			/* Generate the volume photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::EVolumePhotons, m_volumePhotons,
				m_granularity, m_maxDepth-1, m_rrDepth, m_gatherLocally,
				m_autoCancelGathering, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			Log(EDebug, "Volume photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

			ref<PhotonMap> volumePhotonMap = proc->getPhotonMap();
			if (volumePhotonMap->isFull()) {
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

		if (getParent() != NULL && getParent()->getClass()->derivesFrom(MTS_CLASS(SampleIntegrator)))
			m_parentIntegrator = static_cast<SampleIntegrator *>(getParent());
		else
			m_parentIntegrator = this;

		return true;
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
	void wakeup(std::map<std::string, SerializableObject *> &params) {
		if (!m_globalPhotonMap.get() && params.find("globalPhotonMap") != params.end())
			m_globalPhotonMap = static_cast<PhotonMap *>(params["globalPhotonMap"]);
		if (!m_causticPhotonMap.get() && params.find("causticPhotonMap") != params.end())
			m_causticPhotonMap = static_cast<PhotonMap *>(params["causticPhotonMap"]);
		if (!m_bre.get() && params.find("bre") != params.end())
			m_bre = static_cast<BeamRadianceEstimator *>(params["bre"]);

		if (getParent() != NULL && getParent()->getClass()->derivesFrom(MTS_CLASS(SampleIntegrator)))
			m_parentIntegrator = static_cast<SampleIntegrator *>(getParent());
		else
			m_parentIntegrator = this;
	}

	void cancel() {
		SampleIntegrator::cancel();
		if (m_proc)
			Scheduler::getInstance()->cancel(m_proc);
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		Spectrum LiSurf(0.0f), LiMedium(0.0f), transmittance(1.0f);
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
		const Scene *scene = rRec.scene;

		bool cacheQuery = (rRec.extra & RadianceQueryRecord::ECacheQuery);
		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		
		if (rRec.medium) {
			Ray mediumRaySegment(ray, 0, its.t);
			transmittance = rRec.medium->getTransmittance(mediumRaySegment);
			mediumRaySegment.mint = ray.mint;
			if (rRec.type & RadianceQueryRecord::EVolumeRadiance &&
					(rRec.depth < m_maxDepth || m_maxDepth < 0) && m_bre.get() != NULL)
				LiMedium = m_bre->query(mediumRaySegment, rRec.medium);
		}

		if (!its.isValid()) {
			/* If no intersection could be found, possibly return 
			   attenuated radiance from a background luminaire */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				LiSurf = scene->LeBackground(ray);
			return LiSurf * transmittance + LiMedium;
		}

		/* Possibly include emitted radiance if requested */
		if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)) 
			LiSurf += its.Le(-ray.d);

		/* Include radiance from a subsurface integrator if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			LiSurf += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);

		if (rRec.depth >= m_maxDepth && m_maxDepth > 0)
			return LiSurf * transmittance + LiMedium;
		
		if (bsdf == NULL) {
			RadianceQueryRecord rRec2;
			rRec2.recursiveQuery(rRec);

			if (its.isMediumTransition())
				rRec2.medium = its.getTargetMedium(ray.d);

			LiSurf += m_parentIntegrator->Li(RayDifferential(its.p, ray.d, ray.time), rRec2);
			return LiSurf * transmittance + LiMedium;
		}

		unsigned int bsdfType = bsdf->getType() & BSDF::EAll;
		bool isDiffuse = (bsdfType == BSDF::EDiffuseReflection);

		if (isDiffuse || cacheQuery) {
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

		if ((bsdfType & BSDF::EDelta) && (bsdfType & ~BSDF::EDelta) == 0 && rRec.depth < m_maxSpecularDepth && !cacheQuery) {
			if (RadianceQueryRecord::EIndirectSurfaceRadiance) {
				int compCount = bsdf->getComponentCount();
				RadianceQueryRecord rRec2;
				for (int i=0; i<compCount; i++) {
					/* Sample the BSDF and recurse */
					BSDFQueryRecord bRec(its, rRec.sampler, ERadiance);
					bRec.component = i;
					Spectrum bsdfVal = bsdf->sample(bRec, Point2(0.0f));
					if (bsdfVal.isZero())
						continue;
					
					rRec2.recursiveQuery(rRec, RadianceQueryRecord::ERadiance);
					RayDifferential bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);
					if (its.isMediumTransition())
						rRec2.medium = its.getTargetMedium(bsdfRay.d);

					LiSurf += bsdfVal * m_parentIntegrator->Li(bsdfRay, rRec2);
				}
			}
		} else if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) {
			/* Estimate the direct illumination if this is requested */
			Point2 *sampleArray;
			Point2 sample;
			int numLuminaireSamples = m_directSamples,
				   numBSDFSamples;

			Float weightLum, weightBSDF;
	
			if (rRec.depth > 1 || cacheQuery || adaptiveQuery) {
				/* This integrator is used recursively by another integrator.
				   Be less accurate as this sample will not directly be observed. */
				numBSDFSamples = numLuminaireSamples = 1;
				weightLum = weightBSDF = 1.0f;
			} else {
				if (isDiffuse) {
					numBSDFSamples = m_directSamples;
					weightBSDF = weightLum = m_invLuminaireSamples;
				} else {
					numBSDFSamples = m_glossySamples;
					weightLum = m_invLuminaireSamples;
					weightBSDF = m_invGlossySamples;
				}
			}
	
			if (numLuminaireSamples > 1) {
				sampleArray = rRec.sampler->next2DArray(m_directSamples);
			} else {
				sample = rRec.nextSample2D(); sampleArray = &sample;
			}
	
			for (int i=0; i<numLuminaireSamples; ++i) {
				/* Estimate the direct illumination if this is requested */
				if (scene->sampleAttenuatedLuminaire(its, rRec.medium, 
						lRec, sampleArray[i], rRec.sampler)) {
					/* Allocate a record for querying the BSDF */
					BSDFQueryRecord bRec(its, its.toLocal(-lRec.d));
	
					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);
	
					if (!bsdfVal.isZero()) {
						/* Calculate prob. of having sampled that direction
							using BSDF sampling */
						Float bsdfPdf = (lRec.luminaire->isIntersectable() 
								|| lRec.luminaire->isBackgroundLuminaire()) ? 
							bsdf->pdf(bRec) : 0;
	
						/* Weight using the power heuristic */
						const Float weight = miWeight(lRec.pdf * numLuminaireSamples, 
								bsdfPdf * numBSDFSamples) * weightLum;
						LiSurf += lRec.value * bsdfVal * weight;
					}
				}
			}
	
			/* ==================================================================== */
			/*                            BSDF sampling                             */
			/* ==================================================================== */
	
			if (numBSDFSamples > 1) {
				sampleArray = rRec.sampler->next2DArray(
					std::max(m_directSamples, m_glossySamples));
			} else {
				sample = rRec.nextSample2D(); sampleArray = &sample;
			}

			RadianceQueryRecord rRec2;
			Intersection &bsdfIts = rRec2.its;

			for (int i=0; i<numBSDFSamples; ++i) {
				/* Sample BSDF * cos(theta) */
				BSDFQueryRecord bRec(its, rRec.sampler, ERadiance);
				Float bsdfPdf;
				Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
				if (bsdfVal.isZero())
					continue;
	
				/* Trace a ray in this direction */
				RayDifferential bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);

				rRec2.recursiveQuery(rRec, 
					RadianceQueryRecord::ERadianceNoEmission);
					
				if (its.isMediumTransition())
					rRec2.medium = its.getTargetMedium(bsdfRay.d);

				bool indexMatchedMediumTransition = false;
				Spectrum transmittance;
				scene->attenuatedRayIntersect(bsdfRay, rRec2.medium, bsdfIts, 
						indexMatchedMediumTransition, transmittance, rRec.sampler);
				rRec2.type ^= RadianceQueryRecord::EIntersection;

				bool hitLightSource = false;
				if (bsdfIts.isValid()) {
					/* Intersected something - check if it was a luminaire */
					if (bsdfIts.isLuminaire()) {
						lRec = LuminaireSamplingRecord(bsdfIts, -bsdfRay.d);
						lRec.value = bsdfIts.Le(-bsdfRay.d);
						hitLightSource = true;
					}
				} else {
					/* No intersection found. Possibly, there is a background
					   luminaire such as an environment map? */
					if (scene->hasBackgroundLuminaire()) {
						lRec.luminaire = scene->getBackgroundLuminaire();
						lRec.d = -bsdfRay.d;
						lRec.value = lRec.luminaire->Le(bsdfRay);
						hitLightSource = true;
					}
				}

				if (hitLightSource) {
					const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ? 
						scene->pdfLuminaire(its.p, lRec) : 0;
			
					const Float weight = miWeight(bsdfPdf * numBSDFSamples, 
						lumPdf * numLuminaireSamples) * weightBSDF;
					LiSurf += lRec.value * bsdfVal * weight * transmittance;
				}

				if (!isDiffuse && (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance) && !cacheQuery) {
					if (indexMatchedMediumTransition) {
						/* The previous ray intersection code passed through an index-matched
						   medium transition while looking for a luminaire. For the recursion,
						   we need to rewind and account for this transition -- therefore,
						   another ray intersection call is neccessary */
						scene->rayIntersect(bsdfRay, bsdfIts);
					}	

					LiSurf += bsdfVal * m_parentIntegrator->Li(bsdfRay, rRec2) * weightBSDF;
				}
			}
		} else if (!isDiffuse && rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance && !cacheQuery) {
			int numBSDFSamples = (rRec.depth > 1 || adaptiveQuery) ? 1 : m_glossySamples;
			Float weightBSDF;
			Point2 *sampleArray;
			Point2 sample;

			if (numBSDFSamples > 1) {
				sampleArray = rRec.sampler->next2DArray(
					std::max(m_directSamples, m_glossySamples));
				weightBSDF = m_invGlossySamples;
			} else {
				sample = rRec.nextSample2D(); sampleArray = &sample;
				weightBSDF = 1.0f;
			}

			RadianceQueryRecord rRec2;
			for (int i=0; i<numBSDFSamples; ++i) {
				/* Sample BSDF * cos(theta) */
				BSDFQueryRecord bRec(its, rRec.sampler, ERadiance);
				Float bsdfPdf;
				Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
				if (bsdfVal.isZero())
					continue;
				rRec2.recursiveQuery(rRec, 
					RadianceQueryRecord::ERadianceNoEmission);

				RayDifferential bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);
				LiSurf += bsdfVal * m_parentIntegrator->Li(bsdfRay, rRec2) * weightBSDF;
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
	SampleIntegrator *m_parentIntegrator;
	int m_globalPhotonMapID, m_causticPhotonMapID, m_breID;
	size_t m_globalPhotons, m_causticPhotons, m_volumePhotons;
	int m_globalLookupSize, m_causticLookupSize, m_volumeLookupSize;
	Float m_globalLookupRadiusRel, m_globalLookupRadius;
	Float m_causticLookupRadiusRel, m_causticLookupRadius;
	Float m_invLuminaireSamples, m_invGlossySamples;
	int m_granularity, m_directSamples, m_glossySamples;
	int m_rrDepth, m_maxDepth, m_maxSpecularDepth;
	bool m_gatherLocally, m_autoCancelGathering;
};

MTS_IMPLEMENT_CLASS_S(PhotonMapIntegrator, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(PhotonMapIntegrator, "Photon map integrator");
MTS_NAMESPACE_END
