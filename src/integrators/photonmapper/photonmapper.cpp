/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/gatherproc.h>

MTS_NAMESPACE_BEGIN

/**
* Parallel photon gathering, SSE-accelerated lookups, RGBE encoding
* Numbers for caustic/volume are ignored when no specular objects or participating
* media are present.
*/
class PhotonMapIntegrator : public SampleIntegrator {
public:
	PhotonMapIntegrator(const Properties &props) : SampleIntegrator(props) {
		/* Number of luminaire samples for direct illumination */
		m_directSamples = props.getInteger("directSamples", 1);
		/* Number of BSDF samples when intersecting a glossy material */
		m_glossySamples = props.getInteger("glossySamples", 32);
		/* Depth to start using russian roulette when tracing photons */
		m_rrDepth = props.getInteger("rrDepth", 10);
		/* Depth cutoff when tracing photons */
		m_maxDepth = props.getInteger("maxDepth", 40);
		/* Depth cutoff when recursively tracing specular materials */
		m_maxSpecularDepth = props.getInteger("maxSpecularDepth", 6);
		/* Granularity of photon tracing work units (in shot particles) */
		m_granularity = props.getInteger("granularity", 1000);

		/* Number of photons to collect for the global photon map */
		m_globalPhotons = (size_t) props.getLong("globalPhotons", 200000);
		/* Number of photons to collect for the caustic photon map */
		m_causticPhotons = (size_t) props.getLong("causticPhotons", 200000);
		/* Number of photons to collect for the volumetric photon map */
		m_volumePhotons = (size_t) props.getLong("volumePhotons", 200000);
		/* Radius of lookups in the global photon map (relative to the scene size) */
		m_globalLookupRadiusRel = props.getFloat("globalLookupRadius", 0.05f);
		/* Radius of lookups in the caustic photon map (relative to the scene size) */
		m_causticLookupRadiusRel = props.getFloat("causticLookupRadius", 0.0125f);
		/* Radius of lookups in the volumetric photon map (relative to the scene size) */
		m_volumeLookupRadiusRel = props.getFloat("volumeLookupRadius", 0.05f);
		/* Minimum amount of photons to consider a global photon map lookup valid */
		m_globalMinPhotons = props.getInteger("globalMinPhotons", 8);
		/* Minimum amount of photons to consider a caustic photon map lookup valid */
		m_causticMinPhotons = props.getInteger("causticMinPhotons", 100);
		/* Minimum amount of photons to consider a volumetric photon map lookup valid */
		m_volumeMinPhotons = props.getInteger("volumeMinPhotons", 8);
		/* Maximum number of results for global photon map lookups */
		m_globalLookupSize = props.getInteger("globalLookupSize", 200);
		/* Maximum number of results for caustic photon map lookups */
		m_causticLookupSize = props.getInteger("causticLookupSize", 200);
		/* Maximum number of results for volumetric photon map lookups */
		m_volumeLookupSize = props.getInteger("volumeLookupSize", 200);
	}

	/// Unserialize from a binary data stream
	PhotonMapIntegrator(Stream *stream, InstanceManager *manager)
	 : SampleIntegrator(stream, manager) {
		m_directSamples = stream->readInt();
		m_glossySamples = stream->readInt();
		m_maxSpecularDepth = stream->readInt();
		m_globalPhotons = (size_t) stream->readULong();
		m_causticPhotons = (size_t) stream->readULong();
		m_volumePhotons = (size_t) stream->readULong();
		m_globalLookupRadius = stream->readFloat();
		m_causticLookupRadius = stream->readFloat();
		m_volumeLookupRadius = stream->readFloat();
		m_globalMinPhotons = stream->readInt();
		m_causticMinPhotons = stream->readInt();
		m_volumeMinPhotons = stream->readInt();
		m_globalLookupSize = stream->readInt();
		m_causticLookupSize = stream->readInt();
		m_volumeLookupSize = stream->readInt();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SampleIntegrator::serialize(stream, manager);
		stream->writeInt(m_directSamples);
		stream->writeInt(m_glossySamples);
		stream->writeInt(m_maxSpecularDepth);
		stream->writeULong(m_globalPhotons);
		stream->writeULong(m_causticPhotons);
		stream->writeULong(m_volumePhotons);
		stream->writeFloat(m_globalLookupRadius);
		stream->writeFloat(m_causticLookupRadius);
		stream->writeFloat(m_volumeLookupRadius);
		stream->writeInt(m_globalMinPhotons);
		stream->writeInt(m_causticMinPhotons);
		stream->writeInt(m_volumeMinPhotons);
		stream->writeInt(m_globalLookupSize);
		stream->writeInt(m_causticLookupSize);
		stream->writeInt(m_volumeLookupSize);
	}

	/// Configure the sampler for a specified amount of direct illumination samples
	void configureSampler(Sampler *sampler) {
		if (m_directSamples > 1)
			sampler->request2DArray(m_directSamples);
		sampler->request2DArray(m_glossySamples);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, 
			int sceneResID, int cameraResID, int samplerResID) {
		SampleIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);
		/* Create a deterministic sampler for the photon gathering step */
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(Sampler::m_theClass, Properties("halton")));
		int qmcSamplerID = sched->registerResource(sampler);

		/* Don't create a caustic photon map if the scene does not contain specular materials */
		const std::vector<Shape *> &shapes = scene->getShapes();
		bool foundSpecular = false;
		for (size_t i=0; i<shapes.size(); ++i) {
			if (shapes[i]->getBSDF()->getType() & BSDF::EDelta) {
				foundSpecular = true;
				break;
			}
		}
		if (!foundSpecular)
			m_causticPhotons = 0;

		/* Don't create a volumetric photon map if there are no participating media */
		if (!scene->hasMedia())
			m_volumePhotons = 0;

		if (m_globalPhotonMap.get() == NULL && m_globalPhotons > 0) {
			/* Adapt to scene extents */
			m_globalLookupRadius = m_globalLookupRadiusRel * scene->getBSphere().radius;

			/* Generate the global photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::ESurfacePhotons, m_globalPhotons,
				m_granularity, m_maxDepth, m_rrDepth, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("camera", cameraResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			m_globalPhotonMap = proc->getPhotonMap();
			m_globalPhotonMap->setScale(1 / (Float) proc->getShotPhotons());
			m_globalPhotonMap->setMinPhotons(m_globalMinPhotons);
			m_globalPhotonMap->balance();
			Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " photons, excess due to parallelism: " 
				SIZE_T_FMT, proc->getShotPhotons(), proc->getExcess());
			m_globalPhotonMapID = sched->registerResource(m_globalPhotonMap);
		}

		if (m_causticPhotonMap.get() == NULL && m_causticPhotons > 0) {
			/* Adapt to scene extents */
			m_causticLookupRadius = m_causticLookupRadiusRel * scene->getBSphere().radius;

			/* Generate the caustic photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::ECausticPhotons, m_causticPhotons,
				m_granularity, 2, m_rrDepth, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;
	
			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			m_causticPhotonMap = proc->getPhotonMap();
			m_causticPhotonMap->setScale(1 / (Float) proc->getShotPhotons());
			m_causticPhotonMap->setMinPhotons(m_causticMinPhotons);
			m_causticPhotonMap->balance();
			Log(EDebug, "Caustic photon map - excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getExcess());
			m_causticPhotonMapID = sched->registerResource(m_causticPhotonMap);
		}

		if (m_volumePhotonMap.get() == NULL && m_volumePhotons > 0) {
			/* Adapt to scene extents */
			m_volumeLookupRadius = m_volumeLookupRadiusRel * scene->getBSphere().radius;

			/* Generate the volume photon map */
			ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
				GatherPhotonProcess::EVolumePhotons, m_volumePhotons,
				m_granularity, m_maxDepth, m_rrDepth, job);

			proc->bindResource("scene", sceneResID);
			proc->bindResource("sampler", qmcSamplerID);

			m_proc = proc;
			sched->schedule(proc);
			sched->wait(proc);
			m_proc = NULL;

			if (proc->getReturnStatus() != ParallelProcess::ESuccess)
				return false;

			m_volumePhotonMap = proc->getPhotonMap();
			m_volumePhotonMap->setScale(1 / (Float) proc->getShotPhotons());
			m_volumePhotonMap->setMinPhotons(m_volumeMinPhotons);
			m_volumePhotonMap->balance();
			Log(EDebug, "Volumetric photon map - excess photons due to parallelism: " 
				SIZE_T_FMT, proc->getExcess());
			m_volumePhotonMapID = sched->registerResource(m_volumePhotonMap);
		}

		sched->unregisterResource(qmcSamplerID);
		m_parentIntegrator = static_cast<SampleIntegrator *>(getParent());
		return true;
	}

	/// Specify globally shared resources
	void bindUsedResources(ParallelProcess *proc) const {
		if (m_globalPhotonMap.get())
			proc->bindResource("globalPhotonMap", m_globalPhotonMapID);
		if (m_causticPhotonMap.get())
			proc->bindResource("causticPhotonMap", m_causticPhotonMapID);
		if (m_volumePhotonMap.get())
			proc->bindResource("volumePhotonMap", m_volumePhotonMapID);
	}

	/// Connect to globally shared resources
	void wakeup(std::map<std::string, SerializableObject *> &params) {
		if (!m_globalPhotonMap.get() && params.find("globalPhotonMap") != params.end())
			m_globalPhotonMap = static_cast<PhotonMap *>(params["globalPhotonMap"]);
		if (!m_causticPhotonMap.get() && params.find("causticPhotonMap") != params.end())
			m_causticPhotonMap = static_cast<PhotonMap *>(params["causticPhotonMap"]);
		if (!m_volumePhotonMap.get() && params.find("volumetricPhotonMap") != params.end())
			m_volumePhotonMap = static_cast<PhotonMap *>(params["volumetricPhotonMap"]);

		if (getParent() != NULL && getParent()->getClass()->derivesFrom(SampleIntegrator::m_theClass))
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
		Spectrum Li(0.0f), LiVol(0.0f);
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
	
		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);

		if (!its.isValid()) {
			/* If no intersection could be found, possibly return 
			   attenuated radiance from a background luminaire */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				Li += rRec.scene->LeBackground(ray);
			return Li;
		}

		/* Possibly include emitted radiance if requested */
		if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)) 
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface integrator if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(rRec.scene, -ray.d);

		const BSDF *bsdf = its.getBSDF(ray);
		int bsdfType = bsdf->getType();

		Point2 *sampleArray, sample;
		int numDirectSamples = (rRec.depth == 1) ? m_directSamples : 1;
		/* When this integrator is used recursively by another integrator,
			Be less accurate as this sample will not be directly observed. */
		if (numDirectSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numDirectSamples);
		} else {
			sample = rRec.nextSample2D();
			sampleArray = &sample;
		}

		/* Estimate the direct illumination if this is requested */
		if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) {
			Float weight = 1 / (Float) numDirectSamples;

			for (int i=0; i<numDirectSamples; ++i) {
				if (rRec.scene->sampleLuminaire(its, lRec, sampleArray[i])) {
					/* Allocate a record for querying the BSDF */
					const BSDFQueryRecord bRec(rRec, its, its.toLocal(-lRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->fCos(bRec);

					Li += lRec.value * bsdfVal * weight;
				}
			}
		}

//		if (rRec.type & RadianceQueryRecord:EVolumeRadiance) {
			/* Ray marching */
//		}


		if (bsdfType == BSDF::EDiffuseReflection) {
			/* Hit a diffuse material - do a direct photon map visualization. */
			if (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance) 
				Li += m_globalPhotonMap->estimateIrradianceFiltered(its.p, 
					its.shFrame.n, m_globalLookupRadius, m_globalLookupSize)
						* bsdf->getDiffuseReflectance(its) * INV_PI;
			if (rRec.type & RadianceQueryRecord::ECausticRadiance && m_causticPhotonMap.get())
				Li += m_causticPhotonMap->estimateIrradianceFiltered(its.p,
					its.shFrame.n, m_causticLookupRadius, m_causticLookupSize)
							* bsdf->getDiffuseReflectance(its) * INV_PI;
		} else if ((bsdfType & BSDF::EDelta) != 0
				&& (bsdfType & ~BSDF::EDelta) == 0 && !rRec.extra) {
			RadianceQueryRecord rRec2;
			RayDifferential recursiveRay;
			/* Ideal specular material -> recursive ray tracing.
			   Deliberately risk exponential ray growth by spawning
			   several child rays. The impact on the final image is huge
			   and well worth the extra computation. */
			if (rRec.depth+1 < m_maxSpecularDepth) {
				int compCount = bsdf->getComponentCount();
				for (int i=0; i<compCount; i++) {
					/* Sample the BSDF and recurse */
					BSDFQueryRecord bRec(its);
					bRec.component = i;
					Spectrum bsdfVal = bsdf->sampleCos(bRec, Point2(0.0f));
					if (bsdfVal.isZero())
						continue;

					rRec2.recursiveQuery(rRec, RadianceQueryRecord::ERadiance);
					recursiveRay = Ray(its.p, its.toWorld(bRec.wo), ray.time);
					Li += m_parentIntegrator->Li(recursiveRay, rRec2) * bsdfVal;
				}
			}
		} else if (rRec.depth == 1 && (bsdf->getType() & BSDF::EGlossy)) {
			/* Hit a glossy material - MC integration over the hemisphere
			   using BSDF importance sampling */
			sampleArray = rRec.sampler->next2DArray(m_glossySamples);
			RadianceQueryRecord rRec2;
			RayDifferential recursiveRay;
			Float weight = 1 / (Float) m_glossySamples;

			for (int i=0; i<m_glossySamples; ++i) {
				BSDFQueryRecord bRec(rRec, its);
				Spectrum bsdfVal = bsdf->sampleCos(bRec, sampleArray[i]);

				rRec2.recursiveQuery(rRec, RadianceQueryRecord::ERadianceNoEmission);
				recursiveRay = Ray(its.p, its.toWorld(bRec.wo), ray.time);
				Li += m_parentIntegrator->Li(recursiveRay, rRec2) * bsdfVal * weight;
			}
		} else {
			Li += m_globalPhotonMap->estimateRadianceFiltered(its,
				m_globalLookupRadius, m_globalLookupSize);
		}

		return Li * LiVol;
	}
	
	std::string toString() const {
		std::ostringstream oss;
		oss << "PhotonMapIntegrator[" << std::endl
			<< "  directSamples = " << m_directSamples << "," << std::endl
			<< "  glossySamples = " << m_glossySamples << "," << std::endl
			<< "  globalPhotons = " << m_globalPhotons << "," << std::endl
			<< "  causticPhotons = " << m_causticPhotons << "," << std::endl
			<< "  volumePhotons = " << m_volumePhotons << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<PhotonMap> m_globalPhotonMap;
	ref<PhotonMap> m_causticPhotonMap;
	ref<PhotonMap> m_volumePhotonMap;
	ref<ParallelProcess> m_proc;
	SampleIntegrator *m_parentIntegrator;
	int m_globalPhotonMapID, m_causticPhotonMapID, m_volumePhotonMapID;
	size_t m_globalPhotons, m_causticPhotons, m_volumePhotons;
	int m_globalMinPhotons, m_globalLookupSize;
	int m_causticMinPhotons, m_causticLookupSize;
	int m_volumeMinPhotons, m_volumeLookupSize;
	Float m_globalLookupRadiusRel, m_globalLookupRadius;
	Float m_causticLookupRadiusRel, m_causticLookupRadius;
	Float m_volumeLookupRadiusRel, m_volumeLookupRadius;
	int m_granularity;
	int m_directSamples, m_glossySamples;
	int m_rrDepth;
	int m_maxDepth, m_maxSpecularDepth;
};

MTS_IMPLEMENT_CLASS_S(PhotonMapIntegrator, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(PhotonMapIntegrator, "Photon map integrator");
MTS_NAMESPACE_END
