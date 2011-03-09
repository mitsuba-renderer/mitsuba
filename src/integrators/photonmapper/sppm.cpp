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
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/renderqueue.h>
#if !defined(__OSX__) && defined(_OPENMP)
#include <omp.h>
#endif

MTS_NAMESPACE_BEGIN

/**
 * Stochastic progressive photon mapping implementation. Only handles surface
 * interactions. Parallelization is limited to the local cores.
 */
class StochasticProgressivePhotonMapIntegrator : public Integrator {
public:
	struct GatherPoint {
		Intersection its;
		Float radius;
		Spectrum weight;
		Spectrum flux;
		Spectrum emission;
		Float N;
		int depth;
		Point2i pos;

		inline GatherPoint() : weight(0.0f), flux(0.0f), emission(0.0f), N(0.0f) {
		}
	};

	StochasticProgressivePhotonMapIntegrator(const Properties &props) : Integrator(props) {
		/* Initial photon query radius (0 = infer based on scene size and camera resolution) */
		m_initialRadius = props.getFloat("initialRadius", 0);
		/* Alpha parameter from the paper (influences the speed, at which the photon radius is reduced) */
		m_alpha = props.getFloat("alpha", .7);
		/* Number of photons to shoot in each iteration */
		m_photonCount = props.getInteger("photonCount", 250000);
		/* Granularity of the work units used in parallelizing the 
		   particle tracing task (default: 100 samples). */
		m_granularity = props.getInteger("granularity", 500);
		/* Longest visualized path length (<tt>-1</tt>=infinite). When a positive value is
		   specified, it must be greater or equal to <tt>2</tt>, which corresponds to single-bounce
		   (direct-only) illumination */
		m_maxDepth = props.getInteger("maxDepth", 5);
		/* Depth to start using russian roulette */
		m_rrDepth = props.getInteger("rrDepth", 3);
		/* Block size used to parallelize the photon query passes (default: 32x32 pixels). */
		m_blockSize = props.getInteger("blockSize", 32);
		m_mutex = new Mutex();
#if defined(__OSX__)
		Log(EError, "Stochastic progressive photon mapping currently doesn't work "
			"on OSX due to a bug in OpenMP that affects Leopard & Snow Leopard");
#endif
		if (m_maxDepth <= 1) 
			Log(EError, "Maximum depth must be set to \"2\" or higher!");
	}

	StochasticProgressivePhotonMapIntegrator(Stream *stream, InstanceManager *manager)
	 : Integrator(stream, manager) {
	}
		

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		Log(EError, "Network rendering is not supported!");
	}

	void cancel() {
		m_running = false;
	}


	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

		if (m_initialRadius == 0) {
			/* Guess an initial radius if not provided
			  (use scene width / horizontal or vertical pixel count) * 5 */
			Float rad = scene->getBSphere().radius;
			Vector2i filmSize = scene->getCamera()->getFilm()->getSize();

			m_initialRadius = std::min(rad / filmSize.x, rad / filmSize.y) * 5;
		}
		return true;
	}

	bool render(Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID, int unused) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Camera> camera = scene->getCamera();
		ref<Film> film = camera->getFilm();
		size_t nCores = sched->getCoreCount();
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SSE_STR ") ..", 
			film->getCropSize().x, film->getCropSize().y, 
			nCores, nCores == 1 ? "core" : "cores");

		Vector2i cropSize = film->getCropSize();
		Point2i cropOffset = film->getCropOffset();

		m_gatherBlocks.clear();
		m_running = true;
		m_totalEmitted = 0;

		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));

		/* Allocate memory */
		m_bitmap = new Bitmap(film->getSize().x, film->getSize().y, 128);
		m_bitmap->clear();
		for (int yofs=0; yofs<cropSize.y; yofs += m_blockSize) {
			for (int xofs=0; xofs<cropSize.x; xofs += m_blockSize) {
				m_gatherBlocks.push_back(std::vector<GatherPoint>());
				m_offset.push_back(Point2i(cropOffset.x + xofs, cropOffset.y + yofs));
				std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[m_gatherBlocks.size()-1];
				int nPixels = std::min(m_blockSize, cropSize.y-yofs)
							* std::min(m_blockSize, cropSize.x-xofs);
				gatherPoints.resize(nPixels);
				for (int i=0; i<nPixels; ++i)
					gatherPoints[i].radius = m_initialRadius;
			}
		}

		/* Create a sampler instance for every core */
		std::vector<SerializableObject *> samplers(sched->getCoreCount());
		for (size_t i=0; i<sched->getCoreCount(); ++i) {
			ref<Sampler> clonedSampler = sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		int samplerResID = sched->registerManifoldResource(
			static_cast<std::vector<SerializableObject*> &>(samplers)); 

#if !defined(__OSX__) && defined(_OPENMP)
		omp_set_num_threads(nCores);
#endif

#ifdef MTS_DEBUG_FP
		enableFPExceptions();
#endif

		int it=0;
		while (m_running) { 
			distributedRTPass(scene, samplers);
			photonMapPass(++it, queue, job, film, sceneResID, 
					cameraResID, samplerResID);
		}

#ifdef MTS_DEBUG_FP
		disableFPExceptions();
#endif

		for (size_t i=0; i<sched->getCoreCount(); ++i)
			samplers[i]->decRef();
		sched->unregisterResource(samplerResID);
		return true;
	}

	void distributedRTPass(Scene *scene, std::vector<SerializableObject *> &samplers) {
		ref<Camera> camera = scene->getCamera();
		bool needsLensSample = camera->needsLensSample();
		bool needsTimeSample = camera->needsTimeSample();
		ref<Film> film = camera->getFilm();
		Vector2i cropSize = film->getCropSize();
		Point2i cropOffset = film->getCropOffset();

		/* Process the image in parallel using blocks for better memory locality */
		Log(EInfo, "Creating %i gather points", cropSize.x*cropSize.y);
		#pragma omp parallel for schedule(dynamic)
		for (int i=-1; i<(int) m_gatherBlocks.size(); ++i) {
			std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[i];
#if !defined(__OSX__) && defined(_OPENMP)
			Sampler *sampler = static_cast<Sampler *>(samplers[omp_get_thread_num()]);
#else
			Sampler *sampler = static_cast<Sampler *>(samplers[0]);
#endif
			int xofs = m_offset[i].x, yofs = m_offset[i].y;
			int index = 0;
			for (int yofsInt = 0; yofsInt < m_blockSize; ++yofsInt) {
				if (yofsInt + yofs - cropOffset.y >= cropSize.y)
					continue;
				for (int xofsInt = 0; xofsInt < m_blockSize; ++xofsInt) {
					if (xofsInt + xofs - cropOffset.x >= cropSize.x)
						continue;
					Point2 lensSample, sample;
					Float timeSample = 0.0f;
					GatherPoint &gatherPoint = gatherPoints[index++];
					sampler->generate();
					if (needsLensSample)
						lensSample = sampler->next2D();
					if (needsTimeSample)
						timeSample = sampler->next1D();
					gatherPoint.pos = Point2i(xofs + xofsInt, yofs + yofsInt);
					sample = sampler->next2D();
					sample += Vector2((Float) gatherPoint.pos.x, (Float) gatherPoint.pos.y);
					RayDifferential ray;
					camera->generateRayDifferential(sample, lensSample, timeSample, ray);
					Spectrum weight(1.0f);
					int depth = 1;

					while (true) {
						if (depth > m_maxDepth) {
							gatherPoint.depth = -1;
							break;
						}
						if (scene->rayIntersect(ray, gatherPoint.its)) {
							const BSDF *bsdf = gatherPoint.its.shape->getBSDF();
							/* Create hit point if this is a diffuse material or a glossy
							   one, and there has been a previous interaction with
							   a glossy material */
							if (bsdf->getType() == BSDF::EDiffuseReflection || 
								bsdf->getType() == BSDF::EDiffuseTransmission) {
								gatherPoint.weight = weight;
								gatherPoint.depth = depth;
								if (gatherPoint.its.isLuminaire())
									gatherPoint.emission = gatherPoint.its.Le(-ray.d);
								else
									gatherPoint.emission = Spectrum(0.0f);
								break;
							} else {
								/* Recurse for dielectric materials and (specific to SPPM):
								   recursive "final gathering" for glossy materials */
								BSDFQueryRecord bRec(gatherPoint.its);
								weight *= bsdf->sampleCos(bRec, sampler->next2D());
								if (weight.isZero()) {
									gatherPoint.depth = -1;
									break;
								}
								ray = RayDifferential(gatherPoint.its.p, 
									gatherPoint.its.toWorld(bRec.wo), ray.time);
								++depth;
							}
						} else {
							/* Generate an invalid sample */
							gatherPoint.depth = -1;
							break;
						}
					}
					sampler->advance();
				}
			}
		}
	}

	void photonMapPass(int it, RenderQueue *queue, const RenderJob *job,  
			Film *film, int sceneResID, int cameraResID, int samplerResID) {
		Log(EInfo, "Performing a photon mapping pass %i", it);
		ref<Scheduler> sched = Scheduler::getInstance();

		/* Generate the global photon map */
		ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
			GatherPhotonProcess::EAllSurfacePhotons, m_photonCount,
			m_granularity, m_maxDepth-1, m_rrDepth, job);

		proc->bindResource("scene", sceneResID);
		proc->bindResource("camera", cameraResID);
		proc->bindResource("sampler", samplerResID);

		sched->schedule(proc);
		sched->wait(proc);

		ref<PhotonMap> photonMap = proc->getPhotonMap();
		photonMap->balance();
		Log(EInfo, "Global photon map full. Shot " SIZE_T_FMT " photons, excess due to parallelism: " 
			SIZE_T_FMT, proc->getShotPhotons(), proc->getExcess());

		Log(EInfo, "Gathering ..");
		m_totalEmitted += proc->getShotPhotons();
		film->clear();
		#pragma omp parallel for schedule(dynamic)
		for (int blockIdx = 0; blockIdx<(int) m_gatherBlocks.size(); ++blockIdx) {
			std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

			float *target = m_bitmap->getFloatData();
			for (size_t i=0; i<gatherPoints.size(); ++i) {
				GatherPoint &gp = gatherPoints[i];
				Float M, N = gp.N;
				Spectrum flux, contrib;

				if (gp.depth != -1) {
					M = photonMap->estimateRadianceRaw(
						gp.its, gp.radius, flux, m_maxDepth-gp.depth);
				} else {
					M = 0;
					flux = Spectrum(0.0f);
				}

				if (N+M == 0) {
					gp.flux = contrib = Spectrum(0.0f);
				} else {
					Float ratio = (N + m_alpha * M) / (N + M);
					gp.flux = (gp.flux + gp.weight * (flux + 
						gp.emission * proc->getShotPhotons() * M_PI * gp.radius*gp.radius)) * ratio;
					gp.radius = gp.radius * std::sqrt(ratio);
					gp.N = N + m_alpha * M;
					contrib = gp.flux / ((Float) m_totalEmitted * gp.radius*gp.radius * M_PI);
				}

				Float r, g, b;
				contrib.toLinearRGB(r, g, b);
				int pos = (gp.pos.y * m_bitmap->getWidth() + gp.pos.x)*4;
				target[pos + 0] = r; target[pos + 1] = g;
				target[pos + 2] = b; target[pos + 3] = 1;
			}
		}
		film->fromBitmap(m_bitmap);
		queue->signalRefresh(job, NULL);
	}


	std::string toString() const {
		return "StochasticProgressivePhotonMapIntegrator[]";
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<std::vector<GatherPoint> > m_gatherBlocks;
	std::vector<Point2i> m_offset;
	ref<Mutex> m_mutex;
	ref<Bitmap> m_bitmap;
	Float m_initialRadius, m_alpha;
	int m_photonCount, m_granularity;
	int m_maxDepth, m_rrDepth, m_totalEmitted;
	int m_blockSize;
	bool m_running;
};

MTS_IMPLEMENT_CLASS_S(StochasticProgressivePhotonMapIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(StochasticProgressivePhotonMapIntegrator, "Stochastic progressive photon mapper");
MTS_NAMESPACE_END
