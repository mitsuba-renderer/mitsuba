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
#include <mitsuba/render/renderqueue.h>

MTS_NAMESPACE_BEGIN

/**
 * Progressive photon mapping implementation. Only handles surface
 * interactions. Parallelization is limited to the local cores.
 */
class ProgressivePhotonMapIntegrator : public Integrator {
public:
	struct GatherPoint {
		Intersection its;
		Float radius;
		Spectrum weight;
		Spectrum flux;
		Spectrum emission;
		Point2 sample;
		Float N;
		int depth;

		inline GatherPoint() : weight(0.0f), flux(0.0f), emission(0.0f), N(0.0f) {
		}
	};

	ProgressivePhotonMapIntegrator(const Properties &props) : Integrator(props) {
		/* Initial photon query radius (0 = infer based on scene size and camera resolution) */
		m_initialRadius = props.getFloat("initialRadius", 0);
		/* Alpha parameter from the paper (influences the speed, at which the photon radius is reduced) */
		m_alpha = props.getFloat("alpha", .7);
		/* Number of photons to shoot in each iteration */
		m_photonCount = props.getInteger("photonCount", 100000);
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
		/* Indicates if the gathering steps should be canceled if not enough photons are generated. */
		m_autoCancelGathering = props.getBoolean("autoCancelGathering", true);
		m_mutex = new Mutex();
#if MTS_BROKEN_OPENMP == 1
		Log(EError, "Progressive photon mapping currently doesn't work "
			"on OSX due to a bug in OpenMP that affects Leopard & Snow Leopard");
#endif
		if (m_maxDepth <= 1 && m_maxDepth != -1) 
			Log(EError, "Maximum depth must either be set to \"-1\" or \"2\" or higher!");
	}

	ProgressivePhotonMapIntegrator(Stream *stream, InstanceManager *manager)
	 : Integrator(stream, manager) {
	}
		
	virtual ~ProgressivePhotonMapIntegrator() {
		for (size_t i=0; i<m_blocks.size(); ++i)
			m_blocks[i]->decRef();
		m_blocks.clear();
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
			  (scene width / horizontal or vertical pixel count) * 5 */
			Float rad = scene->getBSphere().radius;
			Vector2i filmSize = scene->getCamera()->getFilm()->getSize();

			m_initialRadius = std::min(rad / filmSize.x, rad / filmSize.y) * 5;
		}
		return true;
	}

	bool render(Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Camera> camera = scene->getCamera();
		ref<Film> film = camera->getFilm();
		size_t nCores = sched->getCoreCount();
		Sampler *cameraSampler = (Sampler *) sched->getResource(samplerResID, 0);
	
		size_t sampleCount = cameraSampler->getSampleCount();
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT 
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y, 
			sampleCount, sampleCount == 1 ? "sample" : "samples", nCores, 
			nCores == 1 ? "core" : "cores");

		Vector2i cropSize = film->getCropSize();
		Point2i cropOffset = film->getCropOffset();

		m_gatherPoints.clear();
		m_running = true;
		for (size_t i=0; i<m_blocks.size(); ++i)
			m_blocks[i]->decRef();
		m_blocks.clear();

		m_totalEmitted = 0;
		bool needsLensSample = camera->needsLensSample();
		bool needsTimeSample = camera->needsTimeSample();
		Log(EInfo, "Creating approximately %i gather points", cropSize.x*cropSize.y*sampleCount);
		Point2 lensSample, sample;
		RayDifferential eyeRay;
		Float timeSample = 0;
		m_filter = camera->getFilm()->getTabulatedFilter();
		Vector2 filterSize = m_filter->getFilterSize();
		int borderSize = (int) std::ceil(std::max(filterSize.x, filterSize.y));

		ref<Sampler> independentSampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));

		/* Create a sampler instance for every core */
		std::vector<SerializableObject *> samplers(sched->getCoreCount());
		for (size_t i=0; i<sched->getCoreCount(); ++i) {
			ref<Sampler> clonedSampler = independentSampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		int independentSamplerResID = sched->registerManifoldResource(samplers); 
		for (size_t i=0; i<sched->getCoreCount(); ++i)
			samplers[i]->decRef();

#ifdef MTS_DEBUG_FP
		enableFPExceptions();
#endif

		/* Create gather points in blocks so that gathering can be parallelized later on */
		for (int yofs=0; yofs<cropSize.y; yofs += m_blockSize) {
			for (int xofs=0; xofs<cropSize.x; xofs += m_blockSize) {
				ImageBlock *block = new ImageBlock(Vector2i(m_blockSize, m_blockSize), borderSize, 
					true, true, false, false);
				block->setSize(Vector2i(m_blockSize, m_blockSize));
				block->setOffset(Point2i(cropOffset.x + xofs, cropOffset.y + yofs));
				block->incRef();
				std::vector<GatherPoint> gatherPoints;
				gatherPoints.reserve(m_blockSize*m_blockSize*sampleCount);
				for (int yofsInt = 0; yofsInt < m_blockSize; ++yofsInt) {
					if (yofsInt + yofs >= cropSize.y)
						continue;
					for (int xofsInt = 0; xofsInt < m_blockSize; ++xofsInt) {
						if (xofsInt + xofs >= cropSize.x)
							continue;
						int y = cropOffset.y + yofs + yofsInt;
						int x = cropOffset.x + xofs + xofsInt;
						cameraSampler->generate();
						for (size_t j = 0; j<sampleCount; j++) {
							if (needsLensSample)
								lensSample = cameraSampler->next2D();
							if (needsTimeSample)
								timeSample = cameraSampler->next1D();
							sample = cameraSampler->next2D();
							sample.x += x; sample.y += y;
							camera->generateRayDifferential(sample, 
								lensSample, timeSample, eyeRay);
							size_t offset = gatherPoints.size();
							Float count = (Float) createGatherPoints(scene, eyeRay, sample, 
									cameraSampler, Spectrum(1.0f),
								gatherPoints, 1);
							if (count > 1) { // necessary because of filter weight computation
								for (int i = 0; i<count; ++i)
									gatherPoints[offset+i].weight *= count;
							}

							cameraSampler->advance();
						}
					}
				}
				m_blocks.push_back(block);
				m_gatherPoints.push_back(gatherPoints);
			}
		}

		int it=0;
		while (m_running) 
			photonMapPass(++it, queue, job, film, sceneResID, cameraResID, independentSamplerResID);

#ifdef MTS_DEBUG_FP
		disableFPExceptions();
#endif

		sched->unregisterResource(independentSamplerResID);
		return true;
	}

	int createGatherPoints(Scene *scene, const RayDifferential &ray, 
			const Point2 &sample, Sampler *sampler, const Spectrum &weight, 
			std::vector<GatherPoint> &gatherPoints, int depth) {
		int count = 0;
		if (depth >= m_maxDepth && m_maxDepth != -1)
			return 0;
		GatherPoint p;
		if (scene->rayIntersect(ray, p.its)) {
			const BSDF *bsdf = p.its.shape->getBSDF();
			if (!bsdf) {
				p.radius = 0;
				p.sample = sample;
				gatherPoints.push_back(p);
				++count;
			} else {
				if (bsdf->getType() & BSDF::ESmooth) {
					p.weight = weight;
					p.sample = sample;
					p.radius = m_initialRadius;
					p.depth = depth;
					if (p.its.isLuminaire())
						p.emission = p.its.Le(-ray.d);
					gatherPoints.push_back(p);
					++count;
				}

				if (bsdf->getType() & BSDF::EDelta) {
					int compCount = bsdf->getComponentCount();
					for (int i=0; i<compCount; i++) {
						if ((bsdf->getType(i) & BSDF::EDelta) == 0)
							continue;
						/* Sample the BSDF and recurse */
						BSDFQueryRecord bRec(p.its, sampler);
						bRec.component = i;
						Spectrum bsdfVal = bsdf->sample(bRec, Point2(0.0f));
						if (bsdfVal.isZero())
							continue;
						bsdfVal = bsdf->eval(bRec, EDiscrete);

						const Float rrProb = depth < 4 ? 1 : 0.8f;
						if (sampler->independent1D() < rrProb) {
							RayDifferential recursiveRay(p.its.p, p.its.toWorld(bRec.wo), ray.time);
							count += createGatherPoints(scene, recursiveRay, sample, sampler, 
								weight * bsdfVal / rrProb, gatherPoints, depth+1);
						}
					}
				}
			}
		} else if (depth == 1) {
			/* Generate an invalid sample */
			p.emission = scene->LeBackground(ray);
			p.radius = 0;
			p.sample = sample;
			gatherPoints.push_back(p);
			++count;
		}
		return count;
	}

	void photonMapPass(int it, RenderQueue *queue, const RenderJob *job,  
			Film *film, int sceneResID, int cameraResID, int samplerResID) {
		Log(EInfo, "Performing a photon mapping pass %i", it);
		ref<Scheduler> sched = Scheduler::getInstance();

		/* Generate the global photon map */
		ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
			GatherPhotonProcess::EAllSurfacePhotons, m_photonCount,
			m_granularity, m_maxDepth == -1 ? -1 : (m_maxDepth-1), m_rrDepth, true,
			m_autoCancelGathering, job);

		proc->bindResource("scene", sceneResID);
		proc->bindResource("camera", cameraResID);
		proc->bindResource("sampler", samplerResID);

		sched->schedule(proc);
		sched->wait(proc);

		ref<PhotonMap> photonMap = proc->getPhotonMap();
		photonMap->build();
		Log(EDebug, "Photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: " 
			SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

		Log(EInfo, "Gathering ..");
		m_totalEmitted += proc->getShotParticles();
		film->clear();
		#pragma omp parallel for schedule(dynamic)
		for (int blockIdx = 0; blockIdx<(int) m_blocks.size(); ++blockIdx) {
			ImageBlock *block = m_blocks[blockIdx];
			block->clear();
			std::vector<GatherPoint> &gatherPoints = m_gatherPoints[blockIdx];
			Spectrum flux, contrib;

			for (size_t i=0; i<gatherPoints.size(); ++i) {
				GatherPoint &g = gatherPoints[i];

				if (g.radius == 0) {
					/* Generate a black sample -- necessary for proper sample weight 
					   computation at edges */
					block->putSample(g.sample, g.emission, 1, m_filter);
					continue;
				}

				Float M = (Float) photonMap->estimateRadianceRaw(
					g.its, g.radius, flux, m_maxDepth == -1 ? INT_MAX : (m_maxDepth-g.depth));
				Float N = g.N;

				if (N+M == 0) {
					g.flux = contrib = Spectrum(0.0f);
				} else {
					Float ratio = (N + m_alpha * M) / (N + M);
					g.flux = (g.flux + flux) * ratio;
					g.radius = g.radius * std::sqrt(ratio);
					g.N = N + m_alpha * M;
				}
				contrib = g.flux / ((Float) m_totalEmitted * g.radius*g.radius * M_PI)
					+ g.emission;
				block->putSample(g.sample, contrib * g.weight, 1, m_filter, false);
			}
			m_mutex->lock();
			film->putImageBlock(block);
			m_mutex->unlock();
		}
		queue->signalRefresh(job, NULL);
	}


	std::string toString() const {
		return "ProgressivePhotonMapIntegrator[]";
	}

	MTS_DECLARE_CLASS()
private:
	const TabulatedFilter *m_filter;
	std::vector<ImageBlock *> m_blocks;
	std::vector<std::vector<GatherPoint> > m_gatherPoints;
	ref<Mutex> m_mutex;
	Float m_initialRadius, m_alpha;
	int m_photonCount, m_granularity;
	int m_maxDepth, m_rrDepth;
	size_t m_totalEmitted;
	int m_blockSize;
	bool m_running;
	bool m_autoCancelGathering;
};

MTS_IMPLEMENT_CLASS_S(ProgressivePhotonMapIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(ProgressivePhotonMapIntegrator, "Progressive photon mapper");
MTS_NAMESPACE_END
