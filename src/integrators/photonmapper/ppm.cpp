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
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/renderqueue.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{ppm}{Progressive photon mapping integrator}
 * \order{7}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination,
 *	       and so on. \default{\code{-1}}
 *	   }
 *     \parameter{photonCount}{\Integer}{Number of photons to be shot per iteration\default{250000}}
 *     \parameter{initialRadius}{\Float}{Initial radius of gather points in world space units.
 *         \default{0, i.e. decide automatically}}
 *     \parameter{alpha}{\Float}{Radius reduction parameter \code{alpha} from the paper\default{0.7}}
 *     \parameter{granularity}{\Integer}{
		Granularity of photon tracing work units for the purpose
		of parallelization (in \# of shot particles) \default{0, i.e. decide automatically}
 *     }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	      which the implementation will start to use the ``russian roulette''
 *	      path termination criterion. \default{\code{5}}
 *	   }
 *     \parameter{maxPasses}{\Integer}{Maximum number of passes to render (where \code{-1}
 *        corresponds to rendering until stopped manually). \default{\code{-1}}}
 * }
 * This plugin implements the progressive photon mapping algorithm by Hachisuka et al.
 * \cite{Hachisuka2008Progressive}. Progressive photon mapping is a variant of photon
 * mapping that alternates between photon shooting and gathering passes that involve
 * a relatively small (e.g. 250K) numbers of photons that are subsequently discarded.
 *
 * This is done in a way such that the variance and bias of the resulting output
 * vanish as the number of passes tends to infinity. The progressive nature of this
 * method enables renderings with an effectively arbitrary number of photons
 * without exhausting the available system memory.
 *
 * The desired sample count specified in the sample generator configuration
 * determines how many photon query points are created per pixel. It should not be
 * set too high, since the rendering time is approximately proportional to
 * this number. For good results, use between 2-4 samples along with the
 * \code{ldsampler}. Once started, the rendering process continues indefinitely
 * until it is manually stopped.
 *
 * \remarks{
 *    \item Due to the data dependencies of this algorithm, the parallelization is
 *    limited to the local machine (i.e. cluster-wide renderings are not implemented)
 *    \item This integrator does not handle participating media
 *    \item This integrator does not currently work with subsurface scattering
 *    models.
 * }
 */

class PPMIntegrator : public Integrator {
public:
	/// Represents one individual PPM gather point including relevant statistics
	struct GatherPoint {
		Intersection its;
		Float radius;
		Spectrum weight, flux, emission;
		Point2 sample;
		Float N;
		int depth;

		inline GatherPoint() : weight(0.0f), flux(0.0f), emission(0.0f), N(0.0f) {
		}
	};

	/// Work unit for parallelizaition
	struct PPMWorkUnit {
		ref<ImageBlock> block;
		std::vector<GatherPoint> gatherPoints;
	};

	PPMIntegrator(const Properties &props) : Integrator(props) {
		/* Initial photon query radius (0 = infer based on scene size and sensor resolution) */
		m_initialRadius = props.getFloat("initialRadius", 0);
		/* Alpha parameter from the paper (influences the speed, at which the photon radius is reduced) */
		m_alpha = props.getFloat("alpha", .7);
		/* Number of photons to shoot in each iteration */
		m_photonCount = props.getInteger("photonCount", 250000);
		/* Granularity of the work units used in parallelizing the
		   particle tracing task (default: choose automatically). */
		m_granularity = props.getInteger("granularity", 0);
		/* Longest visualized path length (<tt>-1</tt>=infinite). When a positive value is
		   specified, it must be greater or equal to <tt>2</tt>, which corresponds to single-bounce
		   (direct-only) illumination */
		m_maxDepth = props.getInteger("maxDepth", -1);
		/* Depth to start using russian roulette */
		m_rrDepth = props.getInteger("rrDepth", 3);
		/* Indicates if the gathering steps should be canceled if not enough photons are generated. */
		m_autoCancelGathering = props.getBoolean("autoCancelGathering", true);
        /* Maximum number of passes to render. -1 renders until the process is stopped. */
		m_maxPasses = props.getInteger("maxPasses", -1);

		m_mutex = new Mutex();
		if (m_maxDepth <= 1 && m_maxDepth != -1)
			Log(EError, "Maximum depth must either be set to \"-1\" or \"2\" or higher!");
		if (m_maxPasses <= 0 && m_maxPasses != -1)
			Log(EError, "Maximum number of Passes must either be set to \"-1\" or \"1\" or higher!");
	}

	virtual ~PPMIntegrator() {
		for (size_t i=0; i<m_workUnits.size(); ++i)
			delete m_workUnits[i];
		m_workUnits.clear();
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		/* Prepare the sampler for bucket-based rendering */
		sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
	}

	void cancel() {
		m_running = false;
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

		if (m_initialRadius == 0) {
			/* Guess an initial radius if not provided
			  (scene width / horizontal or vertical pixel count) * 5 */
			Float rad = scene->getBSphere().radius;
			Vector2i filmSize = scene->getSensor()->getFilm()->getSize();

			m_initialRadius = std::min(rad / filmSize.x, rad / filmSize.y) * 5;
		}
		return true;
	}

	bool render(Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sensor> sensor = scene->getSensor();
		ref<Film> film = sensor->getFilm();
		size_t nCores = sched->getCoreCount();
		Sampler *sensorSampler = (Sampler *) sched->getResource(samplerResID, 0);

		size_t sampleCount = sensorSampler->getSampleCount();
		Vector2i cropSize = film->getCropSize();
		Point2i cropOffset = film->getCropOffset();
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
			" %s, " SSE_STR ") ..", cropSize.x, cropSize.y, sampleCount,
			sampleCount == 1 ? "sample" : "samples", nCores, nCores == 1 ? "core" : "cores");

		m_running = true;
		m_totalEmissions = 0;
		m_totalPhotons = 0;
		for (size_t i=0; i<m_workUnits.size(); ++i)
			delete m_workUnits[i];
		m_workUnits.clear();

		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();
		Log(EInfo, "Creating approximately %i gather points ..", cropSize.x*cropSize.y*sampleCount);
		Point2 apertureSample, sample;
		RayDifferential sensorRay;
		Float timeSample = 0;

		ref<Sampler> indepSampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));

		/* Create a sampler instance for every core */
		std::vector<SerializableObject *> samplers(sched->getCoreCount());
		for (size_t i=0; i<sched->getCoreCount(); ++i) {
			ref<Sampler> clonedSampler = indepSampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		int indepSamplerResID = sched->registerMultiResource(samplers);
		for (size_t i=0; i<sched->getCoreCount(); ++i)
			samplers[i]->decRef();

#ifdef MTS_DEBUG_FP
		enableFPExceptions();
#endif

		int blockSize = scene->getBlockSize();
		int blocksW = (cropSize.x + blockSize - 1) / blockSize;
		int blocksH = (cropSize.y + blockSize - 1) / blockSize;
		m_workUnits.resize(blocksW*blocksH);

		/* Create gather points in blocks so that gathering can be parallelized later on */
		for (int yofs=0; yofs<cropSize.y; yofs += blockSize) {
			for (int xofs=0; xofs<cropSize.x; xofs += blockSize) {
				PPMWorkUnit *wu = new PPMWorkUnit();
				m_workUnits[xofs/blockSize+yofs/blockSize*blocksW] = wu;

				wu->block = new ImageBlock(Bitmap::ESpectrumAlphaWeight,
					Vector2i(blockSize), film->getReconstructionFilter());
				wu->block->setOffset(Point2i(cropOffset.x + xofs, cropOffset.y + yofs));
				wu->gatherPoints.clear();
				wu->gatherPoints.reserve(blockSize*blockSize*sampleCount);

				for (int yofsInt = 0; yofsInt < blockSize; ++yofsInt) {
					if (yofsInt + yofs >= cropSize.y)
						continue;
					for (int xofsInt = 0; xofsInt < blockSize; ++xofsInt) {
						if (xofsInt + xofs >= cropSize.x)
							continue;
						int y = cropOffset.y + yofs + yofsInt;
						int x = cropOffset.x + xofs + xofsInt;
						sensorSampler->generate(Point2i(xofs, yofs));

						for (size_t j = 0; j<sampleCount; j++) {
							sample = sensorSampler->next2D();
							if (needsApertureSample)
								apertureSample = sensorSampler->next2D();
							if (needsTimeSample)
								timeSample = sensorSampler->next1D();
							sample.x += x; sample.y += y;
							sensor->sampleRayDifferential(sensorRay, sample,
								apertureSample, timeSample);
							size_t offset = wu->gatherPoints.size();
							int count = createGatherPoints(scene, sensorRay, sample,
								sensorSampler, Spectrum(1.0f), wu->gatherPoints, 1);
							const Float fcount = static_cast<Float>(count);
							for (int i = 0; i<count; ++i)
								wu->gatherPoints[offset+i].weight *= fcount;

							sensorSampler->advance();
						}
					}
				}
			}
		}

		#if defined(MTS_OPENMP)
			Thread::initializeOpenMP(nCores);
		#endif

		int it = 0;
		while (m_running && (m_maxPasses == -1 || it < m_maxPasses)) {
			photonMapPass(++it, queue, job, film, sceneResID, sensorResID, indepSamplerResID);
        }

#ifdef MTS_DEBUG_FP
		disableFPExceptions();
#endif

		sched->unregisterResource(indepSamplerResID);
		return true;
	}

	int createGatherPoints(Scene *scene, const RayDifferential &ray,
			const Point2 &sample, Sampler *sampler, const Spectrum &weight,
			std::vector<GatherPoint> &gatherPoints, int depth) {
		int count = 0;
		if (depth >= m_maxDepth && m_maxDepth != -1)
			return 0;
		GatherPoint p;
		p.emission = Spectrum(0.0f);
		if (scene->rayIntersect(ray, p.its)) {
			const BSDF *bsdf = p.its.getBSDF();
			if (bsdf->getType() & BSDF::ESmooth) {
				p.weight = weight;
				p.sample = sample;
				p.radius = m_initialRadius;
				p.depth = depth;
				if (p.its.isEmitter())
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
					BSDFSamplingRecord bRec(p.its, sampler);
					bRec.component = i;
					Spectrum bsdfVal = bsdf->sample(bRec, Point2(0.0f));
					if (bsdfVal.isZero())
						continue;
					bsdfVal = bsdf->eval(bRec, EDiscrete);

					const Float rrProb = depth < 4 ? 1 : 0.8f;
					if (sampler->next1D() < rrProb) {
						RayDifferential recursiveRay(p.its.p, p.its.toWorld(bRec.wo), ray.time);
						count += createGatherPoints(scene, recursiveRay, sample, sampler,
							weight * bsdfVal / rrProb, gatherPoints, depth+1);
					}
				}
			}
		} else if (depth == 1) {
			/* Generate an invalid sample */
			p.emission = scene->evalEnvironment(ray);
			p.radius = 0;
			p.sample = sample;
			gatherPoints.push_back(p);
			++count;
		}
		return count;
	}

	void photonMapPass(int it, RenderQueue *queue, const RenderJob *job,
			Film *film, int sceneResID, int sensorResID, int samplerResID) {
		Log(EInfo, "Performing a photon mapping pass %i (" SIZE_T_FMT " photons so far)",
				it, m_totalPhotons);
		ref<Scheduler> sched = Scheduler::getInstance();

		/* Generate the global photon map */
		ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
			GatherPhotonProcess::EAllSurfacePhotons, m_photonCount,
			m_granularity, m_maxDepth == -1 ? -1 : (m_maxDepth-1), m_rrDepth, true,
			m_autoCancelGathering, job);

		proc->bindResource("scene", sceneResID);
		proc->bindResource("sensor", sensorResID);
		proc->bindResource("sampler", samplerResID);

		sched->schedule(proc);
		sched->wait(proc);

		ref<PhotonMap> photonMap = proc->getPhotonMap();
		photonMap->build();
		Log(EDebug, "Photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
			SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

		Log(EInfo, "Gathering ..");
		m_totalEmissions += proc->getShotParticles();
		m_totalPhotons += photonMap->size();

		film->clear();
		#if defined(MTS_OPENMP)
			#pragma omp parallel for schedule(dynamic)
		#endif
		for (int wuIdx = 0; wuIdx < (int) m_workUnits.size(); ++wuIdx) {
			PPMWorkUnit *wu = m_workUnits[wuIdx];
			Spectrum flux, contrib;

			wu->block->clear();
			for (size_t i=0; i<wu->gatherPoints.size(); ++i) {
				GatherPoint &g = wu->gatherPoints[i];

				if (g.radius == 0) {
					/* Generate a black sample -- necessary for proper
					   sample weight computation at surface boundaries */
					wu->block->put(g.sample, g.emission, 1);
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
				contrib = g.flux / ((Float) m_totalEmissions * g.radius*g.radius * M_PI)
					+ g.emission;
				wu->block->put(g.sample, contrib * g.weight, 1);
			}
			LockGuard guard(m_mutex);
			film->put(wu->block);
		}
		queue->signalRefresh(job);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SPPMIntegrator[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  initialRadius = " << m_initialRadius << "," << endl
			<< "  alpha = " << m_alpha << "," << endl
			<< "  photonCount = " << m_photonCount << "," << endl
			<< "  granularity = " << m_granularity << "," << endl
			<< "  maxPasses = " << m_maxPasses << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<PPMWorkUnit *> m_workUnits;
	Float m_initialRadius, m_alpha;
	int m_photonCount, m_granularity;
	int m_maxDepth, m_rrDepth;
	size_t m_totalEmissions, m_totalPhotons;
	int m_blockSize;
	bool m_running;
	bool m_autoCancelGathering;
	ref<Mutex> m_mutex;
	int m_maxPasses;
};

MTS_IMPLEMENT_CLASS(PPMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(PPMIntegrator, "Progressive photon mapper");
MTS_NAMESPACE_END
