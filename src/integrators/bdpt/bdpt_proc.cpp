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
#include <mitsuba/core/sfcurve.h>
#include <mitsuba/bidir/util.h>
#include "bdpt_proc.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

class BDPTRenderer : public WorkProcessor {
public:
	BDPTRenderer(const BDPTConfiguration &config) : m_config(config) { }

	BDPTRenderer(Stream *stream, InstanceManager *manager)
		: WorkProcessor(stream, manager), m_config(stream) { }

	virtual ~BDPTRenderer() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		m_config.serialize(stream);
	}

	ref<WorkUnit> createWorkUnit() const {
		return new RectangularWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new BDPTWorkResult(m_config, m_rfilter.get(),
			Vector2i(m_config.blockSize));
	}

	void prepare() {
		Scene *scene = static_cast<Scene *>(getResource("scene"));
		m_scene = new Scene(scene);
		m_sampler = static_cast<Sampler *>(getResource("sampler"));
		m_sensor = static_cast<Sensor *>(getResource("sensor"));
		m_rfilter = m_sensor->getFilm()->getReconstructionFilter();
		m_scene->removeSensor(scene->getSensor());
		m_scene->addSensor(m_sensor);
		m_scene->setSensor(m_sensor);
		m_scene->setSampler(m_sampler);
		m_scene->wakeup(NULL, m_resources);
		m_scene->initializeBidirectional();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
		const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
		BDPTWorkResult *result = static_cast<BDPTWorkResult *>(workResult);
		bool needsTimeSample = m_sensor->needsTimeSample();
		Float time = m_sensor->getShutterOpen();

		result->setOffset(rect->getOffset());
		result->setSize(rect->getSize());
		result->clear();
		m_hilbertCurve.initialize(TVector2<uint8_t>(rect->getSize()));

		#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
		#endif

		Path emitterSubpath;
		Path sensorSubpath;

		/* Determine the necessary random walk depths based on properties of
		   the endpoints */
		int emitterDepth = m_config.maxDepth,
		    sensorDepth = m_config.maxDepth;

		/* Go one extra step if the sensor can be intersected */
		if (!m_scene->hasDegenerateSensor() && emitterDepth != -1)
			++emitterDepth;

		/* Go one extra step if there are emitters that can be intersected */
		if (!m_scene->hasDegenerateEmitters() && sensorDepth != -1)
			++sensorDepth;

		for (size_t i=0; i<m_hilbertCurve.getPointCount(); ++i) {
			Point2i offset = Point2i(m_hilbertCurve[i]) + Vector2i(rect->getOffset());
			m_sampler->generate(offset);

			for (size_t j = 0; j<m_sampler->getSampleCount(); j++) {
				if (stop)
					break;

				if (needsTimeSample)
					time = m_sensor->sampleTime(m_sampler->next1D());

				/* Start new emitter and sensor subpaths */
				emitterSubpath.initialize(m_scene, time, EImportance, m_pool);
				sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

				/* Perform a random walk using alternating steps on each path */
				Path::alternatingRandomWalkFromPixel(m_scene, m_sampler,
					emitterSubpath, emitterDepth, sensorSubpath,
					sensorDepth, offset, m_config.rrDepth, m_pool);

				evaluate(result, emitterSubpath, sensorSubpath);

				emitterSubpath.release(m_pool);
				sensorSubpath.release(m_pool);

				m_sampler->advance();
			}
		}

		#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
		#endif

		/* Make sure that there were no memory leaks */
		Assert(m_pool.unused());
	}

	/// Evaluate the contributions of the given eye and light paths
	void evaluate(BDPTWorkResult *wr,
			Path &emitterSubpath, Path &sensorSubpath) {
		Point2 initialSamplePos = sensorSubpath.vertex(1)->getSamplePosition();
		const Scene *scene = m_scene;
		PathVertex tempEndpoint, tempSample;
		PathEdge tempEdge, connectionEdge;

		/* Compute the combined weights along the two subpaths */
		Spectrum *importanceWeights = (Spectrum *) alloca(emitterSubpath.vertexCount() * sizeof(Spectrum)),
				 *radianceWeights  = (Spectrum *) alloca(sensorSubpath.vertexCount()  * sizeof(Spectrum));

		importanceWeights[0] = radianceWeights[0] = Spectrum(1.0f);
		for (size_t i=1; i<emitterSubpath.vertexCount(); ++i)
			importanceWeights[i] = importanceWeights[i-1] *
				emitterSubpath.vertex(i-1)->weight[EImportance] *
				emitterSubpath.vertex(i-1)->rrWeight *
				emitterSubpath.edge(i-1)->weight[EImportance];

		for (size_t i=1; i<sensorSubpath.vertexCount(); ++i)
			radianceWeights[i] = radianceWeights[i-1] *
				sensorSubpath.vertex(i-1)->weight[ERadiance] *
				sensorSubpath.vertex(i-1)->rrWeight *
				sensorSubpath.edge(i-1)->weight[ERadiance];

		Spectrum sampleValue(0.0f);
		for (int s = (int) emitterSubpath.vertexCount()-1; s >= 0; --s) {
			/* Determine the range of sensor vertices to be traversed,
			   while respecting the specified maximum path length */
			int minT = std::max(2-s, m_config.lightImage ? 0 : 2),
			    maxT = (int) sensorSubpath.vertexCount() - 1;
			if (m_config.maxDepth != -1)
				maxT = std::min(maxT, m_config.maxDepth + 1 - s);

			for (int t = maxT; t >= minT; --t) {
				PathVertex
					*vsPred = emitterSubpath.vertexOrNull(s-1),
					*vtPred = sensorSubpath.vertexOrNull(t-1),
					*vs = emitterSubpath.vertex(s),
					*vt = sensorSubpath.vertex(t);
				PathEdge
					*vsEdge = emitterSubpath.edgeOrNull(s-1),
					*vtEdge = sensorSubpath.edgeOrNull(t-1);

				RestoreMeasureHelper rmh0(vs), rmh1(vt);

				/* Will be set to true if direct sampling was used */
				bool sampleDirect = false;

				/* Stores the pixel position associated with this sample */
				Point2 samplePos = initialSamplePos;

				/* Allowed remaining number of ENull vertices that can
				   be bridged via pathConnect (negative=arbitrarily many) */
				int remaining = m_config.maxDepth - s - t + 1;

				/* Will receive the path weight of the (s, t)-connection */
				Spectrum value;

				/* Account for the terms of the measurement contribution
				   function that are coupled to the connection endpoints */
				if (vs->isEmitterSupernode()) {
					/* If possible, convert 'vt' into an emitter sample */
					if (!vt->cast(scene, PathVertex::EEmitterSample) || vt->isDegenerate())
						continue;

					value = radianceWeights[t] *
						vs->eval(scene, vsPred, vt, EImportance) *
						vt->eval(scene, vtPred, vs, ERadiance);
				} else if (vt->isSensorSupernode()) {
					/* If possible, convert 'vs' into an sensor sample */
					if (!vs->cast(scene, PathVertex::ESensorSample) || vs->isDegenerate())
						continue;

					/* Make note of the changed pixel sample position */
					if (!vs->getSamplePosition(vsPred, samplePos))
						continue;

					value = importanceWeights[s] *
						vs->eval(scene, vsPred, vt, EImportance) *
						vt->eval(scene, vtPred, vs, ERadiance);
				} else if (m_config.sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
					/* s==1/t==1 path: use a direct sampling strategy if requested */
					if (s == 1) {
						if (vt->isDegenerate())
							continue;
						/* Generate a position on an emitter using direct sampling */
						value = radianceWeights[t] * vt->sampleDirect(scene, m_sampler,
							&tempEndpoint, &tempEdge, &tempSample, EImportance);
						if (value.isZero())
							continue;
						vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
						value *= vt->eval(scene, vtPred, vs, ERadiance);
						vt->measure = EArea;
					} else {
						if (vs->isDegenerate())
							continue;
						/* Generate a position on the sensor using direct sampling */
						value = importanceWeights[s] * vs->sampleDirect(scene, m_sampler,
							&tempEndpoint, &tempEdge, &tempSample, ERadiance);
						if (value.isZero())
							continue;
						vt = &tempSample; vtPred = &tempEndpoint; vtEdge = &tempEdge;
						value *= vs->eval(scene, vsPred, vt, EImportance);
						vs->measure = EArea;
					}

					sampleDirect = true;
				} else {
					/* Can't connect degenerate endpoints */
					if (vs->isDegenerate() || vt->isDegenerate())
						continue;

					value = importanceWeights[s] * radianceWeights[t] *
						vs->eval(scene, vsPred, vt, EImportance) *
						vt->eval(scene, vtPred, vs, ERadiance);

					/* Temporarily force vertex measure to EArea. Needed to
					   handle BSDFs with diffuse + specular components */
					vs->measure = vt->measure = EArea;
				}

				/* Attempt to connect the two endpoints, which could result in
				   the creation of additional vertices (index-matched boundaries etc.) */
				int interactions = remaining; // backup
				if (value.isZero() || !connectionEdge.pathConnectAndCollapse(
						scene, vsEdge, vs, vt, vtEdge, interactions))
					continue;

				/* Account for the terms of the measurement contribution
				   function that are coupled to the connection edge */
				if (!sampleDirect)
					value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
				else
					value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
							(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

				if (sampleDirect) {
					/* A direct sampling strategy was used, which generated
					   two new vertices at one of the path ends. Temporarily
					   modify the path to reflect this change */
					if (t == 1)
						sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
					else
						emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
				}

				/* Compute the multiple importance sampling weight */
				Float miWeight = Path::miWeight(scene, emitterSubpath, &connectionEdge,
					sensorSubpath, s, t, m_config.sampleDirect, m_config.lightImage);

				if (sampleDirect) {
					/* Now undo the previous change */
					if (t == 1)
						sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
					else
						emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
				}

				/* Determine the pixel sample position when necessary */
				if (vt->isSensorSample() && !vt->getSamplePosition(vs, samplePos))
					continue;

				#if BDPT_DEBUG == 1
					/* When the debug mode is on, collect samples
					   separately for each sampling strategy. Note: the
					   following piece of code artificially increases the
					   exposure of longer paths */
					Spectrum splatValue = value * (m_config.showWeighted
						? miWeight : 1.0f);// * std::pow(2.0f, s+t-3.0f));
					wr->putDebugSample(s, t, samplePos, splatValue);
				#endif

				if (t >= 2)
					sampleValue += value * miWeight;
				else
					wr->putLightSample(samplePos, value * miWeight);
			}
		}
		wr->putSample(initialSamplePos, sampleValue);
	}

	ref<WorkProcessor> clone() const {
		return new BDPTRenderer(m_config);
	}

	MTS_DECLARE_CLASS()
private:
	ref<Scene> m_scene;
	ref<Sensor> m_sensor;
	ref<Sampler> m_sampler;
	ref<ReconstructionFilter> m_rfilter;
	MemoryPool m_pool;
	BDPTConfiguration m_config;
	HilbertCurve2D<uint8_t> m_hilbertCurve;
};


/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

BDPTProcess::BDPTProcess(const RenderJob *parent, RenderQueue *queue,
		const BDPTConfiguration &config) :
	BlockedRenderProcess(parent, queue, config.blockSize), m_config(config) {
	m_refreshTimer = new Timer();
}

ref<WorkProcessor> BDPTProcess::createWorkProcessor() const {
	return new BDPTRenderer(m_config);
}

void BDPTProcess::develop() {
	if (!m_config.lightImage)
		return;
	LockGuard lock(m_resultMutex);
	const ImageBlock *lightImage = m_result->getLightImage();
	m_film->setBitmap(m_result->getImageBlock()->getBitmap());
	m_film->addBitmap(lightImage->getBitmap(), 1.0f / m_config.sampleCount);
	m_refreshTimer->reset();
	m_queue->signalRefresh(m_parent);
}

void BDPTProcess::processResult(const WorkResult *wr, bool cancelled) {
	if (cancelled)
		return;
	const BDPTWorkResult *result = static_cast<const BDPTWorkResult *>(wr);
	ImageBlock *block = const_cast<ImageBlock *>(result->getImageBlock());
	LockGuard lock(m_resultMutex);
	m_progress->update(++m_resultCount);
	if (m_config.lightImage) {
		const ImageBlock *lightImage = m_result->getLightImage();
		m_result->put(result);
		if (m_parent->isInteractive()) {
			/* Modify the finished image block so that it includes the light image contributions,
			   which creates a more intuitive preview of the rendering process. This is
			   not 100% correct but doesn't matter, as the shown image will be properly re-developed
			   every 2 seconds and once more when the rendering process finishes */

			Float invSampleCount = 1.0f / m_config.sampleCount;
			const Bitmap *sourceBitmap = lightImage->getBitmap();
			Bitmap *destBitmap = block->getBitmap();
			int borderSize = block->getBorderSize();
			Point2i offset = block->getOffset();
			Vector2i size = block->getSize();

			for (int y=0; y<size.y; ++y) {
				const Float *source = sourceBitmap->getFloatData()
					+ (offset.x + (y+offset.y) * sourceBitmap->getWidth()) * SPECTRUM_SAMPLES;
				Float *dest = destBitmap->getFloatData()
					+ (borderSize + (y + borderSize) * destBitmap->getWidth()) * (SPECTRUM_SAMPLES + 2);

				for (int x=0; x<size.x; ++x) {
					Float weight = dest[SPECTRUM_SAMPLES + 1] * invSampleCount;
					for (int k=0; k<SPECTRUM_SAMPLES; ++k)
						*dest++ += *source++ * weight;
					dest += 2;
				}
			}
		}
	}

	m_film->put(block);

	/* Re-develop the entire image every two seconds if partial results are
	   visible (e.g. in a graphical user interface). This only applies when
	   there is a light image. */
	bool developFilm = m_config.lightImage &&
		(m_parent->isInteractive() && m_refreshTimer->getMilliseconds() > 2000);

	m_queue->signalWorkEnd(m_parent, result->getImageBlock(), false);

	if (developFilm)
		develop();
}

void BDPTProcess::bindResource(const std::string &name, int id) {
	BlockedRenderProcess::bindResource(name, id);
	if (name == "sensor" && m_config.lightImage) {
		/* If needed, allocate memory for the light image */
		m_result = new BDPTWorkResult(m_config, NULL, m_film->getCropSize());
		m_result->clear();
	}
}

MTS_IMPLEMENT_CLASS_S(BDPTRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(BDPTProcess, false, BlockedRenderProcess)
MTS_NAMESPACE_END
