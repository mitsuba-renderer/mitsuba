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

#include <mitsuba/bidir/util.h>
#include <mitsuba/bidir/path.h>
#include "pssmlt_proc.h"
#include "pssmlt_sampler.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

StatsCounter largeStepRatio("Primary sample space MLT",
	"Accepted large steps", EPercentage);
StatsCounter smallStepRatio("Primary sample space MLT",
	"Accepted small steps", EPercentage);
StatsCounter acceptanceRate("Primary sample space MLT",
	"Overall acceptance rate", EPercentage);
StatsCounter forcedAcceptance("Primary sample space MLT",
	"Number of forced acceptances");

class PSSMLTRenderer : public WorkProcessor {
public:
	PSSMLTRenderer(const PSSMLTConfiguration &conf)
		: m_config(conf) {
	}

	PSSMLTRenderer(Stream *stream, InstanceManager *manager)
		: WorkProcessor(stream, manager) {
		m_config = PSSMLTConfiguration(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		m_config.serialize(stream);
	}

	ref<WorkUnit> createWorkUnit() const {
		return new SeedWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new ImageBlock(Bitmap::ESpectrum,
			m_film->getCropSize(), m_film->getReconstructionFilter());
	}

	void prepare() {
		Scene *scene = static_cast<Scene *>(getResource("scene"));
		m_origSampler = static_cast<PSSMLTSampler *>(getResource("sampler"));
		m_sensor = static_cast<Sensor *>(getResource("sensor"));
		m_scene = new Scene(scene);
		m_film = m_sensor->getFilm();
		m_scene->setSensor(m_sensor);
		m_scene->setSampler(m_origSampler);
		m_scene->removeSensor(scene->getSensor());
		m_scene->addSensor(m_sensor);
		m_scene->setSensor(m_sensor);
		m_scene->wakeup(NULL, m_resources);
		m_scene->initializeBidirectional();

		m_rplSampler = static_cast<ReplayableSampler*>(
			static_cast<Sampler *>(getResource("rplSampler"))->clone().get());
		m_sensorSampler = new PSSMLTSampler(m_origSampler);
		m_emitterSampler = new PSSMLTSampler(m_origSampler);
		m_directSampler = new PSSMLTSampler(m_origSampler);

		m_pathSampler = new PathSampler(m_config.technique, m_scene,
			m_emitterSampler, m_sensorSampler, m_directSampler, m_config.maxDepth,
			m_config.rrDepth, m_config.separateDirect, m_config.directSampling);
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
		ImageBlock *result = static_cast<ImageBlock *>(workResult);
		const SeedWorkUnit *wu = static_cast<const SeedWorkUnit *>(workUnit);
		const PathSeed &seed = wu->getSeed();
		SplatList *current = new SplatList(), *proposed = new SplatList();

		m_emitterSampler->reset();
		m_sensorSampler->reset();
		m_directSampler->reset();
		m_sensorSampler->setRandom(m_rplSampler->getRandom());
		m_emitterSampler->setRandom(m_rplSampler->getRandom());
		m_directSampler->setRandom(m_rplSampler->getRandom());

		/* Generate the initial sample by replaying the seeding random
		   number stream at the appropriate position. Afterwards, revert
		   back to this worker's own source of random numbers */
		m_rplSampler->setSampleIndex(seed.sampleIndex);

		m_pathSampler->sampleSplats(Point2i(-1), *current);
		result->clear();

		ref<Random> random = m_origSampler->getRandom();
		m_sensorSampler->setRandom(random);
		m_emitterSampler->setRandom(random);
		m_directSampler->setRandom(random);
		m_rplSampler->updateSampleIndex(m_rplSampler->getSampleIndex()
			+ m_sensorSampler->getSampleIndex()
			+ m_emitterSampler->getSampleIndex()
			+ m_directSampler->getSampleIndex());

		m_sensorSampler->accept();
		m_emitterSampler->accept();
		m_directSampler->accept();

		/* Sanity check -- the luminance should match the one from
		   the warmup phase - an error here would indicate inconsistencies
		   regarding the use of random numbers during sample generation */
		if (std::abs((current->luminance - seed.luminance)
				/ seed.luminance) > Epsilon)
			Log(EError, "Error when reconstructing a seed path: luminance "
				"= %f, but expected luminance = %f", current->luminance, seed.luminance);

		ref<Timer> timer = new Timer();

		/* MLT main loop */
		Float cumulativeWeight = 0;
		current->normalize(m_config.importanceMap);
		for (uint64_t mutationCtr=0; mutationCtr<m_config.nMutations && !stop; ++mutationCtr) {
			if (wu->getTimeout() > 0 && (mutationCtr % 8192) == 0
					&& (int) timer->getMilliseconds() > wu->getTimeout())
				break;

			bool largeStep = random->nextFloat() < m_config.pLarge;
			m_sensorSampler->setLargeStep(largeStep);
			m_emitterSampler->setLargeStep(largeStep);
			m_directSampler->setLargeStep(largeStep);

			m_pathSampler->sampleSplats(Point2i(-1), *proposed);
			proposed->normalize(m_config.importanceMap);

			Float a = std::min((Float) 1.0f, proposed->luminance / current->luminance);

			if (std::isnan(proposed->luminance) || proposed->luminance < 0) {
				Log(EWarn, "Encountered a sample with luminance = %f, ignoring!",
						proposed->luminance);
				a = 0;
			}

			bool accept;
			Float currentWeight, proposedWeight;

			if (a > 0) {
				if (m_config.kelemenStyleWeights && !m_config.importanceMap) {
					/* Kelemen-style MLT weights (these don't work for 2-stage MLT) */
					currentWeight = (1 - a) * current->luminance
						/ (current->luminance/m_config.luminance + m_config.pLarge);
					proposedWeight = (a + (largeStep ? 1 : 0)) * proposed->luminance
						/ (proposed->luminance/m_config.luminance + m_config.pLarge);
				} else {
					/* Veach-style use of expectations */
					currentWeight = 1-a;
					proposedWeight = a;
				}
				accept = (a == 1) || (random->nextFloat() < a);
			} else {
				if (m_config.kelemenStyleWeights)
					currentWeight = current->luminance
						/ (current->luminance/m_config.luminance + m_config.pLarge);
				else
					currentWeight = 1;
				proposedWeight = 0;
				accept = false;
			}

			cumulativeWeight += currentWeight;
			if (accept) {
				for (size_t k=0; k<current->size(); ++k) {
					Spectrum value = current->getValue(k) * cumulativeWeight;
					if (!value.isZero())
						result->put(current->getPosition(k), &value[0]);
				}

				cumulativeWeight = proposedWeight;
				std::swap(proposed, current);

				m_sensorSampler->accept();
				m_emitterSampler->accept();
				m_directSampler->accept();
				if (largeStep) {
					largeStepRatio.incrementBase(1);
					++largeStepRatio;
				} else {
					smallStepRatio.incrementBase(1);
					++smallStepRatio;
				}
				acceptanceRate.incrementBase(1);
				++acceptanceRate;
			} else {
				for (size_t k=0; k<proposed->size(); ++k) {
					Spectrum value = proposed->getValue(k) * proposedWeight;
					if (!value.isZero())
						result->put(proposed->getPosition(k), &value[0]);
				}

				m_sensorSampler->reject();
				m_emitterSampler->reject();
				m_directSampler->reject();
				acceptanceRate.incrementBase(1);
				if (largeStep)
					largeStepRatio.incrementBase(1);
				else
					smallStepRatio.incrementBase(1);
			}
		}

		/* Perform the last splat */
		for (size_t k=0; k<current->size(); ++k) {
			Spectrum value = current->getValue(k) * cumulativeWeight;
			if (!value.isZero())
				result->put(current->getPosition(k), &value[0]);
		}


		delete current;
		delete proposed;
	}

	ref<WorkProcessor> clone() const {
		return new PSSMLTRenderer(m_config);
	}

	MTS_DECLARE_CLASS()
private:
	PSSMLTConfiguration m_config;
	ref<Scene> m_scene;
	ref<Sensor> m_sensor;
	ref<Film> m_film;
	ref<PathSampler> m_pathSampler;
	ref<PSSMLTSampler> m_origSampler;
	ref<PSSMLTSampler> m_sensorSampler;
	ref<PSSMLTSampler> m_emitterSampler;
	ref<PSSMLTSampler> m_directSampler;
	ref<ReplayableSampler> m_rplSampler;
};

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

PSSMLTProcess::PSSMLTProcess(const RenderJob *parent, RenderQueue *queue,
	const PSSMLTConfiguration &conf, const Bitmap *directImage,
	const std::vector<PathSeed> &seeds) : m_job(parent), m_queue(queue),
		m_config(conf), m_progress(NULL), m_seeds(seeds) {
	m_directImage = directImage;
	m_timeoutTimer = new Timer();
	m_refreshTimer = new Timer();
	m_resultMutex = new Mutex();
	m_resultCounter = 0;
	m_workCounter = 0;
	m_refreshTimeout = 1;
}

ref<WorkProcessor> PSSMLTProcess::createWorkProcessor() const {
	return new PSSMLTRenderer(m_config);
}

void PSSMLTProcess::develop() {
	LockGuard lock(m_resultMutex);
	size_t pixelCount = m_accum->getBitmap()->getPixelCount();
	const Spectrum *accum = (Spectrum *) m_accum->getBitmap()->getData();
	const Spectrum *direct = m_directImage != NULL ?
		(Spectrum *) m_directImage->getData() : NULL;
	const Float *importanceMap = m_config.importanceMap != NULL ?
			m_config.importanceMap->getFloatData() : NULL;
	Spectrum *target = (Spectrum *) m_developBuffer->getData();

	/* Compute the luminance correction factor */
	Float avgLuminance = 0;
	if (importanceMap) {
		for (size_t i=0; i<pixelCount; ++i)
			avgLuminance += accum[i].getLuminance() * importanceMap[i];
	} else {
		for (size_t i=0; i<pixelCount; ++i)
			avgLuminance += accum[i].getLuminance();
	}

	avgLuminance /= (Float) pixelCount;
	Float luminanceFactor = m_config.luminance / avgLuminance;

	for (size_t i=0; i<pixelCount; ++i) {
		Float correction = luminanceFactor;
		if (importanceMap)
			correction *= importanceMap[i];
		Spectrum value = accum[i] * correction;
		if (direct)
			value += direct[i];
		target[i] = value;
	}
	m_film->setBitmap(m_developBuffer);
	m_refreshTimer->reset();

	m_queue->signalRefresh(m_job);
}

void PSSMLTProcess::processResult(const WorkResult *wr, bool cancelled) {
	LockGuard lock(m_resultMutex);
	const ImageBlock *result = static_cast<const ImageBlock *>(wr);
	m_accum->put(result);
	m_progress->update(++m_resultCounter);
	m_refreshTimeout = std::min(2000U, m_refreshTimeout * 2);

	/* Re-develop the entire image every two seconds if partial results are
	   visible (e.g. in a graphical user interface). */
	if (m_job->isInteractive() && m_refreshTimer->getMilliseconds() > m_refreshTimeout)
		develop();
}

ParallelProcess::EStatus PSSMLTProcess::generateWork(WorkUnit *unit, int worker) {
	int timeout = 0;
	if (m_config.timeout > 0) {
		timeout = static_cast<int>(static_cast<int64_t>(m_config.timeout*1000) -
		          static_cast<int64_t>(m_timeoutTimer->getMilliseconds()));
	}

	if (m_workCounter >= m_config.workUnits || timeout < 0)
		return EFailure;

	SeedWorkUnit *workUnit = static_cast<SeedWorkUnit *>(unit);
	workUnit->setSeed(m_seeds[m_workCounter++]);
	workUnit->setTimeout(timeout);
	return ESuccess;
}

void PSSMLTProcess::bindResource(const std::string &name, int id) {
	ParallelProcess::bindResource(name, id);
	if (name == "sensor") {
		m_film = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id))->getFilm();
		if (m_progress)
			delete m_progress;
		m_progress = new ProgressReporter("Rendering", m_config.workUnits, m_job);
		m_accum = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize());
		m_accum->clear();
		m_developBuffer = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, m_film->getCropSize());
	}
}

MTS_IMPLEMENT_CLASS_S(PSSMLTRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(PSSMLTProcess, false, ParallelProcess)
MTS_IMPLEMENT_CLASS(SeedWorkUnit, false, WorkUnit)

MTS_NAMESPACE_END
