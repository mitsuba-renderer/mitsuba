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

#include "ptracer_proc.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                           Work result impl.                          */
/* ==================================================================== */

void CaptureParticleWorkResult::load(Stream *stream) {
	size_t nEntries = (size_t) m_size.x * (size_t) m_size.y;
	stream->readFloatArray(reinterpret_cast<Float *>(m_bitmap->getFloatData()),
		nEntries * SPECTRUM_SAMPLES);
	m_range->load(stream);
}

void CaptureParticleWorkResult::save(Stream *stream) const {
	size_t nEntries = (size_t) m_size.x * (size_t) m_size.y;
	stream->writeFloatArray(reinterpret_cast<const Float *>(m_bitmap->getFloatData()),
		nEntries * SPECTRUM_SAMPLES);
	m_range->save(stream);
}

/* ==================================================================== */
/*                         Work processor impl.                         */
/* ==================================================================== */

CaptureParticleWorker::CaptureParticleWorker(Stream *stream, InstanceManager *manager)
  : ParticleTracer(stream, manager) {
	  m_maxPathDepth = stream->readInt();
	  m_bruteForce = stream->readBool();
}

void CaptureParticleWorker::serialize(Stream *stream, InstanceManager *manager) const {
	ParticleTracer::serialize(stream, manager);
	stream->writeInt(m_maxPathDepth);
	stream->writeBool(m_bruteForce);
}

void CaptureParticleWorker::prepare() {
	ParticleTracer::prepare();
	m_sensor = static_cast<Sensor *>(getResource("sensor"));
	m_rfilter = m_sensor->getFilm()->getReconstructionFilter();
}

ref<WorkProcessor> CaptureParticleWorker::clone() const {
	return new CaptureParticleWorker(m_maxDepth,
		m_maxPathDepth, m_rrDepth, m_bruteForce);
}

ref<WorkResult> CaptureParticleWorker::createWorkResult() const {
	const Film *film = m_sensor->getFilm();
	return new CaptureParticleWorkResult(film->getCropSize(), m_rfilter.get());
}

void CaptureParticleWorker::process(const WorkUnit *workUnit, WorkResult *workResult,
	const bool &stop) {
	const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
	m_workResult = static_cast<CaptureParticleWorkResult *>(workResult);
	m_workResult->setRangeWorkUnit(range);
	m_workResult->clear();
	ParticleTracer::process(workUnit, workResult, stop);
	m_workResult = NULL;
}

void CaptureParticleWorker::handleEmission(const PositionSamplingRecord &pRec,
		const Medium *medium, const Spectrum &weight) {
	if (m_bruteForce)
		return;

	DirectSamplingRecord dRec(pRec.p, pRec.time);
	int maxInteractions = m_maxPathDepth - 1;

	/* Create a dummy intersection to ensure that sampleAttenuatedSensorDirect()
	   treats the light source vertex as being located on a surface */
	Intersection its;
	its.p = pRec.p;

	Spectrum value = weight * m_scene->sampleAttenuatedSensorDirect(
			dRec, its, medium, maxInteractions, m_sampler->next2D(), m_sampler);

	if (value.isZero())
		return;

	const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
	value *= emitter->evalDirection(DirectionSamplingRecord(dRec.d), pRec);

	/* Splat onto the accumulation buffer */
	m_workResult->put(dRec.uv, (Float *) &value[0]);
}

void CaptureParticleWorker::handleSurfaceInteraction(int depth, int nullInteractions,
		bool caustic, const Intersection &its, const Medium *medium,
		const Spectrum &weight) {

	if (its.isSensor()) {
		if (!m_bruteForce && !caustic)
			return;

		const Sensor *sensor = its.shape->getSensor();
		if (sensor != m_sensor)
			return;

		Vector wi = its.toWorld(its.wi);
		Point2 uv;
		Spectrum value = sensor->eval(its, wi, uv) * weight;
		if (value.isZero())
			return;

		m_workResult->put(uv, (Float *) &value[0]);
		return;
	}

	if (m_bruteForce || (depth >= m_maxPathDepth && m_maxPathDepth > 0))
		return;

	int maxInteractions = m_maxPathDepth - depth - 1;

	DirectSamplingRecord dRec(its);
	Spectrum value = weight * m_scene->sampleAttenuatedSensorDirect(
			dRec, its, medium, maxInteractions,
			m_sampler->next2D(), m_sampler);

	if (value.isZero())
		return;

	const BSDF *bsdf = its.getBSDF();

	Vector wo = dRec.d;
	BSDFSamplingRecord bRec(its, its.toLocal(wo), EImportance);

	/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
	Vector wi = its.toWorld(its.wi);
	Float wiDotGeoN = dot(its.geoFrame.n, wi),
		  woDotGeoN = dot(its.geoFrame.n, wo);
	if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
		woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
		return;

	/* Adjoint BSDF for shading normals -- [Veach, p. 155] */
	Float correction = std::abs(
		(Frame::cosTheta(bRec.wi) * woDotGeoN)/
		(Frame::cosTheta(bRec.wo) * wiDotGeoN));
	value *= bsdf->eval(bRec) * correction;

	/* Splat onto the accumulation buffer */
	m_workResult->put(dRec.uv, (Float *) &value[0]);
}

void CaptureParticleWorker::handleMediumInteraction(int depth, int nullInteractions, bool caustic,
		const MediumSamplingRecord &mRec, const Medium *medium, const Vector &wi,
		const Spectrum &weight) {

	if (m_bruteForce || (depth >= m_maxPathDepth && m_maxPathDepth > 0))
		return;

	DirectSamplingRecord dRec(mRec);

	int maxInteractions = m_maxPathDepth - depth - 1;

	Spectrum value = weight * m_scene->sampleAttenuatedSensorDirect(
		dRec, medium, maxInteractions, m_sampler->next2D(), m_sampler);

	if (value.isZero())
		return;

	/* Evaluate the phase function */
	const PhaseFunction *phase = medium->getPhaseFunction();
	PhaseFunctionSamplingRecord pRec(mRec, wi, dRec.d, EImportance);
	value *= phase->eval(pRec);

	if (value.isZero())
		return;

	/* Splat onto the accumulation buffer */
	m_workResult->put(dRec.uv, (Float *) &value[0]);
}

/* ==================================================================== */
/*                        Parallel process impl.                        */
/* ==================================================================== */

void CaptureParticleProcess::develop() {
	Float weight = (m_accum->getWidth() * m_accum->getHeight())
		/ (Float) m_receivedResultCount;
	m_film->setBitmap(m_accum->getBitmap(), weight);
	m_queue->signalRefresh(m_job);
}

void CaptureParticleProcess::processResult(const WorkResult *wr, bool cancelled) {
	const CaptureParticleWorkResult *result
		= static_cast<const CaptureParticleWorkResult *>(wr);
	const RangeWorkUnit *range = result->getRangeWorkUnit();
	if (cancelled)
		return;

	LockGuard lock(m_resultMutex);
	increaseResultCount(range->getSize());
	m_accum->put(result);
	if (m_job->isInteractive() || m_receivedResultCount == m_workCount)
		develop();
}

void CaptureParticleProcess::bindResource(const std::string &name, int id) {
	if (name == "sensor") {
		Sensor *sensor = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id));
		m_film = sensor->getFilm();
		m_accum = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize(), NULL);
		m_accum->clear();
	}
	ParticleProcess::bindResource(name, id);
}

ref<WorkProcessor> CaptureParticleProcess::createWorkProcessor() const {
	return new CaptureParticleWorker(m_maxDepth, m_maxPathDepth,
			m_rrDepth, m_bruteForce);
}

MTS_IMPLEMENT_CLASS(CaptureParticleProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS(CaptureParticleWorkResult, false, ImageBlock)
MTS_IMPLEMENT_CLASS_S(CaptureParticleWorker, false, ParticleTracer)
MTS_NAMESPACE_END

