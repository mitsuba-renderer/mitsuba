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

#include "ptracer_proc.h"

MTS_NAMESPACE_BEGIN
	
/* ==================================================================== */
/*                           Work result impl.                          */
/* ==================================================================== */

void CaptureParticleWorkResult::load(Stream *stream) {
	Assert(sizeof(Spectrum) == sizeof(Float)*SPECTRUM_SAMPLES);
	size_t nEntries = fullSize.x * fullSize.y;
	stream->readFloatArray(reinterpret_cast<Float *>(pixels), nEntries*SPECTRUM_SAMPLES);
	m_range->load(stream);
}

void CaptureParticleWorkResult::save(Stream *stream) const {
	Assert(sizeof(Spectrum) == sizeof(Float)*SPECTRUM_SAMPLES);
	size_t nEntries = fullSize.x * fullSize.y;
	stream->writeFloatArray(reinterpret_cast<Float *>(pixels), nEntries*SPECTRUM_SAMPLES);
	m_range->save(stream);
}

/* ==================================================================== */
/*                         Work processor impl.                         */
/* ==================================================================== */

void CaptureParticleWorker::prepare() {
	ParticleTracer::prepare();
	m_camera = static_cast<Camera *>(getResource("camera"));
	m_isPerspectiveCamera = m_camera->getClass()->derivesFrom(MTS_CLASS(PerspectiveCamera));
	m_filter = m_camera->getFilm()->getTabulatedFilter();
}

ref<WorkProcessor> CaptureParticleWorker::clone() const {
	return new CaptureParticleWorker(m_maxDepth, 
		m_rrDepth);
}

ref<WorkResult> CaptureParticleWorker::createWorkResult() const {
	const Film *film = m_camera->getFilm();
	const int border = (int) std::ceil(std::max(m_filter->getFilterSize().x,
		m_filter->getFilterSize().y) - 0.5f);
	return new CaptureParticleWorkResult(film->getCropOffset(), film->getCropSize(), border);
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

void CaptureParticleWorker::handleSurfaceInteraction(int depth,
		bool caustic, const Intersection &its, const Medium *medium,
		const Spectrum &weight) {
	const ProjectiveCamera *camera = static_cast<const ProjectiveCamera *>(m_camera.get());
	Point2 screenSample;

	if (camera->positionToSample(its.p, screenSample)) {
		Point cameraPosition = camera->getPosition(screenSample);
	
		Float t = dot(camera->getImagePlaneNormal(), its.p-cameraPosition);
		if (t < camera->getNearClip() || t > camera->getFarClip())
			return;

		if (its.isMediumTransition()) 
			medium = its.getTargetMedium(cameraPosition - its.p);

		Spectrum transmittance = m_scene->getTransmittance(its.p,
				cameraPosition, its.time, medium);

		if (transmittance.isZero())
			return;

		const BSDF *bsdf = its.shape->getBSDF();
		Vector wo = cameraPosition - its.p;
		Float dist = wo.length(); wo /= dist;

		BSDFQueryRecord bRec(its, its.toLocal(wo), EImportance);

		Float importance; 
		if (m_isPerspectiveCamera)
			importance = ((const PerspectiveCamera *) camera)->importance(screenSample) / (dist * dist);
		else
			importance = 1/camera->areaDensity(screenSample);

		Vector wi = its.toWorld(its.wi);

		/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
		Float wiDotGeoN = dot(its.geoFrame.n, wi),
			  woDotGeoN = dot(its.geoFrame.n, wo);
		if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 || 
			woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
			return;

		/* Adjoint BSDF for shading normals -- [Veach, p. 155] */
		Float correction = std::abs(
			(Frame::cosTheta(bRec.wi) * woDotGeoN)/
			(Frame::cosTheta(bRec.wo) * wiDotGeoN));

		/* Splat onto the accumulation buffer */
		Ray ray(its.p, wo, 0, dist, its.time);
		Spectrum sampleVal = weight * bsdf->eval(bRec) 
			* transmittance * (importance * correction);

		m_workResult->splat(screenSample, sampleVal, m_filter);
	}
}

void CaptureParticleWorker::handleMediumInteraction(int depth, bool caustic,
		const MediumSamplingRecord &mRec, const Medium *medium,
		Float time, const Vector &wi, const Spectrum &weight) {
	const ProjectiveCamera *camera = static_cast<const ProjectiveCamera *>(m_camera.get());
	Point2 screenSample;

	if (camera->positionToSample(mRec.p, screenSample)) {
		Point cameraPosition = camera->getPosition(screenSample);

		Float t = dot(camera->getImagePlaneNormal(), mRec.p-cameraPosition);
		if (t < camera->getNearClip() || t > camera->getFarClip())
			return;

		Spectrum transmittance = m_scene->getTransmittance(mRec.p,
			cameraPosition, time, medium);

		if (transmittance.isZero())
			return;

		Vector wo = cameraPosition - mRec.p;
		Float dist = wo.length(); wo /= dist;

		Float importance; 
		if (m_isPerspectiveCamera)
			importance = ((const PerspectiveCamera *) camera)->importance(screenSample) / (dist * dist);
		else
			importance = 1/camera->areaDensity(screenSample);

		/* Splat onto the accumulation buffer */
		Ray ray(mRec.p, wo, 0, dist, time);

		Spectrum sampleVal = weight * medium->getPhaseFunction()->eval(
			  PhaseFunctionQueryRecord(mRec, wi, wo, EImportance)) 
			  * transmittance * importance;

		m_workResult->splat(screenSample, sampleVal, m_filter);
	}
}

/* ==================================================================== */
/*                        Parallel process impl.                        */
/* ==================================================================== */

void CaptureParticleProcess::develop() {
	float *accumImageData = m_accumBitmap->getFloatData();
	float *finalImageData = m_finalBitmap->getFloatData();
	size_t size = m_accumBitmap->getWidth() * m_accumBitmap->getHeight() * 4;
	Float weight = (m_accumBitmap->getWidth() * m_accumBitmap->getHeight()) 
		/ (Float) m_receivedResultCount;
	for (size_t i=0; i<size; i+=4) {
		for (int j=0; j<3; ++j) 
			finalImageData[i+j] = accumImageData[i+j] * weight;
		finalImageData[i+3] = 1.0f;
	}
	m_film->fromBitmap(m_finalBitmap);

	m_queue->signalRefresh(m_job, m_finalBitmap);
}

void CaptureParticleProcess::processResult(const WorkResult *wr, bool cancelled) {
	const CaptureParticleWorkResult *result 
		= static_cast<const CaptureParticleWorkResult *>(wr);
	const RangeWorkUnit *range = result->getRangeWorkUnit();
	if (cancelled) 
		return;

	m_resultMutex->lock();
	increaseResultCount(range->getSize());

	/* Accumulate the received pixel data */
	float *imageData = m_accumBitmap->getFloatData();
	size_t pixelIndex, imagePixelIndex = 0;
	Float r, g, b;
	Vector2i start(result->getBorder(), result->getBorder());
	Vector2i end(result->getFullSize().x - result->getBorder(), 
		result->getFullSize().y - result->getBorder());

	for (int y=start.y; y<end.y; ++y) {
		pixelIndex = y*result->getFullSize().x + start.x;
		for (int x=start.x; x<end.x; ++x) {
			Spectrum spec(result->getPixel(pixelIndex));
			spec.toLinearRGB(r,g,b);
			imageData[imagePixelIndex++] += r;
			imageData[imagePixelIndex++] += g;
			imageData[imagePixelIndex++] += b;
			++imagePixelIndex;
			++pixelIndex;
		}
	}

	develop();

	m_resultMutex->unlock();
}

void CaptureParticleProcess::bindResource(const std::string &name, int id) {
	if (name == "camera") {
		Camera *camera = static_cast<Camera *>(Scheduler::getInstance()->getResource(id));
		m_film = camera->getFilm();
		const Vector2i res(m_film->getCropSize());
		m_accumBitmap = new Bitmap(res.x, res.y, 128);
		m_finalBitmap = new Bitmap(res.x, res.y, 128);
		m_accumBitmap->clear();
	}
	ParticleProcess::bindResource(name, id);
}

ref<WorkProcessor> CaptureParticleProcess::createWorkProcessor() const {
	return new CaptureParticleWorker(m_maxDepth-1, m_rrDepth);
}

MTS_IMPLEMENT_CLASS(CaptureParticleProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS(CaptureParticleWorkResult, false, ImageBlock)
MTS_IMPLEMENT_CLASS_S(CaptureParticleWorker, false, ParticleTracer)
MTS_NAMESPACE_END

