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

#include <mitsuba/core/statistics.h>
#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/range.h>

MTS_NAMESPACE_BEGIN

ParticleProcess::ParticleProcess(EMode mode, size_t workCount, 
	size_t granularity, const std::string &progressText,
	const void *progressReporterPayload) : 
	m_mode(mode), m_workCount(workCount), m_numGenerated(0), 
	m_granularity(granularity), m_receivedResultCount(0) {
	m_progress = new ProgressReporter(progressText, workCount, 
		progressReporterPayload);
	m_resultMutex = new Mutex();
}

ParticleProcess::~ParticleProcess() {
	delete m_progress;
}

ParallelProcess::EStatus ParticleProcess::generateWork(WorkUnit *unit, int worker) {
	if (m_mode == ETrace) {
		if (m_numGenerated == m_workCount)
			return EFailure;
		/* Reserve a sequence of at most 'granularity' particles */
		size_t workSize = std::min(m_granularity, m_workCount - m_numGenerated);
		RangeWorkUnit *range = static_cast<RangeWorkUnit *>(unit);
		range->setRange(m_numGenerated, m_numGenerated + workSize - 1);
		m_numGenerated += workSize;
	} else {
		if (m_receivedResultCount >= m_workCount)
			return EFailure;
		/* Reserve a sequence of exactly 'granularity' particles */
		RangeWorkUnit *range = static_cast<RangeWorkUnit *>(unit);
		range->setRange(m_numGenerated, m_numGenerated + m_granularity - 1);
		m_numGenerated += m_granularity;
	}

	return ESuccess;
}

void ParticleProcess::increaseResultCount(size_t resultCount) {
	m_resultMutex->lock();
	m_receivedResultCount += resultCount;
	m_progress->update(m_receivedResultCount);
	m_resultMutex->unlock();
}

ParticleTracer::ParticleTracer(Stream *stream, InstanceManager *manager)
	: WorkProcessor(stream, manager) {

	m_maxDepth = stream->readInt();	
	m_multipleScattering = stream->readBool();	
	m_rrDepth = stream->readInt();	
}

void ParticleTracer::serialize(Stream *stream, InstanceManager *manager) const {
	stream->writeInt(m_maxDepth);
	stream->writeBool(m_multipleScattering);
	stream->writeInt(m_rrDepth);
}

ref<WorkUnit> ParticleTracer::createWorkUnit() const {
	return new RangeWorkUnit();
}

void ParticleTracer::prepare() {
	m_scene = static_cast<Scene *>(getResource("scene"));
	m_sampler = static_cast<Sampler *>(getResource("sampler"));
}

void ParticleTracer::process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop) {	
	const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
	EmissionRecord eRec;
	MediumSamplingRecord mRec;
	Ray ray;
	Intersection its;
	Spectrum weight, bsdfVal;
	int depth;
	bool caustic;
	ref<Camera> camera = m_scene->getCamera();
	Float shutterOpen     = camera->getShutterOpen(), 
		  shutterOpenTime = camera->getShutterOpenTime();

	m_sampler->generate();
	for (size_t index = range->getRangeStart(); index <= range->getRangeEnd() && !stop; ++index) {
		m_sampler->setSampleIndex(index);
		Point2 areaSample = m_sampler->next2D(), 
		       dirSample  = m_sampler->next2D();	   

		/* Sample an emitted particle */
		m_scene->sampleEmission(eRec, areaSample, dirSample);

		ray = Ray(eRec.sRec.p, eRec.d, shutterOpen + shutterOpenTime * m_sampler->next1D());
		weight = eRec.P;
		depth = 1;
		caustic = true;

		while (!weight.isBlack() && (depth <= m_maxDepth || m_maxDepth < 0)) {
			m_scene->rayIntersect(ray, its); 

            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */
			if (m_scene->sampleDistance(ray, its.t, mRec, m_sampler)) {
				Vector wo;
				/* Sample the integral
				  \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/

				weight *= mRec.sigmaS * mRec.attenuation / mRec.pdf;
				handleMediumInteraction(depth, caustic, mRec, ray.time, -ray.d, weight);
	
				if (!m_multipleScattering)
					break;

				PhaseFunction::ESampledType sampledType; // ignored
				weight *= mRec.medium->getPhaseFunction()->sample(mRec, -ray.d, 
					wo, sampledType, m_sampler->next2D());
				caustic = false;

				/* Russian roulette */
				if (depth >= m_rrDepth) {
					if (m_sampler->next1D() > mRec.albedo)
						break;
					else
						weight /= mRec.albedo;
				}

				ray = Ray(mRec.p, wo, ray.time);
			} else if (its.t == std::numeric_limits<Float>::infinity()) {
				/* There is no surface in this direction */
				break;
			} else {
				ray.mint = 0; ray.maxt = its.t;

				/* Sample 
					tau(x, y) * (Surface integral). This happens with probability mRec.pdf
					Divide this out and multiply with the proper per color channel attenuation.
				*/
				weight *= m_scene->getAttenuation(ray) / mRec.pdf;
				handleSurfaceInteraction(depth, caustic, its, weight);

				const BSDF *bsdf = its.shape->getBSDF();
				BSDFQueryRecord bRec(its, m_sampler->next2D());
				bRec.quantity = EImportance;
				bsdfVal = bsdf->sampleCos(bRec);
				if (bsdfVal.isBlack())
					break;

				/* Russian roulette */
				if (depth >= m_rrDepth) {
					/* Assuming that BSDF importance sampling is perfect,
					   the following should equal the maximum albedo
					   over all spectral samples */
					Float approxAlbedo = std::min((Float) 1, bsdfVal.max());
					if (m_sampler->next1D() > approxAlbedo)
						break;
					else
						weight /= approxAlbedo;
				}

				weight *= bsdfVal;
				Vector wi = -ray.d, wo = its.toWorld(bRec.wo);
				ray = Ray(its.p, wo, ray.time);

				/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
				Float wiDotGeoN = dot(its.geoFrame.n, wi),
						woDotGeoN = dot(its.geoFrame.n, wo);
				if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 || 
					woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
					break;

				/* Adjoint BSDF for shading normals -- [Veach, p. 155] */
				weight *= std::abs(
					(Frame::cosTheta(bRec.wi) * woDotGeoN)/
					(Frame::cosTheta(bRec.wo) * wiDotGeoN));

				caustic &= (bRec.sampledType & BSDF::EDelta) ? true : false;

			}
			++depth;
		}
	}
}

void ParticleTracer::handleSurfaceInteraction(int depth, bool caustic,
	const Intersection &its, const Spectrum &weight) {
}

void ParticleTracer::handleMediumInteraction(int depth, bool caustic,
	const MediumSamplingRecord &mRec, Float time, const Vector &wi, 
	const Spectrum &weight) {
}

MTS_IMPLEMENT_CLASS(RangeWorkUnit, false, WorkUnit)
MTS_IMPLEMENT_CLASS(ParticleProcess, true, ParallelProcess)
MTS_IMPLEMENT_CLASS(ParticleTracer, true, WorkProcessor)
MTS_NAMESPACE_END
