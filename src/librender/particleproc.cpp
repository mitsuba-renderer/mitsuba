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
#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/range.h>

MTS_NAMESPACE_BEGIN

ParticleProcess::ParticleProcess(EMode mode, size_t workCount, size_t granularity,
        const std::string &progressText, const void *progressReporterPayload)
    : m_mode(mode), m_workCount(workCount), m_numGenerated(0),
      m_granularity(granularity), m_receivedResultCount(0) {

    /* Choose a suitable work unit granularity if none was specified */
    if (m_granularity == 0)
        m_granularity = std::max((size_t) 1, workCount /
            (16 * Scheduler::getInstance()->getWorkerCount()));

    /* Create a visual progress reporter */
    m_progress = new ProgressReporter(progressText, workCount,
        progressReporterPayload);
    m_resultMutex = new Mutex();
}

ParticleProcess::~ParticleProcess() {
    delete m_progress;
}

ParallelProcess::EStatus ParticleProcess::generateWork(WorkUnit *unit, int worker) {
    RangeWorkUnit *range = static_cast<RangeWorkUnit *>(unit);
    size_t workUnitSize;

    if (m_mode == ETrace) {
        if (m_numGenerated == m_workCount)
            return EFailure; // There is no more work

        workUnitSize = std::min(m_granularity, m_workCount - m_numGenerated);
    } else {
        if (m_receivedResultCount >= m_workCount)
            return EFailure; // There is no more work

        workUnitSize = m_granularity;
    }

    range->setRange(m_numGenerated, m_numGenerated + workUnitSize - 1);
    m_numGenerated += workUnitSize;

    return ESuccess;
}

void ParticleProcess::increaseResultCount(size_t resultCount) {
    LockGuard lock(m_resultMutex);
    m_receivedResultCount += resultCount;
    m_progress->update(m_receivedResultCount);
}

ParticleTracer::ParticleTracer(int maxDepth, int rrDepth, bool emissionEvents)
    : m_maxDepth(maxDepth), m_rrDepth(rrDepth), m_emissionEvents(emissionEvents) { }

ParticleTracer::ParticleTracer(Stream *stream, InstanceManager *manager)
    : WorkProcessor(stream, manager) {

    m_maxDepth = stream->readInt();
    m_rrDepth = stream->readInt();
    m_emissionEvents = stream->readBool();
}

void ParticleTracer::serialize(Stream *stream, InstanceManager *manager) const {
    stream->writeInt(m_maxDepth);
    stream->writeInt(m_rrDepth);
    stream->writeBool(m_emissionEvents);
}

ref<WorkUnit> ParticleTracer::createWorkUnit() const {
    return new RangeWorkUnit();
}

void ParticleTracer::prepare() {
    Scene *scene = static_cast<Scene *>(getResource("scene"));
    m_scene = new Scene(scene);
    m_sampler = static_cast<Sampler *>(getResource("sampler"));
    Sensor *newSensor = static_cast<Sensor *>(getResource("sensor"));
    m_scene->removeSensor(scene->getSensor());
    m_scene->addSensor(newSensor);
    m_scene->setSensor(newSensor);
    m_scene->initializeBidirectional();
}

void ParticleTracer::process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop) {
    const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
    MediumSamplingRecord mRec;
    Intersection its;
    ref<Sensor> sensor    = m_scene->getSensor();
    bool needsTimeSample  = sensor->needsTimeSample();
    PositionSamplingRecord pRec(sensor->getShutterOpen()
        + 0.5f * sensor->getShutterOpenTime());
    Ray ray;

    m_sampler->generate(Point2i(0));

    for (size_t index = range->getRangeStart(); index <= range->getRangeEnd() && !stop; ++index) {
        m_sampler->setSampleIndex(index);

        /* Sample an emission */
        if (needsTimeSample)
            pRec.time = sensor->sampleTime(m_sampler->next1D());

        const Emitter *emitter = NULL;
        const Medium *medium;

        Spectrum power;
        Ray ray;

        if (m_emissionEvents) {
            /* Sample the position and direction component separately to
               generate emission events */
            power = m_scene->sampleEmitterPosition(pRec, m_sampler->next2D());
            emitter = static_cast<const Emitter *>(pRec.object);
            medium = emitter->getMedium();

            /* Forward the sampling event to the attached handler */
            handleEmission(pRec, medium, power);

            DirectionSamplingRecord dRec;
            power *= emitter->sampleDirection(dRec, pRec,
                    emitter->needsDirectionSample() ? m_sampler->next2D() : Point2(0.5f));
            ray.setTime(pRec.time);
            ray.setOrigin(pRec.p);
            ray.setDirection(dRec.d);
        } else {
            /* Sample both components together, which is potentially
               faster / uses a better sampling strategy */

            power = m_scene->sampleEmitterRay(ray, emitter,
                m_sampler->next2D(), m_sampler->next2D(), pRec.time);
            medium = emitter->getMedium();
            handleNewParticle();
        }

        int depth = 1, nullInteractions = 0;
        bool delta = false;

        Spectrum throughput(1.0f); // unitless path throughput (used for russian roulette)
        while (!throughput.isZero() && (depth <= m_maxDepth || m_maxDepth < 0)) {
            m_scene->rayIntersectAll(ray, its);

            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */
            if (medium && medium->sampleDistance(Ray(ray, 0, its.t), mRec, m_sampler)) {
                /* Sample the integral
                  \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
                */

                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                /* Forward the medium scattering event to the attached handler */
                handleMediumInteraction(depth, nullInteractions,
                        delta, mRec, medium, -ray.d, throughput*power);

                PhaseFunctionSamplingRecord pRec(mRec, -ray.d, EImportance);

                throughput *= medium->getPhaseFunction()->sample(pRec, m_sampler);
                delta = false;

                ray = Ray(mRec.p, pRec.wo, ray.time);
                ray.mint = 0;
            } else if (its.t == std::numeric_limits<Float>::infinity()) {
                /* There is no surface in this direction */
                break;
            } else {
                /* Sample
                    tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
                    Account for this and multiply by the proper per-color-channel transmittance.
                */
                if (medium)
                    throughput *= mRec.transmittance / mRec.pdfFailure;

                const BSDF *bsdf = its.getBSDF();

                /* Forward the surface scattering event to the attached handler */
                handleSurfaceInteraction(depth, nullInteractions, delta, its, medium, throughput*power);

                BSDFSamplingRecord bRec(its, m_sampler, EImportance);
                Spectrum bsdfWeight = bsdf->sample(bRec, m_sampler->next2D());
                if (bsdfWeight.isZero())
                    break;

                /* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
                Vector wi = -ray.d, wo = its.toWorld(bRec.wo);
                Float wiDotGeoN = dot(its.geoFrame.n, wi),
                      woDotGeoN = dot(its.geoFrame.n, wo);
                if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                    woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    break;

                /* Keep track of the weight, medium and relative
                   refractive index along the path */
                throughput *= bsdfWeight;
                if (its.isMediumTransition())
                    medium = its.getTargetMedium(woDotGeoN);

                if (bRec.sampledType & BSDF::ENull)
                    ++nullInteractions;
                else
                    delta = bRec.sampledType & BSDF::EDelta;

#if 0
                /* This is somewhat unfortunate: for accuracy, we'd really want the
                   correction factor below to match the path tracing interpretation
                   of a scene with shading normals. However, this factor can become
                   extremely large, which adds unacceptable variance to output
                   renderings.

                   So for now, it is disabled. The adjoint particle tracer and the
                   photon mapping variants still use this factor for the last
                   bounce -- just not for the intermediate ones, which introduces
                   a small (though in practice not noticeable) amount of error. This
                   is also what the implementation of SPPM by Toshiya Hachisuka does.

                   Ultimately, we'll need better adjoint BSDF sampling strategies
                   that incorporate these extra terms */

                /* Adjoint BSDF for shading normals -- [Veach, p. 155] */
                throughput *= std::abs(
                    (Frame::cosTheta(bRec.wi) * woDotGeoN)/
                    (Frame::cosTheta(bRec.wo) * wiDotGeoN));
#endif

                ray.setOrigin(its.p);
                ray.setDirection(wo);
                ray.mint = Epsilon;
            }

            if (depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max(), (Float) 0.95f);
                if (m_sampler->next1D() >= q)
                    break;
                throughput /= q;
            }
        }
    }
}

void ParticleTracer::handleEmission(const PositionSamplingRecord &pRec,
        const Medium *medium, const Spectrum &weight) { }

void ParticleTracer::handleNewParticle() { }

void ParticleTracer::handleSurfaceInteraction(int depth, int nullInteractions,
    bool delta, const Intersection &its, const Medium *medium,
    const Spectrum &weight) { }

void ParticleTracer::handleMediumInteraction(int depth, int nullInteractions,
    bool delta, const MediumSamplingRecord &mRec, const Medium *medium,
    const Vector &wi, const Spectrum &weight) { }

MTS_IMPLEMENT_CLASS(RangeWorkUnit, false, WorkUnit)
MTS_IMPLEMENT_CLASS(ParticleProcess, true, ParallelProcess)
MTS_IMPLEMENT_CLASS(ParticleTracer, true, WorkProcessor)
MTS_NAMESPACE_END
