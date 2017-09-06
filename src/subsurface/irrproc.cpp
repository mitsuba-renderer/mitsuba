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
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/range.h>
#include "irrproc.h"

MTS_NAMESPACE_BEGIN

/* Parallel irradiance sampling implementation (worker) */
class IrradianceSamplingWorker : public WorkProcessor {
public:
    IrradianceSamplingWorker(int irrSamples, bool irrIndirect, Float time)
        : m_irrSamples(irrSamples), m_irrIndirect(irrIndirect), m_time(time) {
    }

    IrradianceSamplingWorker(Stream *stream, InstanceManager *manager) {
        m_irrSamples = stream->readInt();
        m_irrIndirect = stream->readBool();
        m_time = stream->readFloat();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        stream->writeInt(m_irrSamples);
        stream->writeBool(m_irrIndirect);
        stream->writeFloat(m_time);
    }

    ref<WorkUnit> createWorkUnit() const {
        return new PositionSampleVector();
    }

    ref<WorkResult> createWorkResult() const {
        return new IrradianceSampleVector();
    }

    void prepare() {
        m_scene = static_cast<Scene *>(getResource("scene"));
        m_sampler = static_cast<Sampler *>(getResource("sampler"));
        m_integrator = static_cast<SamplingIntegrator *>(getResource("integrator"));
        m_scene->wakeup(NULL, m_resources);
        m_integrator->wakeup(NULL, m_resources);
    }

    void process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop) {
        const PositionSampleVector &positions = *static_cast<const PositionSampleVector *>(workUnit);
        IrradianceSampleVector *result = static_cast<IrradianceSampleVector *>(workResult);
        const SamplingIntegrator *integrator = m_integrator.get();

        result->clear();

        for (size_t i=0; i<positions.size(); ++i) {
            /* Create a fake intersection record */
            const PositionSample &sample = positions[i];
            Intersection its;
            its.p = sample.p;
            its.shFrame = Frame(sample.n);
            its.shape = m_scene->getShapes()[sample.shapeIndex].get();
            its.time = m_time;
            its.hasUVPartials = false;

            result->put(IrradianceSample(
                its.p,
                integrator->E(m_scene.get(), its, its.shape->getExteriorMedium(), m_sampler,
                    m_irrSamples, m_irrIndirect)
            ));
        }
    }

    ref<WorkProcessor> clone() const {
        return new IrradianceSamplingWorker(m_irrSamples, m_irrIndirect, m_time);
    }

    MTS_DECLARE_CLASS()
protected:
    virtual ~IrradianceSamplingWorker() { }
private:
    ref<Scene> m_scene;
    ref<Sampler> m_sampler;
    ref<SamplingIntegrator> m_integrator;
    int m_irrSamples;
    bool m_irrIndirect;
    Float m_time;
};

void PositionSampleVector::load(Stream *stream) {
    clear();
    size_t count = stream->readSize();
    m_samples.resize(count);
    for (size_t i=0; i<count; ++i)
        m_samples[i] = PositionSample(stream);
}

void PositionSampleVector::save(Stream *stream) const {
    stream->writeSize(m_samples.size());
    for (size_t i=0; i<m_samples.size(); ++i)
        m_samples[i].serialize(stream);
}

void PositionSampleVector::set(const WorkUnit *workUnit) {
    m_samples = ((PositionSampleVector *) workUnit)->m_samples;
}

std::string PositionSampleVector::toString() const {
    std::ostringstream oss;
    oss << "PositionSampleVector[size="
        << m_samples.size() << "]";
    return oss.str();
}

void IrradianceSampleVector::load(Stream *stream) {
    clear();
    size_t count = stream->readSize();
    m_samples.resize(count);
    for (size_t i=0; i<count; ++i)
        m_samples[i] = IrradianceSample(stream);
}

void IrradianceSampleVector::save(Stream *stream) const {
    stream->writeSize(m_samples.size());
    for (size_t i=0; i<m_samples.size(); ++i)
        m_samples[i].serialize(stream);
}

std::string IrradianceSampleVector::toString() const {
    std::ostringstream oss;
    oss << "IrradianceSampleVector[size="
        << m_samples.size() << "]";
    return oss.str();
}

IrradianceSamplingProcess::IrradianceSamplingProcess(PositionSampleVector *positions,
        size_t granularity, int irrSamples, bool irrIndirect, Float time,
        const void *data)
    : m_positionSamples(positions), m_granularity(granularity),
      m_irrSamples(irrSamples), m_irrIndirect(irrIndirect), m_time(time) {
    m_resultMutex = new Mutex();
    m_irradianceSamples = new IrradianceSampleVector();
    m_irradianceSamples->reserve(positions->size());
    m_samplesRequested = 0;
    m_progress = new ProgressReporter("Sampling irradiance", positions->size(), data);
}

IrradianceSamplingProcess::~IrradianceSamplingProcess() {
    if (m_progress)
        delete m_progress;
}

ref<WorkProcessor> IrradianceSamplingProcess::createWorkProcessor() const {
    return new IrradianceSamplingWorker(m_irrSamples, m_irrIndirect, m_time);
}

ParallelProcess::EStatus IrradianceSamplingProcess::generateWork(WorkUnit *unit, int worker) {
    if (m_samplesRequested == m_positionSamples->size())
        return EFailure;

    /* Reserve a sequence of at most 'granularity' samples */
    size_t workSize = std::min(m_granularity, m_positionSamples->size() - m_samplesRequested);

    std::vector<PositionSample> &samples = static_cast<PositionSampleVector *>(unit)->get();
    const std::vector<PositionSample> &source = m_positionSamples->get();

    samples.clear();
    samples.insert(samples.begin(),
            source.begin() + m_samplesRequested,
            source.begin() + m_samplesRequested + workSize);
    m_samplesRequested += workSize;

    return ESuccess;
}

void IrradianceSamplingProcess::processResult(const WorkResult *wr, bool cancelled) {
    const IrradianceSampleVector *result = static_cast<const IrradianceSampleVector *>(wr);
    LockGuard lock(m_resultMutex);
    for (size_t i=0; i<result->size(); ++i)
        m_irradianceSamples->put((*result)[i]);
    m_progress->update(m_irradianceSamples->size());
}

MTS_IMPLEMENT_CLASS(PositionSampleVector, false, WorkUnit);
MTS_IMPLEMENT_CLASS(IrradianceSampleVector, false, WorkResult);
MTS_IMPLEMENT_CLASS_S(IrradianceSamplingWorker, false, WorkProcessor);
MTS_IMPLEMENT_CLASS(IrradianceSamplingProcess, false, ParallelProcess);
MTS_NAMESPACE_END
