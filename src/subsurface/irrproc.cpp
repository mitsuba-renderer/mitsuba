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
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/range.h>
#include "irrproc.h"

MTS_NAMESPACE_BEGIN

/* Parallel irradiance sampling implementation (worker) */
class IrradianceSamplingWorker : public WorkProcessor {
public:
	IrradianceSamplingWorker(size_t sampleCount, int ssIndex, int irrSamples, bool irrIndirect) 
		: m_sampleCount(sampleCount), m_ssIndex(ssIndex), 
		m_irrSamples(irrSamples), m_irrIndirect(irrIndirect) {
	}

	IrradianceSamplingWorker(Stream *stream, InstanceManager *manager) {
		m_sampleCount = (size_t) stream->readLong();
		m_ssIndex = stream->readInt();
		m_irrSamples = stream->readInt();
		m_irrIndirect = stream->readBool();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		stream->writeLong(m_sampleCount);
		stream->writeInt(m_ssIndex);
		stream->writeInt(m_irrSamples);
		stream->writeBool(m_irrIndirect);
	}

	ref<WorkUnit> createWorkUnit() const {
		return new RangeWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new IrradianceRecordVector();
	}

	void prepare() {
		m_scene = static_cast<Scene *>(getResource("scene"));
		m_integrator = static_cast<SampleIntegrator *>(m_scene->getIntegrator());
		Properties props;
		props.setLong("sampleCount", m_sampleCount);
		props.setPluginName("hammersley");
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		props.setPluginName("independent");
		m_independentSampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), props));
		m_scene->wakeup(m_resources);
		const Subsurface *ss = m_scene->getSubsurfaceIntegrators()[m_ssIndex];
		m_shapes = ss->getShapes();
		for (size_t i=0; i<m_shapes.size(); ++i)
			m_areaPDF.put(m_shapes[i]->getSurfaceArea());
		m_areaPDF.build();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop) {
		const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
		IrradianceRecordVector *result = static_cast<IrradianceRecordVector *>(workResult);
		const SampleIntegrator *integrator = m_integrator.get();
		ref<Camera> camera = m_scene->getCamera();

		result->clear();
		for (size_t i=range->getRangeStart(); i<=range->getRangeEnd(); ++i) {
			m_sampler->setSampleIndex(i);
			Point2 sample = m_sampler->next2D();

			Float expSamples;
			unsigned int index = m_areaPDF.sampleReuse(sample.x, expSamples);
			expSamples *= m_sampleCount;
			ShapeSamplingRecord sRec;
			Float pdf = m_shapes[index]->sampleArea(sRec, sample) * expSamples;
			Float time = camera->getShutterOpen() + m_sampler->next1D() * camera->getShutterOpenTime();

			result->put(IrradianceSample(
				sRec.p,
				integrator->E(m_scene.get(), sRec.p, sRec.n, time, 
					camera->getMedium(), m_independentSampler,
					m_irrSamples, m_irrIndirect), 1/pdf
			));
		}
	}

	ref<WorkProcessor> clone() const {
		return new IrradianceSamplingWorker(m_sampleCount, m_ssIndex, m_irrSamples, m_irrIndirect);
	}

	MTS_DECLARE_CLASS()
protected:
	virtual ~IrradianceSamplingWorker() { }
private:
	ref<Scene> m_scene;
	ref<Camera> m_camera;
	ref<Sampler> m_sampler, m_independentSampler;
	ref<SampleIntegrator> m_integrator;
	DiscretePDF m_areaPDF;
	std::vector<Shape *> m_shapes;
	size_t m_sampleCount;
	int m_ssIndex;
	int m_irrSamples;
	bool m_irrIndirect;
};

void IrradianceRecordVector::load(Stream *stream) {
	clear();
	size_t count = stream->readSize();
	m_samples.resize(count);
	for (size_t i=0; i<count; ++i)
		m_samples[i] = IrradianceSample(stream);
}

void IrradianceRecordVector::save(Stream *stream) const {
	stream->writeSize(m_samples.size());
	for (size_t i=0; i<m_samples.size(); ++i)
		m_samples[i].serialize(stream);
}

std::string IrradianceRecordVector::toString() const {
	std::ostringstream oss;
	oss << "IrradianceRecordVector[size="
		<< m_samples.size() << "]";
	return oss.str();
}

IrradianceSamplingProcess::IrradianceSamplingProcess(size_t sampleCount, 
	size_t granularity, int ssIndex, int irrSamples, bool irrIndirect,
	const void *progressReporterPayload) 
	: m_sampleCount(sampleCount), m_granularity(granularity), m_ssIndex(ssIndex), 
	m_irrSamples(irrSamples), m_irrIndirect(irrIndirect) {
	m_resultCount = 0;
	m_resultMutex = new Mutex();
	m_samples = new IrradianceRecordVector();
	m_samplesRequested = 0;
	m_progress = new ProgressReporter("Sampling irradiance", sampleCount,
		progressReporterPayload); 
}

IrradianceSamplingProcess::~IrradianceSamplingProcess() {
	if (m_progress)
		delete m_progress;
}

ref<WorkProcessor> IrradianceSamplingProcess::createWorkProcessor() const {
	return new IrradianceSamplingWorker(m_sampleCount, m_ssIndex, 
			m_irrSamples, m_irrIndirect);
}

ParallelProcess::EStatus IrradianceSamplingProcess::generateWork(WorkUnit *unit, int worker) {
	if (m_samplesRequested == m_sampleCount)
		return EFailure;

	/* Reserve a sequence of at most 'granularity' samples */
	size_t workSize = std::min(m_granularity, m_sampleCount - m_samplesRequested);
	RangeWorkUnit *range = static_cast<RangeWorkUnit *>(unit);
	range->setRange(m_samplesRequested, m_samplesRequested + workSize - 1);
	m_samplesRequested += workSize;

	return ESuccess;
}

void IrradianceSamplingProcess::processResult(const WorkResult *wr, bool cancelled) {
	const IrradianceRecordVector *result = static_cast<const IrradianceRecordVector *>(wr);
	m_resultMutex->lock();
	for (size_t i=0; i<result->size(); ++i)
		m_samples->put((*result)[i]);
	m_resultCount += result->size();
	m_progress->update(m_resultCount);
	m_resultMutex->unlock();
}

MTS_IMPLEMENT_CLASS(IrradianceRecordVector, false, WorkResult);
MTS_IMPLEMENT_CLASS_S(IrradianceSamplingWorker, false, WorkProcessor);
MTS_IMPLEMENT_CLASS(IrradianceSamplingProcess, false, ParallelProcess);
MTS_NAMESPACE_END
