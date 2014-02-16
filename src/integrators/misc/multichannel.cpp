/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{direct}{Multi-channel integrator}
 * \order{16}
 *
 */

class MultiChannelIntegrator : public SamplingIntegrator {
public:
	MultiChannelIntegrator(const Properties &props) : SamplingIntegrator(props) { }

	MultiChannelIntegrator(Stream *stream, InstanceManager *manager)
	 : SamplingIntegrator(stream, manager) {
		m_integrators.resize(stream->readSize());
		for (size_t i=0; i<m_integrators.size(); ++i)
			m_integrators[i] = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);

		stream->writeSize(m_integrators.size());
		for (size_t i=0; i<m_integrators.size(); ++i)
			manager->serialize(stream, m_integrators[i].get());
	}

	bool preprocess(const Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID,
		int samplerResID) {
		if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID,
				sensorResID, samplerResID))
			return false;
		for (size_t i=0; i<m_integrators.size(); ++i) {
			if (!m_integrators[i]->preprocess(scene, queue, job, sceneResID,
					sensorResID, samplerResID))
				return false;
		}
		return true;
	}

	bool render(Scene *scene,
			RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
		ref<Film> film = sensor->getFilm();

		size_t nCores = sched->getCoreCount();
		const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
		size_t sampleCount = sampler->getSampleCount();

		if (m_integrators.empty())
			Log(EError, "No sub-integrators were supplied to the multi-channel integrator!");

		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
			sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
			nCores == 1 ? "core" : "cores");

		/* This is a sampling-based integrator - parallelize */
		ref<BlockedRenderProcess> proc = new BlockedRenderProcess(job,
			queue, scene->getBlockSize());

		proc->setPixelFormat(Bitmap::EMultiSpectrumAlphaWeight, (int) (m_integrators.size() * SPECTRUM_SAMPLES + 2), false);

		int integratorResID = sched->registerResource(this);
		proc->bindResource("integrator", integratorResID);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("sensor", sensorResID);
		proc->bindResource("sampler", samplerResID);
		scene->bindUsedResources(proc);
		bindUsedResources(proc);
		sched->schedule(proc);

		m_process = proc;
		sched->wait(proc);
		m_process = NULL;
		sched->unregisterResource(integratorResID);

		return proc->getReturnStatus() == ParallelProcess::ESuccess;
	}

	void renderBlock(const Scene *scene,
			const Sensor *sensor, Sampler *sampler, ImageBlock *block,
			const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {

		Float diffScaleFactor = 1.0f /
			std::sqrt((Float) sampler->getSampleCount());

		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();

		RadianceQueryRecord rRec(scene, sampler);
		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;
		RayDifferential sensorRay;

		block->clear();

		uint32_t queryType = RadianceQueryRecord::ESensorRay;
		Float *temp = (Float *) alloca(sizeof(Float) * (m_integrators.size() * SPECTRUM_SAMPLES + 1));

		for (size_t i = 0; i<points.size(); ++i) {
			Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
			if (stop)
				break;

			sampler->generate(offset);

			for (size_t j = 0; j<sampler->getSampleCount(); j++) {
				rRec.newQuery(queryType, sensor->getMedium());
				Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

				if (needsApertureSample)
					apertureSample = rRec.nextSample2D();

				if (needsTimeSample)
					timeSample = rRec.nextSample1D();

				Spectrum spec = sensor->sampleRayDifferential(
					sensorRay, samplePos, apertureSample, timeSample);

				sensorRay.scaleDifferential(diffScaleFactor);
				rRec.rayIntersect(sensorRay);

				int offset = 0;
				for (size_t k = 0; k<m_integrators.size(); ++k) {
					RadianceQueryRecord rRec2(rRec);
					rRec2.its = rRec.its;
					Spectrum result = spec * m_integrators[k]->Li(sensorRay, rRec2);
					for (int l = 0; l<SPECTRUM_SAMPLES; ++l)
						temp[offset++] = result[l];
				}
				temp[offset++] = rRec.alpha;
				temp[offset] = 1.0f;
				block->put(samplePos, temp);
				sampler->advance();
			}
		}
	}

	void bindUsedResources(ParallelProcess *proc) const {
		SamplingIntegrator::bindUsedResources(proc);

		for (size_t i=0; i<m_integrators.size(); ++i)
			m_integrators[i]->bindUsedResources(proc);
	}

	void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
		SamplingIntegrator::wakeup(parent, params);

		for (size_t i=0; i<m_integrators.size(); ++i)
			m_integrators[i]->wakeup(parent, params);
	}


	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);
		for (size_t i=0; i<m_integrators.size(); ++i)
			m_integrators[i]->configureSampler(scene, sampler);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator))) {
			SamplingIntegrator *integrator = static_cast<SamplingIntegrator *>(child);
			m_integrators.push_back(integrator);
			integrator->incRef();
		} else {
			SamplingIntegrator::addChild(name, child);
		}
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		NotImplementedError("Li");
	}

	const Integrator *getSubIntegrator(int idx) const {
		if (idx < 0 || idx >= (int) m_integrators.size())
			return NULL;
		return m_integrators[idx].get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MultiChannelIntegrator[" << endl
			<< "  integrators = {" << endl;
		for (size_t i=0; i<m_integrators.size(); ++i)
			oss << "    " << indent(m_integrators[i]->toString(), 2) << "," << endl;
		oss << "  }" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref_vector<SamplingIntegrator> m_integrators;
};

MTS_IMPLEMENT_CLASS_S(MultiChannelIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MultiChannelIntegrator, "Multi-channel integrator");
MTS_NAMESPACE_END
