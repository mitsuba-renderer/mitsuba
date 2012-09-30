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

#include <mitsuba/bidir/util.h>
#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include "erpt_proc.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{erpt}{Energy redistribution path tracing}
 * \order{11}
 * \cite{Cline2005Energy}
 */
class EnergyRedistributionPathTracing : public Integrator {
public:
	EnergyRedistributionPathTracing(const Properties &props) : Integrator(props) {
		m_config.maxDepth = props.getInteger("maxDepth", -1);

		/* Specifies the number of Markov Chains that, on average, are
		   started per pixel */
		m_config.numChains = props.getFloat("numChains", 1.0f);
		m_config.maxChains = props.getInteger("maxChains", 0);

		/* Specifies the number of mutations to be performed in each
		   Markov Chain */
		m_config.chainLength = props.getInteger("chainLength", 100);

		m_config.directSamples = props.getInteger("directSamples", 16);
		m_config.separateDirect = m_config.directSamples >= 0;

		/* Number of samples used to estimate the average contribution of a
		   single sample. Usually, this parameter can be left untouched. */
		m_config.luminanceSamples = props.getInteger("luminanceSamples", 15000);
	
		/* Selectively enable/disable the bidirectional mutation */
		m_config.bidirectionalMutation = props.getBoolean("bidirectionalMutation", false);

		/* Selectively enable/disable the lens perturbation */
		m_config.lensPerturbation = props.getBoolean("lensPerturbation", true);

		/* Selectively enable/disable the caustic perturbation */
		m_config.causticPerturbation = props.getBoolean("causticPerturbation", true);

		/* Selectively enable/disable the multi-chain perturbation */
		m_config.multiChainPerturbation = props.getBoolean("multiChainPerturbation", true);

		/* Selectively enable/disable the manifold perturbation */ 
		m_config.manifoldPerturbation = props.getBoolean("manifoldPerturbation", false);
		m_config.probFactor = props.getFloat("probFactor", 50);
		m_config.avgAngleChangeSurface = props.getFloat("avgAngleChangeSurface", 0);
		m_config.avgAngleChangeMedium = props.getFloat("avgAngleChangeMedium", 0);

		if (m_config.maxDepth <= 0 && m_config.maxDepth != -1)
			Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
	}

	/// Unserialize from a binary data stream
	EnergyRedistributionPathTracing(Stream *stream, InstanceManager *manager)
	 : Integrator(stream, manager) {
		m_config = ERPTConfiguration(stream);
	}

	virtual ~EnergyRedistributionPathTracing() { }

	void configureSampler(const Scene *scene, Sampler *sampler) { 
		/* Prepare the sampler for tile-based rendering */
		sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		m_config.serialize(stream);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, 
			const RenderJob *job, int sceneResID, int sensorResID,
			int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID,
				sensorResID, samplerResID);

		if (scene->getSubsurfaceIntegrators().size() > 0)
			Log(EError, "Subsurface integrators are not supported by ERPT!");

		return true;
	}

	void cancel() {
		ref<RenderJob> nested = m_nestedJob;
		if (nested)
			nested->cancel();
		Scheduler::getInstance()->cancel(m_process);
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sensor> sensor = scene->getSensor();
		ref<Sampler> sampler = sensor->getSampler();
		const Film *film = sensor->getFilm();
		size_t nCores = sched->getCoreCount();
		size_t sampleCount = sampler->getSampleCount();

		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT
			" %s, " SSE_STR ", " SIZE_T_FMT " samples/pixel) ..", 
			film->getCropSize().x, film->getCropSize().y, nCores, 
			nCores == 1 ? "core" : "cores", sampleCount);

		ref<Sampler> indepSampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Sampler), Properties("independent")));
		indepSampler->configure();

		ref<PathSampler> pathSampler = new PathSampler(PathSampler::EBidirectional, scene, 
			indepSampler, indepSampler, indepSampler, m_config.maxDepth, 10,
			m_config.separateDirect, true, true);

		m_config.luminance = pathSampler->computeAverageLuminance(
				m_config.luminanceSamples);
		m_config.blockSize = scene->getBlockSize();

		m_config.dump();

		ref<Bitmap> directImage;
		if (m_config.separateDirect && m_config.directSamples > 0) {
			directImage = BidirectionalUtils::renderDirectComponent(scene, 
				sceneResID, sensorResID, queue, job, m_config.directSamples);
			if (directImage == NULL)
				return false;
		}

		ref<ERPTProcess> process = new ERPTProcess(job, queue, m_config, directImage);

		/* Create an independent sampler for use by the MLT chains */
		std::vector<SerializableObject *> indepSamplers(sched->getCoreCount());
		for (size_t i=0; i<sched->getCoreCount(); ++i) {
			ref<Sampler> clonedSampler = indepSampler->clone();
			clonedSampler->incRef();
			indepSamplers[i] = clonedSampler.get();
		}
		int indepSamplerResID = sched->registerMultiResource(indepSamplers); 
		for (size_t i=0; i<indepSamplers.size(); ++i)
			indepSamplers[i]->decRef();
	
		process->bindResource("scene", sceneResID);
		process->bindResource("sensor", sensorResID);
		process->bindResource("sampler", samplerResID);
		process->bindResource("indepSampler", indepSamplerResID);

		m_process = process;
		sched->schedule(process);
		sched->wait(process);
		m_process = NULL;

		sched->unregisterResource(indepSamplerResID);

		return process->getReturnStatus() == ParallelProcess::ESuccess;
	}

	MTS_DECLARE_CLASS()
private:
	ref<ParallelProcess> m_process;
	ref<RenderJob> m_nestedJob;
	ERPTConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(EnergyRedistributionPathTracing, false, Integrator)
MTS_EXPORT_PLUGIN(EnergyRedistributionPathTracing, "Energy redistribution path tracing");
MTS_NAMESPACE_END
