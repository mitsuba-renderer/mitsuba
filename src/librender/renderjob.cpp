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

#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/renderproc.h>
#include <boost/filesystem.hpp>

MTS_NAMESPACE_BEGIN

RenderJob::RenderJob(const std::string &threadName,
	Scene *scene, RenderQueue *queue, int sceneResID, int sensorResID,
	int samplerResID, bool threadIsCritical, bool interactive)
	: Thread(threadName), m_scene(scene), m_queue(queue), m_interactive(interactive) {

	/* Optional: bring the process down when this thread crashes */
	setCritical(threadIsCritical);

	m_queue->addJob(this);
	ref<Scheduler> sched = Scheduler::getInstance();

	ref<Sensor> sensor = m_scene->getSensor();
	ref<Sampler> sampler = m_scene->getSampler();

	/* Register the scene with the scheduler if needed */
	if (sceneResID == -1) {
		m_sceneResID = sched->registerResource(m_scene);
		m_ownsSceneResource = true;
	} else {
		m_sceneResID = sceneResID;
		m_ownsSceneResource = false;
	}

	/* Register the sensor with the scheduler if needed */
	if (sensorResID == -1) {
		m_sensorResID = sched->registerResource(sensor);
		m_ownsSensorResource = true;
	} else {
		m_sensorResID = sensorResID;
		m_ownsSensorResource = false;
	}

	/* Register the sampler with the scheduler if needed */
	if (samplerResID == -1) {
		/* Create a sampler instance for every core */
		std::vector<SerializableObject *> samplers(sched->getCoreCount());
		for (size_t i=0; i<sched->getCoreCount(); ++i) {
			ref<Sampler> clonedSampler = sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}
		m_samplerResID = sched->registerMultiResource(samplers);
		for (size_t i=0; i<sched->getCoreCount(); ++i)
			samplers[i]->decRef();
		m_ownsSamplerResource = true;
	} else {
		m_samplerResID = samplerResID;
		m_ownsSamplerResource = false;
	}
	m_cancelled = false;
}

RenderJob::~RenderJob() {
	Scheduler *sched = Scheduler::getInstance();
	if (m_ownsSceneResource)
		sched->unregisterResource(m_sceneResID);
	if (m_ownsSamplerResource)
		sched->unregisterResource(m_samplerResID);
	if (m_ownsSensorResource)
		sched->unregisterResource(m_sensorResID);
}

void RenderJob::run() {
	ref<Film> film = m_scene->getFilm();
	ref<Sampler> sampler = m_scene->getSampler();
	m_cancelled = false;

	try {
		m_scene->getFilm()->setDestinationFile(m_scene->getDestinationFile(),
			m_scene->getBlockSize());

		if (!m_scene->preprocess(m_queue, this, m_sceneResID, m_sensorResID, m_samplerResID)) {
			m_cancelled = true;
			Log(EWarn, "Preprocessing of scene \"%s\" did not complete successfully!",
				m_scene->getSourceFile().filename().string().c_str());
		}

		if (!m_cancelled) {
			if (!m_scene->render(m_queue, this, m_sceneResID, m_sensorResID, m_samplerResID)) {
				m_cancelled = true;
				Log(EWarn, "Rendering of scene \"%s\" did not complete successfully!",
					m_scene->getSourceFile().filename().string().c_str());
			}
			Log(EInfo, "Render time: %s", timeString(m_queue->getRenderTime(this), true).c_str());
			m_scene->postprocess(m_queue, this, m_sceneResID, m_sensorResID, m_samplerResID);
		}
	} catch (const std::exception &ex) {
		Log(EWarn, "Rendering of scene \"%s\" did not complete successfully, caught exception: %s",
			m_scene->getSourceFile().filename().string().c_str(), ex.what());
		m_cancelled = true;
	}

	m_queue->removeJob(this, m_cancelled);
}

MTS_IMPLEMENT_CLASS(RenderJob, false, Thread)
MTS_NAMESPACE_END
