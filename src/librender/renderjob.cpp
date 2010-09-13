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

#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/renderproc.h>

#ifdef MTS_DEBUG_FP
#include <fenv.h>
#endif

MTS_NAMESPACE_BEGIN
	
RenderJob::RenderJob(const std::string &threadName, 
	Scene *scene, RenderQueue *queue, TestSupervisor *testSupervisor,
	int sceneResID, int cameraResID, int samplerResID, bool threadIsCritical) : Thread(threadName),
	m_scene(scene), m_queue(queue), m_testSupervisor(testSupervisor) {

	/* Optional: bring the process down when this thread crashes */
	setCritical(threadIsCritical); 

	m_queue->addJob(this);
	ref<Scheduler> sched = Scheduler::getInstance();

	ref<Camera> camera = m_scene->getCamera();
	ref<Sampler> sampler = m_scene->getSampler();

	/* Register the scene with the scheduler if needed */
	if (sceneResID == -1) {
		m_sceneResID = sched->registerResource(m_scene);
		m_ownsSceneResource = true; 
	} else {
		m_sceneResID = sceneResID;
		m_ownsSceneResource = false; 
	}

	/* Register the camera with the scheduler if needed */
	if (cameraResID == -1) {
		m_cameraResID = sched->registerResource(camera);
		m_ownsCameraResource = true;
	} else {
		m_cameraResID = cameraResID;
		m_ownsCameraResource = false;
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
		m_samplerResID = sched->registerManifoldResource(samplers); 
		for (size_t i=0; i<sched->getCoreCount(); ++i)
			samplers[i]->decRef();
		m_ownsSamplerResource = true;
	} else {
		m_samplerResID = samplerResID;
		m_ownsSamplerResource = false;
	}
}

RenderJob::~RenderJob() {
	Scheduler *sched = Scheduler::getInstance();
	if (m_ownsSceneResource)
		sched->unregisterResource(m_sceneResID);
	if (m_ownsSamplerResource)
		sched->unregisterResource(m_samplerResID);
	if (m_ownsCameraResource)
		sched->unregisterResource(m_cameraResID);
}

void RenderJob::run() {
	ref<Film> film = m_scene->getFilm();
	ref<Sampler> sampler = m_scene->getSampler();

	if (m_testSupervisor.get()) {
		if (film->getClass()->getName() != "MFilm")
			Log(EError, "Only the MATLAB M-file film is supported when "
				"running in test case mode!");
		if (m_scene->getTestType() == Scene::ETTest) {
			if (film->getTabulatedFilter()->getName() != "BoxFilter")
				Log(EError, "Only the box reconstruction filter is supported when "
					"performing a t-test in test case mode!");
			if (!m_scene->getIntegrator()->getClass()->derivesFrom(SampleIntegrator::m_theClass)) 
				Log(EError, "Only sampling-based integrators are supported when "
					"performing a t-test in test case mode!");
		}
	}

	bool cancelled = false;
	m_scene->initialize();
	try {
		m_scene->preprocess(m_queue, this, m_sceneResID, m_cameraResID, m_samplerResID);

		if (!m_scene->render(m_queue, this, m_sceneResID, m_cameraResID, m_samplerResID)) {
			cancelled = true;
			Log(EWarn, "Rendering of scene \"%s\" did not complete successfully!",
				m_scene->getSourceFile().leaf().c_str());
		}
		m_scene->postprocess(m_queue, this, m_sceneResID, m_cameraResID, m_samplerResID);

		if (m_testSupervisor.get()) 
			m_testSupervisor->analyze(m_scene);
	} catch (const std::exception &ex) {
		Log(EWarn, "Rendering of scene \"%s\" did not complete successfully, caught exception: %s",
			m_scene->getSourceFile().leaf().c_str(), ex.what());
		cancelled = true;
	}

	m_queue->removeJob(this, cancelled);
}

MTS_IMPLEMENT_CLASS(RenderJob, false, Thread)
MTS_NAMESPACE_END
