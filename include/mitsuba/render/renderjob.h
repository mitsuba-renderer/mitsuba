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

#if !defined(__RENDERJOB_H)
#define __RENDERJOB_H

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/render/testcase.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Coordinates the process of rendering a single image. 
 *
 * Implemented as a thread so that multiple jobs can 
 * be executed concurrently.
 *
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER RenderJob : public Thread {
public:
	/**
	 * \brief Create a new render job for the given scene. 
	 *
	 * When the Resource ID parameters (\c sceneResID, \c cameraResID, ..) are
	 * set to \c -1, the implementation will automatically register the
	 * associated objects (scene, camera, sampler) with the scheduler and 
	 * forward copies to all involved network rendering workers. When some 
	 * of these resources have already been registered with
	 * the scheduler, their IDs can be provided to avoid this extra
	 * communication cost.
	 *
	 * \param threadName
	 *     Thread name identifier for this render job
	 * \param scene
	 *     Scene to be rendered
	 * \param queue
	 *     Pointer to a queue, to which this job should be added
	 * \param sceneResID
	 *     Resource ID of \c scene (or \c -1)
	 * \param cameraResID
	 *     Resource ID of \c scene->getCamera() (or \c -1)
	 * \param samplerResID
	 *     Resource ID of the sample generator (or \c -1)
	 * \param threadIsCritical
	 *     When set to \c true, the entire program will terminate 
	 *     if this thread fails unexpectedly.
	 * \param testSupervisor
	 *     When this image is being rendered as part of a test suite,
	 *     this parameter points to the associated \ref TestSupervisor
	 *     instance.
	 */
	RenderJob(const std::string &threadName, 
		Scene *scene, RenderQueue *queue,
		int sceneResID = -1,
		int cameraResID = -1,
		int samplerResID = -1,
		bool threadIsCritical = true,
		TestSupervisor *testSupervisor = NULL);

	/// Write out the current (partially rendered) image
	inline void flush() { m_scene->flush(); }

	/// Cancel a running render job
	inline void cancel() { m_scene->cancel(); }

	/// Wait for the job to finish and return whether it was successful
	inline bool wait() { join(); return !m_cancelled; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~RenderJob();
	/// Run method
	void run();
private:
	ref<Scene> m_scene;
	ref<RenderQueue> m_queue;
	ref<TestSupervisor> m_testSupervisor;
	int m_sceneResID, m_samplerResID, m_cameraResID;
	bool m_ownsSceneResource;
	bool m_ownsCameraResource;
	bool m_ownsSamplerResource;
	bool m_cancelled;
};

MTS_NAMESPACE_END

#endif /* __RENDERJOB_H */
