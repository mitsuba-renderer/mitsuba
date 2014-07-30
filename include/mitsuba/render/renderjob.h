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

#pragma once
#if !defined(__MITSUBA_RENDER_RENDERJOB_H_)
#define __MITSUBA_RENDER_RENDERJOB_H_

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderqueue.h>

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
	 * When the Resource ID parameters (\c sceneResID, \c sensorResID, ..) are
	 * set to \c -1, the implementation will automatically register the
	 * associated objects (scene, sensor, sampler) with the scheduler and
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
	 * \param sensorResID
	 *     Resource ID of \c scene->getSensor() (or \c -1)
	 * \param samplerResID
	 *     Resource ID of the sample generator (or \c -1)
	 * \param threadIsCritical
	 *     When set to \c true, the entire program will terminate
	 *     if this thread fails unexpectedly.
	 * \param interactive
	 *     Are partial results of the rendering process visible, e.g. in
	 *     a graphical user interface?
	 */
	RenderJob(const std::string &threadName,
		Scene *scene, RenderQueue *queue,
		int sceneResID = -1,
		int sensorResID = -1,
		int samplerResID = -1,
		bool threadIsCritical = true,
		bool interactive = false);

	/// Write out the current (partially rendered) image
	inline void flush() { m_scene->flush(m_queue, this); }

	/// Cancel a running render job
	inline void cancel() { m_scene->cancel(); }

	/// Wait for the job to finish and return whether it was successful
	inline bool wait() { join(); return !m_cancelled; }

	/**
	 * \brief Are partial results of the rendering process visible, e.g. in
	 * a graphical user interface?
	 *
	 * Some integrators may choose to invest more time on generating high-quality
	 * intermediate results in this case.
	 */
	inline bool isInteractive() const { return m_interactive; }

	/// Define whether or not this is an interactive job
	inline void setInteractive(bool interactive) { m_interactive = interactive; }

	/// Get a pointer to the underlying scene
	inline Scene *getScene() { return m_scene.get(); }

	/// Get a pointer to the underlying scene (const version)
	inline const Scene *getScene() const { return m_scene.get(); }

	/// Get a pointer to the underlying render queue
	inline RenderQueue *getRenderQueue() { return m_queue.get(); }

	/// Get a pointer to the underlying render queue (const version)
	inline const RenderQueue *getRenderQueue() const { return m_queue.get(); }

	/// Return the amount of time spent rendering the given job (in seconds)
	inline Float getRenderTime() const { return m_queue->getRenderTime(this); }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~RenderJob();
	/// Run method
	void run();
private:
	ref<Scene> m_scene;
	ref<RenderQueue> m_queue;
	int m_sceneResID, m_samplerResID, m_sensorResID;
	bool m_ownsSceneResource;
	bool m_ownsSensorResource;
	bool m_ownsSamplerResource;
	bool m_cancelled;
	bool m_interactive;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_RENDERJOB_H_ */
