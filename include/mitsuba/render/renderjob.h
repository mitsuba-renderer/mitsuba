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

#if !defined(__RENDERJOB_H)
#define __RENDERJOB_H

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/render/testcase.h>

MTS_NAMESPACE_BEGIN

/**
 * Render job - coordinates the process of rendering a single
 * image. Implemented as a thread so that multiple jobs can
 * be executed concurrently.
 */
class MTS_EXPORT_RENDER RenderJob : public Thread {
public:
	/**
	 * Create a new render job for the given scene. When the
	 * scene, sampler or camera are already registered with the scheduler, 
	 * the last parameters can optionally be specified (that way 
	 * they do not have to be re-sent to network rendering servers).
	 */
	RenderJob(const std::string &threadName, 
		Scene *scene, RenderQueue *queue,
		TestSupervisor *testSupervisor,
		int sceneResID = -1,
		int cameraResID = -1,
		int samplerResID = -1,
		bool threadIsCritical = true);

	/// Write out the current (partially rendered) image
	inline void flush() { m_scene->flush(); }

	/// Cancel a running render job
	inline void cancel() { m_scene->cancel(); }

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
};

MTS_NAMESPACE_END

#endif /* __RENDERJOB_H */
