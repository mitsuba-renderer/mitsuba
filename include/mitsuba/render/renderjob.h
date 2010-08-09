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
	ref<FileResolver> m_fileResolver;
	int m_sceneResID, m_samplerResID, m_cameraResID;
	bool m_ownsSceneResource;
	bool m_ownsCameraResource;
	bool m_ownsSamplerResource;
};

MTS_NAMESPACE_END

#endif /* __RENDERJOB_H */
