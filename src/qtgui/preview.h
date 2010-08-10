#if !defined(__PREVIEW_H)
#define __PREVIEW_H

#include <mitsuba/core/platform.h>
#include <QtGui>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/render/preview.h>
#include "common.h"
#include "preview_proc.h"

using namespace mitsuba;

class SceneContext;

/**
 * Asynchronous preview rendering thread
 */
class PreviewThread : public QObject, public Thread {
	Q_OBJECT
public:
	PreviewThread(Device *parentDevice, Renderer *parentRenderer);

	/**
	 * Change the scene context. 
	 */
	void setSceneContext(SceneContext *context, bool swapContext, bool motion);

	/**
	 * Resume rendering 
	 */
	void resume();

	/**
	 * Wait until the thread has started
	 */
	inline void waitUntilStarted() { m_started->wait(); }

	/**
	 * Acquire the best current approximation of the rendered scene.
	 * If that takes longer than 'ms' milliseconds, a failure occurs
	 */
	PreviewQueueEntry acquireBuffer(int ms = -1);

	/// Return the buffer to the renderer
	void releaseBuffer(PreviewQueueEntry &entry);
	
	/// Terminate the preview thread
	void quit();
signals:
	void caughtException(const QString &what);
	void statusMessage(const QString &status);
protected:
	/// Preview thread main loop
	virtual void run();
	/// Virtual destructure
	virtual ~PreviewThread();
	/// Render a single VPL using OpenGL
	void oglRenderVPL(PreviewQueueEntry &target, const VPL &vpl);
	/// Render a single VPL using real-time coherent ray tracing
	void rtrtRenderVPL(PreviewQueueEntry &target, const VPL &vpl);
private:
	ref<Session> m_session;
	ref<Device> m_device, m_parentDevice;
	ref<Renderer> m_renderer, m_parentRenderer;
	ref<VPLShaderManager> m_shaderManager;
	ref<GPUTexture> m_framebuffer;
	ref<GPUProgram> m_accumProgram;
	ref<PreviewProcess> m_previewProc;
	int m_accumProgramParam_source1;
	int m_accumProgramParam_source2;
	const GPUTexture *m_accumBuffer;
	ref<Mutex> m_mutex;
	ref<FileResolver> m_fileResolver;
	ref<Random> m_random;
	ref<ConditionVariable> m_queueCV;
	ref<Timer> m_timer;
	ref<WaitFlag> m_started;
	std::list<PreviewQueueEntry> m_readyQueue, m_recycleQueue;
	SceneContext *m_context;
	int m_vplSampleOffset, m_minVPLs, m_vplCount;
	int m_vplsPerSecond, m_raysPerSecond;
	int m_bufferCount, m_queueEntryIndex;
	std::deque<VPL> m_vpls;
	std::vector<GPUTexture *> m_releaseList;
	Point m_camPos;
	Transform m_camViewTransform;
	bool m_quit, m_sleep, m_motion, m_useSync;
	bool m_refreshScene;
};

#endif /* __PREVIEW_H */
