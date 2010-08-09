#if !defined(__RENDERPROC_H)
#define __RENDERPROC_H

#include <mitsuba/render/scene.h>
#include <mitsuba/render/imageproc.h>
#include <mitsuba/render/renderqueue.h>

MTS_NAMESPACE_BEGIN

class RenderJob;

/**
 * Parallel process for rendering with sampling-based integrators. Splits
 * an image into square pixel regions, which can be processed independently.
 */
class MTS_EXPORT_RENDER BlockedRenderProcess : public BlockedImageProcess {
public:
	BlockedRenderProcess(const RenderJob *parent, RenderQueue *queue, 
		int blockSize);

	/* ParallelProcess interface */
	ref<WorkProcessor> createWorkProcessor() const;
	void processResult(const WorkResult *result, bool cancelled);
	void bindResource(const std::string &name, int id);
	EStatus generateWork(WorkUnit *unit, int worker);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~BlockedRenderProcess();
protected:
	ref<RenderQueue> m_queue;
	ref<Scene> m_scene;
	ref<Film> m_film;
	const RenderJob *m_parent;
	int m_resultCount;
	ref<Mutex> m_resultMutex;
	ProgressReporter *m_progress;
	int m_borderSize;
};

MTS_NAMESPACE_END

#endif /* __RENDERPROC_H */
