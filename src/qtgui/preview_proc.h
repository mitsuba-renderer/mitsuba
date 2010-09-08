#if !defined(__PREVIEW_PROC_H)
#define __PREVIEW_PROC_H

#include <mitsuba/render/scene.h>
#include <mitsuba/render/imageproc.h>
#include <mitsuba/render/preview.h>
#include <mitsuba/core/bitmap.h>

using namespace mitsuba;

/**
 * Parallel process for rendering a preview frame using 
 * realtime coherent ray tracing.
 */
class PreviewProcess : public BlockedImageProcess {
public:
	PreviewProcess(const Scene *scene, int sceneResID, int blockSize);

	void configure(const VPL &vpl, Float minDist, const Point2 &jitter, 
		const Bitmap *source, Bitmap *target, bool coherent,
		bool diffuseSources, bool diffuseReceivers); 
	inline int getRayCount() const { return m_numRays; }
	inline const Scene *getScene() const { return m_scene; }

	/* ParallelProcess interface */
	ref<WorkProcessor> createWorkProcessor() const;
	void processResult(const WorkResult *result, bool cancelled);
	bool isLocal() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~PreviewProcess() { }
private:
	const Bitmap *m_source;
	Bitmap *m_target;
	const Scene *m_scene;
	const Film *m_film;
	Point m_cameraO;
	Vector m_cameraTL, m_cameraDx, m_cameraDy;
	int m_numRays;
	ref<Mutex> m_mutex;
	const VPL *m_vpl;
	Float m_minDist;
	bool m_coherent;
	bool m_diffuseSources;
	bool m_diffuseReceivers;
};

#endif /* __PREVIEW_PROC_H */
