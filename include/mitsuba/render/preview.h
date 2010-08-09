#if !defined(__PREVIEW_H)
#define __PREVIEW_H

#include <mitsuba/render/vpl.h>

MTS_NAMESPACE_BEGIN

/**
 * Preview worker - can be used to render a quick preview of a scene
 * (illuminated by a VPL). The implementation uses coherent ray tracing 
 * when compiled in single precision and SSE is available.
 */
class MTS_EXPORT_RENDER PreviewWorker : public WorkProcessor {
public:
	inline PreviewWorker(int blockSize, Point cameraO, Vector cameraTL, 
		Vector cameraDx, Vector cameraDy, const VPL &vpl, Float minDist, bool coherent) 
		: m_blockSize(blockSize), m_cameraO(cameraO), m_cameraTL(cameraTL),
		m_cameraDx(cameraDx), m_cameraDy(cameraDy), m_vpl(vpl), 
		m_minDist(minDist), m_coherent(coherent) {
	}

	void processIncoherent(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop);

	void processCoherent(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop);

	/* WorkProcessor interface */
	ref<WorkUnit> createWorkUnit() const;
	ref<WorkResult> createWorkResult() const;
	void prepare();
	void process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop);
	void serialize(Stream *stream, InstanceManager *manager) const;
	ref<WorkProcessor> clone() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~PreviewWorker() { }
private:
	ref<Scene> m_scene;
	ref<Camera> m_camera;
	ref<KDTree> m_kdtree;
	int m_blockSize;
	Point m_cameraO;
	Vector m_cameraTL;
	Vector m_cameraDx;
	Vector m_cameraDy;
	const VPL &m_vpl;
	Float m_minDist;
	const std::vector<const Shape *> *m_shapes;
	bool m_coherent;
};

MTS_NAMESPACE_END

#endif /* __PREVIEW_H */
