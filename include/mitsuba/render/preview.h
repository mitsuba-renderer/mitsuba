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

#if !defined(__PREVIEW_H)
#define __PREVIEW_H

#include <mitsuba/render/vpl.h>

MTS_NAMESPACE_BEGIN

/**
 * Preview worker - can be used to render a quick preview of a scene
 * (illuminated by a VPL). The implementation uses coherent ray tracing 
 * when compiled in single precision and SSE is available.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER PreviewWorker : public WorkProcessor {
public:
	inline PreviewWorker(int blockSize, Point cameraO, Vector cameraTL, 
		Vector cameraDx, Vector cameraDy, const VPL &vpl, Float minDist, bool coherent,
		bool diffuseSources, bool diffuseReceivers, Float backgroundScale) 
		: m_blockSize(blockSize), m_cameraO(cameraO), m_cameraTL(cameraTL),
		m_cameraDx(cameraDx), m_cameraDy(cameraDy), m_vpl(vpl), 
		m_minDist(minDist), m_coherent(coherent), m_diffuseSources(diffuseSources),
		m_diffuseReceivers(diffuseReceivers), m_backgroundScale(backgroundScale) {
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
	ref<ShapeKDTree> m_kdtree;
	int m_blockSize;
	Point m_cameraO;
	Vector m_cameraTL;
	Vector m_cameraDx;
	Vector m_cameraDy;
	const VPL &m_vpl;
	Float m_minDist;
	const std::vector<const Shape *> *m_shapes;
	bool m_coherent;
	bool m_diffuseSources, m_diffuseReceivers;
	Float m_backgroundScale;
};

MTS_NAMESPACE_END

#endif /* __PREVIEW_H */
