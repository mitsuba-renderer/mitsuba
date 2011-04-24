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
		bool diffuseSources, bool diffuseReceivers, Float backgroundScale); 
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
	Float m_backgroundScale;
};

#endif /* __PREVIEW_PROC_H */
