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

#if !defined(__PTRACER_PROC_H)
#define __PTRACER_PROC_H

#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                             Work result                              */
/* ==================================================================== */

/**
 * Packages the result of a particle tracing work unit. Contains
 * the range of traced particles plus a snapshot of the camera film.
 */
class CaptureParticleWorkResult : public ImageBlock {
public:
	inline CaptureParticleWorkResult(const Point2i &offset, const Vector2i &res, int border) 
	 : ImageBlock(res, border, false, false, false, false) {
		setOffset(offset);
		setSize(res);
		m_range = new RangeWorkUnit();
	}

	inline const RangeWorkUnit *getRangeWorkUnit() const {
		return m_range.get();
	}

	inline void setRangeWorkUnit(const RangeWorkUnit *range) {
		m_range->set(range);
	}

	/* Work unit implementation */
	void load(Stream *stream);
	void save(Stream *stream) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CaptureParticleWorkResult() { }
protected:
	ref<RangeWorkUnit> m_range;
};


/* ==================================================================== */
/*                             Work processor                           */
/* ==================================================================== */

/**
 * Particle tracing worker - looks for volume and surface interactions
 * and tries to accumulate the resulting information at the image plane.
 */
class CaptureParticleWorker : public ParticleTracer {
public:
	inline CaptureParticleWorker(int maxDepth, int rrDepth)
		: ParticleTracer(maxDepth, rrDepth) { }

	inline CaptureParticleWorker(Stream *stream, InstanceManager *manager) 
	 : ParticleTracer(stream, manager) { }

	/* ParticleTracer impl. */
	void prepare();
	ref<WorkProcessor> clone() const;
	ref<WorkResult> createWorkResult() const;
	void process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop);

	/**
	 * Handles particles interacting with a surface - if a reflection to the
	 * camera is possible, compute the importance and accumulate in the proper 
	 * pixel of the accumulation buffer.
	 */
	void handleSurfaceInteraction(int depth, bool caustic,
		const Intersection &its, const Medium *medium,
		const Spectrum &weight);

	/**
	 * Handles particles interacting with a medium - if a reflection to the
	 * camera is possible, compute the importance and accumulate in the proper 
	 * pixel of the accumulation buffer.
	 */
	void handleMediumInteraction(int depth, bool caustic,
			const MediumSamplingRecord &mRec, const Medium *medium,
			Float time, const Vector &wi, const Spectrum &weight);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CaptureParticleWorker() { }
private:
	ref<const Camera> m_camera;
	ref<const TabulatedFilter> m_filter;
	ref<CaptureParticleWorkResult> m_workResult;
	bool m_isPerspectiveCamera;
};

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */
/**
 * Parallel particle tracing process - used to run this over
 * a group of machines
 */
class CaptureParticleProcess : public ParticleProcess {
public:
	CaptureParticleProcess(const RenderJob *job, RenderQueue *queue, 
		size_t sampleCount, size_t granularity, int maxDepth, int rrDepth)
		: ParticleProcess(ParticleProcess::ETrace, sampleCount, 
		  granularity, "Rendering", job), m_job(job), m_queue(queue), 
		  m_maxDepth(maxDepth), m_rrDepth(rrDepth) {
	}

	void develop();

	/* ParallelProcess impl. */
	void processResult(const WorkResult *wr, bool cancelled);
	void bindResource(const std::string &name, int id);
	ref<WorkProcessor> createWorkProcessor() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CaptureParticleProcess() { }
private:
	ref<const RenderJob> m_job;
	ref<RenderQueue> m_queue;
	ref<Film> m_film;
	ref<Bitmap> m_accumBitmap;
	ref<Bitmap> m_finalBitmap;
	int m_maxDepth;
	int m_rrDepth;
};

MTS_NAMESPACE_END

#endif /* __PTRACER_PROC_H */
