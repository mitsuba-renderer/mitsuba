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
 * \brief Packages the result of a particle tracing work unit. Contains
 * the range of traced particles plus a snapshot of the sensor film.
 */
class CaptureParticleWorkResult : public ImageBlock {
public:
	inline CaptureParticleWorkResult(const Vector2i &res, const ReconstructionFilter *filter)
	 : ImageBlock(Bitmap::ESpectrum, res, filter) {
		setOffset(Point2i(0, 0));
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
 * \brief Particle tracing worker -- looks for volume and surface interactions
 * and tries to accumulate the resulting information at the image plane.
 */
class CaptureParticleWorker : public ParticleTracer {
public:
	inline CaptureParticleWorker(int maxDepth, int maxPathDepth,
		int rrDepth, bool bruteForce) : ParticleTracer(maxDepth, rrDepth, true),
		m_maxPathDepth(maxPathDepth), m_bruteForce(bruteForce) { }

	CaptureParticleWorker(Stream *stream, InstanceManager *manager);

	void serialize(Stream *stream, InstanceManager *manager) const;

	void prepare();
	ref<WorkProcessor> clone() const;
	ref<WorkResult> createWorkResult() const;
	void process(const WorkUnit *workUnit, WorkResult *workResult,
		const bool &stop);

	/**
	 * \brief Handles particles emitted by a light source
	 *
	 * If a connection to the sensor is possible, compute the importance
	 * and accumulate in the proper pixel of the accumulation buffer.
	 */
	void handleEmission(const PositionSamplingRecord &pRec,
			const Medium *medium, const Spectrum &weight);

	/**
	 * \brief Handles particles interacting with a surface
	 *
	 * If a connection to the sensor is possible, compute the importance
	 * and accumulate in the proper pixel of the accumulation buffer.
	 */
	void handleSurfaceInteraction(int depth, int nullInteractions, bool caustic,
		const Intersection &its, const Medium *medium,
		const Spectrum &weight);

	/**
	 * \brief Handles particles interacting with a medium
	 *
	 * If a connection to the sensor is possible, compute the importance
	 * and accumulate in the proper pixel of the accumulation buffer.
	 */
	void handleMediumInteraction(int depth, int nullInteractions, bool caustic,
			const MediumSamplingRecord &mRec, const Medium *medium,
			const Vector &wi, const Spectrum &weight);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CaptureParticleWorker() { }
private:
	ref<const Sensor> m_sensor;
	ref<const ReconstructionFilter> m_rfilter;
	ref<CaptureParticleWorkResult> m_workResult;
	int m_maxPathDepth;
	bool m_bruteForce;
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
			size_t sampleCount, size_t granularity, int maxDepth,
			int maxPathDepth, int rrDepth, bool bruteForce)
		: ParticleProcess(ParticleProcess::ETrace, sampleCount,
		  granularity, "Rendering", job), m_job(job), m_queue(queue),
		  m_maxDepth(maxDepth), m_maxPathDepth(maxPathDepth),
		  m_rrDepth(rrDepth), m_bruteForce(bruteForce) {
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
	ref<ImageBlock> m_accum;
	int m_maxDepth;
	int m_maxPathDepth;
	int m_rrDepth;
	bool m_bruteForce;
};

MTS_NAMESPACE_END

#endif /* __PTRACER_PROC_H */
