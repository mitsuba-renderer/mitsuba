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

#if !defined(__PARTICLEPROC_H)
#define __PARTICLEPROC_H

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Abstract parallel particle tracing process - distributes the work of 
 * particle tracing and gathers results as computed by a subclass of 
 * <tt>ParticleTracer</tt>
 */
class MTS_EXPORT_RENDER ParticleProcess : public ParallelProcess {
public:
	enum EMode {
		/// Trace a specified amount of photons and then stop
		ETrace = 0,

		/**
		 * Trace an unspecified amount of photons until a certain 
		 * number of events has been recorded. Due to asynchronous
		 * nature introduced by parallelism, this comes at the
		 * risk of distributing a bit too much work.
		 */
		EGather
	};

	/* ParallelProcess interface */
	virtual EStatus generateWork(WorkUnit *unit, int worker);

	MTS_DECLARE_CLASS()
protected:

	/**
	 * Create a new particle process 
	 *
	 * @param mode
	 *    Particle tracing mode - see above
	 * @param workCount
	 *    Total # of particles to trace / # events to record
	 * @param granularity
	 *    Number of particles in each work unit
	 * @param progressText
	 *    Title of the progress bar
	 * @param progressReporterPayload
	 *    Custom pointer payload to be delivered with progress messages
	 */
	ParticleProcess(EMode mode, size_t workCount,
		size_t granularity, const std::string &progressText,
		const void *progressReporterPayload);

	void increaseResultCount(size_t resultCount);

	/// Virtual destructor
	virtual ~ParticleProcess();
protected:
	EMode m_mode;
	ProgressReporter *m_progress;
	size_t m_workCount;
	size_t m_numGenerated;
	size_t m_granularity;
	ref<Mutex> m_resultMutex;
	size_t m_receivedResultCount;
};

/**
 * Abstract work processor - traces particles and performs a
 * customizable action every time a surface or volume interaction occurs.
 */
class MTS_EXPORT_RENDER ParticleTracer : public WorkProcessor {
public:
	/* WorkProcessor interface */
	virtual ref<WorkUnit> createWorkUnit() const;
	virtual void prepare();
	virtual void process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop);
	void serialize(Stream *stream, InstanceManager *manager) const;

	/**
     * Abstract method - handle a surface interation, 
	 * which occurred while tracing particles.
     */
	virtual void handleSurfaceInteraction(int depth, bool caustic,
		const Intersection &its, const Spectrum &weight);

	/**
     * Abstract method - handle an interation with a participating medium, 
	 * which occurred while tracing particles.
     */
	virtual void handleMediumInteraction(int depth, bool caustic,
		const MediumSamplingRecord &mRec, Float time, const Vector &wi, 
		const Spectrum &weight);

	MTS_DECLARE_CLASS()
protected:
	/// Protected constructor
	inline ParticleTracer(int maxDepth, bool multipleScattering, int rrDepth) 
		: m_maxDepth(maxDepth), m_multipleScattering(multipleScattering), 
		m_rrDepth(rrDepth) {
	}
	/// Protected constructor
	ParticleTracer(Stream *stream, InstanceManager *manager);
	/// Virtual destructor
	virtual ~ParticleTracer() { }
protected:
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	int m_maxDepth;
	bool m_multipleScattering;
	int m_rrDepth;
};

MTS_NAMESPACE_END

#endif /* __PARTICLEPROC_H */
