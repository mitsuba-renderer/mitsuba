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

#if !defined(__PARTICLEPROC_H)
#define __PARTICLEPROC_H

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract parallel particle tracing process
 *
 * This class implements a particle tracer similar to what is
 * described in appendix 4.A of Eric Veach's PhD thesis. Particles
 * are emitted from the light source and subsequently perform a random
 * walk that includes both surface and medium scattering events. The 
 * work is spread out over multiple cores/machines. For every such
 * event,a custom routine is invoked.
 *
 * To actually use this class, you must extend it and implement the
 * function \ref createWorkProcessor(), which should return a subclass
 * of \ref ParticleTracer with overridden functions
 * \ref ParticleTracer::handleSurfaceInteraction and
 * \ref ParticleTracer::handleMediumInteraction.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER ParticleProcess : public ParallelProcess {
public:
	/// The particle tracer supports two principal modes of operation
	enum EMode {
		/**
		 * \brief Trace a fixed number of \a particles
		 * 
		 * In this mode, a specified number of particles will be emitted, and
		 * a customizable action is performed for every scattering event.
		 * Note that the number of resulting \a events will generally be 
		 * different from the number of traced \a particles.
		 *
		 * This mode is used for instance by the \c ptracer plugin.
		 */
		ETrace = 0,

		/**
		 * \brief `Gather' a fixed number of scattering \a events
		 *
		 * In this mode, the number of particles to be emitted is
		 * unknown ahead of time. Instead, the implementation traces
		 * particles until a a certain number of scattering events
		 * have been recorded.
		 * 
		 * This mode is used to create photon maps. See
		 * \ref GatherPhotonProcess for an implementation.
		 */
		EGather
	};

	// =============================================================
	//! @{ \name Implementation of the ParallelProcess interface
	// =============================================================

	virtual EStatus generateWork(WorkUnit *unit, int worker);
	
	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:

	/**
	 * Create a new particle process 
	 *
	 * \param mode
	 *    Particle tracing mode - see above
	 * \param workCount
	 *    Total # of particles to trace / # events to record
	 * \param granularity
	 *    Number of particles in each work unit. When set to zero,
	 *    a suitable number will be automatically chosen.
	 * \param progressText
	 *    Title of the progress bar
	 * \param progressReporterPayload
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
 * \brief Abstract particle tracer implementation
 *
 * Traces particles and performs a customizable action every time a
 * surface or volume interaction occurs.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER ParticleTracer : public WorkProcessor {
public:
	// =============================================================
	//! @{ \name Implementation of the WorkProcessor interface
	// =============================================================

	virtual ref<WorkUnit> createWorkUnit() const;
	virtual void prepare();
	virtual void process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop);
	void serialize(Stream *stream, InstanceManager *manager) const;

	//! @}
	// =============================================================
	
	/**
	 * \brief Handle a particle emission event
	 *
	 * To be overridden in a subclass. The default implementation
	 * does nothing
	 *
	 * \param eRec
	 *    Emission sampling record
	 * \param medium
	 *    Pointer to the current medium
	 * \param time
	 *    Time value associated with the particle
	 */
	virtual void handleEmission(const EmissionRecord &eRec,
			const Medium *medium, Float time);

	/**
     * \brief Handle a surface interaction event
	 *
	 * To be overridden in a subclass. The default implementation
	 * does nothing
	 *
	 * \param depth 
	 *    Depth of the interaction in path space (with 1
	 *    corresponding to the first bounce) 
	 * \param caustic
	 *    Is this a caustic path? (This flag is \c false when there
	 *    was any non-specular interaction on the path so far)
	 * \param its
	 *    Associated intersection record
	 * \param medium
	 *    Pointer to the current medium \a before the surface
	 *    interaction
	 * \param weight
	 *    Particle weight along the preceding subpath
     */
	virtual void handleSurfaceInteraction(int depth, bool caustic,
		const Intersection &its, const Medium *medium,
		const Spectrum &weight);

    /**
     * \brief Handle a medium interaction event
	 *
	 * To be overridden in a subclass. The default implementation
	 * does nothing
	 *
	 * \param depth 
	 *    Depth of the interaction in path space (with 1
	 *    corresponding to the first bounce) 
	 * \param caustic
	 *    Is this a caustic path? (This flag is \c false when there
	 *    was any non-specular interaction on the path so far)
	 * \param mRec
	 *    Associated medium sampling record
	 * \param medium
	 *    Pointer to the current medium
	 * \param time
	 *    Time value associated with the particle
	 * \param weight
	 *    Particle weight along the preceding subpath
     */
	virtual void handleMediumInteraction(int depth, bool caustic,
		const MediumSamplingRecord &mRec, const Medium *medium,
		Float time, const Vector &wi, const Spectrum &weight);

	MTS_DECLARE_CLASS()
protected:
	/// Protected constructor
	inline ParticleTracer(int maxDepth, int rrDepth) 
		: m_maxDepth(maxDepth), m_rrDepth(rrDepth) { }
	/// Protected constructor
	ParticleTracer(Stream *stream, InstanceManager *manager);
	/// Virtual destructor
	virtual ~ParticleTracer() { }
protected:
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	int m_maxDepth;
	int m_rrDepth;
};

MTS_NAMESPACE_END

#endif /* __PARTICLEPROC_H */
