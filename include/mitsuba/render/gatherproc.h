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

#pragma once
#if !defined(__MITSUBA_RENDER_GATHERPROC_H_)
#define __MITSUBA_RENDER_GATHERPROC_H_

#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/photonmap.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Process for parallel photon map construction
 *
 * Given a number and type (surface/caustic/volume) of photons, this
 * class distributes the work over an arbitrary number of machines.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER GatherPhotonProcess : public ParticleProcess {
public:
	enum EGatherType {
		/// Surface photons (indirect on diffuse surfaces, last bounce was not through a delta BSDF)
		ESurfacePhotons,
		/// Caustic photons (indirect on diffuse surfaces, last bounce was through a delta BSDF)
		ECausticPhotons,
		/// Surface photons (all of them, even direct illumination)
		EAllSurfacePhotons,
		/// Volumetric photons
		EVolumePhotons
	};

	/**
	 * Create a new process for parallel photon gathering
	 * \param type
	 *     Specifies the type of requested photons (surface/caustic/volume)
	 * \param photonCount
	 *     Specifies the number of requested photons
	 * \param granularity
	 *     Size of the internally used work units (in photons)
	 * \param isLocal
	 *     Should the parallel process only be executed locally? (sending
	 *     photons over the network may be unnecessary and wasteful)
	 * \param autoCancel
	 *     Indicates if the gathering process should be canceled if there
	 *     are not enough photons generated
	 * \param progressReporterPayload
	 *    Custom pointer payload to be delivered with progress messages
	 */
	GatherPhotonProcess(EGatherType type, size_t photonCount,
		size_t granularity, int maxDepth, int rrDepth, bool isLocal,
		bool autoCancel, const void *progressReporterPayload);

	/**
	 * Once the process has finished, this returns a reference
	 * to the (still unbalanced) photon map
	 */
	inline PhotonMap *getPhotonMap() { return m_photonMap; }

	/**
	 * \brief Return the number of discarded photons
	 *
	 * Due to asynchronous processing, some excess photons
	 * will generally be produced. This function returns the number
	 * of excess photons that had to be discarded. If this is too
	 * high, the granularity should be decreased.
	 */
	inline size_t getExcessPhotons() const { return m_excess; }

	/**
	 * \brief Lists the nuber of particles that had to be shot
	 * in order to fill the photon map.
	 */
	inline size_t getShotParticles() const { return m_numShot; }

	// ======================================================================
	/// @{ \name ParallelProcess implementation
	// ======================================================================

	bool isLocal() const;
	ref<WorkProcessor> createWorkProcessor() const;
	void processResult(const WorkResult *wr, bool cancelled);
	EStatus generateWork(WorkUnit *unit, int worker);

	/// @}
	// ======================================================================

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GatherPhotonProcess() { }

	/**
	 * \brief Checks if the configuration of needed, generated and shot
	 * photons indicates an unsuccessful progress of the gathering. This
	 * check is taken from PBRT.
	 */
	inline bool unsuccessful(size_t needed, size_t gen, size_t shot) {
		return (gen < needed && (gen == 0 || gen < shot/1024));
	}
protected:
	EGatherType m_type;
	ref<PhotonMap> m_photonMap;
	size_t m_photonCount;
	int m_maxDepth;
	int m_rrDepth;
	bool m_isLocal;
	bool m_autoCancel;
	size_t m_excess, m_numShot;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_GATHERPROC_H_ */
