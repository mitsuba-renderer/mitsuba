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
#if !defined(__MITSUBA_RENDER_IRRCACHE_H_)
#define __MITSUBA_RENDER_IRRCACHE_H_

#include <mitsuba/render/scene.h>
#include <mitsuba/core/octree.h>

MTS_NAMESPACE_BEGIN

/* 3 (X, Y, and Z) components for each spectral sample */
typedef Spectrum RotationalGradient[3];
typedef Spectrum TranslationalGradient[3];

/**
 * \brief Utility data structure for hemispherical sampling and
 * translational/rotational gradient computation.
 *
 * Uses the improved translational gradients proposed in
 * the paper "Improved radiance gradient computation" by
 * Krivanek J., Gautron P., Bouatouch K., Pattanaik S.
 * (Proceedings of SCCG 2005)
 *
 * \author Wenzel Jakob
 * \ingroup librender
 */
class MTS_EXPORT_RENDER HemisphereSampler : public Object {
public:
	struct SampleEntry {
		Vector d;
		Spectrum L;
		Float dist;
		Float cosTheta;
		Float sinTheta;
	};

	/**
	 * Allocate storage for a hemispherical sample with M
	 * elevation and N azimuthal samples.
	 */
	HemisphereSampler(uint32_t M, uint32_t N);

	/// Return the elevational resolution
	inline uint32_t getM() const { return m_M; }

	/// Return the azimuthal resolution
	inline uint32_t getN() const { return m_N; }

	/// Access a cell by index
	inline SampleEntry &operator() (uint32_t j, uint32_t k) {
		return m_entries[j*m_N + k];
	}

	/// Access a cell by index (const version)
	inline const SampleEntry &operator() (uint32_t j, uint32_t k) const {
		return m_entries[j*m_N + k];
	}

	/// Generate a set of projected solid angle-distributed directions
	void generateDirections(const Intersection &its);

	/// Compute mean distance, gradients etc.
	void process(const Intersection &its);

	/// Return the irradiance gradient with respect to rotation
	inline const RotationalGradient &getRotationalGradient() const {
		return m_rGrad;
	}

	/// Return the irradiance gradient with respect to translation
	inline const RotationalGradient &getTranslationalGradient() const {
		return m_tGrad;
	}

	/// Return the average distance over all cells (harmonic mean)
	inline Float getHarmonicMeanDistance() const {
		return m_hMean;
	}

	/// Return the minimum distance over all cells
	inline Float getMinimumDistance() const {
		return m_hMin;
	}

	/// Return the minimum distance over all cells (>10 deg in elevation)
	inline Float getMinimumDistanceRestricted() const {
		return m_hMinRestricted;
	}

	/// Return the computed irradiance
	inline const Spectrum &getIrradiance() const {
		return m_E;
	}

	MTS_DECLARE_CLASS()
protected:
	/// Free all memory
	virtual ~HemisphereSampler();
private:
	uint32_t m_M, m_N;
	SampleEntry *m_entries;
	Vector *m_uk, *m_vk, *m_vkMinus;
	Spectrum m_E;
	RotationalGradient m_rGrad;
	TranslationalGradient m_tGrad;
	Float m_hMean, m_hMin, m_hMinRestricted;
	ref<Random> m_random;
};


/** \brief Irradiance cache data structure based on "A Ray Tracing Solution
 * for Diffuse Interreflection" by Greg J. Ward, Francis M. Rubinstein and
 * Robert D. Clear (Computer Graphics, Volume 22, Number 4, August 1988)
 *
 * with extensions from
 * "Irradiance Gradients" by Greg J. Ward and Paul S. Heckbert
 * (1992 Eurographics Workshop on Rendering).
 *
 * Includes optimizations proposed by Jaroslav Krivanek et al. in
 * "Making Radiance and Irradiance Caching Practical: Adaptive Caching and
 * Neighbor Clamping" (in EGSR 2006),
 *
 * and
 *
 * "An Approximate Global Illumination System for Computer Generated Films"
 * by E. Tabellion and A. Lamorlette (SIGGRAPH 2004)
 *
 * \author Wenzel Jakob
 * \ingroup librender
 */
class MTS_EXPORT_RENDER IrradianceCache : public SerializableObject {
public:
	struct Record;

    /* ===================================================================== */
    /*                        Public access methods                          */
    /* ===================================================================== */

	/**
	 * Create an empty irradiance of the given size
	 */
	IrradianceCache(const AABB &aabb);

	/**
	 * Unserialize an irradiance cache from a binary data stream
	 */
	IrradianceCache(Stream *stream, InstanceManager *manager);

	/**
	 * Set the quality parameter \kappa from the
	 * Tabellion and Lamorlette paper. Once samples
	 * have been stored, this should only be decreased.
	 */
	inline void setQuality(Float quality) { m_kappa = quality; }

	/**
	 * Set the influence region cutoff values of samples in the
	 * cache. Will be multiplied by the scene size
	 */
	void clampInfluence(Float min, Float max);

	/**
	 * Minimal influence region falloff with increasing distance
	 */
	inline void clampScreen(bool active) { m_clampScreen = active; }

	/**
	 * Enable neighbor clamping?
	 */
	inline void clampNeighbor(bool active) { m_clampNeighbor = active; }

	/**
	 * Enable/disable irradiance gradients
	 */
	inline void useGradients(bool active) { m_useGradients = active; }

	/**
	 * Add a sample to the irradiance cache
	 *
	 * \param ray
	 * 		Ray differentials (if they exist)
	 * \param its
	 * 		The position/normal of the surface in question
	 * \param sample
	 *      Record containing all hemispherical samples and
	 *      derived gradient information
	 */
	Record *put(const RayDifferential &ray, const Intersection &its,
		const HemisphereSampler &hs);

	/**
	 * Use the irradiance cache to interpolate/extrapolate an
	 * irradiance value for the given position and surface normal
	 * Returns false on a cache miss
	 */
	bool get(const Intersection &its, Spectrum &E) const;

	/// Manually insert an irradiance record
	void insert(Record *rec);

	/**
	 * Serialize an irradiance cache to a binary data stream
	 */
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return a string representation
	std::string toString() const;

    /* ===================================================================== */
    /*                       Internal data structures                        */
    /* ===================================================================== */

	struct Record {
		/* Sample position */
		Point p;
		/* Normal vector of the associated surface */
		Normal n;
		/* Minimum intersection distance */
		Float R0;
		/* Unclamped distance - must be stored for
		   neighbor clamping */
		Float originalR0;
		/* Minimum/Maximum distance - must be
		   stored for neighbor clamping. */
		Float R0_min, R0_max;
		/* Irradiance value */
		Spectrum E;
		/* Rotational gradient for improved interpolation */
		RotationalGradient rGrad;
		/* Translational gradient for improved interpolation */
		TranslationalGradient tGrad;

		/// Dummy constructor
		inline Record() { }

		/// Copy constructor
		inline Record(const Record *rec)
		  : p(rec->p), n(rec->n), R0(rec->R0), originalR0(rec->originalR0),
		    R0_min(rec->R0_min), R0_max(rec->R0_max), E(rec->E) {
			for (int i=0; i<3; ++i) {
				tGrad[i] = rec->tGrad[i];
				rGrad[i] = rec->rGrad[i];
			}
		}

		/// Unserialize from a binary data stream
		inline Record(Stream *stream) {
			p = Point(stream);
			n = Normal(stream);
			R0 = stream->readFloat();
			originalR0 = stream->readFloat();
			R0_min = stream->readFloat();
			R0_max = stream->readFloat();
			E = Spectrum(stream);
			for (int i=0; i<3; ++i)
				rGrad[i] = Spectrum(stream);
			for (int i=0; i<3; ++i)
				tGrad[i] = Spectrum(stream);
		}

		/// Serialize to a binary data stream
		inline void serialize(Stream *stream) const {
			p.serialize(stream);
			n.serialize(stream);
			stream->writeFloat(R0);
			stream->writeFloat(originalR0);
			stream->writeFloat(R0_min);
			stream->writeFloat(R0_max);
			E.serialize(stream);
			for (int i=0; i<3; ++i)
				rGrad[i].serialize(stream);
			for (int i=0; i<3; ++i)
				tGrad[i].serialize(stream);
		}

		/**
		 * Calculate contribution of this sample if it were used
		 * to calculate an interpolated value at the given position
		 */
		inline Float getWeight(const Point &p2, const Normal &n2, Float kappa) const {
			Float dp = dot(n, n2);

			/* Quickly discard opposite-facing samples */
			if (dp < 0.0f)
				return 0.0f;
			else if (dp > 1.0f)
				dp = 1.0f;

			/* Reject illuminance values 'in front' of P2 */
			if (dot(p2 - p, n + n2) < -0.05f)
				return 0.0f;

			/* Ad-hoc weight function (Tabellion & Lamorlette) */
			Float ePI = (p-p2).length() / (.5f * R0);
			Float eNI = std::sqrt(1.0f - std::abs(dp))
				/ 0.12326f;
			Float weight = 1 - kappa*std::max(ePI, eNI);

			if (weight < 0)
				return 0.0f;

			return weight;
		}
	};

	MTS_DECLARE_CLASS()
protected:
	/// Release all memory
	virtual ~IrradianceCache();
protected:
    /* ===================================================================== */
    /*                        Protected attributes                           */
    /* ===================================================================== */

	DynamicOctree<Record *> m_octree;
	std::vector<Record *> m_records;
	Float m_kappa;
	Float m_sceneSize;
	Float m_minDist, m_maxDist;
	bool m_clampScreen, m_clampNeighbor, m_useGradients;
	ref<Mutex> m_mutex;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_IRRCACHE_H_ */
