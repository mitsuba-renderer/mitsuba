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

#if !defined(__PHOTON_H)
#define __PHOTON_H

#include <mitsuba/core/serialization.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/** \brief Memory-efficient photon representation
 *
 * Requires 24 bytes when Mitsuba is compiled with single precision
 * and RGB-based color spectra.
 */
struct MTS_EXPORT_RENDER Photon {
	friend class PhotonMap;
public:
	// ======================================================================
	/// @{ \name Photon attributes
	// ======================================================================

	float pos[3];			//!< Photon position in single precision
#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
	Spectrum power;         //!< Accurate spectral photon power representation
#else
	uint8_t power[4];		//!< Photon power stored in Greg Ward's RGBE format
#endif
	uint8_t theta;			//!< Discretized photon direction (\a theta component)
	uint8_t phi;			//!< Discretized photon direction (\a phi component)
	uint8_t thetaN;			//!< Discretized surface normal (\a theta component)
	uint8_t phiN;			//!< Discretized surface normal (\a phi component)
	uint16_t depth;			//!< Photon depth (number of preceding interactions)
	uint8_t axis;			//!< Split axis in the associated KD-tree
	uint8_t unused;			//!< Unused 8-bit field (needed for alignment)

	/// @}
	// ======================================================================

	/// Dummy constructor
	inline Photon() { }

	/// Construct from a photon interaction 
	Photon(const Point &pos, const Normal &normal,
			const Vector &dir, const Spectrum &power,
			uint16_t depth);
	
	/// Unserialize from a binary data stream
	Photon(Stream *stream);

	/// @}
	// ======================================================================

	/// Return the depth (in # of interactions)
	inline int getDepth() const {
		return depth;
	}

	/// Compute the squared distance between this photon and some point.
	inline float distSquared(const float *q) const {
		float dist1 = pos[0]-q[0], dist2 = pos[1]-q[1],
		      dist3 = pos[2]-q[2];
		return dist1*dist1 + dist2*dist2 + dist3*dist3;
	}

	/**
	 * Convert the photon direction from quantized spherical coordinates
	 * to a floating point vector value. Precomputation idea based on 
	 * Jensen's implementation.
	 */
	inline Vector getDirection() const {
		return Vector(
			m_cosPhi[phi] * m_sinTheta[theta],
			m_sinPhi[phi] * m_sinTheta[theta],
			m_cosTheta[theta]
		);
	}

	/**
	 * Convert the normal direction from quantized spherical coordinates
	 * to a floating point vector value.
	 */
	inline Normal getNormal() const {
		return Normal(
			m_cosPhi[phiN] * m_sinTheta[thetaN],
			m_sinPhi[phiN] * m_sinTheta[thetaN],
			m_cosTheta[thetaN]
		);
	}

	/// Return the photon position as a vector
	inline Point getPosition() const {
		return Point(pos[0], pos[1], pos[2]);
	}

	/// Convert the photon power from RGBE to floating point
	inline Spectrum getPower() const {
#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
		return power;
#else
		Spectrum result;
		result.fromRGBE(power);
		return result;
#endif
	}

	/// Serialize to a binary data stream
	inline void serialize(Stream *stream) const {
		stream->writeSingleArray(pos, 3);
		#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
			power.serialize(stream);
			stream->writeUChar(phi);
			stream->writeUChar(theta);
			stream->writeUChar(phiN);
			stream->writeUChar(thetaN);
		#else
			stream->write(power, 8);
		#endif
		stream->writeUShort(depth);
		stream->writeUChar(axis);
	}

	/// Return a string representation (for debugging)
	std::string toString() const {
		std::ostringstream oss;
		oss << "Photon[pos = [" << pos[0] << ", "
			<< pos[1] << ", " << pos[2] << "]"
			<< ", power = " << getPower().toString()
			<< ", direction = " << getDirection().toString()
			<< ", normal = " << getNormal().toString()
			<< ", axis = " << axis
			<< ", depth = " << depth
			<< "]";
		return oss.str();
	}
	
protected:
	// ======================================================================
	/// @{ \name Precomputed lookup tables
	// ======================================================================

	static Float m_cosTheta[256];
	static Float m_sinTheta[256];
	static Float m_cosPhi[256];
	static Float m_sinPhi[256];
	static Float m_expTable[256];
	static bool m_precompTableReady;

	/// @}
	// ======================================================================

	/// Initialize the precomputed lookup tables
	static bool createPrecompTables();
};

#if defined(SINGLE_PRECISION) && SPECTRUM_SAMPLES == 3
/* Compiler sanity check */
BOOST_STATIC_ASSERT(sizeof(Photon) == 24);
#endif

MTS_NAMESPACE_END

#endif /* __PHOTON_H */
