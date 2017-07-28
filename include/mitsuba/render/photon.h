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
#if !defined(__MITSUBA_RENDER_PHOTON_H_)
#define __MITSUBA_RENDER_PHOTON_H_

#include <mitsuba/core/serialization.h>
#include <mitsuba/core/kdtree.h>

/**
 * \brief Should Mitsuba use a left-balanced photon map?
 *
 * This saves some memory, but at a noticeable cost in query
 * performance. The default is to build an unbalanced
 * photon map using the sliding midpoint rule.
 */
#define MTS_PHOTONMAP_LEFT_BALANCED 0

MTS_NAMESPACE_BEGIN

/// Internal data record used by \ref Photon
struct PhotonData {
#if defined(SINGLE_PRECISION) && SPECTRUM_SAMPLES == 3
    uint8_t power[4];       //!< Photon power stored in Greg Ward's RGBE format
#else
    Spectrum power;         //!< Accurate spectral photon power representation
#endif
    uint8_t theta;          //!< Discretized photon direction (\a theta component)
    uint8_t phi;            //!< Discretized photon direction (\a phi component)
    uint8_t thetaN;         //!< Discretized surface normal (\a theta component)
    uint8_t phiN;           //!< Discretized surface normal (\a phi component)
    uint16_t depth;         //!< Photon depth (number of preceding interactions)
};

/** \brief Memory-efficient photon representation for use with
 * \ref PointKDTree
 *
 * \ingroup librender
 * \sa PhotonMap
 */
struct MTS_EXPORT_RENDER Photon :
#if MTS_PHOTONMAP_LEFT_BALANCED == 1
    public LeftBalancedKDNode<Point, PhotonData> {
#else
    public SimpleKDNode<Point, PhotonData> {
#endif
    friend class PhotonMap;
public:
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
        return data.depth;
    }

    /**
     * Convert the photon direction from quantized spherical coordinates
     * to a floating point vector value. Precomputation idea based on
     * Jensen's implementation.
     */
    inline Vector getDirection() const {
        return Vector(
            m_cosPhi[data.phi] * m_sinTheta[data.theta],
            m_sinPhi[data.phi] * m_sinTheta[data.theta],
            m_cosTheta[data.theta]
        );
    }

    /**
     * Convert the normal direction from quantized spherical coordinates
     * to a floating point vector value.
     */
    inline Normal getNormal() const {
        return Normal(
            m_cosPhi[data.phiN] * m_sinTheta[data.thetaN],
            m_sinPhi[data.phiN] * m_sinTheta[data.thetaN],
            m_cosTheta[data.thetaN]
        );
    }

    /// Convert the photon power from RGBE to floating point
    inline Spectrum getPower() const {
#if defined(SINGLE_PRECISION) && SPECTRUM_SAMPLES == 3
        Spectrum result;
        result.fromRGBE(data.power);
        return result;
#else
        return data.power;
#endif
    }

    /// Serialize to a binary data stream
    void serialize(Stream *stream) const;

    /// Return a string representation (for debugging)
    std::string toString() const;
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

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_PHOTON_H_ */
