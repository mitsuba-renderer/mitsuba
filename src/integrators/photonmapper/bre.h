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

#if !defined(__BEAM_RADIANCE_ESTIMATOR_H)
#define __BEAM_RADIANCE_ESTIMATOR_H

#include <mitsuba/render/photonmap.h>

MTS_NAMESPACE_BEGIN

/**
 * Implements the beam radiance estimate described in
 * "The Beam Radiance Estimator for Volumetric Photon Mapping"
 * by Wojciech Jarosz, Matthias Zwicker, and Henrik Wann Jensen.
 */

class BeamRadianceEstimator : public SerializableObject {
public:
	/**
	 * \brief Create a BRE acceleration data structure from
	 * an existing volumetric photon map
	 */
	BeamRadianceEstimator(const PhotonMap *pmap, size_t lookupSize);

	/**
	 * \brief Unserialize a BRE acceleration data structure from
	 * a binary data stream
	 */
	BeamRadianceEstimator(Stream *stream, InstanceManager *manager);

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Compute the beam radiance estimate for the given ray segment and medium
	Spectrum query(const Ray &ray, const Medium *medium) const;

	MTS_DECLARE_CLASS()
protected:
	/// Release all memory
	virtual ~BeamRadianceEstimator();

	/// Fit a hierarchy of bounding boxes to the stored photons
	AABB buildHierarchy(size_t index);

	/// Heap convenience routines
	inline uint32_t leftChild(uint32_t index) const { return 2*index; }
	inline uint32_t rightChild(uint32_t index) const { return 2*index + 1; }
	inline bool isInnerNode(uint32_t index) const { return index <= m_lastInnerNode; }
	inline bool hasRightChild(uint32_t index) const { return index <= m_lastRChildNode; }
private:
	struct BRENode {
		AABB aabb;
		Photon photon;
		Float radius;
	};

	BRENode *m_nodes;
	size_t m_photonCount;
	size_t m_lastInnerNode;
	size_t m_lastRChildNode;
	int m_depth;
	Float m_scaleFactor;
};

MTS_NAMESPACE_END

#endif /* __BEAM_RADIANCE_ESTIMATOR_H */
