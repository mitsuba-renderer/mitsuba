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

#if !defined(__PHOTONMAP_H)
#define __PHOTONMAP_H

#include <mitsuba/render/photon.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Hard-coded size limit - permits to use a recursion-free
 * nearest neighbor search.
 */
#define MAX_PHOTONMAP_DEPTH 30 // corresponds to 1,073,741,824 photons

/** \brief Generic photon map implementation based on Henrik Wann Jensen's book
 * "Realistic Image Synthesis Using Photon Mapping". Uses 64-bit addressing
 * for really really large photon maps. STL-ified implementation with some 
 * comments and explanations of the individial algorithms.
 *
 * \remark The reference implementation in the book suffers from an indexing 
 *  bug, which results in a small amount of photons being lost in the final 
 *  estimate. The bug has been fixed in this re-implementation.
 *
 * \author Wenzel Jakob
 */
class MTS_EXPORT_RENDER PhotonMap : public SerializableObject {
public:
    /* ===================================================================== */
    /*                        Public access methods                          */
    /* ===================================================================== */

	/**
	 * Create an empty photon map and reserve memory for a specified
	 * number of photons.
	 */
	PhotonMap(size_t maxPhotons);

	/**
	 * Unserialize a photon map from a binary data stream
	 */
	PhotonMap(Stream *stream, InstanceManager *manager);

	/**
	 * Store a photon in the (still unbalanced) photon map
	 */
	bool storePhoton(const Point &pos, const Normal &normal,
					 const Vector &dir, const Spectrum &power,
					 uint16_t depth);

	/**
	 * Store a photon in the (still unbalanced) photon map
	 */
	bool storePhoton(const Photon &photon);

	/**
	 * Scale all photon power values contained in this photon map
	 */
	void setScale(Float value);

	/**
	 * Recursively build a left-balanced kd-tree. This has to be
	 * done once after all photons have been stored, but prior to
	 * executing any queries.
	 */
	void balance();

	/**
	 * Using the photon map, estimate the irradiance on a surface (Unfiltered)
	 *
	 * \param p
	 * 		The surface position for the estimate
	 * \param n
	 * 		Normal vector of the surface in question
	 * \param searchRadius
	 * 		Size of the spherical photon search region
	 * \param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 */
	Spectrum estimateIrradiance(
		const Point &p, const Normal &n, 
		Float searchRadius,
		size_t maxPhotons) const;

	/**
	 * Using the photon map, estimate the irradiance on a surface (filtered
	 * using Simpson's kernel)
	 *
	 * \param p
	 * 		The surface position for the estimate
	 * \param n
	 * 		Normal vector of the surface in question
	 * \param searchRadius
	 * 		Size of the spherical photon search region
	 * \param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 */
	Spectrum estimateIrradianceFiltered(
		const Point &p, const Normal &n, 
		Float searchRadius,
		size_t maxPhotons) const;

	/**
	 * Using the photon map and a surface intersection, estimate the
	 * radiant exitance into a given direction
	 *
	 * \param its
	 * 		A surface interaction
	 * \param searchRadius
	 * 		Size of the spherical photon search region
	 * \param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 */
	Spectrum estimateRadianceFiltered(
		const Intersection &its,
		Float searchRadius,
		size_t maxPhotons) const;

	/**
	 * Compute scattered contributions from all photons within
	 * the specified radius. Does no weighting/filtering/dynamic
	 * search radius reduction and simply sums over all photons. Only
	 * considers photons with a depth value less than or equal
	 * to the 'maxDepth' parameter. This function is meant to be 
	 * used with progressive photon mapping. 
	 */
	size_t estimateRadianceRaw(const Intersection &its,
		Float searchRadius, Spectrum &result, int maxDepth) const;

	/**
	 * Using the photon map and an outgoing ray, estimate the
	 * scattered volume radiance into a given direction
	 *
	 * \param ray 
	 * 		The outgoing ray
	 * \param searchRadius
	 * 		Size of the spherical photon search region
	 * \param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 * \param medium
	 *      The associated participating medium
	 */
	Spectrum estimateVolumeRadiance(
		const MediumSamplingRecord &mRec,
		const Ray &ray,
		Float searchRadius,
		size_t maxPhotons,
		const Medium *medium) const;

	/// Determine if the photon map is completely filled
	inline bool isFull() const {
		return m_photonCount >= m_maxPhotons;
	}

	/// Return the currently stored amount of photons
	inline size_t getPhotonCount() const {
		return m_photonCount;
	}

	/// Has the photon map been balanced?
	inline bool isBalanced() const {
		return m_balanced;
	}

	/// Return the position of a photon in the photon map (for debugging)
	inline Point getPhotonPosition(size_t i) const {
		return m_photons[i].getPosition();
	}

	/// Set the minimum amount of photons to consider an estimate valid
	void setMinPhotons(int minPhotons);

	/**
	 * Serialize a photon map to a binary data stream
	 */
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Dump the photons to an OBJ file to analyze their spatial distribution
	void dumpOBJ(const std::string &filename);
	
	/// Initialize the static lookup tables (executed once at program startup)
	static bool createPrecompTables();

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
    /* ===================================================================== */
    /*                      Protected data structures                        */
    /* ===================================================================== */

	typedef Photon *photon_ptr;
	typedef const Photon *const_photon_ptr;

	typedef std::vector<photon_ptr>::iterator photon_iterator;
	typedef std::pair<float, const_photon_ptr> search_result;

	/// \cond
	/* Photon position comparison functor (< comparison). */
	struct comparePhotonLess : public std::unary_function<photon_ptr, bool> {
	public:
		inline comparePhotonLess(int splitAxis, float splitValue) 
			: axis(splitAxis), value(splitValue) {
		}

		inline bool operator()(photon_ptr photon) const {
			return photon->pos[axis] < value;
		}
	private:
		int axis;
		float value;
	};

	/* Photon position comparison functor (> comparison). */
	struct comparePhotonGreater : public 
		std::unary_function<photon_ptr, bool> {
	public:
		inline comparePhotonGreater(int splitAxis, float splitValue) 
			: axis(splitAxis), value(splitValue) {
		}

		inline bool operator()(photon_ptr photon) const {
			return photon->pos[axis] > value;
		}
	private:
		int axis;
		float value;
	};

	/// Distance metric used to build MAX-heaps for nearest-neighbor searches
	struct distanceMetric : public 
		std::binary_function<search_result, search_result, bool> {
	public:
		inline bool operator()(search_result &a, search_result &b) const {
			return a.first < b.first;
		}
	};

	/// \endcond
protected:
    /* ===================================================================== */
    /*                         Protected subroutines                         */
    /* ===================================================================== */

	/// Virtual destructor
	virtual ~PhotonMap();

	/**
	 * Perform a new nearest-neighbor search query
	 *
	 * \param pos
	 * 		Nearest-neighbor search position
	 * \param searchRadiusSquared
	 * 		Squared search radius - is updated should the search
	 *      radius be decreased
	 * \param maxSize
	 * 		Maximum number of photons, which will be returned. If
	 * 		this value is ever exceeded, the search uses a priority
	 * 		queue to only return the closest entries
	 * \param results
	 *      Pre-allocated search result data structure. Should have
	 *      one extra entry for internal use. 
	 * \return
	 *      The number of results
	 */
	size_t nnSearch(const Point &p, Float &searchRadiusSquared, 
		size_t maxSize, search_result *results) const;

	/**
	 * Partition a list of photons so that there is an ordering between
	 * the photon at index 'pivot' and all other photons. See implementation
	 * for further documentation.
	 */
	void quickPartition(photon_iterator left, 
		photon_iterator right, photon_iterator pivot, int axis) const;

	/**
 	 * Given a number of entries, this method calculates the maximum amount of
	 * nodes on the left subtree of a left-balanced tree. See implementation
	 * for further documentation
	 */
	size_t leftSubtreeSize(size_t treeSize) const;

	/**
	 * Recursive balancing algorithm. See implementation for
	 * further documentation.
	 *
	 * \param sortStart
	 * 		The beginning of the list to be balanced
	 * \param sortEnd
	 * 		The end of the list to be balanced
	 * \param heapPointers
	 * 		An array, which will hold the photon tree as it is constructed
	 * \param heapIndex
	 * 		The heap array index of the next pivot
	 */
	void balanceRecursive(
		photon_iterator basePtr,
		photon_iterator sortStart,
		photon_iterator sortEnd,
		std::vector<size_t> &heapPermutation,
		AABB &aabb, size_t heapIndex) const;

	/// Heap access routines
	inline size_t leftChild(size_t index) const { return 2*index; }
	inline size_t rightChild(size_t index) const { return 2*index + 1; }
	inline bool isInnerNode(size_t index) const { return index <= m_lastInnerNode; }
	inline bool hasRightChild(size_t index) const { return index <= m_lastRChildNode; }
private:
    /* ===================================================================== */
    /*                        Protected attributes                           */
    /* ===================================================================== */
	Photon *m_photons;
	AABB m_aabb;
	size_t m_photonCount;
	size_t m_maxPhotons;
	size_t m_minPhotons;
	size_t m_lastInnerNode;
	size_t m_lastRChildNode;
	bool m_balanced;
	Float m_scale;
};

MTS_NAMESPACE_END

#endif /* __PHOTONMAP_H */
