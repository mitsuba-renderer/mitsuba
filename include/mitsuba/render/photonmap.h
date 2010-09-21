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

#include <mitsuba/core/serialization.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/**
 * Hard-coded size limit - permits to use a recursion-free
 * nearest neighbor search.
 */
#define MAX_PHOTONMAP_DEPTH 30 // corresponds to 1,073,741,824 photons

/** \brief Generic photon map implementation based on Henrik Wann Jensen's book
 * "Realistic Image Synthesis Using Photon Mapping". Uses 64-bit addressing
 * for really really large photon maps. STL-ified implementation with some 
 * comments and explanations of the individial algorithms.
 *
 * @remark The reference implementation in the book suffers from an indexing 
 *  bug, which results in a small amount of photons being lost in the final 
 *  estimate. The bug has been fixed in this re-implementation.
 *
 * @author Wenzel Jakob
 */
class MTS_EXPORT_RENDER PhotonMap : public SerializableObject {
public:
	struct Photon;

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
	 * If the photon map is to be filled by several processors, 
	 * it must be prepared using this method
	 */
	void prepareSMP(int numThreads);

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
	 * Store a photon in the (still unbalanced) photon map.
	 * SMP version - additionally requires a thread ID
	 */
	bool storePhotonSMP(int thread, const Point &pos, const Normal &normal,
		const Vector &dir, const Spectrum &power, uint16_t depth);

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
	 * @param p
	 * 		The surface position for the estimate
	 * @param n
	 * 		Normal vector of the surface in question
	 * @param searchRadius
	 * 		Size of the spherical photon search region
	 * @param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 */
	Spectrum estimateIrradiance(
		const Point &p, const Normal &n, 
		Float searchRadius,
		unsigned int maxPhotons) const;

	/**
	 * Using the photon map, estimate the irradiance on a surface (filtered
	 * using Simpson's kernel)
	 *
	 * @param p
	 * 		The surface position for the estimate
	 * @param n
	 * 		Normal vector of the surface in question
	 * @param searchRadius
	 * 		Size of the spherical photon search region
	 * @param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 */
	Spectrum estimateIrradianceFiltered(
		const Point &p, const Normal &n, 
		Float searchRadius,
		unsigned int maxPhotons) const;

	/**
	 * Using the photon map and a surface intersection, estimate the
	 * radiant exitance into a given direction
	 *
	 * @param its
	 * 		A surface interaction
	 * @param searchRadius
	 * 		Size of the spherical photon search region
	 * @param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 */
	Spectrum estimateRadianceFiltered(
		const Intersection &its,
		Float searchRadius,
		unsigned int maxPhotons) const;

	/**
	 * Compute scattered contributions from all photons within
	 * the specified radius. Does no weighting/filtering/dynamic
	 * search radius reduction and simply sums over all photons. Only
	 * considers photons with a depth value less than or equal
	 * to the 'maxDepth' parameter. This function is meant to be 
	 * used with progressive photon mapping. 
	 */
	int estimateRadianceRaw(const Intersection &its,
		Float searchRadius, Spectrum &result, int maxDepth) const;

	/**
	 * Using the photon map and an outgoing ray, estimate the
	 * scattered volume radiance into a given direction
	 *
	 * @param ray 
	 * 		The outgoing ray
	 * @param searchRadius
	 * 		Size of the spherical photon search region
	 * @param maxPhotons
	 * 		How many photon should (at most) be used in the estimate?
	 * @param medium
	 *      The associated participating medium
	 */
	Spectrum estimateVolumeRadiance(
		const MediumSamplingRecord &mRec,
		const Ray &ray,
		Float searchRadius,
		unsigned int maxPhotons,
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
	inline Point getPhotonPosition(int i) const {
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
public:
    /* ===================================================================== */
    /*                        Shared data structures                         */
    /* ===================================================================== */

	/** \brief Memory-efficient photon data structure as suggested in
		"Realistic Image Synthesis Using Photon Mapping" by Henrik Wann Jensen.
		The Photon Map will always use single-precision floating point
		(even if mitsuba is configured to use double precision).
	*/
	struct Photon { // 4*3 + 4 + 2 + 2 = 20 bytes
		float pos[3];			/* Photon position in single precision */
#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
		Spectrum power;         /* Super-accurate photon map :) */
#else
		uint8_t power[4];		/* Photon power in Ward's RGBE format */
#endif
		uint8_t phi, theta;		/* Approximate photon direction */
		uint8_t phiN, thetaN;	/* Approximate associated normal direction */
		uint16_t depth;			/* Photon depth (in # of interactions) */
		uint8_t axis;			/* Split axis */
		uint8_t unused;			/* Word alignment :( */

		/// Dummy constructor
		inline Photon() { }

		/// Construct from a photon interaction 
		Photon(const Point &pos, const Normal &normal,
			   const Vector &dir, const Spectrum &power,
			   uint16_t depth);

		/// Unserialize from a binary data stream
		inline Photon(Stream *stream) {
			stream->readSingleArray(pos, 3);
#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
			power = Spectrum(stream);
			phi = stream->readUChar();
			theta = stream->readUChar();
			phiN = stream->readUChar();
			thetaN = stream->readUChar();
#else
			stream->read(power, 8);
#endif
			depth = stream->readUShort();
			axis = stream->readUChar();
			unused = 0;
		}

		/// Return the depth (in # of interactions)
		inline int getDepth() const {
			return depth;
		}

		/// Compute the squared distance between this photon and some point.
		inline float distSquared(const float *q) const {
			float dist1 = pos[0]-q[0],
				  dist2 = pos[1]-q[1],
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
	};
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
			return a.first <= b.first;
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
	 * @param pos
	 * 		Nearest-neighbor search position
	 * @param searchRadiusSquared
	 * 		Squared search radius - is updated should the search
	 *      radius be decreased
	 * @param maxSize
	 * 		Maximum number of photons, which will be returned. If
	 * 		this value is ever exceeded, the search uses a priority
	 * 		queue to only return the closest entries
	 * @param results
	 *      Pre-allocated search result data structure. Should have
	 *      one extra entry for internal use. 
	 * @return
	 *      The number of results
	 */
	unsigned int nnSearch(const Point &p, Float &searchRadiusSquared, 
		unsigned int maxSize, search_result *results) const;

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
	 * @param sortStart
	 * 		The beginning of the list to be balanced
	 * @param sortEnd
	 * 		The end of the list to be balanced
	 * @param heapPointers
	 * 		An array, which will hold the photon tree as it is constructed
	 * 	@param heapIndex
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
	struct ThreadContext {
		AABB aabb;
		int photonOffset;
		int photonCount;
		int maxPhotons;
		uint8_t unused[128]; // Avoid false sharing
	};

	/* Precomputed lookup tables */
	static Float m_cosTheta[256];
	static Float m_sinTheta[256];
	static Float m_cosPhi[256];
	static Float m_sinPhi[256];
	static Float m_expTable[256];
	static bool m_precompTableReady;

	Photon *m_photons;
	AABB m_aabb;
	size_t m_photonCount;
	size_t m_maxPhotons;
	size_t m_minPhotons;
	size_t m_lastInnerNode;
	size_t m_lastRChildNode;
	int m_numThreads;
	bool m_balanced;
	Float m_scale;
	ThreadContext *m_context;
};

MTS_NAMESPACE_END

#endif /* __PHOTONMAP_H */
