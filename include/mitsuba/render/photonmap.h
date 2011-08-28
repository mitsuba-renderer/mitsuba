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

#if !defined(__PHOTONMAP_H)
#define __PHOTONMAP_H

#include <mitsuba/render/photon.h>

MTS_NAMESPACE_BEGIN

/** \brief Implementation of the photon map data structure
 *
 * Based on Henrik Wann Jensen's book "Realistic Image Synthesis 
 * Using Photon Mapping". 
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER PhotonMap : public SerializableObject {
public:
	typedef PointKDTree<Photon> PhotonTree;

    /* ===================================================================== */
    /*                        Public access methods                          */
    /* ===================================================================== */

	/**
	 * \brief Create an empty photon map and reserve memory 
	 * for a specified number of photons.
	 */
	PhotonMap(size_t photonCount = 0);

	/**
	 * \brief Unserialize a photon map from a binary data stream
	 */
	PhotonMap(Stream *stream, InstanceManager *manager);

	// =============================================================
	//! @{ \name \c stl::vector-like interface
	// =============================================================
	/// Clear the kd-tree array
	inline void clear() { m_kdtree.clear(); }
	/// Resize the kd-tree array
	inline void resize(size_t size) { m_kdtree.resize(size); }
	/// Reserve a certain amount of memory for the kd-tree array
	inline void reserve(size_t size) { m_kdtree.reserve(size); }
	/// Return the size of the kd-tree
	inline size_t size() const { return m_kdtree.size(); }
	/// Return the capacity of the kd-tree
	inline size_t capacity() const { return m_kdtree.capacity(); }
	/// Append a kd-tree photon to the photon array
	inline void push_back(const Photon &photon) { m_kdtree.push_back(photon); }
	/// Return one of the photons by index
	inline Photon &operator[](size_t idx) { return m_kdtree[idx]; }
	/// Return one of the photons by index (const version)
	inline const Photon &operator[](size_t idx) const { return m_kdtree[idx]; }
	//! @}
	// =============================================================

	/**
	 * \brief Try to append a photon to the photon map
	 *
	 * \return \c false If the photon map is full
	 */
	inline bool tryAppend(const Photon &photon) {
		if (size() < capacity()) {
			push_back(photon);
			return true;
		} else {
			return false;
		}
	}

	/// Scale all photon power values contained in this photon map
	inline void setScaleFactor(Float value) { m_scale = value; }

	/// Return the power scale factor of this photon map
	inline Float getScaleFactor() const { return m_scale; }

	/**
	 * \brief Build a photon map over the supplied photons. 
	 *
	 * This has to be done once after all photons have been stored, 
	 * but prior to executing any queries.
	 */
	inline void build() { m_kdtree.build(); }

	/// Determine if the photon map is completely filled
	inline bool isFull() const {
		return capacity() == size();
	}

	/**
	 * Serialize a photon map to a binary data stream
	 */
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Dump the photons to an OBJ file to analyze their spatial distribution
	void dumpOBJ(const std::string &filename);
	
	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~PhotonMap();
protected:
	PhotonTree m_kdtree;
	Float m_scale;
};

MTS_NAMESPACE_END

#endif /* __PHOTONMAP_H */
