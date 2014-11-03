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
#if !defined(__MITSUBA_BIDIR_EDGE_H_)
#define __MITSUBA_BIDIR_EDGE_H_

#include <mitsuba/bidir/common.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Bidirectional path edge data structure
 *
 * The path edge data structure is responsible for representing the transport of
 * light between pairs of scattering or emission events.
 * Amongst other things, it keeps track of the medium that fills the space between
 * adjacent vertices of a \ref Path. Furthermore, it can be used to evaluate and
 * sample the visibility and transmittance functions of the scene.
 *
 * Although they do not correspond to any real transport, this implementation
 * also places edges next to "supernode" vertices (see \ref PathVertex for a
 * description), which simplifies the implementations of various rendering
 * algorithms that make use of this framework.
 *
 * \sa PathVertex
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
struct MTS_EXPORT_BIDIR PathEdge {
	/* ==================================================================== */
	//! @{ \name             Enumerations and Fields
	/* ==================================================================== */

	/// Pointer to the medium that contains this edge (where \c NULL is vacuum)
	const Medium *medium;

	/**
	 * \brief Normalized direction vector associated with this edge
	 *
	 * The direction always points along the light path (from the light)
	 */
	Vector d;

	/**
	 * \brief Length of this edge in world-space distance units
	 *
	 * Note that edges adjacent to supernodes have length zero to
	 * mark them as such.
	 */
	Float length;

	/**
	 * \brief Measurement contribution weight
	 *
	 * This field stores the terms of the path-space measurement contribution
	 * function that are coupled to this specific edge divided by the
	 * associated density function.
	 *
	 * More specifically, it stores the transmittance of the medium across
	 * this edge divided by the density per unit length of the adjacent
	 * vertices int the radiance and importance transport directions (hence,
	 * it is an array with two entries).
	 *
	 * Note that this field only accounts for medium-related terms. The
	 * interactions with vertices are captured by \ref PathVertex::weight.
	 */
	Spectrum weight[ETransportModes];

	/**
	 * \brief Medium sampling density of the adjacent vertices
	 *
	 * This field stores the probability of sampling the preceding and
	 * successive path vertices using the sampling technique implemented by
	 * the function \ref PathEdge::sampleNext(). Depending on whether or not
	 * they are medium interactions, this eintries either store a density per
	 * unit length or a discrete probability.
	 */
	Float pdf[ETransportModes];

	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name              Sampling-related functions
	/* ==================================================================== */

	/**
	 * \brief Given a ray \c ray, sample a distance in this direction and
	 * fill the edge data structure, as well as its target vertex with content.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param pred
	 *     Pointer to the preceding vertex
	 * \param ray
	 *     Specifies the direction and origin associated with one
	 *     endpoint of the edge. The sampling routine will then determine
	 *     the other endpoint.
	 * \param succ
	 *     Pointer to an unused vertex data structure, which will be filled
	 *     with information about the successor vertex
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \return \c true on success
	 */
	bool sampleNext(const Scene *scene, Sampler *sampler,
			const PathVertex *pred, const Ray &ray, PathVertex *next,
			ETransportMode mode);

	/**
	 * \brief Create a perturbed successor vertex and edge
	 *
	 * This function behaves similar to \ref sampleNext() in that it
	 * generates a successor edge and vertex.
	 *
	 * The main difference is that the desired direction, distance, and
	 * type of the successor vertex are all specified, which makes the
	 * sampling process completely deterministic. This is useful for
	 * implementing path-space perturbation strategies.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex
	 * \param ray
	 *     Specifies the direction and origin associated with one
	 *     endpoint of the edge. The sampling routine will then determine
	 *     the other endpoint.
	 * \param dist
	 *     Specifies the desired distance between the current vertex and \c succ
	 *     (this only applies when <tt>desiredType=EMediumInteraction</tt>)
	 * \param desiredType
	 *     Specifies the desired vertex type of \c succ.
	 * \param succ
	 *     Pointer to an unused vertex data structure, which will be filled
	 *     with information about the successor vertex
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \return \c true on success
	 */
	bool perturbDirection(const Scene *scene, const PathVertex *pred,
			const Ray &ray, Float dist, PathVertex::EVertexType desiredType,
			PathVertex *next, ETransportMode mode);

	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                    Query functions
	/* ==================================================================== */

	enum ECachedValues {
		EValueImp                 = 0x01,
		EValueRad                 = 0x02,
		ECosineImp                = 0x04,
		ECosineRad                = 0x08,
		EValue                    = EValueImp | EValueRad,
		ECosine                   = ECosineImp | ECosineRad,
		EValueCosineImp           = EValueImp | ECosineImp,
		EValueCosineRad           = EValueRad | ECosineRad,
		EInverseSquareFalloff     = 0x10,
		ETransmittance            = 0x20,
		EGeometricTerm            = ECosine | EInverseSquareFalloff,
		EGeneralizedGeometricTerm = EGeometricTerm | ETransmittance,
		EEverything               = EValue | EGeneralizedGeometricTerm
	};

	/**
	 * \brief Evaluate cached quantities associated with this edge
	 *
	 * This function computes the product of certain terms that are cached
	 * in this edge and its adjacent vertices. The \c what parameter specifies
	 * the terms to be included; it must be a combination of the flags
	 * in the enumeration \ref ECachedValues.
	 *
	 * \remark This function assumes that \c pred and \c succ are the
	 * vertices associated with this edge, and that they have not been
	 * modified since the edge was created.
	 *
	 * \param pred
	 *     The predecessor vertex of this edge
	 * \param succ
	 *     The successor vertex of this edge
	 */
	Spectrum evalCached(const PathVertex *pred, const PathVertex *succ,
			unsigned int what) const;

	/**
	 * \brief Compute the density of a successor node
	 *
	 * This function computes the hypothetical transport-related sampling density
	 * of a given successor node conditioned on a specified predecessor when
	 * using the sampling technique implemented by \ref sampleNext(). Depending
	 * on whether or not the successor node is a medium interaction, the returned
	 * value is either a density per unit length or a discrete probability.
	 *
	 * Note: this function only computes terms associated with the transport
	 * between vertices -- to account for the vertices themselves,
	 * refer to  \ref PathEdge::evalPdf.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex
	 * \param succ
	 *     Pointer to the successor vertex
	 *
	 * \return The computed probability density
	 */
	Float evalPdf(const PathVertex *pred, const PathVertex *succ) const;

	/**
	 * \brief Compute the transmittance between an arbitrary
	 * pair of vertices
	 *
	 * This function queries the medium associated with this edge for the
	 * transmittance between an arbitrary pair of nodes, \c pred and \c succ.
	 *
	 * \param pred
	 *     Pointer to the preceding vertex
	 * \param succ
	 *     Pointer to the successor vertex
	 *
	 * \return A spectrally varying transmittance value
	 */
	Spectrum evalTransmittance(const PathVertex *pred, const PathVertex *succ) const;

	/**
	 * \brief Return the transmittance value associated with this edge
	 *
	 * \param pred
	 *     Pointer to the preceding vertex
	 * \param succ
	 *     Pointer to the successor vertex
	 *
	 * \return A spectrally varying transmittance value
	 */
	Spectrum evalTransmittance() const;

	/**
	 * \brief Evaluate the one of the cosine factors associated with the geometric
	 * term over an edge
	 */
	Float evalCosine(const PathVertex *pred, const PathVertex *succ, const PathVertex *base) const;


	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                    Miscellaneous
  	/* ==================================================================== */

	/**
	 * \brief Verify the cached values stored in this path vertex
	 * for consistency
	 *
	 * This function re-evaluates a series of quantities associated with
	 * this edge and compares them to locally cached values.
	 * If any mismatch is found, the function sends debug output to a
	 * specified output stream and returns \c false.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the vertex adjacent in the emitter direction or \c NULL
	 * \param succ
	 *     Pointer to the vertex adjacent in the sensor direction or \c NULL
	 * \param mode
	 *     Transport mode -- disambiguates the meaning of \c pred and \c succ.
	 * \param os
	 *     Target output stream for error messages
	 */
	bool verify(const Scene *scene, const PathVertex *adjL,
		const PathVertex *adjE, ETransportMode mode, std::ostream &os) const;

	/**
	 * \brief Create a connection edge between two vertices
	 *
	 * This function can be used to create an edge data structure
	 * when connecting two separately sampled path vertices. This
	 * involves checking that they are mutually visible and computing
	 * the attenuation caused by the medium in between (if any).
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param predEdge
	 *     Pointer to an edge between \c vs and its predecessor
	 *     (which is not needed by this function)
	 * \param vs
	 *     First path vertex to be connected.
	 * \param vt
	 *     Second path vertex to be connected.
	 * \param succEdge
	 *     Pointer to an edge between \c vt and its successor
	 *     (which is not needed by this function)
	 * \return \c true upon success, \c false when there is no
	 *     throughput or an inconsistency has been detected.
	 */
	bool connect(const Scene *scene, const PathEdge *predEdge,
		const PathVertex *vs, const PathVertex *vt,
		const PathEdge *succEdge);

	/**
	 * \brief Create a connection path between two vertices
	 *
	 * This function is conceptually similar to \ref connect().
	 * However, instead of a single edge, it potentially generates
	 * an entire connection path, where intermediate vertices are
	 * either index-matched medium transitions or other surface
	 * scattering events of type \ref BSDF::ENull.
	 *
	 * This is important to support efficient direct illumination sampling
	 * through such surfaces (e.g. a heterogeneous medium or a leaf with
	 * textured alpha transparency).
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param predEdge
	 *     Pointer to an edge between \c vs and its predecessor
	 *     (which is not needed by this function)
	 * \param vs
	 *     First path vertex to be connected.
	 * \param result
	 *     A path data structure that will be filled with the
	 *     created vertices and edges
	 * \param vt
	 *     Second path vertex to be connected.
	 * \param succEdge
	 *     Pointer to an edge between \c vt and its successor
	 *     (which is not needed by this function)
	 * \param maxInteractions
	 *     Specifies the maximum permissible number of intermediate
	 *     vertices (-1 == arbitrarily many)
	 * \param pool
	 *     Reference to memory pool that will be used to allocate
	 *     edges and vertices.
	 * \return \c true upon success, \c false when there is no
	 *     throughput or an inconsistency has been detected.
	 */
	static bool pathConnect(const Scene *scene, const PathEdge *predEdge,
		const PathVertex *vs, Path &result, const PathVertex *vt,
		const PathEdge *succEdge, int maxInteractions, MemoryPool &pool);

	/**
	 * \brief Create a connection path between two vertices and
	 * collapse it into a single edge that summarizes its properties
	 *
	 * This function can be thought of as being half-way in between
	 * \c connect() and \c pathConnect(). Like \c pathConnect(), it
	 * potentially generates an entire connection path between the
	 * specified endpoints, where intermediate vertices are
	 * either index-matched medium transitions or other surface
	 * scattering events of type \ref BSDF::ENull.
	 *
	 * This is important to support efficient direct illumination sampling
	 * through such surfaces (e.g. a heterogeneous medium or a leaf with
	 * textured alpha transparency).
	 *
	 * However, this variant does not return the intermediate vertices
	 * and edges -- instead, everything is collapsed into a single
	 * edge that captures the aggregate weight and probability densities.
	 *
	 * This function is used by bidirectional path tracing, since it creates
	 * connections through index-matched boundaries but does not require
	 * explicit knowledge about the associated path vertices.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param predEdge
	 *     Pointer to an edge between \c vs and its predecessor
	 *     (which is not needed by this function)
	 * \param vs
	 *     First path vertex to be connected.
	 * \param vt
	 *     Second path vertex to be connected.
	 * \param succEdge
	 *     Pointer to an edge between \c vt and its successor
	 *     (which is not needed by this function)
	 * \param interactions
	 *    Specifies the maximum permissible number of index-matched medium
	 *    transitions or \ref BSDF::ENull scattering events on the way
	 *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
	 *    When the function is successful, this parameter will
	 *    additionally be used to return the actual number of intermediate
	 *    interactions.
	 * \return \c true upon success, \c false when there is no
	 *     throughput or an inconsistency has been detected.
	 */
	bool pathConnectAndCollapse(const Scene *scene, const PathEdge *predEdge,
		const PathVertex *vs, const PathVertex *vt,
		const PathEdge *succEdge, int &interactions);

	/// Create a deep copy of this edge
	PathEdge *clone(MemoryPool &pool) const;

	/// Compare this edge against another edge
	bool operator==(const PathEdge &edge) const;

	/// Compare this edge against another edge
	inline bool operator!=(const PathEdge &edge) const {
		return !operator==(edge);
	}

	/// Return a string representation of the information stored in this vertex
	std::string toString() const;

	//! @}
	/* ==================================================================== */
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_EDGE_H_ */
