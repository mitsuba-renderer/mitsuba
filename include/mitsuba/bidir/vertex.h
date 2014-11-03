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
#if !defined(__MITSUBA_BIDIR_VERTEX_H_)
#define __MITSUBA_BIDIR_VERTEX_H_

#include <mitsuba/bidir/common.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Bidirectional path vertex data structure
 *
 * The path vertex data structure represents a basic interaction within the
 * path-space light transport framework. It abstracts away the specifics
 * of the underlying interaction, which simplifies the implementation of
 * bidirectional rendering techniques such as Bidirectional Path Tracing
 * or Veach-style Metropolis Light Transport.
 *
 * This data structure can describe a several different types of interactions,
 * including surface and medium scattering events, as well as emission events
 * by a light source or a sensor (in the bidirectional framework, the response
 * of a sensor is treated as an emitted quantity)
 *
 * A path vertex only describes what happens at the location of a scattering
 * or emission event -- what happens <em>in between</em> such interactions is
 * captured in a separate data structure named \ref PathEdge.
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
struct MTS_EXPORT_BIDIR PathVertex {
	/* ==================================================================== */
	//! @{ \name             Enumerations and Fields
	/* ==================================================================== */

	/// Denotes the size of the auxiliary data section associated with each node
	enum {
		EDataSize = MAX(MAX(sizeof(Intersection),
			sizeof(MediumSamplingRecord)), MAX(
			sizeof(PositionSamplingRecord), sizeof(EndpointRecord)))
	};

	/**
	 * \brief What kind of vertex is this (e.g. medium, surface, emitter)?
	 *
	 * The two special 'supernode' types are used as path endpoints, which greatly
	 * simplifies the control flow of various functions in this data structure and the
	 * MLT/BDPT implementations.
	 */
	enum EVertexType {
		/// Invalid/uninitialized vertex
		EInvalid = 0,
		/// Special sensor 'supernode' -- just before the sensor sample
		ESensorSupernode = 1,
		/// Special emitter 'supernode' -- just before the emitter sample
		EEmitterSupernode = 2,
		/// Sampled position on the surface of a sensor
		ESensorSample = 4,
		/// Sampled position on the surface of an emitter
		EEmitterSample = 8,
		/// Interaction with a scene surface
		ESurfaceInteraction = 16,
		/// Interaction with a participating medium
		EMediumInteraction = 32,
		/// Union of the two supernode types
		ESupernode = ESensorSupernode | EEmitterSupernode,
		/// Union of all other vertex types
		ENormal = ESensorSample | EEmitterSample
			| ESurfaceInteraction | EMediumInteraction
	};

	/**
	 * \brief Specifies one of several possible path vertex types
	 *
	 * \sa EVertexType
	 */
	EVertexType type : 7;

	/**
	 * \brief Denotes whether this vertex \a only supports
	 * sampling from degenerate distributions.
	 *
	 * It's useful to cache this information, since it allows
	 * to quickly determine when certain pairs of vertices
	 * cannot be deterministically connected.
	 *
	 * Examples of degenerate vertices are surface interactions with
	 * dielectric boundaries and certain special cases. For instance, the
	 * sensor supernode can be degenerate when the used sensor does not
	 * cover any area (e.g. because it is a point camera).
	 */
	bool degenerate : 1;

	/**
	 * \brief Denotes the measure associated with the
	 * probability densities stored in \ref pdf.
	 *
	 * Certain vertices use sampling methods, whose probability
	 * mass lies on domains that have different associated
	 * measures. An example would be a surface interaction with
	 * a smooth plastic material, where the specular reflection
	 * component is degenerate and the glossy reflection component
	 * is non-degenerate.
	 *
	 * This attribute stores the actual measure (an enumeration
	 * item of type \ref EMeasure) and therefore clarifies whether
	 * or not a degenerate component was sampled. When no sampling
	 * event was generated yet, it is undefined.
	 *
	 * Currently, only three values are permissible:
	 * \ref EArea, \ref EDiscrete, or \ref EInvalidMeasure.
	 *
	 * \sa EMeasure
	 */
	EMeasure measure : 8;

	/**
	 * \brief When the current vertex supports sampling
	 * from several components (this currently only applies
	 * to %BSDF lobes), this attribute records the type
	 * sampled component.
	 *
	 * Otherwise, it will be set to zero.
	 */
	uint16_t componentType;

	/**
	 * \brief Measurement contribution weight
	 *
	 * This field stores the terms of the path-space measurement contribution
	 * function that are coupled to this specific vertex, divided by the
	 * density of the adjacent vertices in the radiance and importance transport
	 * directions (hence, it is an array with two entries).
	 *
	 * More precisely, it stores
	 * \f[
	 * f(\mathbf{x}_{i-1}\to \mathbf{x}_i\to\mathbf{x}_{i+1})\,
	 * G(\mathbf{x}_i\leftrightarrow\mathbf{x}_{i+1})\,p_A(\mathbf{x}_{i-1}\to
	 * \mathbf{x}_i\to\mathbf{x}_{i+1})^{-1}
	 * \f]
	 *
	 * and
	 *
	 * \f[
	 * f(\mathbf{x}_{i+1}\to \mathbf{x}_i\to\mathbf{x}_{i-1})\,
	 * G(\mathbf{x}_i\leftrightarrow\mathbf{x}_{i-1})\,p_A(\mathbf{x}_{i+1}\to
	 * \mathbf{x}_i\to\mathbf{x}_{i-1})^{-1}
	 * \f]
	 *
	 * Where \f$G\f$ is the geometric term, \f$p_A\f$ is an area density, and
	 * increasing indices are closer to the sensor. Generally much cancellation
	 * will occur in the above expressions. For instance, for a surface
	 * iteractions, this is equal to the BRDF times a cosine foreshortening
	 * factor divided by the solid angle density of the default sampling
	 * method.
	 *
	 * Note that this field does not account for medium-related terms. These
	 * can be found in \ref PathEdge::weight
	 */
	Spectrum weight[ETransportModes];

	/**
	 * \brief Area density of the two adjacent vertices
	 *
	 * This field stores the density of the predecessor and sucessor nodes
	 * with respect of the sampling technique implemented by \ref sampleNext().
	 * The measure of this value is specified by the \ref measure field
	 * (generally, it is the density per unit area).
	 *
	 * When one of the adjacent vertices is a medium interaction (i.e. it is
	 * not located on a surface), the stored probability will specify the density
	 * on a hypothetical surface oriented perpendicularly to the transport
	 * direction.
	 *
	 * Note that this field does not account for medium-related terms. When
	 * an adjacent vertex is a medium interaction, its volume density can be
	 * recovered by computing the product of \c pdf and \ref PathEdge::pdf
	 * of the associated transport edge.
	 */
	Float pdf[ETransportModes];

	/// \brief Termination weight due to russian roulette (used by BDPT)
	Float rrWeight;

	/**
	 * \brief Auxilary node-depependent data associated with each vertex
	 *
	 * This "payload" field is large enough to describe any possible kind of
	 * vertex. Currently, it is cast to the desired type of data structure
	 * which wreaks havoc with strict aliasing (and hence it is disabled
	 * for all bidirectional code).
	 *
	 * Once C++11 is widely supported across all target platforms, this
	 * should be replaced by an unrestricted union.
	 */
	uint8_t data[EDataSize];

	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name              Sampling-related functions
	/* ==================================================================== */

	/**
	 * \brief Generate a path endpoint that can be used to start
	 * a random walk
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param mode
	 *     Specifies the desired mode of transport, i.e. radiance or
	 *     importance transport
	 * \param time
	 *     Denotes the time value that will be associated with this
	 *     endpoint and any paths generated starting from there.
	 */
	void makeEndpoint(const Scene *scene, Float time, ETransportMode mode);

	/**
	 * \brief Sample the next vertex in a random walk using the default
	 * sampling technique implemented by the current vertex.
	 *
	 * Given a vertex, its predecessor, as well as the edge in between them,
	 * this function samples a new successor edge and vertex.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param pred
	 *     Pointer to the preceding vertex (if any) and \c NULL otherwise
	 * \param predEdge
	 *     Edge to the preceding edge (if any) and \c NULL otherwise
	 * \param succEdge
	 *     Pointer to an unused edge data structure, which will be filled
	 *     with information about the edge between the current vertex and
	 *     the newly sampled successor.
	 * \param succ
	 *     Pointer to an unused vertex data structure, which will be filled
	 *     with information about the successor vertex
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \param russianRoulette
	 *     Should russian roulette be used while sampling the successor?
	 *     Note that the effects of this are only captured in the \c rrWeight
	 *     field -- the \ref weight and \ref pdf fields intentionally remain unchanged.
	 * \param throughput
	 *     If russian roulette is active, this parameter should point to a
	 *     spectrum value that is used to record the aggregate path weight
	 *     thus far. It will be updated automatically to account for the current
	 *     interaction.
	 * \return \c true on success
	 */
	bool sampleNext(const Scene *scene, Sampler *sampler,
		const PathVertex *pred, const PathEdge *predEdge,
		PathEdge *succEdge, PathVertex *succ,
		ETransportMode mode, bool russianRoulette = false,
		Spectrum *throughput = NULL);

	/**
	 * \brief \a Direct sampling: given the current vertex as a reference
	 * sample an emitter (or sensor) position that has a nonzero emission
	 * (or response) towards it.
	 *
	 * This can be seen as a generalization of direct illumination sampling
	 * that can be used for both emitter and sensor endpoints.
	 *
	 * Ideally, the implementation should importance sample the product of
	 * the emission or response profile and the geometry term between the
	 * reference point and the sampled position. In practice, one of these
	 * usually has to be sacrificed.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param endpoint
	 *     Unused vertex data structure, which will be configured as
	 *     the endpoint associated with \c sample.
	 * \param edge
	 *     Unused edge data structure, which will be configured
	 *     as the edge between \c endpoint and \c sample.
	 * \param sample
	 *     Unused vertex data structure, which will hold the
	 *     sampled sensor or emitter position
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \return The emitted radiance or importance divided by the
	 *     sample probability per unit area per unit solid angle.
	 */
	Spectrum sampleDirect(const Scene *scene, Sampler *sampler,
		PathVertex *endpoint, PathEdge *edge, PathVertex *sample,
		ETransportMode mode) const;

	/**
	 * \brief Sample the first vertices on a sensor subpath such that
	 * they contribute to a specified pixel in the output image
	 *
	 * This function samples the spatial and directional components of the sensor
	 * and is similar to calling \ref sampleNext() two times in sequence starting
	 * from a sensor supernode. The main difference is that the resulting subpath
	 * passes through a specified pixel position, which is important to implement
	 * algorithms that parallelize rendering of images by processing it in separate
	 * blocks. If this function is called once for every pixel in the output
	 * image, the resulting path distribution is identical to what would have
	 * been obtained via \ref sampleNext().
	 *
	 * The function throws an exception when the current vertex is not
	 * a sensor supernode.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param pixelPosition
	 *     Specifies the desired pixel position
	 * \param e0
	 *     Pointer to the first edge on the sensor subpath
	 * \param e0
	 *     Pointer to the first edge on the sensor subpath
	 * \param v1
	 *     Pointer to the second vertex on the sensor subpath
	 * \param v2
	 *     Pointer to the third vertex on the sensor subpath
	 * \return The number of successful sampling operations (i.e.
	 *    \c 0 if sensor sampling failed, \c 1 if no surface/medium
	 *    interaction was encountered, or \c 2 if all data structures
	 *    were successfully filled)
	 */
	int sampleSensor(const Scene *scene, Sampler *sampler, const Point2i &pixelPosition,
		PathEdge *e0, PathVertex *v1, PathEdge *e1, PathVertex *v2);

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
	 * This function only applies to non-supernode vertices.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex (if any) and \c NULL otherwise
	 * \param predEdge
	 *     Edge to the preceding edge (if any) and \c NULL otherwise
	 * \param succEdge
	 *     Pointer to an unused edge data structure, which will be filled
	 *     with information about the edge between the current vertex and
	 *     the new successor.
	 * \param succ
	 *     Pointer to an unused vertex data structure, which will be filled
	 *     with information about the successor vertex
	 * \param d
	 *     Specifies the desired outgoing direction at the current vertex
	 * \param dist
	 *     Specifies the desired distance between the current vertex and \c succ
	 *     (this only applies when <tt>desiredType=EMediumInteraction</tt>)
	 * \param desiredType
	 *     Specifies the desired vertex type of \c succ.
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \return \c true on success
	 */
	bool perturbDirection(const Scene *scene, const PathVertex *pred,
		const PathEdge *predEdge, PathEdge *succEdge, PathVertex *succ,
		const Vector &d, Float dist, EVertexType desiredType, ETransportMode mode);

	/**
	 *
	 */
	bool perturbPosition(const Scene *scene, Sampler *sampler, Float stddev);
	Float perturbPositionPdf(const PathVertex *target, Float stddev) const;

	/**
	 * \brief Propagate a perturbation through an ideally specular interaction
	 *
	 * This function behaves similar to \ref sampleNext() and \ref
	 * perturbDirection() in that it generates a successor edge and vertex.
	 *
	 * The main difference is that it only works for specular interactions,
	 * where the requested type of interaction (reflection/refraction) is
	 * additionally specified, which makese the sampling process completely
	 * deterministic. This is useful for implementing path-space
	 * perturbation strategies. For now, it is only used by the perturbations
	 * of Veach and Guibas.
	 *
	 * This function only applies to surface interaction vertices.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex (if any) and \c NULL otherwise
	 * \param predEdge
	 *     Edge to the preceding edge (if any) and \c NULL otherwise
	 * \param succEdge
	 *     Pointer to an unused edge data structure, which will be filled
	 *     with information about the edge between the current vertex and
	 *     the new successor.
	 * \param dist
	 *     Specifies the desired distance between the current vertex and \c succ
	 *     (this only applies when <tt>desiredType=EMediumInteraction</tt>)
	 * \param desiredType
	 *     Specifies the desired vertex type of \c succ.
	 * \param componentType
	 *     Specifies the desired type of scattering interaction (equivalent
	 *     to \ref BSDF::typeMask)
	 * \param succ
	 *     Pointer to an unused vertex data structure, which will be filled
	 *     with information about the successor vertex
	 * \param dist
	 *     Specifies the desired distance between the current vertex and
	 *     \c succ (this only applies to medium interactions)
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \return \c true on success
	 */
	bool propagatePerturbation(const Scene *scene, const PathVertex *pred,
		const PathEdge *predEdge, PathEdge *succEdge, PathVertex *succ,
		unsigned int componentType, Float dist, EVertexType desiredType,
		ETransportMode mode);

	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                    Query functions
	/* ==================================================================== */

	/**
	 * \brief Evaluate the terms of the measurement contribution function
	 * that are associated with this vertex.
	 *
	 * Compute a single term of the contribution weighting function associated
	 * with the current node given a predecessor and successor node.
	 *
	 * Note: this function only accounts for the factor associated with this
	 * specific vertex -- to account for an adjacent edge, refer to
	 * \ref PathEdge::eval.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex (if any) and \c NULL otherwise
	 * \param succ
	 *     Pointer to the successor vertex (if any) and \c NULL otherwise
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \param measure
	 *     Specifies the measure of the queried component. This is necessary
	 *     to handle mixture scattering functions, whose components are
	 *     defined on spaces with different measures.
	 * \return The contribution weighting factor
	 */
	Spectrum eval(const Scene *scene, const PathVertex *pred,
		const PathVertex *succ, ETransportMode mode, EMeasure measure = EArea) const;

	/**
	 * \brief Compute the density of a successor node
	 *
	 * This function computes the hypothetical scattering-related sampling density
	 * of a given successor node when using the sampling technique implemented by
	 * \ref sampleNext(). Since this technique conditions on a predecessor vertex,
	 * it must also be provided here. The desired measure (e.g. area/solid angle/discrete)
	 * can be provided as an extra parameter.
	 *
	 * Note: this function only computes probability associated with the
	 * vertices -- to account for edges, refer to \ref PathEdge::evalPdf.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex (if any) and \c NULL otherwise
	 * \param succ
	 *     Pointer to the successor vertex (if any) and \c NULL otherwise
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \param measure
	 *     Specifies the measure of the queried component. This is necessary
	 *     to handle mixture scattering functions, whose components are
	 *     defined on spaces with different measures.
	 * \return The computed probability density
	 */
	Float evalPdf(const Scene *scene, const PathVertex *pred,
		const PathVertex *succ, ETransportMode mode, EMeasure measure = EArea) const;

	/**
	 * \brief Compute the area density of a provided emitter or sensor
	 * sample with respect the \a direct sampling technique implemented in
	 * \ref sampleDirect().
	 *
	 * The current vertex is taken to be the reference point of the direct
	 * sampling technique.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param mode
	 *     Specifies whether radiance or importance is being transported
	 * \param sample
	 *     An emitter or sensor sample
	 * \param measure
	 *     Specifies the measure of the queried component. This is necessary
	 *     to handle scattering functions, whose components are
	 *     defined on spaces with different measures.
	 * \return The density of \c sample conditioned on this vertex.
	 */
	Float evalPdfDirect(const Scene *scene, const PathVertex *sample,
		ETransportMode mode, EMeasure measure = EArea) const;

	/**
	 * \brief Determine the medium that fills the space between
	 * the current vertex and a specified successor
	 *
	 * This function assumes that there is no surface between \c this
	 * and \c succ, hence it does not account for intermediate medium
	 * changes.
	 *
	 * \param predEdge
	 *    Pointer to an edge that connects the current node to its predecessor
	 * \param succ
	 *    Pointer to the successor vertex in question
	 * \return
	 *    The medium between \c this and \c succ
	 */
	const Medium *getTargetMedium(const PathEdge *predEdge,
			const PathVertex *succ) const;

	/**
	 * \brief Determine the medium that fills the space in direction \c d
	 *
	 * \param predEdge
	 *    Pointer to an edge that connects the current node to its predecessor
	 * \param d
	 *    A world-space direction
	 * \return
	 *    The medium containing the ray <tt>(this->getPosition(), d)</tt>
	 */
	const Medium *getTargetMedium(const PathEdge *predEdge,
			const Vector &d) const;

	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                     Accessors
    /* ==================================================================== */

	/// Return the type associated with this vertex
	inline EVertexType getType() const { return (EVertexType) type; }

	/**
	 * \brief Return the position associated with this vertex
	 *
	 * Throws an exception when called on a supernode
	 */
	Point getPosition() const;

	/// Check if the vertex lies on a surface
	inline bool isOnSurface() const {
		return (type == ESurfaceInteraction) || ((type == EEmitterSample ||
			type == ESensorSample) && static_cast<const AbstractEmitter *>(
			getPositionSamplingRecord().object)->getType() & AbstractEmitter::EOnSurface);
	}

	/**
	 * \brief Returns whether or not this vertex describes a "null" scattering interaction.
	 *
	 * A null interaction is a degenerate scattering event with a Dirac delta peak in the
	 * forward direction. Apart from a potential influence on their weight, particles
	 * will pass through such an interface unchanged.
 	 */
	inline bool isNullInteraction() const {
		return type == ESurfaceInteraction && componentType == BSDF::ENull;
	}

	/**
	 * \brief Returns whether or not this vertex describes a diffuse surface
	 * scattering interaction.
 	 */
	inline bool isDiffuseInteraction() const {
		return type == ESurfaceInteraction &&
			(componentType == BSDF::EDiffuseReflection || componentType == BSDF::EDiffuseTransmission);
	}

	/**
	 * \brief Returns whether or not this vertex describes a 100% absorbing surface
	 *
	 * Such is the case on emitters/sensors that don't have an explicit BSDF assigned
	 * to them. It is useful to be able to query this to avoid some useless connection
	 * attempts involving these vertices.
 	 */
	inline bool isAbsorbing() const {
		return (type == ESurfaceInteraction &&
			(getIntersection().shape->getBSDF()->getType() & BSDF::EAll) == 0);
	}

	/**
	 * \brief Return the component type associated with this vertex
	 *
	 * This currently only applies to surface interactions. The returned
	 * result will consist of flags in \ref BSDF::EBSDFType.
	 */
	inline unsigned int getComponentType() const { return (unsigned int) componentType; }

	/**
	 * \brief Return the geometric surface normal associated with this vertex
	 *
	 * Throws an exception when called on a supernode or a medium interaction
	 */
	Normal getGeometricNormal() const;

	/**
	 * \brief Return the shading surface normal associated with this vertex
	 *
	 * Throws an exception when called on a supernode or a medium interaction
	 */
	Normal getShadingNormal() const;

	/// Return the time value associated with this node
	Float getTime() const;

	/// Is this vertex a supernode?
	inline bool isSupernode() const { return type & ESupernode; }
	/// Is this vertex a sensor super-node?
	inline bool isSensorSupernode() const { return type == ESensorSupernode; }
	/// Is this vertex a emitter super-node?
	inline bool isEmitterSupernode() const { return type == EEmitterSupernode; }
	/// Is this vertex an emitter sample?
	inline bool isEmitterSample() const { return type == EEmitterSample; }
	/// Is this vertex a lens sample?
	inline bool isSensorSample() const { return type == ESensorSample; }
	/// Is this vertex a surface interaction?
	inline bool isSurfaceInteraction() const { return type == ESurfaceInteraction; }
	/// Is this vertex a medium interaction?
	inline bool isMediumInteraction() const { return type == EMediumInteraction; }

	/// Return the endpoint record associated with this node
	inline EndpointRecord &getEndpointRecord() {
		EndpointRecord* ptr MTS_MAY_ALIAS =
			reinterpret_cast<EndpointRecord*>(&data[0]);
		return *ptr;
	}

	/// Return the endpoint record associated with this node (const)
	inline const EndpointRecord &getEndpointRecord() const {
		const EndpointRecord* ptr MTS_MAY_ALIAS =
			reinterpret_cast<const EndpointRecord*>(&data[0]);
		return *ptr;
	}

	/// Return the position sampling record associated with this node
	inline PositionSamplingRecord &getPositionSamplingRecord() {
		PositionSamplingRecord* ptr MTS_MAY_ALIAS =
			reinterpret_cast<PositionSamplingRecord*>(&data[0]);
		return *ptr;
	}

	/// Return the position sampling record associated with this node (const)
	inline const PositionSamplingRecord &getPositionSamplingRecord() const {
		const PositionSamplingRecord* ptr MTS_MAY_ALIAS =
			reinterpret_cast<const PositionSamplingRecord*>(&data[0]);
		return *ptr;
	}

	/// Return the intersection record associated with this node
	inline Intersection &getIntersection() {
		Intersection* ptr MTS_MAY_ALIAS =
			reinterpret_cast<Intersection*>(&data[0]);
		return *ptr;
	}

	/// Return the intersection record associated with this node (const)
	inline const Intersection &getIntersection() const {
		const Intersection* ptr MTS_MAY_ALIAS =
			reinterpret_cast<const Intersection*>(&data[0]);
		return *ptr;
	}

	/// Return the medium sampling record associated with this node
	inline MediumSamplingRecord &getMediumSamplingRecord() {
		MediumSamplingRecord* ptr MTS_MAY_ALIAS =
			reinterpret_cast<MediumSamplingRecord*>(&data[0]);
		return *ptr;
	}

	/// Return the medium sampling record associated with this node
	inline const MediumSamplingRecord &getMediumSamplingRecord() const {
		const MediumSamplingRecord* ptr MTS_MAY_ALIAS =
			reinterpret_cast<const MediumSamplingRecord*>(&data[0]);
		return *ptr;
	}

	/// Return the fractional pixel position associated with a sensor sample
	inline const Point2 &getSamplePosition() const {
		return getPositionSamplingRecord().uv;
	}

	/// Return the abstract emitter associated with a sensor/emitter sample
	inline const AbstractEmitter *getAbstractEmitter() const {
		return static_cast<const AbstractEmitter *>(
				getPositionSamplingRecord().object);
	}

	/**
	 * \brief Returns whether or not this vertex is degenerate, i.e.
	 * its distribution has measure zero
	 */
	inline bool isDegenerate() const { return degenerate; }

	/**
	 * \brief Returns whether or not this vertex can be deterministically
	 * connected to other vertices.
	 *
	 * This is the case when <tt>\ref degenerate == false</tt> and
	 * <tt>\ref measure != EDiscrete</tt>.
	 */
	inline bool isConnectable() const { return !degenerate && measure != EDiscrete; }

	/**
	 * \brief Special routine for sensor sample vertices: given a successor
	 * \c succ, update the fractional pixel position stored in the vertex.
	 *
	 * \param v
	 *    Pointer to the target vertex
	 *
	 * \return \c true upon success (i.e. when the point is in the
	 *     sensor's field of view)
	 */
	bool updateSamplePosition(const PathVertex *succ);

	/**
	 * \brief Special routine for sensor sample vertices: given a successor
	 * \c succ, return its associated fractional pixel position.
	 *
	 * \param v
	 *    Pointer to the target vertex
	 *
	 * \param result
	 *    Reference to a 2D point that will be set to the fractional
	 *    pixel coordinates associated with \c succ.
	 *
	 * \return \c true upon success (i.e. when the point is in the
	 *     sensor's field of view)
	 */
	bool getSamplePosition(const PathVertex *succ, Point2 &result) const;

	//! @}
	/* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                    Miscellaneous
 	/* ==================================================================== */

	/**
	 * \brief Cast this vertex into an equivalent from having a different type
	 *
	 * Sometimes it is necessary to cast a vertex into a different type.
	 * An example when this occurs is when a surface interaction vertex
	 * lies on the surface of an emitter or a sensor. In such a situation,
	 * it may be useful to retroactively turn it into an emitter or sensor
	 * sample vertex located at the same position. This function allows
	 * to do precisely that.
	 *
	 * \param scene
	 *     A pointer to the underlying scene
	 * \param desired
	 *     Desired type after the cast
	 * \return \c true on success. When returning \c false, the
	 *     function did not make any changes.
	 */
	bool cast(const Scene *scene, EVertexType desired);

	/**
	 * \brief Verify the cached values stored in this path vertex
	 * for consistency
	 *
	 * This function re-evaluates a series of quantities associated with
	 * this vertex and compares them to locally cached values including
	 * \ref pdf, \ref value, and \ref degenerate. If any mismatch
	 * is found, the function sends debug output to a specified output
	 * stream and returns \c false.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the preceding vertex (if any) and \c NULL otherwise
	 * \param succ
	 *     Pointer to the successor vertex (if any) and \c NULL otherwise
	 * \param mode
	 *     Transport mode -- disambiguates the meaning of \c pred and \c succ.
	 * \param os
	 *     Target output stream for error messages
	 */
	bool verify(const Scene *scene, const PathVertex *adjL,
		const PathVertex *adjE, ETransportMode mode, std::ostream &os) const;

	/**
	 * \brief Given the specified predecessor and successor, update
	 * the cached values stored in this vertex
	 *
	 * \param pred
	 *     Pointer to the predecessor vertex (if any) and \c NULL otherwise
	 * \param succ
	 *     Pointer to the successor vertex (if any) and \c NULL otherwise
	 * \param mode
	 *     Specifies the direction of light transport
	 * \return \c false when there is no throughput
	 */
	bool update(const Scene *scene, const PathVertex *pred,
		const PathVertex *succ, ETransportMode mode, EMeasure measure = EArea);

	/**
	 * \brief Create a connection between two disconnected subpaths
	 *
	 * This function can be used to connect two seperately created emitter
	 * and sensor subpaths so that they can be merged into a \ref Path data
	 * structure. The function checks that the vertices \c vs and \c vt are
	 * mutually visible, and that there is a nonzero throughput between them.
	 * If that is the case, it updates the cached values stored in \c vs,
	 * \c edge, and \c vt.
	 *
	 * The expected order of the parameters in path-space is
	 * <pre>
	 * (pred) -> predEdge -> (vs) -> edge -> (vt) -> succEdge -> (succ)
	 * </pre>
	 * where entries in parentheses denote vertices,
	 * \c pred is the closer to the light source, and \c succ
	 * is the closer to the sensor.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param pred
	 *     Pointer to the predecessor vertex of \c vs (towards the emitter)
	 * \param predEdge
	 *     Pointer to an edge between \c pred and \c vs.
	 * \param vs
	 *     Last vertex of the emitter subpath to be connected. The
	 *     cached values of this vertex will be updated should the
	 *     connection attempt succeed.
	 * \param edge
	 *     Pointer to an unused edge data structure, which will be
	 *     annotated with information about the medium-related transport
	 *     between \c vs and \c vt.
	 * \param vt
	 *     Last vertex of the sensor subpath to be connected. The
	 *     cached values of this vertex will be updated should the
	 *     connection attempt succeed.
	 * \param succEdge
	 *     Pointer to an edge between \c vt and \c succ.
	 * \param succ
	 *     Pointer to the successor vertex of \c vt (towards the sensor)
	 * \return \c true upon success, \c false when there is no
	 *     throughput or an inconsistency has been detected.
	 */
	static bool connect(const Scene *scene,
			const PathVertex *pred, const PathEdge *predEdge,
			PathVertex *vs, PathEdge *edge, PathVertex *vt,
			const PathEdge *succEdge, const PathVertex *succ);

	/// Like the above, but can be used to connect delta endpoints
	static bool connect(const Scene *scene,
			const PathVertex *pred, const PathEdge *predEdge,
			PathVertex *vs, PathEdge *edge, PathVertex *vt,
			const PathEdge *succEdge, const PathVertex *succ,
			EMeasure vsMeasure, EMeasure vtMeasure);

	/// Create a deep copy of this vertex
	PathVertex *clone(MemoryPool &pool) const;

	/// Return a string representation of the information stored in this vertex
	std::string toString() const;

	/// Compare this vertex against another vertex
	bool operator==(const PathVertex &vertex) const;

	/// Compare this vertex against another vertex
	inline bool operator!=(const PathVertex &vertex) const {
		return !operator==(vertex);
	}

	//! @}
	/* ==================================================================== */
};

/// \cond
extern MTS_EXPORT_BIDIR
	std::ostream &operator<<(std::ostream &os, PathVertex::EVertexType type);
/// \endcond

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_VERTEX_H_ */
