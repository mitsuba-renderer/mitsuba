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
#if !defined(__MITSUBA_BIDIR_PATH_H_)
#define __MITSUBA_BIDIR_PATH_H_

#include <mitsuba/bidir/vertex.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/bidir/mempool.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Bidirectional path data structure
 *
 * In the path-space light transport framework, a path is represented as a
 * linear sequence of interactions (expressed using vertices) and transport
 * (expressed using edges).
 *
 * The \ref Path data structure is responsible for the storage of this
 * information. It also contains useful utility functions, for instance
 * to perform a random walk, or to splice and connect path segments.
 *
 * \sa PathVertex
 * \sa PathEdge
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
struct MTS_EXPORT_BIDIR Path {
public:
	typedef PathEdge *          PathEdgePtr;
	typedef PathVertex *        PathVertexPtr;

	/* ==================================================================== */
	//! @{ \name              Path construction
	/* ==================================================================== */

	/// Create a new, empty path
	inline Path() { }

	/// Create a new path of the specified size
	inline Path(size_t size) : m_vertices(size) { }

	/// Copy constructor
	inline Path(const Path &path) : m_vertices(path.m_vertices),
		m_edges(path.m_edges) { }

	/**
	 * \brief Initialize the path with an endpoint vertex
	 *
	 * This function clears the path and initializes it with a single
	 * endpoint vertex of the type implied by the \c mode parameter.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param time
	 *     Specifies the time value to be associated with the path.
	 * \param mode
	 *     Specifies the desired endpoint type
	 * \param pool
	 *     Reference to a memory pool that will be used to release
	 *     and allocate edges and vertices.
	 */
	void initialize(const Scene *scene, Float time,
		ETransportMode mode, MemoryPool &pool);

	/**
	 * \brief Perform/continue a random walk starting from the current endpoint
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param nSteps
	 *     Desired number of random walk steps (<tt>-1</tt>=infinite)
	 * \param rrStart
	 *     Depth to start using russian roulette
	 *     (<tt>-1</tt>=never, <tt>0</tt>=starting at the first
	 *     bounce, and so on)
	 * \param mode
	 *     Denotes whether radiance or importance are being transported
	 * \param pool
	 *     Reference to a memory pool that will be used to allocate
	 *     edges and vertices.
	 * \return The number of successful steps performed by the random walk.
	 */
	int randomWalk(const Scene *scene, Sampler *sampler,
		int nSteps, int rrStart, ETransportMode mode,
		MemoryPool &pool);

	/**
	 * \brief Perform a random walk starting at a specified
	 * pixel on the sensor
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param nSteps
	 *     Desired number of random walk steps (<tt>-1</tt>=infinite)
	 * \param pixelPosition
	 *     Pixel position associated with the newly created
	 *     sensor subpath
	 * \param rrStart
	 *     Depth to start using russian roulette
	 *     (<tt>-1</tt>=never, <tt>0</tt>=starting at the first
	 *     bounce, and so on)
	 * \param pool
	 *     Reference to a memory pool that will be used to allocate
	 *     edges and vertices.
	 * \return The number of successful steps performed by the random walk.
	 */
	int randomWalkFromPixel(const Scene *scene, Sampler *sampler,
		int nSteps, const Point2i &pixelPosition, int rrStart,
		MemoryPool &pool);

	/**
	 * \brief Perform two random walks on an emitter and sensor subpath
	 *
	 * This function is almost identical to calling \ref randomWalk() twice
	 * in sequence. The main difference is that it performs the random
	 * walk steps in a staggered order (i.e. one step on the emitter subpath,
	 * one step on the sensor subpath, and so on..), which is important for
	 * obtaining good results with QMC random number sequences.
	 * Additinally, it ensures that the sensor path passes through a specified
	 * pixel instead of sampling it uniformly.
	 *
	 * Used by bidirectional path tracing.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param sampler
	 *     Pointer to a sample generator
	 * \param emitterPath
	 *     Reference to the emitter subpath to be filled
	 * \param nEmitterSteps
	 *     Desired number of random walk steps on the emitter subpath
	 *     (<tt>-1</tt>=infinite)
	 * \param sensorPath
	 *     Reference to the sensor subpath to be filled
	 * \param nSensorSteps
	 *     Desired number of random walk steps on the sensor subpath
	 *     (<tt>-1</tt>=infinite)
	 * \param pixelPosition
	 *     Pixel position associated with the newly created
	 *     sensor subpath
	 * \param rrStart
	 *     Depth to start using russian roulette
	 *     (<tt>-1</tt>=never, <tt>0</tt>=starting at the first
	 *     bounce, and so on)
	 * \param pool
	 *     Reference to a memory pool that will be used to allocate
	 *     edges and vertices.
	 * \return The number of successful steps performed by the random walk
	 *         on the emitter and sensor subpath, respectively.
	 */
	static std::pair<int, int> alternatingRandomWalkFromPixel(const Scene *scene,
		Sampler *sampler, Path &emitterPath, int nEmitterSteps,
		Path &sensorPath, int nSensorSteps, const Point2i &pixelPosition,
		int rrStart, MemoryPool &pool);

	/**
	 * \brief Verify the cached values stored in this path
	 *
	 * This function re-evaluates a series of quantities associated with
	 * each vertex and edge and compares them to locally cached values.
	 * If any mismatch is found, the function sends debug output to a
	 * specified output stream and returns \c false.
	 *
	 * \param scene
	 *     Pointer to the underlying scene
	 * \param mode
	 *     Disambiguates the order of the vertices in this path
	 * \param os
	 *     Target output stream for error messages
	 */
	bool verify(const Scene *scene, ETransportMode mode, std::ostream &os) const;

	//! @}
    /* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                     Accessors
    /* ==================================================================== */

	/**
	 * \brief Return the number of vertices stored in this path
	 *
	 * For a nonempty path, the number of vertices is always equal
	 * To \ref edgeCount()+1.
	 */
	inline size_t vertexCount() const {
		return m_vertices.size();
	}

	/**
	 * \brief Return the number of edges stored in this path
	 *
	 * For a nonempty path, the number of vertices is always equal
	 * To \ref vertexCount()-1.
	 */
	inline size_t edgeCount() const {
		return m_edges.size();
	}

	/**
	 * \brief Return the length of the path. This is just the
	 * number of edges.
	 */
	inline int length() const {
		return (int) m_edges.size();
	}

	/// Return an vertex by its index
	inline PathVertexPtr &vertex(size_t index) {
		#if MTS_BD_DEBUG == 1
			if (index >= m_vertices.size())
				SLog(EError, "Path vertex index " SIZE_T_FMT
					" is out of bounds, array size: " SIZE_T_FMT,
					index, m_vertices.size());
		#endif
		return m_vertices[index];
	}

	/// Return an vertex by its index (const version)
	inline const PathVertexPtr &vertex(size_t index) const {
		#if MTS_BD_DEBUG == 1
			if (index >= m_vertices.size())
				SLog(EError, "Path vertex index " SIZE_T_FMT
					" is out of bounds, array size: " SIZE_T_FMT,
					index, m_vertices.size());
		#endif
		return m_vertices[index];
	}

	/// Return an vertex by its index (or NULL if out of bounds)
	inline PathVertexPtr vertexOrNull(size_t index) {
		if (index >= m_vertices.size())
			return NULL;
		return m_vertices[index];
	}

	/// Return an vertex by its index (or NULL if out of bounds, const version)
	inline PathVertexPtr vertexOrNull(size_t index) const {
		if (index >= m_vertices.size())
			return NULL;
		return m_vertices[index];
	}

	/// Return an edge by its index
	inline PathEdgePtr &edge(size_t index) {
		#if MTS_BD_DEBUG == 1
			if (index >= m_edges.size())
				SLog(EError, "Path edge index " SIZE_T_FMT
					" is out of bounds, array size: " SIZE_T_FMT,
					index, m_edges.size());
		#endif
		return m_edges[index];
	}

	/// Return an edge by its index (const version)
	inline const PathEdgePtr &edge(size_t index) const {
		#if MTS_BD_DEBUG == 1
			if (index >= m_edges.size())
				SLog(EError, "Path edge index " SIZE_T_FMT
					" is out of bounds, array size: " SIZE_T_FMT,
					index, m_edges.size());
		#endif
		return m_edges[index];
	}

	/// Return an edge by its index (or \c NULL if out of bounds)
	inline PathEdgePtr edgeOrNull(size_t index) {
		if (index >= m_edges.size())
			return NULL;
		return m_edges[index];
	}

	/// Return an edge by its index (or \c NULL if out of bounds, const version)
	inline PathEdgePtr edgeOrNull(size_t index) const {
		if (index >= m_edges.size())
			return NULL;
		return m_edges[index];
	}


	//! @}
    /* ==================================================================== */

	/* ==================================================================== */
	//! @{ \name                    Miscellaneous
 	/* ==================================================================== */

	/**
	 * \brief Return the spectrally varying path weight that corresponds to
	 * a prefix of length \c l and a suffix of length \c m.
	 *
	 * This operation is used by path mutation strategies that operate
	 * by cutting out a pice of a path and replacing it with something else.
	 * In that case, it is necessary know about the effects of vertices
	 * and edges \a outside of the modified range.
	 */
	inline Spectrum getPrefixSuffixWeight(int l, int m) const {
		Spectrum weight(1.0f);

		for (int s=0; s<l; ++s)
			weight *= m_vertices[s]->weight[EImportance]
				* m_edges[s]->weight[EImportance];

		for (int t=length(); t>m; --t)
			weight *= m_vertices[t]->weight[ERadiance]
				* m_edges[t-1]->weight[ERadiance];

		return weight;
	}

	/**
	 * \brief Determine whether another path of the same length
	 * matches the configuration of this path
	 *
	 * More specifically, this function verifies that the length matches,
	 * that all vertices are of the same type, and that the paths have
	 * identical values of \ref isConnectable() for each vertex.
	 */
	inline bool matchesConfiguration(const Path &p) const {
		if (p.length() != length())
			return false;

		for (size_t i=0; i<m_vertices.size(); ++i) {
			if (m_vertices[i]->type != p.vertex(i)->type ||
				m_vertices[i]->isConnectable() != p.vertex(i)->isConnectable())
				return false;
		}

		return true;
	}

	/**
	 * \brief Return the relative spectrally varying path weight associated
	 * with this path
	 *
	 * This function computes the product of all the importance weights
	 * (see the \ref weight field) along the path, and the result is
	 * normalized so that it has luminance 1. This quantity is required
	 * by the sample splatting implementation in Veach-MLT.
	 */
	inline Spectrum getRelativeWeight() const {
		Spectrum weight(1.0f);
		int k = length();

		for (int s=0; s<k-1; ++s)
			weight *= m_vertices[s]->weight[EImportance]
				* m_edges[s]->weight[EImportance];

		weight = weight
			* m_vertices[k]->weight[ERadiance]
			* m_vertices[k-1]->weight[ERadiance]
			* m_edges[k-1]->weight[ERadiance];

		Float lum = weight.getLuminance();
		return lum != 0.0f ? (weight / lum) : Spectrum(0.0f);
	}

	/**
	 * \brief Compute the multiple importance sampling weight of the <tt>(s,t)</tt>
	 * sampling strategy in BDPT.
	 *
	 * This implementation uses the power heuristic with exponent 2 and
	 * repeatedly evaluates equation (10.9) from Eric Veach's PhD thesis to
	 * compute the weight in an efficient and numerically stable manner.
	 *
	 * The function completely ignores the effects of russian roulette, since
	 * this allows for a more efficient implementation. The resulting estimator
	 * is still unbiased despite this apparent inaccuracy.
	 *
	 * \param scene
	 *    Pointer to the underlying scene
	 * \param emitterSubpath
	 *    Reference to the emitter subpath
	 * \param connectionEdge
	 *    Pointer to an edge data structure associated with the
	 *    transport between <tt>emitterSubpath[s]</tt> and
	 *    <tt>sensorSubpath[t]</tt>.
	 * \param sensorSubpath
	 *    Reference to the sensor subpath
	 * \param s
	 *    Number of steps to take along the emitter subpath
	 * \param t
	 *    Number of steps to take along the sensor subpath
	 * \param direct
	 *    When the parameter \c direct is set to \c true, the implementation
	 *    accounts for the fact that specialized direct sampling strategies
	 *    are used for paths with <tt>s==1</tt> and <tt>t==1</tt>.
	 * \param lightImage
	 *    Denotes whether or not rendering strategies that require a 'light image'
	 *    (specifically, those with <tt>t==0</tt> or <tt>t==1</tt>) are included
	 *    in the rendering process.
	 */
	static Float miWeight(const Scene *scene,
			const Path &emitterSubpath,
			const PathEdge *connectionEdge,
			const Path &sensorSubpath, int s, int t,
			bool direct, bool lightImage);

	/**
	 * \brief Collapse a path into an entire edge that summarizes the aggregate
	 * transport and sampling densities
	 *
	 * This is only allowed when the path only consists of \ref BSDF::ENull
	 * surface scattering events (e.g. index-matched medium transitions)
	 */
	void collapseTo(PathEdge &edge) const;

	/// Reverse a path
	void reverse();

	/// Swap the two endpoint vertices of a path with the provided values
	inline void swapEndpoints(PathVertexPtr &supernode, PathEdgePtr &edge, PathVertexPtr &sample) {
		std::swap(supernode, m_vertices[0]);
		std::swap(sample, m_vertices[1]);
		std::swap(edge, m_edges[0]);
	}

	/// \brief Append a vertex to this path
	inline void append(PathVertex *vertex) { m_vertices.push_back(vertex); }

	/// \brief Append an edge to this path
	inline void append(PathEdge *edge) { m_edges.push_back(edge); }

	/**
	 * \brief Append an edge and a vertex to this path
	 *
	 * \param edge
	 *     And edge from <tt>vertex(vertexCount()-1)</tt> to \c vertex.
	 * \param vertex
	 *     A vertex to be appended at the end of the path
	 */
	inline void append(PathEdge *edge, PathVertex *vertex) {
		m_edges.push_back(edge);
		m_vertices.push_back(vertex);
	}

	/// Append an entire path to this path
	void append(const Path &path);

	/**
	 * \brief Append the vertex range <tt>[start, end)</tt> of \c path
	 * (and all intermediate edges) to the current path.
	 *
	 * \param path
	 *     Source of vertices and edges to be added
	 * \param l
	 *     Specifies the start of the range <tt>[start, end)</tt>.
	 * \param m
	 *     Specifies the end of the range <tt>[start, end)</tt>.
	 * \param reverse
	 *     Should the vertices and edges be added in \a reverse order?
	 */
	void append(const Path &path, size_t start, size_t end, bool reverse = false);

	/// Clear the path and release all elements to the memory pool
	void release(MemoryPool &pool);

	/// Release a certain subpath [start, end) to the memory pool
	void release(size_t start, size_t end, MemoryPool &pool);

	/// Return the sample position associated with the path
	inline const Point2 &getSamplePosition() const {
		return m_vertices[length()-1]->getSamplePosition();
	}

	/// Clear the path
	inline void clear() {
		m_vertices.clear();
		m_edges.clear();
	}

	/// Compare this path against another path
	bool operator==(const Path &path) const;

	/// Compare this path against another path
	inline bool operator!=(const Path &path) const {
		return !operator==(path);
	}

	/// Create a deep copy of this path
	void clone(Path &target, MemoryPool &pool) const;

	/// Return a string representation of the path
	std::string toString() const;

	/// Return a basic string summary of the path
	std::string summarize() const;

	//! @}
    /* ==================================================================== */
private:
	std::vector<PathVertexPtr> m_vertices;
	std::vector<PathEdgePtr>   m_edges;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_PATH_H_ */
