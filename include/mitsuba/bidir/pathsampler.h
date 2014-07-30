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
#if !defined(__MITSUBA_BIDIR_PATHSAMPLER_H_)
#define __MITSUBA_BIDIR_PATHSAMPLER_H_

#include <mitsuba/bidir/path.h>
#include <boost/function.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Implements a sampling strategy that is able to produce paths using
 * bidirectional path tracing or unidirectional volumetric path tracing.
 *
 * This versatile class does the heavy lifting under the hood of Mitsuba's
 * PSSMLT implementation. It is also used to provide the Veach-MLT
 * implementation with a luminance estimate and seed paths.
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
class MTS_EXPORT_BIDIR PathSampler : public Object {
public:
	/**
	 * \brief Callback type for use with \ref samplePaths()
	 *
	 * The arguments are (s, t, weight, path) where \c s denotes the number
	 * of steps from the emitter, \c is the number of steps from the sensor,
	 * and \c weight contains the importance weight associated with the sample.
	 */
	typedef boost::function<void (int, int, Float, Path &)> PathCallback;

	/// Specifies the sampling algorithm that is internally used
	enum ETechnique {
		/// Bidirectional path tracing
		EBidirectional,

		/// Unidirectional path tracing (via the 'volpath' plugin)
		EUnidirectional
	};

	/**
	 * Construct a new path sampler
	 *
	 * \param technique
	 *     What path generation technique should be used (unidirectional
	 *     or bidirectional path tracing?)
	 *
	 * \param scene
	 *     \ref A pointer to the underlying scene
	 *
	 * \param emitterSampler
	 *     Sample generator that should be used for the random walk
	 *     from the emitter direction
	 *
	 * \param sensorSampler
	 *     Sample generator that should be used for the random walk
	 *     from the sensor direction
	 *
	 * \param directSampler
	 *     Sample generator that should be used for direct sampling
	 *     strategies (or \c NULL, when \c sampleDirect=\c false)
	 *
	 * \param maxDepth
	 *     Maximum path depth to be visualized (-1==infinite)
	 *
	 * \param rrDepth
	 *     Depth to begin using russian roulette
	 *
	 * \param excludeDirectIllum
	 *     If set to true, the direct illumination component will
	 *     be ignored. Note that this parameter is unrelated
	 *     to the next one (\a sampleDirect) although they are
	 *     named similarly.
	 *
	 * \param sampleDirect
	 *     When this parameter is set to true, specialized direct
	 *     sampling strategies are used for s=1 and t=1 paths.
	 *
	 * \param lightImage
	 *    Denotes whether or not rendering strategies that require a 'light image'
	 *    (specifically, those with <tt>t==0</tt> or <tt>t==1</tt>) are included
	 *    in the rendering process.
	 */
	PathSampler(ETechnique technique, const Scene *scene, Sampler *emitterSampler,
		Sampler *sensorSampler, Sampler *directSampler, int maxDepth, int rrDepth,
		bool excludeDirectIllum, bool sampleDirect, bool lightImage = true);

	/**
	 * \brief Generate a sample using the configured sampling strategy
	 *
	 * The result is stored as a series of screen-space splats (pixel position
	 * and spectral value pairs) within the parameter \c list. These can be
	 * used to implement algorithms like Bidirectional Path Tracing or Primary
	 * Sample Space MLT.
	 *
	 * \param offset
	 *    Specifies the desired integer pixel position of the sample. The special
	 *    value <tt>Point2i(-1)</tt> results in uniform sampling in screen space.
	 *
	 * \param list
	 *    Output parameter that will receive a list of splats
	 */
	void sampleSplats(const Point2i &offset, SplatList &list);

	/**
	 * \brief Sample a series of paths and invoke the specified callback
	 * function for each one.
	 *
	 * This function is similar to \ref sampleSplats(), but instead of
	 * returning only the contribution of the samples paths in the form of
	 * screen-space "splats", it returns the actual paths by invoking a
	 * specified callback function multiple times.
	 *
	 * This function is currently only implemented for the bidirectional
	 * sampling strategy -- i.e. it cannot be used with the unidirectional
	 * path tracer.
	 *
	 * \param offset
	 *    Specifies the desired integer pixel position of the sample. The special
	 *    value <tt>Point2i(-1)</tt> results in uniform sampling in screen space.
	 *
	 * \param pathCallback
	 *    A callback function that will be invoked once for each
	 *    path generated by the BDPT sampling strategy. The first argument
	 *    specifies the path importance weight.
	 */
	void samplePaths(const Point2i &offset, PathCallback &callback);

	/**
	 * \brief Generates a sequence of seeds that are suitable for
	 * starting a MLT Markov Chain
	 *
	 * This function additionally computes the average luminance
	 * over the image plane.
	 *
	 * \param sampleCount
	 *     The number of luminance samples that will be taken
	 * \param seedCount
	 *     The desired number of MLT seeds (must be > \c sampleCount)
	 * \param fineGrained
	 *     This parameter only matters when the technique is set to
	 *     \ref EBidirectional. It specifies whether to generate \ref PathSeed
	 *     records at the granularity of entire sensor/emitter subpaths or at
	 *     the granularity of their constituent sampling strategies.
	 * \param seeds
	 *     A vector of resulting MLT seeds
	 * \return The average luminance over the image plane
	 */
	Float generateSeeds(size_t sampleCount, size_t seedCount,
			bool fineGrained, const Bitmap *importanceMap,
			std::vector<PathSeed> &seeds);

	/**
	 * \brief Compute the average luminance over the image plane
	 * \param sampleCount
	 *     The number of luminance samples that will be taken
	 */
	Float computeAverageLuminance(size_t sampleCount);

	/**
	 * \brief Reconstruct a path from a \ref PathSeed record
	 *
	 * Given a \ref PathSeed data structure, this function rewinds
	 * the random number stream of the underlying \ref ReplayableSampler
	 * to the indicated position and recreates the associated path.
	 */
	void reconstructPath(const PathSeed &seed, 
		const Bitmap *importanceMap, Path &result);

	/// Return the underlying memory pool
	inline MemoryPool &getMemoryPool() { return m_pool; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~PathSampler();
protected:
	ETechnique m_technique;
	ref<const Scene> m_scene;
	ref<SamplingIntegrator> m_integrator;
	ref<Sampler> m_emitterSampler;
	ref<Sampler> m_sensorSampler;
	ref<Sampler> m_directSampler;
	int m_maxDepth;
	int m_rrDepth;
	bool m_excludeDirectIllum;
	bool m_sampleDirect;
	bool m_lightImage;
	int m_emitterDepth, m_sensorDepth;
	Path m_emitterSubpath, m_sensorSubpath;
	Path m_connectionSubpath, m_fullPath;
	MemoryPool m_pool;
};

/**
 * \brief Stores information required to re-create a seed path (e.g. for MLT)
 *
 * This class makes it possible to transmit a path over the network or store
 * it locally, while requiring very little storage to do so. This is done by
 * describing a path using an index into a random number stream, which allows
 * to generate it cheaply when needed.
 */
struct PathSeed {
	size_t sampleIndex; ///< Index into a rewindable random number stream
	Float luminance;    ///< Luminance value of the path (for sanity checks)
	int s;              ///< Number of steps from the luminaire
	int t;              ///< Number of steps from the eye

	inline PathSeed() { }

	inline PathSeed(size_t sampleIndex, Float luminance, int s = 0, int t = 0)
		: sampleIndex(sampleIndex), luminance(luminance), s(s), t(t) { }

	inline PathSeed(Stream *stream) {
		sampleIndex = stream->readSize();
		luminance = stream->readFloat();
		s = stream->readInt();
		t = stream->readInt();
	}

	void serialize(Stream *stream) const {
		stream->writeSize(sampleIndex);
		stream->writeFloat(luminance);
		stream->writeInt(s);
		stream->writeInt(t);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PathSeed[" << endl
			<< "  sampleIndex = " << sampleIndex << "," << endl
			<< "  luminance = " << luminance << "," << endl
			<< "  s = " << s << "," << endl
			<< "  t = " << t << endl
			<< "]";
		return oss.str();
	}
};

/**
 * MLT work unit -- wraps a \ref PathSeed into a
 * \ref WorkUnit instance.
 */
class SeedWorkUnit : public WorkUnit {
public:
	inline void set(const WorkUnit *wu) {
		m_seed = static_cast<const SeedWorkUnit *>(wu)->m_seed;
		m_timeout = static_cast<const SeedWorkUnit *>(wu)->m_timeout;
	}

	inline const PathSeed &getSeed() const {
		return m_seed;
	}

	inline void setSeed(const PathSeed &seed) {
		m_seed = seed;
	}

	inline int getTimeout() const {
		return m_timeout;
	}

	inline void setTimeout(int timeout) {
		m_timeout = timeout;
	}

	inline void load(Stream *stream) {
		m_seed = PathSeed(stream);
		m_timeout = stream->readInt();
	}

	inline void save(Stream *stream) const {
		m_seed.serialize(stream);
		stream->writeInt(m_timeout);
	}

	inline std::string toString() const {
		return "SeedWorkUnit[]";
	}

	MTS_DECLARE_CLASS()
private:
	PathSeed m_seed;
	int m_timeout;
};

/**
 * \brief List storage for the image-space contributions ("splats") of a
 * path sample generated using \ref PathSampler.
 */
struct MTS_EXPORT_BIDIR SplatList {
	/// Represents a screen-space splat produced by a path sampling technique
	typedef std::pair<Point2, Spectrum> Splat;

	/// A series of splats associated with the current sample
	std::vector<Splat> splats;
	/// Combined luminance of all splats in this sample
	Float luminance;
	/// Total number of samples in the splat list
	int nSamples;

	inline SplatList() : luminance(0.0f), nSamples(0) { }

	/// Appends a splat entry to the list
	inline void append(const Point2 &samplePos, const Spectrum &value) {
		splats.push_back(std::make_pair(samplePos, value));
		luminance += value.getLuminance();
		++nSamples;
	}

	/// Increases the contribution of an existing splat
	inline void accum(size_t i, const Spectrum &value) {
		splats[i].second += value;
		luminance += value.getLuminance();
		++nSamples;
	}

	/// Returns the number of contributions
	inline size_t size() const {
		return splats.size();
	}

	/// Clear the splat list
	inline void clear() {
		luminance = 0;
		nSamples = 0;
		splats.clear();
	}

	/// Return the position associated with a splat in the list
	inline const Point2 &getPosition(size_t i) const { return splats[i].first; }

	/// Return the spectral contribution associated with a splat in the list
	inline const Spectrum &getValue(size_t i) const { return splats[i].second; }

	/**
	 * \brief Normalize the splat list
	 *
	 * This function divides all splats so that they have unit
	 * luminance (though it leaves the \c luminance field untouched).
	 * When given an optional importance map in 2-stage MLT approaches,
	 * it divides the splat values by the associated importance
	 * map values
	 */
	void normalize(const Bitmap *importanceMap = NULL);

	/// Return a string representation
	std::string toString() const;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_PATHSAMPLER_H_ */
