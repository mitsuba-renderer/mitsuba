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

#if !defined(__SAMPLER_H)
#define __SAMPLER_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Base class of all sample generators. 
 *
 * For each sample in a pixel, a sample generator produces a (hypothetical)
 * point in the infinite dimensional random number cube. A rendering 
 * algorithm can then request subsequent 1D or 2D components of this point 
 * using the \ref next1D() and \ref next2D() functions. Some implementations
 * make certain guarantees about the stratification of the first n components
 * with respect to the other points that are sampled within a pixel. (the 
 * low-discrepancy and stratified samplers do this for instance). 
 *
 * The general interaction between a sampler and a rendering algorithm is as 
 * follows: Before beginning to render a pixel, the rendering algorithm calls 
 * \ref generate(). The first pixel sample can now be computed, after which
 * \ref advance() needs to be invoked. This repeats until all pixel samples have
 * been generated. Note that some implementations need to be configured for a 
 * certain number of pixel samples, and exceeding these will lead to an 
 * exception being thrown. While computing a pixel sample, the rendering 
 * algorithm usually requests batches of (pseudo-) random numbers using 
 * the \ref next1D(), \ref next2D(), \ref next1DArray() and 
 * \ref next2DArray() functions.
 *
 * The difference between calling \ref next1D(), \ref next2D() a number of
 * times versus using the array variants \ref next1DArray() and 
 * \ref next2DArray() is that the latter can provide stratification 
 * not only with respect to random numbers obtained within another pixel
 * sample, but also within the array itself. This is useful e.g. in the
 * direct illumination strategy, which spawns several rays after hitting
 * a surface when computing a pixel sample. Since this is done over and 
 * over again for each one of the pixel samples, it makes sense to 
 * stratify over all of the rays that are ultimately generated, and the
 * \ref next1DArray() and \ref next2DArray() methods allow to do this.
 * See the file \c direct.cpp for an example.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Sampler : public ConfigurableObject {
public:
	/**
	 * Create a clone of this sampler. The clone is allowed to be different
	 * to some extent, e.g. a pseudorandom generator should be based on a
	 * different random seed compared to the original. All other parameters,
	 * are copied exactly.
	 */
	virtual ref<Sampler> clone() = 0;

	/**
	 * Generate new samples - called initially and every time the generated 
	 * samples have been exhausted. When used in conjunction with a 
	 * SampleIntegrator, this will be called before starting to render 
	 * each pixel.
	 */
	virtual void generate();

	/// Advance to the next sample
	virtual void advance();

	/// Manually set the current sample index
	virtual void setSampleIndex(size_t sampleIndex);

	/// Retrieve the next component value from the current sample
	virtual Float next1D() = 0;

	/// Retrieve the next two component values from the current sample
	virtual Point2 next2D() = 0;

	/**
	 * Retrieve the next 2D array of values from the current sample.
	 * Note that this is different from just calling <tt>next2D()</tt>
	 * repeatedly - this function will generally return a set of 2D vectors,
	 * which are not only well-laid out over all samples at the current pixel,
	 * but also with respect to each other. Note that this 2D array has to be
	 * requested initially using <tt>request2DArray</tt> and later, they have 
	 * to be retrieved in the same same order and size configuration as the 
	 * requests. An exception is thrown when a mismatch is detected.
	 *
	 * This function is useful to support things such as a direct illumination
	 * rendering technique with "n" pixel samples and "m" shading samples,
	 * while ensuring that the "n*m" sampled positions on an area light source
	 * are all well-stratified with respect to each other.
	 */
	Point2 *next2DArray(unsigned int size);

	/// Same as above, but 1D
	Float *next1DArray(unsigned int size);

	/**
	 * Request that a 2D array will be made available for 
	 * later consumption by next2DArray(). This must be called
	 * before generate(). See 'next2DArray' for a description
	 * of this feature.
	 */
	virtual void request2DArray(unsigned int size);

	/// Same as above, but 1D
	virtual void request1DArray(unsigned int size);

	/**
	 * Return an uniformly distributed number on [0, 1).
	 * Throws an error when the underlying implementation is 
	 * fully deterministic (e.g. QMC).
	 */
	virtual Float independent1D() = 0;

	/**
	 * Return an uniformly distributed 2D vector on [0, 1)x[0, 1).
	 * Throws an error when the underlying implementation is 
	 * fully deterministic (e.g. QMC).
	 */
	virtual Point2 independent2D() = 0;

	/// Return total number of samples
	inline size_t getSampleCount() const { return m_sampleCount; }
	
	/// Return the current sample index
	inline size_t getSampleIndex() const { return m_sampleIndex; }

	/// Serialize this sampler to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return the properties of this sampler
	inline const Properties &getProperties() const { return m_properties; }

	/// Return a string description
	virtual std::string toString() const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Construct a new sampler
	Sampler(const Properties &props);

	/// Unserialize a sampler
	Sampler(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Sampler();
protected:
	size_t m_sampleCount;
	size_t m_sampleIndex;
	std::vector<unsigned int> m_req1D, m_req2D;
	std::vector<Float *> m_sampleArrays1D;
	std::vector<Point2 *> m_sampleArrays2D;
	int m_sampleDepth1DArray, m_sampleDepth2DArray;
	Properties m_properties;
};

MTS_NAMESPACE_END

#endif /* __SAMPLER_H */
