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
#if !defined(__MITSUBA_CORE_RFILTER_H_)
#define __MITSUBA_CORE_RFILTER_H_

#include <mitsuba/core/cobject.h>

MTS_NAMESPACE_BEGIN

/// Reconstruction filters will be tabulated at this resolution
#define MTS_FILTER_RESOLUTION 31

/**
 * \brief Generic interface to separable image reconstruction filters
 *
 * When resampling bitmaps or adding radiance-valued samples to a rendering in
 * progress, Mitsuba first convolves them with a so-called image reconstruction
 * filter. Various kinds are implemented as subclasses of this interface.
 *
 * Because image filters are generally too expensive to evaluate for
 * each sample, the implementation of this class internally precomputes
 * an discrete representation (resolution given by \ref MTS_FILTER_RESOLUTION)
 *
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_CORE ReconstructionFilter : public ConfigurableObject {
public:
    /**
     * \brief When resampling data to a different resolution using
     * \ref Resampler::resample(), this enumeration specifies how lookups
     * <em>outside</em> of the input domain are handled.
     *
     * \see Resampler
     */
    enum EBoundaryCondition {
        /// Clamp to the outermost sample position
        EClamp = 0,
        /// Assume that the input repeats in a periodic fashion
        ERepeat,
        /// Assume that the input is mirrored along the boundary
        EMirror,
        /// Assume that the input function is zero outside of the defined domain
        EZero,
        /// Assume that the input function is equal to one outside of the defined domain
        EOne
    };

    /// Return the filter's width
    inline Float getRadius() const { return m_radius; }

    /// Return the block border size required when rendering with this filter
    inline int getBorderSize() const { return m_borderSize; }

    /// Evaluate the filter function
    virtual Float eval(Float x) const = 0;

    /// Perform a lookup into the discretized version
    inline Float evalDiscretized(Float x) const { return m_values[
        std::min((int) std::abs(x * m_scaleFactor), MTS_FILTER_RESOLUTION)]; }

    /// Serialize the filter to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    /** \brief Configure the object (called \a once after construction) */
    void configure();

    MTS_DECLARE_CLASS()
protected:
    /// Create a new reconstruction filter
    ReconstructionFilter(const Properties &props);

    /// Unserialize a filter
    ReconstructionFilter(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~ReconstructionFilter();
protected:
    Float m_radius, m_scaleFactor;
    Float m_values[MTS_FILTER_RESOLUTION+1];
    int m_borderSize;
};

/**
 * \brief Utility class for efficiently resampling discrete datasets to different resolutions
 * \tparam Scalar
 *      Denotes the underlying floating point data type (i.e. <tt>half</tt>, <tt>float</tt>,
 *      or <tt>double</tt>)
 */
template <typename Scalar> struct Resampler {
    /**
     * \brief Create a new Resampler object that transforms between the specified resolutions
     *
     * This constructor precomputes all information needed to efficiently perform the
     * desired resampling operation. For that reason, it is most efficient if it can
     * be used over and over again (e.g. to resample the equal-sized rows of a bitmap)
     *
     * \param sourceRes
     *      Source resolution
     * \param targetRes
     *      Desired target resolution
     * \param bc
     *      Boundary conditions that should be observed when looking up samples
     *      outside of the defined input domain.
     */
    Resampler(const ReconstructionFilter *rfilter, ReconstructionFilter::EBoundaryCondition bc,
            int sourceRes, int targetRes) : m_bc(bc), m_sourceRes(sourceRes), m_targetRes(targetRes),
            m_start(NULL), m_weights(NULL) {
        SAssert(sourceRes > 0 && targetRes > 0);
        Float filterRadius = rfilter->getRadius(), scale = (Float)1.0, invScale = (Float)1.0;

        /* Low-pass filter: scale reconstruction filters when downsampling */
        if (targetRes < sourceRes) {
            scale = (Float) sourceRes / (Float) targetRes;
            invScale = 1 / scale;
            filterRadius *= scale;
        }

        m_taps = math::ceilToInt(filterRadius * 2);
        if (sourceRes == targetRes && (m_taps % 2) != 1)
            --m_taps;
        m_halfTaps = m_taps / 2;

        if (sourceRes != targetRes) { /* Resampling mode */
            m_start = new int[targetRes];
            m_weights = new Scalar[m_taps * targetRes];
            m_fastStart = 0;
            m_fastEnd = m_targetRes;

            for (int i=0; i<targetRes; i++) {
                /* Compute the fractional coordinates of the new sample i in the original coordinates */
                Float center = (i + (Float) 0.5f) / targetRes * sourceRes;

                /* Determine the index of the first original sample that might contribute */
                m_start[i] = math::floorToInt(center - filterRadius + (Float) 0.5f);

                /* Determine the size of center region, on which to run fast non condition-aware code */
                if (m_start[i] < 0)
                    m_fastStart = std::max(m_fastStart, i + 1);
                else if (m_start[i] + m_taps - 1 >= m_sourceRes)
                    m_fastEnd = std::min(m_fastEnd, i - 1);

                Float sum = 0;
                for (int j=0; j<m_taps; j++) {
                    /* Compute the the position where the filter should be evaluated */
                    Float pos = m_start[i] + j + (Float) 0.5f - center;

                    /* Perform the evaluation and record the weight */
                    Float weight = rfilter->eval(pos * invScale);
                    m_weights[i * m_taps + j] = (Scalar) weight;
                    sum += weight;
                }

                /* Normalize the contribution of each sample */
                Float normalization = 1.0f / sum;
                for (int j=0; j<m_taps; j++) {
                    Scalar &value = m_weights[i * m_taps + j];
                    value = (Scalar) ((Float) value * normalization);
                }
            }
        } else { /* Filtering mode */
            m_weights = new Scalar[m_taps];
            Float sum = 0;
            for (int i=0; i<m_taps; i++) {
                Scalar weight = (Scalar) rfilter->eval((Float) (i-m_halfTaps));
                m_weights[i] = weight;
                sum += (Float) weight;
            }
            Float normalization = 1.0f / sum;
            for (int i=0; i<m_taps; i++) {
                Scalar &value = m_weights[i];
                value = (Scalar) ((Float) value * normalization);
            }
            m_fastStart = std::min(m_halfTaps, m_targetRes-1);
            m_fastEnd = std::max(m_targetRes-m_halfTaps-1, 0);
        }

        /* Don't have overlapping fast start/end intervals when
           the target image is very small compared to the source image */
        m_fastStart = std::min(m_fastStart, m_fastEnd);
    }

    /// Release all memory
    ~Resampler() {
        if (m_start)
            delete[] m_start;
        if (m_weights)
            delete[] m_weights;
    }

    /**
     * \brief Resample a multi-channel array and clamp the results
     * to a specified valid range
     *
     * This function is preferred if too large positive/negative values
     * due to ringing are unacceptable.
     *
     * \param source
     *     Source array of samples
     * \param target
     *     Target array of samples
     * \param sourceStride
     *     Stride of samples in the source array. A value
     *     of '1' implies that they are densely packed.
     * \param targetStride
     *     Stride of samples in the source array. A value
     *     of '1' implies that they are densely packed.
     * \param channels
     *     Number of channels to be resampled
     * \param min
     *     Minumum sample value after resampling
     * \param max
     *     Maximum sample value after resampling
     */
    void resampleAndClamp(const Scalar *source, size_t sourceStride,
            Scalar *target, size_t targetStride, int channels,
            Scalar min = (Scalar) 0, Scalar max = (Scalar) 1) {
        const int taps = m_taps, halfTaps = m_halfTaps;
        targetStride = channels * (targetStride - 1);
        sourceStride *= channels;

        if (m_start) {
            /* Resample the left border region, while accounting for the boundary conditions */
            for (int i=0; i<m_fastStart; ++i) {
                int start = m_start[i];

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[i * taps + j];
                    *target++ = std::min(max, std::max(min, result));
                }

                target += targetStride;
            }

            /* Use a faster, vectorizable loop for resampling the main portion */
            for (int i=m_fastStart; i<m_fastEnd; ++i) {
                int start = m_start[i];

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += source[sourceStride * (start + j) + ch] * m_weights[i * taps + j];
                    *target++ = std::min(max, std::max(min, result));
                }

                target += targetStride;
            }

            /* Resample the right border region, while accounting for the boundary conditions */
            for (int i=m_fastEnd; i<m_targetRes; ++i) {
                int start = m_start[i];

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[i * taps + j];
                    *target++ = std::min(max, std::max(min, result));
                }

                target += targetStride;
            }
        } else {
            /* Filter the left border region, while accounting for the boundary conditions */
            for (int i=0; i<m_fastStart; ++i) {
                int start = i - halfTaps;

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[j];
                    *target++ = std::min(max, std::max(min, result));
                }

                target += targetStride;
            }

            /* Use a faster, vectorizable loop for filtering the main portion */
            for (int i=m_fastStart; i<m_fastEnd; ++i) {
                int start = i - halfTaps;

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += source[sourceStride * (start + j) + ch] * m_weights[j];
                    *target++ = std::min(max, std::max(min, result));
                }

                target += targetStride;
            }

            /* Filter the right border region, while accounting for the boundary conditions */
            for (int i=m_fastEnd; i<m_targetRes; ++i) {
                int start = i - halfTaps;

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[j];
                    *target++ = std::min(max, std::max(min, result));
                }

                target += targetStride;
            }
        }
    }

    /**
     * \brief Resample a multi-channel array
     *
     * \param source
     *     Source array of samples
     * \param target
     *     Target array of samples
     * \param sourceStride
     *     Stride of samples in the source array. A value
     *     of '1' implies that they are densely packed.
     * \param targetStride
     *     Stride of samples in the source array. A value
     *     of '1' implies that they are densely packed.
     * \param channels
     *     Number of channels to be resampled
     */
    void resample(const Scalar *source, size_t sourceStride,
            Scalar *target, size_t targetStride, int channels) {
        const int taps = m_taps, halfTaps = m_halfTaps;

        targetStride = channels * (targetStride - 1);
        sourceStride *= channels;

        if (m_start) {
            /* Resample the left border region, while accounting for the boundary conditions */
            for (int i=0; i<m_fastStart; ++i) {
                int start = m_start[i];

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[i * taps + j];
                    *target++ = result;
                }

                target += targetStride;
            }

            /* Use a faster, vectorizable loop for resampling the main portion */
            for (int i=m_fastStart; i<m_fastEnd; ++i) {
                int start = m_start[i];

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += source[sourceStride * (start + j) + ch] * m_weights[i * taps + j];
                    *target++ = result;
                }

                target += targetStride;
            }

            /* Resample the right border region, while accounting for the boundary conditions */
            for (int i=m_fastEnd; i<m_targetRes; ++i) {
                int start = m_start[i];

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[i * taps + j];
                    *target++ = result;
                }

                target += targetStride;
            }
        } else {
            /* Filter the left border region, while accounting for the boundary conditions */
            for (int i=0; i<m_fastStart; ++i) {
                int start = i - halfTaps;

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[j];
                    *target++ = result;
                }

                target += targetStride;
            }

            /* Use a faster, vectorizable loop for filtering the main portion */
            for (int i=m_fastStart; i<m_fastEnd; ++i) {
                int start = i - halfTaps;

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += source[sourceStride * (start + j) + ch] * m_weights[j];
                    *target++ = result;
                }

                target += targetStride;
            }

            /* Filter the right border region, while accounting for the boundary conditions */
            for (int i=m_fastEnd; i<m_targetRes; ++i) {
                int start = i - halfTaps;

                for (int ch=0; ch<channels; ++ch) {
                    Scalar result = 0;
                    for (int j=0; j<taps; ++j)
                        result += lookup(source, start + j, sourceStride, ch) * m_weights[j];
                    *target++ = result;
                }

                target += targetStride;
            }
        }
    }

private:
    FINLINE Scalar lookup(const Scalar *source, int pos, size_t stride, int offset) const {
        if (EXPECT_NOT_TAKEN(pos < 0 || pos >= m_sourceRes)) {
            switch (m_bc) {
                case ReconstructionFilter::EClamp:
                    pos = math::clamp(pos, 0, m_sourceRes - 1);
                    break;
                case ReconstructionFilter::ERepeat:
                    pos = math::modulo(pos, m_sourceRes);
                    break;
                case ReconstructionFilter::EMirror:
                    pos = math::modulo(pos, 2*m_sourceRes);
                    if (pos >= m_sourceRes)
                        pos = 2*m_sourceRes - pos - 1;
                    break;
                case ReconstructionFilter::EZero:
                    return (Scalar) 0;
                case ReconstructionFilter::EOne:
                    return (Scalar) 1;
            }
        }
        return source[stride * pos + offset];
    }

private:
    ReconstructionFilter::EBoundaryCondition m_bc;
    int m_sourceRes;
    int m_targetRes;
    int *m_start;
    Scalar *m_weights;
    int m_fastStart, m_fastEnd;
    int m_taps, m_halfTaps;
};

extern MTS_EXPORT_CORE std::ostream &operator<<(std::ostream &os, const ReconstructionFilter::EBoundaryCondition &value);

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_RFILTER_H_ */
