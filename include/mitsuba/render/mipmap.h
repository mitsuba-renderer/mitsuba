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
#if !defined(__MITSUBA_RENDER_MIPMAP_H_)
#define __MITSUBA_RENDER_MIPMAP_H_

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/barray.h>
#include <mitsuba/core/mmap.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/statistics.h>
#include <boost/filesystem/fstream.hpp>

MTS_NAMESPACE_BEGIN

/// Use a blocked array to store MIP map data (slightly faster)
#define MTS_MIPMAP_BLOCKED 1

/// Look-up table size for a tabulated Gaussian filter
#define MTS_MIPMAP_LUT_SIZE 64

/// MIP map cache file version
#define MTS_MIPMAP_CACHE_VERSION 0x01

/// Make sure that the actual cache contents start on a cache line
#define MTS_MIPMAP_CACHE_ALIGNMENT 64

/* Some statistics counters */
namespace stats {
    extern MTS_EXPORT_RENDER StatsCounter avgEWASamples;
    extern MTS_EXPORT_RENDER StatsCounter clampedAnisotropy;
    extern MTS_EXPORT_RENDER StatsCounter mipStorage;
    extern MTS_EXPORT_RENDER StatsCounter filteredLookups;
};

/// Specifies the desired antialiasing filter
enum EMIPFilterType {
    /// No filtering, nearest neighbor lookups
    ENearest = 0,
    /// No filtering, only bilinear interpolation
    EBilinear = 1,
    /// Basic trilinear filtering
    ETrilinear = 2,
    /// Elliptically weighted average
    EEWA = 3
};

/**
 * \brief MIP map class with support for elliptically weighted averages
 *
 * This class stores a precomputed collection of images that provide
 * a hierarchy of resolution levels of an input image. These are used
 * to prefilter texture lookups at render time, reducing aliasing
 * artifacts. The implementation here supports non-power-of two images,
 * different data types and quantizations thereof, as well as the
 * computation of elliptically weighted averages that account for the
 * anisotropy of texture lookups in UV space.
 *
 * Generating good mip maps is costly, and therefore this class provides
 * the means to cache them on disk if desired.
 *
 * \tparam Value
 *    This class can be parameterized to yield MIP map classes for
 *    RGB values, color spectra, or just plain floats. This parameter
 *    specifies the underlying type.
 *
 * \tparam QuantizedValue
 *    If desired, the MIP map can store its contents using a quantized
 *    representation (e.g. half precision). In that case, this type
 *    denotes the associated data type.
 *
 * \ingroup librender
 */
template <typename Value, typename QuantizedValue> class TMIPMap : public Object {
public:
#if MTS_MIPMAP_BLOCKED == 1
    /// Use a blocked array to store MIP map data
    typedef BlockedArray<QuantizedValue> Array2DType;
#else
    /// Use a non-blocked array to store MIP map data
    typedef LinearArray<QuantizedValue> Array2DType;
#endif

    /// Shortcut
    typedef ReconstructionFilter::EBoundaryCondition EBoundaryCondition;

    /**
     * \brief Construct a new MIP map from the given bitmap
     *
     * \param bitmap
     *    An arbitrary input bitmap that will (if necessary) be
     *    transformed into a representation that is compatible
     *    with the desired MIP map pixel format.
     *
     * \ref pixelFormat
     *    A pixel format that is compatible with the
     *    \c Value and \c QuantizedValue types
     *
     * \ref componentFormat
     *    A component format that is compatible with the
     *    \c Value type.
     *
     * \param rfilter
     *    An image reconstruction filter that is used to create
     *    progressively lower-resolution versions of the input image.
     *
     * \param bcu
     *    Specifies how to handle texture lookups outside of the
     *    horizontal range [0, 1]
     *
     * \param bcv
     *    Specifies how to handle texture lookups outside of the
     *    vertical range [0, 1]
     *
     * \param filterType
     *    Denotes the desired filter type to use for lookups;
     *    the default is to use elliptically weighted averages.
     *
     * \param maxAnisotropy
     *    Denotes the highest tolerated anisotropy of the lookup
     *    kernel. This is necessary to bound the computational
     *    cost of filtered lookups.
     *
     * \param cacheFilename
     *    Optional filename of a memory-mapped cache file that is used to keep
     *    MIP map data out of core, and to avoid having to load and
     *    downsample textures over and over again in subsequent Mitsuba runs.
     *
     * \param maxValue
     *    Maximum image value. This is used to clamp out-of-range values
     *    that occur during the image resampling process
     *
     * \param intent
     *    When an RGB image is transformed into a spectral representation,
     *    this parameter specifies what conversion method should be used.
     *    See \ref Spectrum::EConversionIntent for further details.
     */
    TMIPMap(Bitmap *bitmap_,
            Bitmap::EPixelFormat pixelFormat,
            Bitmap::EComponentFormat componentFormat,
            const ReconstructionFilter *rfilter,
            EBoundaryCondition bcu = ReconstructionFilter::ERepeat,
            EBoundaryCondition bcv = ReconstructionFilter::ERepeat,
            EMIPFilterType filterType = EEWA,
            Float maxAnisotropy = 20.0f,
            fs::path cacheFilename = fs::path(),
            uint64_t timestamp = 0,
            Float maxValue = 1.0f,
            Spectrum::EConversionIntent intent = Spectrum::EReflectance)
        : m_pixelFormat(pixelFormat), m_bcu(bcu), m_bcv(bcv), m_filterType(filterType),
          m_weightLut(NULL), m_maxAnisotropy(maxAnisotropy) {

        /* Keep track of time */
        ref<Timer> timer = new Timer();

        /* Compute the size of the MIP map cache file (for now
           assuming that one should be created) */
        size_t padding = sizeof(MIPMapHeader) % MTS_MIPMAP_CACHE_ALIGNMENT;
        if (padding)
            padding = MTS_MIPMAP_CACHE_ALIGNMENT - padding;
        size_t cacheSize = sizeof(MIPMapHeader) + padding +
            Array2DType::bufferSize(bitmap_->getSize());

        /* 1. Determine the number of MIP levels. The following
              code also handles non-power-of-2 input. */
        m_levels = 1;
        if (m_filterType != ENearest && m_filterType != EBilinear) {
            Vector2i size = bitmap_->getSize();
            while (size.x > 1 || size.y > 1) {
                size.x = std::max(1, (size.x + 1) / 2);
                size.y = std::max(1, (size.y + 1) / 2);
                cacheSize += Array2DType::bufferSize(size);
                ++m_levels;
            }
        }

        stats::mipStorage += cacheSize;

        /* Potentially create a MIP map cache file */
        uint8_t *mmapData = NULL, *mmapPtr = NULL;
        if (!cacheFilename.empty()) {
            Log(EInfo, "Generating MIP map cache file \"%s\" ..", cacheFilename.string().c_str());
            try {
                m_mmap = new MemoryMappedFile(cacheFilename, cacheSize);
            } catch (std::runtime_error &e) {
                Log(EWarn, "Unable to create MIP map cache file \"%s\" -- "
                    "retrying with a temporary file. Error message was: %s",
                    cacheFilename.string().c_str(), e.what());
                m_mmap = MemoryMappedFile::createTemporary(cacheSize);
            }
            mmapData = mmapPtr = (uint8_t *) m_mmap->getData();
        }

        /* 2. Store the base image in a suitable memory layout */
        m_pyramid = new Array2DType[m_levels];
        m_sizeRatio = new Vector2[m_levels];

        /* Allocate memory for the first MIP map level */
        if (mmapPtr) {
            mmapPtr += sizeof(MIPMapHeader) + padding;
            m_pyramid[0].map(mmapPtr, bitmap_->getSize());
            mmapPtr += m_pyramid[0].getBufferSize();
        } else {
            m_pyramid[0].alloc(bitmap_->getSize());
        }

        /* Initialize the first mip map level and extract some general
           information (i.e. the minimum, maximum, and average texture value) */
        ref<Bitmap> bitmap = bitmap_->expand()->convert(pixelFormat,
            componentFormat, 1.0f, 1.0f, intent);

        m_pyramid[0].cleanup();
        m_pyramid[0].init((Value *) bitmap->getData(), m_minimum, m_maximum, m_average);

        if (m_minimum.min() < 0) {
            Log(EWarn, "The texture contains negative pixel values! These will be clamped!");
            Value *value = (Value *) bitmap->getData();

            for (size_t i=0, count=bitmap->getPixelCount(); i<count; ++i)
                (*value++).clampNegative();

            m_pyramid[0].init((Value *) bitmap->getData(), m_minimum, m_maximum, m_average);
        }

        m_sizeRatio[0] = Vector2(1, 1);

        /* 3. Progressively downsample until only a 1x1 image is left */
        if (m_filterType != ENearest && m_filterType != EBilinear) {
            Vector2i size = bitmap_->getSize();
            m_levels = 1;
            while (size.x > 1 || size.y > 1) {
                /* Compute the size of the next downsampled layer */
                size.x = std::max(1, (size.x + 1) / 2);
                size.y = std::max(1, (size.y + 1) / 2);

                /* Either allocate memory or index into the memory map file */
                if (mmapPtr) {
                    m_pyramid[m_levels].map(mmapPtr, size);
                    mmapPtr += m_pyramid[m_levels].getBufferSize();
                } else {
                    m_pyramid[m_levels].alloc(size);
                }

                bitmap = bitmap->resample(rfilter, bcu, bcv, size, 0.0f, maxValue);
                m_pyramid[m_levels].cleanup();
                m_pyramid[m_levels].init((Value *) bitmap->getData());
                m_sizeRatio[m_levels] = Vector2(
                    (Float) size.x / (Float) m_pyramid[0].getWidth(),
                    (Float) size.y / (Float) m_pyramid[0].getHeight());

                ++m_levels;
            }
        }

        if (mmapData) {
            /* If a cache file was requested, create a header that
               describes the current MIP map configuration */
            MIPMapHeader header;
            memcpy(header.identifier, "MIP", 3);
            header.version = MTS_MIPMAP_CACHE_VERSION;
            header.pixelFormat = (uint8_t) m_pixelFormat;
            header.levels = (uint8_t) m_levels;
            header.bcu = (uint8_t) bcu;
            header.bcv = (uint8_t) bcv;
            header.filterType = (uint8_t) m_filterType;
            header.gamma = (float) bitmap_->getGamma();
            header.width = bitmap_->getWidth();
            header.height = bitmap_->getHeight();
            header.timestamp = timestamp;
            header.minimum = m_minimum;
            header.maximum = m_maximum;
            header.average = m_average;
            memcpy(mmapData, &header, sizeof(MIPMapHeader));
        }

        Log(EDebug, "Created %s of MIP maps in %i ms", memString(
            getBufferSize()).c_str(), timer->getMilliseconds());

        if (m_filterType == EEWA) {
            m_weightLut = static_cast<Float *>(allocAligned(sizeof(Float) * MTS_MIPMAP_LUT_SIZE));
            for (int i=0; i<MTS_MIPMAP_LUT_SIZE; ++i) {
                Float r2 = (Float) i / (Float) (MTS_MIPMAP_LUT_SIZE-1);
                m_weightLut[i] = math::fastexp(-2.0f * r2) - math::fastexp(-2.0f);
            }
        }
    }

    /**
     * \brief Construct a new MIP map from a previously created cache file
     *
     * \param cacheFilename
     *    Filename of a memory-mapped cache file that is used to keep
     *    MIP map data out of core, and to avoid having to load and
     *    downsample textures over and over again in subsequent Mitsuba runs.
     *
     * \param maxAnisotropy
     *    Denotes the highest tolerated anisotropy of the lookup
     *    kernel. This is necessary to bound the computational
     *    cost of filtered lookups. This parameter is independent of the
     *    cache file that was previously created.
     */
    TMIPMap(fs::path cacheFilename, Float maxAnisotropy = 20.0f)
            : m_weightLut(NULL), m_maxAnisotropy(maxAnisotropy) {
        m_mmap = new MemoryMappedFile(cacheFilename);
        uint8_t *mmapPtr = (uint8_t *) m_mmap->getData();
        Log(EInfo, "Mapped MIP map cache file \"%s\" into memory (%s).", cacheFilename.string().c_str(),
            memString(m_mmap->getSize()).c_str());

        stats::mipStorage += m_mmap->getSize();

        /* Load the file header, and run some santity checks */
        MIPMapHeader header;
        memcpy(&header, mmapPtr, sizeof(MIPMapHeader));
        Assert(header.identifier[0] == 'M' && header.identifier[1] == 'I'
            && header.identifier[2] == 'P' && header.version == MTS_MIPMAP_CACHE_VERSION);
        m_pixelFormat = (Bitmap::EPixelFormat) header.pixelFormat;
        m_levels = (int) header.levels;
        m_bcu = (EBoundaryCondition) header.bcu;
        m_bcv = (EBoundaryCondition) header.bcv;
        m_filterType = (EMIPFilterType) header.filterType;
        m_minimum = header.minimum;
        m_maximum = header.maximum;
        m_average = header.average;

        /* Move the pointer to the beginning of the MIP map data */
        size_t padding = sizeof(MIPMapHeader) % MTS_MIPMAP_CACHE_ALIGNMENT;
        if (padding)
            padding = MTS_MIPMAP_CACHE_ALIGNMENT - padding;
        mmapPtr += sizeof(MIPMapHeader) + padding;

        /* Map the highest resolution level */
        m_pyramid = new Array2DType[m_levels];
        m_sizeRatio = new Vector2[m_levels];
        Vector2i size(header.width, header.height);
        m_pyramid[0].map(mmapPtr, size);
        mmapPtr += m_pyramid[0].getBufferSize();
        m_sizeRatio[0] = Vector2(1, 1);

        if (m_filterType != ENearest && m_filterType != EBilinear) {
            /* Map the remainder of the image pyramid */
            int level = 1;
            while (size.x > 1 || size.y > 1) {
                size.x = std::max(1, (size.x + 1) / 2);
                size.y = std::max(1, (size.y + 1) / 2);
                m_pyramid[level].map(mmapPtr, size);
                m_sizeRatio[level] = Vector2(
                    (Float) size.x / (Float) m_pyramid[0].getWidth(),
                    (Float) size.y / (Float) m_pyramid[0].getHeight());
                mmapPtr += m_pyramid[level++].getBufferSize();
            }
            Assert(level == m_levels);
        }

        if (m_filterType == EEWA) {
            m_weightLut = static_cast<Float *>(allocAligned(sizeof(Float) * MTS_MIPMAP_LUT_SIZE));
            for (int i=0; i<MTS_MIPMAP_LUT_SIZE; ++i) {
                Float r2 = (Float) i / (Float) (MTS_MIPMAP_LUT_SIZE-1);
                m_weightLut[i] = math::fastexp(-2.0f * r2) - math::fastexp(-2.0f);
            }
        }
    }

    /// Release all memory
    ~TMIPMap() {
        delete[] m_pyramid;
        delete[] m_sizeRatio;
        if (m_weightLut)
            freeAligned(m_weightLut);
    }

    /**
     * \brief Check if a MIP map cache is up-to-date and matches the
     * desired configuration
     *
     * \param path
     *    File system path of the MIP map cache
     * \param timestamp
     *    Timestamp of the original texture file
     * \param bcu
     *    Horizontal boundary condition used when creating the MIP map pyramid
     * \param bcv
     *    Vertical boundary condition used when creating the MIP map pyramid
     * \param gamma
     *    If nonzero, it is verified that the provided gamma value
     *    matches that of the cache file.
     * \return \c true if the texture file is good for use
     */
    static bool validateCacheFile(const fs::path &path, uint64_t timestamp,
            Bitmap::EPixelFormat pixelFormat, EBoundaryCondition bcu,
            EBoundaryCondition bcv, EMIPFilterType filterType, Float gamma) {
        fs::ifstream is(path);
        if (!is.good())
            return false;

        MIPMapHeader header;
        is.read((char *) &header, sizeof(MIPMapHeader));
        if (is.fail())
            return false;

        if (header.identifier[0] != 'M' || header.identifier[1] != 'I'
            || header.identifier[2] != 'P' || header.version != MTS_MIPMAP_CACHE_VERSION
            || header.timestamp != timestamp
            || header.bcu != (uint8_t) bcu || header.bcv != (uint8_t) bcv
            || header.pixelFormat != (uint8_t) pixelFormat
            || header.filterType != (uint8_t) filterType)
            return false;

        if (gamma != 0 && (float) gamma != header.gamma)
            return false;

        /* Sanity check on the expected file size */
        size_t padding = sizeof(MIPMapHeader) % MTS_MIPMAP_CACHE_ALIGNMENT;
        if (padding)
            padding = MTS_MIPMAP_CACHE_ALIGNMENT - padding;

        Vector2i size(header.width, header.height);
        size_t expectedFileSize = sizeof(MIPMapHeader) + padding
            + Array2DType::bufferSize(size);

        if (filterType != ENearest && filterType != EBilinear) {
            while (size.x > 1 || size.y > 1) {
                size.x = std::max(1, (size.x + 1) / 2);
                size.y = std::max(1, (size.y + 1) / 2);
                expectedFileSize += Array2DType::bufferSize(size);
            }
        }

        return fs::file_size(path) == expectedFileSize;
    }

    /// Return the size of all buffers
    size_t getBufferSize() const {
        size_t size = 0;
        for (int i=0; i<m_levels; ++i)
            size += m_pyramid[i].getBufferSize();
        return size;
    }

    /// Return the size of the underlying full resolution texture
    inline const Vector2i &getSize() const { return m_pyramid[0].getSize(); }

    /// Return the width of the represented texture
    inline int getWidth() const { return getSize().x; }

    /// Return the height of the represented texture
    inline int getHeight() const { return getSize().y; }

    /// Return the number of mip-map levels
    inline int getLevels() const { return m_levels; }

    /// Return the filter type that is used to pre-filter lookups
    inline EMIPFilterType getFilterType() const { return m_filterType; }

    /// Get the component-wise maximum at the zero level
    inline const Value &getMinimum() const { return m_minimum; }

    /// Get the component-wise minimum
    inline const Value &getMaximum() const { return m_maximum; }

    /// Get the component-wise average
    inline const Value &getAverage() const { return m_average; }

    /// Return the blocked array used to store a given MIP level
    inline const Array2DType &getArray(int level = 0) const {
        return m_pyramid[level];
    }

    /// Return a bitmap representation of the given level
    ref<Bitmap> toBitmap(int level = 0) const {
        const Array2DType &array = m_pyramid[level];
        ref<Bitmap> result = new Bitmap(
            m_pixelFormat,
            Bitmap::componentFormat<typename QuantizedValue::Scalar>(),
            array.getSize()
        );

        array.copyTo((QuantizedValue *) result->getData());

        return result;
    }

    /**
     * \brief Return the texture value at a texel specified using integer
     * coordinates, while accounting for boundary conditions
     */
    inline Value evalTexel(int level, int x, int y) const {
        const Vector2i &size = m_pyramid[level].getSize();

        if (x < 0 || x >= size.x) {
            /* Encountered an out of bounds access -- determine what to do */
            switch (m_bcu) {
                case ReconstructionFilter::ERepeat:
                    // Assume that the input repeats in a periodic fashion
                    x = math::modulo(x, size.x);
                    break;
                case ReconstructionFilter::EClamp:
                    // Clamp to the outermost sample position
                    x = math::clamp(x, 0, size.x - 1);
                    break;
                case ReconstructionFilter::EMirror:
                    // Assume that the input is mirrored along the boundary
                    x = math::modulo(x, 2*size.x);
                    if (x >= size.x)
                        x = 2*size.x - x - 1;
                    break;
                case ReconstructionFilter::EZero:
                    // Assume that the input function is zero
                    // outside of the defined domain
                    return Value(0.0f);
                case ReconstructionFilter::EOne:
                    // Assume that the input function is equal
                    // to one outside of the defined domain
                    return Value(1.0f);
            }
        }

        if (y < 0 || y >= size.y) {
            /* Encountered an out of bounds access -- determine what to do */
            switch (m_bcv) {
                case ReconstructionFilter::ERepeat:
                    // Assume that the input repeats in a periodic fashion
                    y = math::modulo(y, size.y);
                    break;
                case ReconstructionFilter::EClamp:
                    // Clamp to the outermost sample position
                    y = math::clamp(y, 0, size.y - 1);
                    break;
                case ReconstructionFilter::EMirror:
                    // Assume that the input is mirrored along the boundary
                    y = math::modulo(y, 2*size.y);
                    if (y >= size.y)
                        y = 2*size.y - y - 1;
                    break;
                case ReconstructionFilter::EZero:
                    // Assume that the input function is zero
                    // outside of the defined domain
                    return Value(0.0f);
                case ReconstructionFilter::EOne:
                    // Assume that the input function is equal
                    // to one outside of the defined domain
                    return Value(1.0f);
            }
        }

        return Value(m_pyramid[level](x, y));
    }

    /// Evaluate the texture at the given resolution using a box filter
    inline Value evalBox(int level, const Point2 &uv) const {
        const Vector2i &size = m_pyramid[level].getSize();
        return evalTexel(level, math::floorToInt(uv.x*size.x), math::floorToInt(uv.y*size.y));
    }

    /**
     * \brief Evaluate the texture using fractional texture coordinates and
     * bilinear interpolation. The desired MIP level must be specified
     */
    inline Value evalBilinear(int level, const Point2 &uv) const {
        if (EXPECT_NOT_TAKEN(!std::isfinite(uv.x) || !std::isfinite(uv.y))) {
            Log(EWarn, "evalBilinear(): encountered a NaN!");
            return Value(0.0f);
        } else if (EXPECT_NOT_TAKEN(level >= m_levels)) {
            /* The lookup is larger than the entire texture */
            return evalBox(m_levels-1, uv);
        }

        /* Convert to fractional pixel coordinates on the specified level */
        const Vector2i &size = m_pyramid[level].getSize();
        Float u = uv.x * size.x - 0.5f, v = uv.y * size.y - 0.5f;

        int xPos = math::floorToInt(u), yPos = math::floorToInt(v);
        Float dx1 = u - xPos, dx2 = 1.0f - dx1,
              dy1 = v - yPos, dy2 = 1.0f - dy1;

        return evalTexel(level, xPos, yPos) * dx2 * dy2
             + evalTexel(level, xPos, yPos + 1) * dx2 * dy1
             + evalTexel(level, xPos + 1, yPos) * dx1 * dy2
             + evalTexel(level, xPos + 1, yPos + 1) * dx1 * dy1;
    }

    /**
     * \brief Evaluate the gradient of the texture at the given MIP level
     */
    inline void evalGradientBilinear(int level, const Point2 &uv, Value *gradient) const {
        if (EXPECT_NOT_TAKEN(!std::isfinite(uv.x) || !std::isfinite(uv.y))) {
            Log(EWarn, "evalGradientBilinear(): encountered a NaN!");
            gradient[0] = gradient[1] = Value(0.0f);
            return;
        } else if (EXPECT_NOT_TAKEN(level >= m_levels)) {
            evalGradientBilinear(m_levels-1, uv, gradient);
            return;
        }

        /* Convert to fractional pixel coordinates on the specified level */
        const Vector2i &size = m_pyramid[level].getSize();
        Float u = uv.x * size.x - 0.5f, v = uv.y * size.y - 0.5f;

        int xPos = math::floorToInt(u), yPos = math::floorToInt(v);
        Float dx = u - xPos, dy = v - yPos;

        const Value p00 = evalTexel(level, xPos,   yPos);
        const Value p10 = evalTexel(level, xPos+1, yPos);
        const Value p01 = evalTexel(level, xPos,   yPos+1);
        const Value p11 = evalTexel(level, xPos+1, yPos+1);
        Value tmp = p01 + p10 - p11;

        gradient[0] = (p10 + p00*(dy-1) - tmp*dy) * static_cast<Float> (size.x);
        gradient[1] = (p01 + p00*(dx-1) - tmp*dx) * static_cast<Float> (size.y);
    }

    /// \brief Perform a filtered texture lookup using the configured method
    Value eval(const Point2 &uv, const Vector2 &d0, const Vector2 &d1) const {
        if (m_filterType == ENearest)
            return evalBox(0, uv);
        else if (m_filterType == EBilinear)
            return evalBilinear(0, uv);

        /* Convert into texel coordinates */
        const Vector2i &size = m_pyramid[0].getSize();
        Float du0 = d0.x * size.x, dv0 = d0.y * size.y,
              du1 = d1.x * size.x, dv1 = d1.y * size.y;

        /* Turn the texture-space Jacobian into the coefficients of an
           implicitly defined ellipse. */
        Float A = dv0*dv0 + dv1*dv1,
              B = -2.0f * (du0*dv0 + du1*dv1),
              C = du0*du0 + du1*du1,
              F = A*C - B*B*0.25f;

        /* Compute the major and minor radii */
        Float root = math::hypot2(A-C, B),
              Aprime = 0.5f * (A + C - root),
              Cprime = 0.5f * (A + C + root),
              majorRadius = Aprime != 0 ? std::sqrt(F / Aprime) : 0,
              minorRadius = Cprime != 0 ? std::sqrt(F / Cprime) : 0;

        if (m_filterType == ETrilinear || !(minorRadius > 0) || !(majorRadius > 0) || F < 0) {
            /* Determine a suitable mip map level, while preferring
               blurring over aliasing */
            Float level = math::log2(std::max(majorRadius, Epsilon));
            int ilevel = math::floorToInt(level);

            if (ilevel < 0) {
                /* Bilinear interpolation (lookup is smaller than 1 pixel) */
                return evalBilinear(0, uv);
            } else {
                /* Trilinear interpolation between two mipmap levels */
                Float a = level - ilevel;
                return evalBilinear(ilevel,   uv) * (1.0f - a)
                     + evalBilinear(ilevel+1, uv) * a;
            }
        } else {
            /* Artificially enlarge ellipses that are too skinny
               to avoid having to compute excessively many samples */
            if (minorRadius * m_maxAnisotropy < majorRadius) {
                /* Enlarge the minor radius, which will cause extra blurring */
                minorRadius = majorRadius / m_maxAnisotropy;

                /* We need to find the coefficients of the adjusted ellipse, which
                   unfortunately involves expensive trig and arctrig functions.
                   Fortunately, this is somewhat of a corner case and won't
                   happen overly often in practice. */
                Float theta = 0.5f * std::atan(B / (A-C)), sinTheta, cosTheta;
                math::sincos(theta, &sinTheta, &cosTheta);

                Float a2 = majorRadius*majorRadius,
                      b2 = minorRadius*minorRadius,
                      sinTheta2 = sinTheta*sinTheta,
                      cosTheta2 = cosTheta*cosTheta,
                      sin2Theta = 2*sinTheta*cosTheta;

                A = a2*cosTheta2 + b2*sinTheta2;
                B = (a2-b2) * sin2Theta;
                C = a2*sinTheta2 + b2*cosTheta2;
                F = a2*b2;

                ++stats::clampedAnisotropy;
            }
            stats::clampedAnisotropy.incrementBase();

            /* Switch to normalized coefficients */
            Float scale = 1.0f / F;
            A *= scale; B *= scale; C *= scale;

            /* Determine a suitable MIP map level, such that the filter
               covers a reasonable amount of pixels */
            Float level = std::max((Float) 0.0f, math::log2(minorRadius));
            int ilevel = (int) level;
            Float a = level - ilevel;

            /* Switch to bilinear interpolation, be wary of round-off errors */
            if (majorRadius < 1 || !(A > 0 && C > 0))
                return evalBilinear(ilevel, uv);
            else
                return evalEWA(ilevel,   uv, A, B, C) * (1.0f-a) +
                       evalEWA(ilevel+1, uv, A, B, C) * a;
        }
    }

    /// Return a human-readable string representation
    std::string toString() const {
        std::ostringstream oss;
        oss << "TMIPMap[" << endl
            << "   pixelFormat = " << m_pixelFormat << "," << endl
            << "   size = " << memString(getBufferSize()) << "," << endl
            << "   levels = " << m_levels << "," << endl
            << "   cached = " << (m_mmap.get() ? "yes" : "no") << "," << endl
            << "   filterType = ";

        switch (m_filterType) {
            case ENearest: oss << "nearest," << endl; break;
            case EBilinear: oss << "bilinear," << endl; break;
            case ETrilinear: oss << "trilinear," << endl; break;
            case EEWA: oss << "ewa," << endl; break;
        }

        oss << "   bc = [" << m_bcu << ", " << m_bcv << "]," << endl
            << "   minimum = " << m_minimum.toString() << "," << endl
            << "   maximum = " << m_maximum.toString() << "," << endl
            << "   average = " << m_average.toString() << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    /// Header file for MIP map cache files
    struct MIPMapHeader {
        char identifier[3];
        uint8_t version;
        uint8_t pixelFormat;
        uint8_t levels;
        uint8_t bcu:4;
        uint8_t bcv:4;
        uint8_t filterType;
        float gamma;
        int width;
        int height;
        uint64_t timestamp;
        Value minimum;
        Value maximum;
        Value average;
    };


    /// Calculate the elliptically weighted average of a sample and associated Jacobian
    Value evalEWA(int level, const Point2 &uv, Float A, Float B, Float C) const {
        Assert(A > 0);
        if (EXPECT_NOT_TAKEN(!std::isfinite(A+B+C+uv.x+uv.y))) {
            Log(EWarn, "evalEWA(): encountered a NaN!");
            return Value(0.0f);
        } else if (EXPECT_NOT_TAKEN(level >= m_levels)) {
            /* The lookup is larger than the entire texture */
            return evalBox(m_levels-1, uv);
        }

        /* Convert to fractional pixel coordinates on the specified level */
        const Vector2i &size = m_pyramid[level].getSize();
        Float u = uv.x * size.x - 0.5f;
        Float v = uv.y * size.y - 0.5f;

        /* Do the same to the ellipse coefficients */
        const Vector2 &ratio = m_sizeRatio[level];
        A /= ratio.x * ratio.x;
        B /= ratio.x * ratio.y;
        C /= ratio.y * ratio.y;

        /* Compute the ellipse's bounding box in texture space */
        Float invDet = 1.0f / (-B*B + 4.0f*A*C),
              deltaU = 2.0f * std::sqrt(C * invDet),
              deltaV = 2.0f * std::sqrt(A * invDet);
        int u0 = math::ceilToInt(u - deltaU), u1 = math::floorToInt(u + deltaU);
        int v0 = math::ceilToInt(v - deltaV), v1 = math::floorToInt(v + deltaV);

        /* Scale the coefficients by the size of the Gaussian lookup table */
        Float As = A * MTS_MIPMAP_LUT_SIZE,
              Bs = B * MTS_MIPMAP_LUT_SIZE,
              Cs = C * MTS_MIPMAP_LUT_SIZE;

        Value result(0.0f);
        Float denominator = 0.0f;
        Float ddq = 2*As, uu0 = (Float) u0 - u;
        int nSamples = 0;

        for (int vt = v0; vt <= v1; ++vt) {
            const Float vv = (Float) vt - v;

            Float q  = As*uu0*uu0 + (Bs*uu0 + Cs*vv)*vv;
            Float dq = As*(2*uu0 + 1) + Bs*vv;

            for (int ut = u0; ut <= u1; ++ut) {
                if (q < (Float) MTS_MIPMAP_LUT_SIZE) {
                    uint32_t qi = (uint32_t) q;
                    if (qi < MTS_MIPMAP_LUT_SIZE) {
                        const Float weight = m_weightLut[(int) q];
                        result += evalTexel(level, ut, vt) * weight;
                        denominator += weight;
                        ++nSamples;
                    }
                }

                q += dq;
                dq += ddq;
            }
        }

        if (denominator == 0) {
            /* The filter did not cover any samples..
               Revert to bilinear interpolation */
            return evalBilinear(level, uv);
        } else {
            stats::avgEWASamples += nSamples;
            stats::avgEWASamples.incrementBase();
        }

        return result / denominator;
    }
private:
    ref<MemoryMappedFile> m_mmap;
    Bitmap::EPixelFormat m_pixelFormat;
    EBoundaryCondition m_bcu, m_bcv;
    EMIPFilterType m_filterType;
    Float *m_weightLut;
    Float m_maxAnisotropy;
    Vector2 *m_sizeRatio;
    Array2DType *m_pyramid;
    int m_levels;
    Value m_minimum;
    Value m_maximum;
    Value m_average;
};

template <typename Value, typename QuantizedValue>
    Class *TMIPMap<Value, QuantizedValue>::m_theClass
        = new Class("MIPMap", false, "Object");

template <typename Value, typename QuantizedValue>
    const Class *TMIPMap<Value, QuantizedValue>::getClass() const {
    return m_theClass;
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_MIPMAP_H_ */
