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
#if !defined(__MITSUBA_CORE_BITMAP_H_)
#define __MITSUBA_CORE_BITMAP_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/half.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief General-purpose bitmap class with read and write support
 * for several common file formats.
 *
 * This class handles loading of PNG, JPEG, BMP, TGA, as well as
 * OpenEXR files, and it supports writing of PNG, JPEG and OpenEXR files.
 *
 * PNG and OpenEXR files are optionally annotated with string-valued
 * metadata, and the gamma setting can be stored as well. Please see
 * the class methods and enumerations for further detail.
 *
 * The Python version of this API contains thee additional member functions:
 * <tt>fromByteArray</tt> and <tt>toByteArray</tt> copy image data
 * between the Bitmap instance and a Python <tt>bytearray</tt>. The function
 * <tt>buffer</tt> returns a <tt>memoryview</tt>-compatible buffer
 * object that can also interoperate with <tt>numpy</tt>'s <tt>ndarray</tt>.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Bitmap : public Object {
public:
    // ======================================================================
    //! @{ \name Constructors and Enumerations
    // ======================================================================

    /**
     * This enumeration lists all pixel format types supported
     * by the \ref Bitmap class. This both determines the
     * number of channels, and how they should be interpreted
     */
    enum EPixelFormat {
        /// Single-channel luminance bitmap
        ELuminance                = 0x00,

        /// Two-channel luminance + alpha bitmap
        ELuminanceAlpha           = 0x01,

        /// RGB bitmap
        ERGB                      = 0x02,

        /// RGB bitmap + alpha channel
        ERGBA                     = 0x03,

        /// XYZ tristimulus bitmap
        EXYZ                      = 0x04,

        /// XYZ tristimulus + alpha channel
        EXYZA                     = 0x05,

        /// Spectrum bitmap
        ESpectrum                 = 0x06,

        /// Spectrum bitmap + alpha channel
        ESpectrumAlpha            = 0x07,

        /// Spectrum bitmap + alpha + weight channel (Mitsuba's internal render bucket representation)
        ESpectrumAlphaWeight      = 0x08,

        /// Bitmap with multiple spectra + alpha + weight channel (render buckets used by the 'multichannel' plugin)
        EMultiSpectrumAlphaWeight = 0x09,

        /// Arbitrary multi-channel bitmap without a fixed interpretation
        EMultiChannel             = 0x10
    };

    /// Supported per-component data formats
    enum EComponentFormat {
        /**
         * \brief 1-bit (mask) component encoding.
         *
         * Note that you can use this class to load, save, and access bitmasks -- however,
         * many of the manipulation operations (i.e. \ref crop(), \ref accumulate(), etc.)
         * do not currently handle them.
         *
         * Default gamma value: linear (1.0)
         */
        EBitmask = 0,

        /**
         * \brief 8-bit unsigned integer (\c uint8_t) component encoding
         *
         * Default gamma value: sRGB (approx. 2.2)
         */
        EUInt8,

        /**
         * \brief 16-bit unsigned integer (\c uint16_t) component encoding
         *
         * Default gamma value: linear (1.0)
         */
        EUInt16,

        /**
         * \brief 32-bit unsigned integer (\c uint32_t) component encoding
         *
         * Default gamma value: linear (1.0)
         */
        EUInt32,

        /**
         * 16-bit floating point (\c half) HDR component encoding
         *
         * Default gamma value: linear (1.0)
         */
        EFloat16,

        /**
         * \brief 32-bit floating point (\c float) HDR component encoding
         *
         * Default gamma value: linear (1.0)
         */
        EFloat32,

        /**
         * \brief 64-bit floating point (\c double) HDR component encoding
         *
         * Default gamma value: linear (1.0)
         */
        EFloat64,

        /**
         * \brief Invalid component format (used to report error conditions)
         */
        EInvalid,

        /**
         * \brief Floating point (\c float or \c double depending on the
         * compilation settings) HDR component encoding
         *
         * Default gamma value: linear (1.0)
         */
#if defined(SINGLE_PRECISION)
        EFloat = EFloat32
#else
        EFloat = EFloat64
#endif
    };

    /// Supported file formats
    enum EFileFormat {
        /**
         * \brief Portable network graphics
         *
         * The following is supported:
         * <ul>
         * <li> Loading and saving of 8/16-bit per component bitmaps for
         *   all pixel formats (ELuminance, ELuminanceAlpha, ERGB, ERGBA)</li>
         * <li> Loading and saving of 1-bit per component mask bitmaps</li>
         * <li> Loading and saving of string-valued metadata fields</li>
         * </ul>
         */
        EPNG = 0,

        /**
         * \brief OpenEXR high dynamic range file format developed by
         * Industrial Light & Magic (ILM)
         *
         * The following is supported:
         * <ul>
         *   <li>Loading and saving of \ref Eloat16 / \ref EFloat32/ \ref
         *   EUInt32 bitmaps with all supported RGB/Luminance/Alpha combinations</li>
         *   <li>Loading and saving of spectral bitmaps</tt>
         *   <li>Loading and saving of XYZ tristimulus bitmaps</tt>
         *   <li>Loading and saving of string-valued metadata fields</li>
         * </ul>
         *
         * The following is <em>not</em> supported:
         * <ul>
         *   <li>Saving of tiled images, tile-based read access</li>
         *   <li>Display windows that are different than the data window</li>
         *   <li>Loading of spectrum-valued bitmaps</li>
         * </ul>
         */
        EOpenEXR,

        /**
         * \brief RGBE image format by Greg Ward
         *
         * The following is supported
         * <ul>
         *   <li>Loading and saving of \ref EFloat32 - based RGB bitmaps</li>
         * </ul>
         */
        ERGBE,

        /**
         * \brief PFM (Portable Float Map) image format
         *
         * The following is supported
         * <ul>
         *   <li>Loading and saving of \ref EFloat32 - based Luminance or RGB bitmaps</li>
         * </ul>
         */
        EPFM,

        /**
         * \brief PPM (Portable Pixel Map) image format
         *
         * The following is supported
         * <ul>
         *   <li>Loading and saving of \ref EUInt8 and \ref EUInt16 - based RGB bitmaps</li>
         * </ul>
         */
        EPPM,

        /**
         * \brief Joint Photographic Experts Group file format
         *
         * The following is supported:
         * <ul><li>
         * Loading and saving of 8 bit per component RGB and
         * luminance bitmaps
         * </li></ul>
         */
        EJPEG,

        /**
         * \brief Truevision Advanced Raster Graphics Array file format
         *
         * The following is supported:
         * <ul><li>Loading of uncompressed 8-bit RGB/RGBA files</ul></li>
         */
        ETGA,

        /**
         * \brief Windows Bitmap file format
         *
         * The following is supported:
         * <ul><li>Loading of uncompressed 8-bit luminance and RGBA bitmaps</ul></li>
         */
        EBMP,

        /**
         * \brief Automatically detect the file format
         *
         * Note: this flag only applies when loading a file. In this case,
         * the source stream must support the \c seek() operation.
         */
        EAuto
    };

    /// List of different rotation/flip types that can be passed to \ref rotateFlip()
    enum ERotateFlipType {
        ERotateNoneFlipNone = 0,
        ERotate180FlipXY    = ERotateNoneFlipNone,
        ERotate90FlipNone   = 1,
        ERotate270FlipXY    = ERotate90FlipNone,
        ERotate180FlipNone  = 2,
        ERotateNoneFlipXY   = ERotate180FlipNone,
        ERotate270FlipNone  = 3,
        ERotate90FlipXY     = ERotate270FlipNone,
        ERotateNoneFlipX    = 4,
        ERotate180FlipY     = ERotateNoneFlipX,
        ERotate90FlipX      = 5,
        ERotate270FlipY     = ERotate90FlipX,
        ERotate180FlipX     = 6,
        ERotateNoneFlipY    = ERotate180FlipX,
        ERotate270FlipX     = 7,
        ERotate90FlipY      = ERotate270FlipX
    };

    /// The four basic arithmetic operations for use with \ref arithmeticOperation()
    enum EArithmeticOperation {
        EAddition = 0,
        ESubtraction,
        EMultiplication,
        EDivision
    };

    /**
     * \brief Describes a sub-layer of a multilayer bitmap (e.g. OpenEXR)
     *
     * A layer is defined as a named collection of bitmap channels and a pixel format
     *
     * This data structure is used by \ref Bitmap::getLayers()
     */
    struct Layer {
        /// Descriptive name of the bitmap layer
        std::string name;
        /// Pixel format that should be associated with underlying combination of bitmap channels
        EPixelFormat format;
        /// References to bitmap channels
        std::vector<int> channels;
    };

    /**
     * \brief Create a bitmap of the specified type and allocate
     * the necessary amount of memory
     *
     * \param pFmt
     *    Specifies the pixel format (e.g. RGBA or Luminance-only)
     *
     * \param cFmt
     *    Specifies how the per-pixel components are encoded
     *    (e.g. unsigned 8 bit integers or 32-bit floating point values)
     *
     * \param size
     *    Specifies the horizontal and vertical bitmap size in pixels
     *
     * \param channelCount
     *    Channel count of the image. This parameter is only required when
     *    \c pFmt = \ref EMultiChannel or \ref EMultiSpectrumAlphaWeight
     *
     * \param data
     *    External pointer to the image data. If set to \c NULL, the
     *    implementation will allocate memory itself.
     */
    Bitmap(EPixelFormat pFmt, EComponentFormat cFmt, const Vector2i &size,
        uint8_t channelCount = 0, uint8_t *data = NULL);

    /**
     * \brief Load a bitmap from an arbitrary stream data source
     *
     * \param format
     *    File format to be read (PNG/EXR/Auto-detect ...)
     *
     * \param stream
     *    Pointer to an arbitrary stream data source
     *
     * \param prefix
     *    Only consider image layers whose identifier begins with \c prefix.
     *    This is currently only supported by the OpenEXR format loader.
     */
    Bitmap(EFileFormat format, Stream *stream, const std::string &prefix = "");

    /**
     * \brief Load a bitmap from a file on disk
     *
     * \param path
     *    Path to a bitmap file (format will be auto-detected)
     *
     * \param prefix
     *    Only consider image layers whose identifier begins with \c prefix.
     *    This is currently only supported by the OpenEXR format loader.
     */
    Bitmap(const fs::path &path, const std::string &prefix = "");

    //! @}
    // ======================================================================

    // ======================================================================
    //! @{ \name Miscellaneous
    // ======================================================================

    /// Return the bitmap component format for a specified type
    template <typename T> inline static EComponentFormat componentFormat();

    /// Return the pixel format of this bitmap
    inline EPixelFormat getPixelFormat() const { return m_pixelFormat; }

    /// Return the component format of this bitmap
    inline EComponentFormat getComponentFormat() const { return m_componentFormat; }

    /// Return the bitmap's resolution in pixels
    inline const Vector2i &getSize() const { return m_size; }

    /// Return the total number of pixels
    inline size_t getPixelCount() const { return (size_t) m_size.x * (size_t) m_size.y; }

    /// Return the bitmap's width in pixels
    inline int getWidth() const { return m_size.x; }

    /// Return the bitmap's height in pixels
    inline int getHeight() const { return m_size.y; }

    /// Return the number of channels used by this bitmap
    inline int getChannelCount() const { return (int) m_channelCount; }

    /// Return whether this image has matching width and height
    inline bool isSquare() const { return m_size.x == m_size.y; }

    /// Return a string representation of the name of a channel
    std::string getChannelName(int channelIndex) const;

    /// Return whether this image has an alpha channel
    inline bool hasAlpha() const {
        return
            m_pixelFormat == ELuminanceAlpha ||
            m_pixelFormat == ERGBA ||
            m_pixelFormat == EXYZA ||
            m_pixelFormat == ESpectrumAlpha;
    }

    /// Return whether this image has a weight channel
    inline bool hasWeight() const {
        return m_pixelFormat == ESpectrumAlphaWeight ||
               m_pixelFormat == EMultiSpectrumAlphaWeight;
    }

    /// Return whether this is a generic multi-channel image
    inline bool isMultiChannel() const {
        return m_pixelFormat == EMultiChannel ||
               m_pixelFormat == EMultiSpectrumAlphaWeight;
    }

    /**
     * \brief Return the number bits per component
     *
     * When the component format is set to \c EBitmask,
     * this function returns 1; \c EUInt8 maps to 8, etc.
     */
    int getBitsPerComponent() const;

    /**
     * \brief Return the number bytes per component
     *
     * When the component format is set to \c EBitmask,
     * this function throws an exception.
     */
    int getBytesPerComponent() const;

    /**
     * \brief Return the number of bytes per pixel
     *
     * When the component format is set to \c EBitmask,
     * this function throws an exception.
     */
    inline int getBytesPerPixel() const {
        return getBytesPerComponent() * getChannelCount(); }

    /// Return the bitmap size in bytes
    size_t getBufferSize() const;

    /**
     * \brief Return the spectral color value associated with
     * a given pixel.
     *
     * When this is not a spectrum image, conversions are applied
     * as appropriate.
     *
     * \remark This function is provided here for convenience and
     * debugging purposes. Its use is discouraged, particularly for
     * performance critical code (direct buffer access will be
     * faster by at least an order of magnitude)
     */
    Spectrum getPixel(const Point2i &pos) const;

    /**
     * \brief Set the spectral color value associated with
     * a given pixel.
     *
     * When this is not a spectrum image, conversions are applied
     * as appropriate.
     *
     * \remark This function is provided here for convenience and
     * debugging purposes. Its use is discouraged, particularly for
     * performance critical code (direct buffer access will be
     * faster by at least order of magnitude)
     */
    void setPixel(const Point2i &pos, const Spectrum &value);

    /**
     * \brief Draw a horizontal line at the specified position
     *
     * Draws from x1<x2 up to and including x2.
     */
    void drawHLine(int y, int x1, int x2, const Spectrum &value);

    /**
     * \brief Draw a vertical line at the specified position
     *
     * Draws from y1<y2 up to and including y2.
     */
    void drawVLine(int x, int y1, int y2, const Spectrum &value);

    /// Draw a rectangle with the specified position and size
    void drawRect(const Point2i &offset, const Vector2i &size, const Spectrum &value);

    /// Draw a filled rectangle with the specified position and size
    void fillRect(Point2i offset, Vector2i size, const Spectrum &value);

    /**
     * \brief Convenience function to visually indicate that a thread is
     * working on a certain part of an image
     */
    void drawWorkUnit(const Point2i &offset, const Vector2i &size, int worker);

    /// Bitmap equality operator (useful for unit-tests etc.)
    bool operator==(const Bitmap &bitmap) const;

    /// Bitmap inequality operator (useful for unit-tests etc.)
    inline bool operator!=(const Bitmap &bitmap) const
        { return !operator==(bitmap); }

    /// Create an identical copy
    ref<Bitmap> clone() const;

    /// Clear the bitmap to zero
    void clear();

    /**
     * Write an encoded form of the bitmap to a stream using the specified file format
     *
     * \param format
     *    Target file format (\ref EOpenEXR, \ref EPNG, or \ref EOpenEXR)
     *
     * \param stream
     *    Target stream that will receive the encoded output
     *
     * \param compression
     *    For PNG images, this parameter can be used to control how
     *    strongly libpng will try to compress the output (with 1 being
     *    the lowest and 9 denoting the highest compression). Note that
     *    saving files with the highest compression will be very slow.
     *    For JPEG files, this denotes the desired quality (between 0 and 100,
     *    the latter being best). The default argument (-1) uses compression
     *    5 for PNG and 100 for JPEG files.
     */
    void write(EFileFormat format, Stream *stream, int compression = -1) const;

    /**
     * Write an encoded form of the bitmap to a file using the specified file format
     *
     * \param format
     *    Target file format (\ref EOpenEXR, \ref EPNG, or \ref EOpenEXR)
     *
     * \param stream
     *    Target stream that will receive the encoded output
     *
     * \param compression
     *    For PNG images, this parameter can be used to control how
     *    strongly libpng will try to compress the output (with 1 being
     *    the lowest and 9 denoting the highest compression). Note that
     *    saving files with the highest compression will be very slow.
     *    For JPEG files, this denotes the desired quality (between 0 and 100,
     *    the latter being best). The default argument (-1) uses compression
     *    5 for PNG and 100 for JPEG files.
     */
    void write(EFileFormat format, const fs::path &filename, int compression = -1) const;

    /**
     * Write an encoded form of the bitmap to a file (auto-detecting the file format)
     *
     * \param stream
     *    Target stream that will receive the encoded output
     *
     * \param compression
     *    For PNG images, this parameter can be used to control how
     *    strongly libpng will try to compress the output (with 1 being
     *    the lowest and 9 denoting the highest compression). Note that
     *    saving files with the highest compression will be very slow.
     *    For JPEG files, this denotes the desired quality (between 0 and 100,
     *    the latter being best). The default argument (-1) uses compression
     *    5 for PNG and 100 for JPEG files.
     */
    void write(const fs::path &filename, int compression = -1) const;

    //! @}
    // ======================================================================


    // ======================================================================
    //! @{ \name Conversions and transformations
    // ======================================================================

    /**
     * \brief Convert the bitmap into another pixel and/or component format
     *
     * This helper function can be used to efficiently convert a bitmap
     * between different underlying representations. For instance, it can
     * translate a 24-bit sRGB bitmap to a linear spectrum-valued representation
     * based on half-, single- or double-precision floating point-backed storage.
     *
     * This function roughly does the following:
     *
     * <ul>
     * <li>For each pixel and channel, it converts the associated value
     *   into a normalized linear-space form (any gamma of the source
     *   bitmap is removed)</li>
     * <li>The multiplier and gamma correction specified in
     *      \c targetGamma is applied</li>
     * <li>The corrected value is clamped against the representable range
     *   of the desired component format.</li>
     * <li>The clamped gamma-corrected value is then written to
     *   the new bitmap</li>
     *
     * If the pixel formats differ, this function will also perform basic
     * conversions (e.g. spectrum to rgb, luminance to uniform spectrum
     * values, etc.)
     *
     * Note that the alpha channel is assumed to be linear in both
     * the source and target bitmap, hence it won't be affected by
     * any gamma-related transformations.
     *
     * \remark This <tt>convert()</tt> variant takes a pointer to an existing
     * bitmap, which must have a matching size. It then copies all data, while
     * performing the necessary pixel/component format and gamma conversions.
     *
     * \ref target
     *      Pointer to the target bitmap, which should be of the
     *      same size.
     * \ref multiplier
     *      An optional multiplicative factor that can be
     *      applied to all color/luminance values in linear
     *      space (alpha will not be affected).
     * \ref intent
     *      When converting from RGB to spectral color values, this flag
     *      specifies how ambiguities in this highly under-constrained problem
     *      should be resolved.
     */
    void convert(Bitmap *target, Float multiplier = 1.0f,
        Spectrum::EConversionIntent intent = Spectrum::EReflectance) const;

    /**
     * \brief Convert the bitmap into another pixel and/or component format
     *
     * This helper function can be used to efficiently convert a bitmap
     * between different underlying representations. For instance, it can
     * translate a 24-bit sRGB bitmap to a linear spectrum-valued representation
     * based on half-, single- or double-precision floating point-backed storage.
     *
     * This function roughly does the following:
     *
     * <ul>
     * <li>For each pixel and channel, it converts the associated value
     *   into a normalized linear-space form (any gamma of the source
     *   bitmap is removed)</li>
     * <li>The multiplier and gamma correction specified in
     *      \c targetGamma is applied</li>
     * <li>The corrected value is clamped against the representable range
     *   of the desired component format.</li>
     * <li>The clamped gamma-corrected value is then written to
     *   the new bitmap</li>
     *
     * If the pixel formats differ, this function will also perform basic
     * conversions (e.g. spectrum to rgb, luminance to uniform spectrum
     * values, etc.)
     *
     * Note that the alpha channel is assumed to be linear in both
     * the source and target bitmap, hence it won't be affected by
     * any gamma-related transformations.
     *
     * \remark This <tt>convert()</tt> variant usually returns a new
     * bitmap instance. When the conversion would just involve copying
     * the original bitmap, the function becomes a no-op and returns
     * the current instance.
     *
     * \ref pixelFormat
     *      Specifies the desired pixel format
     * \ref componentFormat
     *      Specifies the desired component format
     * \ref gamma
     *      Specifies the desired gamma value.
     *      Special values: \c 1.0 denotes a linear space, and
     *      \c -1.0 corresponds to sRGB.
     * \ref multiplier
     *      An optional multiplicative factor that can be
     *      applied to all color/luminance values in linear
     *      space (alpha will not be affected).
     * \ref intent
     *      When converting from RGB to spectral color values, this flag
     *      specifies how ambiguities in this highly under-constrained problem
     *      should be resolved.
     */
    ref<Bitmap> convert(EPixelFormat pixelFormat, EComponentFormat componentFormat,
            Float gamma = 1.0f, Float multiplier = 1.0f,
            Spectrum::EConversionIntent intent = Spectrum::EReflectance);

    /**
     * \brief Convert the bitmap into another pixel and/or component format
     *
     * This helper function can be used to efficiently convert a bitmap
     * between different underlying representations. For instance, it can
     * translate a 24-bit sRGB bitmap to a linear spectrum-valued representation
     * based on half-, single- or double-precision floating point-backed storage.
     *
     * This function roughly does the following:
     *
     * <ul>
     * <li>For each pixel and channel, it converts the associated value
     *   into a normalized linear-space form (any gamma of the source
     *   bitmap is removed)</li>
     * <li>The multiplier and gamma correction specified in
     *      \c targetGamma is applied</li>
     * <li>The corrected value is clamped against the representable range
     *   of the desired component format.</li>
     * <li>The clamped gamma-corrected value is then written to
     *   the new bitmap</li>
     *
     * If the pixel formats differ, this function will also perform basic
     * conversions (e.g. spectrum to rgb, luminance to uniform spectrum
     * values, etc.)
     *
     * Note that the alpha channel is assumed to be linear in both
     * the source and target bitmap, hence it won't be affected by
     * any gamma-related transformations.
     *
     * \remark This <tt>convert()</tt> variant writes to a raw bitmap
     * buffer provided as a pointer.
     *
     * \ref target
     *      A pointer to a "raw" image pixel buffer
     * \ref pixelFormat
     *      Specifies the desired pixel format
     * \ref componentFormat
     *      Specifies the desired component format
     * \ref gamma
     *      Specifies the desired gamma value.
     *      Special values: \c 1.0 denotes a linear space, and
     *      \c -1.0 corresponds to sRGB.
     * \ref multiplier
     *      An optional multiplicative factor that can be
     *      applied to all color/luminance values in linear
     *      space (alpha will not be affected).
     * \ref intent
     *      When converting from RGB to spectral color values, this flag
     *      specifies how ambiguities in this highly under-constrained problem
     *      should be resolved.
     */
    void convert(void *target,
            EPixelFormat pixelFormat, EComponentFormat componentFormat,
            Float gamma = 1.0f, Float multiplier = 1.0f,
            Spectrum::EConversionIntent intent = Spectrum::EReflectance) const;

    /**
     * \brief Specialized conversion method for multi-channel HDR images
     *
     * The bitmap class pixel format \ref EMultiSpectrumAlphaWeight is used by the
     * 'multichannel' plugin to render multiple related images at the same time.
     * This function implements a conversion function analogous to \ref convert()
     * that adjusts each of the sub-images so that it has a desired pixel format.
     * The returned bitmap has the combined pixel format \ref EMultiChannel and the
     * specified component format. Names for each of the resulting channels should
     * be provided via the \c channelNames parameters.
     *
     * This feature is currently used by the \c hdrfilm and \c tiledhdrfilm plugins.
     */
    ref<Bitmap> convertMultiSpectrumAlphaWeight(const std::vector<EPixelFormat> &pixelFormats,
            EComponentFormat componentFormat, const std::vector<std::string> &channelNames) const;

    /// Similar to the above, but writes to an already existing image
    static void convertMultiSpectrumAlphaWeight(const Bitmap *source,
        const uint8_t *sourcePtr, const Bitmap *target, uint8_t *targetPtr,
        const std::vector<EPixelFormat> &pixelFormats,
        EComponentFormat componentFormat, size_t count);

    /**
     * \brief Apply Reinhard et al's tonemapper in chromaticity space
     *
     * This function
     * <ol>
     *  <li>Computes and stores the maximum and log-average luminance of the image
     *  If the \c logAvgLuminace and \c maxLuminance parameters are nonzero, it is
     *  assumed that these have already been computed, and this step is omitted.
     *  Explicitly specifying them can be useful to prevent flickering when tonemapping
     *  multiple frames of an animation.</li>
     *
     *  <li>Converts the image to xyY, applies Reinhard's photographic tonemapper in this
     *  space, and then transforms the result back to its original
     *  representation (i.e. RGB or XYZ). The global version of the
     *  tonemapping algorithm is used (i.e. the one without automatic
     *  dodging and burning).</li>
     * </ol>
     *
     * The <tt>key</tt> parameter specifies whether a low-key or high-key image is
     * desired. The value must be in the interval <tt>[0,1]</tt>. A good starting
     * point is a middle-grey, i.e. <tt>key=0.18</tt>.
     *
     * When <tt>burn=0</tt>, the entire image is re-mapped into the displayable range. This
     * may not always be desireable -- for instance, one might choose to let a highlight
     * burn out, and this is achieved by choosing burn in (0, 1].
     *
     * For reference, see
     * "Photographic Tone Reproduction for Digital Images" by Erik Reinhard, Michael Stark,
     * Peter Shirley, and James Fewerda, in ACM Transactions on Graphics 2002, Vol. 21, 3
     *
     * \remark The implementation assumes that the image has a RGB(A), XYZ(A), or Luminance(Alpha)
     * pixel format, that <tt>gamma=1</tt>, and that it uses a EFloat16/EFloat32/EFloat64
     * component format. The conversion process is destructive in the sense that it overwrites
     * the original image.
     *
     * \remark In the Python bindings, the signature of this function is:
     * <tt>(logAvgLuminance, maxLuminance) = tonemapReinhard(logAvgLuminance, maxLuminance, key, burn)</tt>
     */
    void tonemapReinhard(Float &logAvgLuminance, Float &maxLuminance,
            Float key, Float burn);

    /**
     * \brief Expand bitmask images
     *
     * When this image is a bitmask, the implementation returns an expanded
     * version in \ref EUInt8 format. Otherwise, it returns a pointer to
     * the current instance.
     */
    ref<Bitmap> expand();

    /**
     * \brief Extract one of the color channels in the
     * bitmap and return it as a separate luminance bitmap
     *
     * When this is already a single-channel bitmap, the function
     * returns a pointer to the current instance.
     */
    ref<Bitmap> extractChannel(int channelIndex) const;

    /**
     * \brief Extract several color channels of a multi-channel
     * bitmap and return them as a bitmap with the given pixel format
     */
    ref<Bitmap> extractChannels(Bitmap::EPixelFormat fmt, const std::vector<int> &channels) const;

    /**
     * \brief Extracts layer information from the bitmap
     *
     * This is a convenience function which analyzes the bitmap's
     * channel names and extracts groups matching a name
     * convention specified by the OpenEXR standard.
     *
     * A series of channel names referred to as
     *  'diffuse.R', 'diffuse.B', and 'diffuse.G'
     * will be identified as a single \ref Bitmap::ERGB "layer".
     * This works for RGB[A]/XYZ[A], and Luminance[A] images
     */
    std::vector<Layer> getLayers() const;

    /**
     * \brief Split an multi-channel image buffer (e.g. from an OpenEXR image
     * with lots of AOVs) into its constituent layers
     *
     * This operation internally calls \ref getLayers() to extract bitmap
     * layers followed by one or more calls to \ref extractChannels()
     */
    std::map<std::string, Bitmap *> split() const;

    /**
     * \brief Merges multiple bitmaps of the same type and resolution
     * into one with a larger number of channels
     *
     * \param fmt The desired pixel format of the combined bitmap
     *
     * \param sourceBitmaps An array containing the input bitmaps
     * They must all be of the same size and component format.
     *
     * \return A newly created bitmap containing all input images
     * as color channels
     */
    static ref<Bitmap> join(EPixelFormat fmt,
            const std::vector<Bitmap *> &sourceBitmaps);

    /// Crop the bitmap to the specified rectangle
    ref<Bitmap> crop(const Point2i &offset, const Vector2i &size) const;

    /// Vertically flip the image contents
    void flipVertically();

    /// Compute the average value of the bitmap
    Spectrum average() const;

    /// Perform the specified rotatation & flip operation
    ref<Bitmap> rotateFlip(ERotateFlipType type) const;

    /**
     * \brief Scale the entire image by a certain value
     *
     * Skips the image's alpha channel, if it has one. When the image uses
     * a fixed point representation and a pixel value overflows during the
     * scale operation, it is clamped to the representable range.
     */
    void scale(Float value);

    /**
     * \brief Raise the entire image to a certain value
     *
     * Skips the image's alpha channel, if it has one. When the image uses
     * a fixed point representation and a pixel value overflows during the
     * scale operation, it is clamped to the representable range.
     */
    void pow(Float value);

    /**
     * \brief Color balancing: apply the given scale factors to the
     * red, green, and blue channels of the image
     *
     * When the image is not an \c EFloat16, \c EFloat32, or
     * \c EFloat64-based RGB/RGBA image, the function throws an exception
     */
    void colorBalance(Float r, Float g, Float b);

    /**
     * Apply a color transformation matrix to the contents of the bitmap
     *
     * The implementation assumes that the contents have the
     * RGB, RGBA, XYZ, or XYZA pixel format and a floating point
     * component format.
     */
    void applyMatrix(Float matrix[3][3]);

    /**
     * \brief Copy the contents of another bitmap into the
     * region with the specified offset
     *
     * Out-of-bounds regions are ignored. It is assumed that
     * <tt>bitmap != this</tt>.
     *
     * \remark This function throws an exception when the bitmaps
     * use different component formats or channels, or when the
     * component format is \ref EBitmask.
     */
    void copyFrom(const Bitmap *bitmap, Point2i sourceOffset,
            Point2i targetOffset, Vector2i size);

    /**
     * \brief Copy the contents of another bitmap into the
     * region with the specified offset
     *
     * This convenience function calls the main <tt>copyFrom()</tt>
     * implementation with <tt>size</tt> set to <tt>bitmap->getSize()</tt>
     * and <tt>sourceOffset</tt> set to zero. Out-of-bounds regions are
     * ignored. It is assumed that <tt>bitmap != this</tt>.
     *
     * \remark This function throws an exception when the bitmaps
     * use different component formats or channels, or when the
     * component format is \ref EBitmask.
     */
    inline void copyFrom(const Bitmap *bitmap, Point2i targetOffset) {
        copyFrom(bitmap, Point2i(0), targetOffset, bitmap->getSize());
    }

    /**
     * \brief Copy the contents of another bitmap into the
     * region with the specified offset
     *
     * This convenience function calls the main <tt>copyFrom()</tt>
     * implementation with <tt>size</tt> set to <tt>bitmap->getSize()</tt>
     * and <tt>sourceOffset</tt> and <tt>targetOffset</tt>tt> set to zero.
     * Out-of-bounds regions are ignored. It is assumed
     * that <tt>bitmap != this</tt>.
     *
     * \remark This function throws an exception when the bitmaps
     * use different component formats or channels, or when the
     * component format is \ref EBitmask.
     */
    inline void copyFrom(const Bitmap *bitmap) {
        copyFrom(bitmap, Point2i(0), Point2i(0), bitmap->getSize());
    }

    /**
     * \brief Accumulate the contents of another bitmap into the
     * region with the specified offset
     *
     * Out-of-bounds regions are ignored. It is assumed that
     * <tt>bitmap != this</tt>.
     *
     * \remark This function throws an exception when the bitmaps
     * use different component formats or channels, or when the
     * component format is \ref EBitmask.
     */
    void accumulate(const Bitmap *bitmap, Point2i sourceOffset,
            Point2i targetOffset, Vector2i size);

    /**
     * \brief Accumulate the contents of another bitmap into the
     * region with the specified offset
     *
     * This convenience function calls the main <tt>accumulate()</tt>
     * implementation with <tt>size</tt> set to <tt>bitmap->getSize()</tt>
     * and <tt>sourceOffset</tt> set to zero. Out-of-bounds regions are
     * ignored. It is assumed that <tt>bitmap != this</tt>.
     *
     * \remark This function throws an exception when the bitmaps
     * use different component formats or channels, or when the
     * component format is \ref EBitmask.
     */
    inline void accumulate(const Bitmap *bitmap, Point2i targetOffset) {
        accumulate(bitmap, Point2i(0), targetOffset, bitmap->getSize());
    }

    /**
     * \brief Accumulate the contents of another bitmap into the
     * region with the specified offset
     *
     * This convenience function calls the main <tt>accumulate()</tt>
     * implementation with <tt>size</tt> set to <tt>bitmap->getSize()</tt>
     * and <tt>sourceOffset</tt> and <tt>targetOffset</tt>tt> set to zero.
     * Out-of-bounds regions are ignored. It is assumed
     * that <tt>bitmap != this</tt>.
     *
     * \remark This function throws an exception when the bitmaps
     * use different component formats or channels, or when the
     * component format is \ref EBitmask.
     */
    inline void accumulate(const Bitmap *bitmap) {
        accumulate(bitmap, Point2i(0), Point2i(0), bitmap->getSize());
    }

    /**
     * \brief Convolve the image with a (centered) convolution kernel
     *
     * When compiled with FFTW, Mitsuba will do the convolution
     * in frequency space using three FFT operations. Otherwise,
     * it falls back to a brute force method with quadratic
     * complexity.
     *
     * The image can have any resolution; the kernel should be
     * square and of odd resolution. Both images must be of the
     * same floating point-valued component format. The kernel can
     * either have one color channel or as many color channels as
     * the image to be convolved.
     *
     * The convolution is always performed in double precision.
     * irrespective of the precision of the underlying data.
     */
    void convolve(const Bitmap *kernel);

    /**
     * \brief Perform an arithmetic operation using two images
     *
     * This function can add, subtract, multiply, or divide arbitrary
     * images. If the input images have different sizes or component
     * and pixel formats, the implementation first resamples and
     * converts them into the most "expressive" format that subsumes
     * both input images (at the cost of some temporary dynamic
     * memory allocations).
     *
     * To keep the implementation simple, there is currently no
     * special treatment of integer over/underflows if the component
     * format is \ref EUInt8, \ref EUInt16, or \ref EUInt32.
     */
    static ref<Bitmap> arithmeticOperation(EArithmeticOperation operation,
            const Bitmap *bitmap1, const Bitmap *bitmap2);

    /**
     * \brief Up- or down-sample this image to a different resolution
     *
     * Uses the provided reconstruction filter and observes the specified
     * horizontal and vertical boundary conditions when looking up data
     * outside of the input domain.
     *
     * A maximum and maximum image value can be specified to prevent to
     * prevent out-of-range values that are created by the resampling process.
     *
     * The optional 'temp' parameter can be used to pass an image of
     * resolution <tt>Vector2i(target->getWidth(), this->getHeight())</tt>
     * to avoid intermediate memory allocations.
     */
    void resample(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        Bitmap *target, Bitmap *temp = NULL,
        Float minValue = 0.0f, Float maxValue = 1.0f) const;

    /**
     * \brief Up- or down-sample this image to a different resolution
     *
     * Uses the provided reconstruction filter and observes the specified
     * horizontal and vertical boundary conditions when looking up data
     * outside of the input domain.
     *
     * A maximum and maximum image value can be specified to prevent to prevent
     * out-of-range values that are created by the resampling process.
     *
     * This function allocates a new output image and returns it.
     */
    ref<Bitmap> resample(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        const Vector2i &size, Float minValue = 0.0f,
        Float maxValue = 1.0f) const;

    /**
     * \brief Apply a separable convolution filter to the image
     *
     * Applies the provided filter while observing the specified
     * horizontal and vertical boundary conditions when looking up data
     * outside of the input domain.
     *
     * In comparison to \ref convolve(), this function operates
     * in the primal domain.
     *
     * A maximum and maximum image value can be specified to prevent to
     * prevent out-of-range values that are created by the filtering process.
     *
     * The optional 'temp' parameter can be used to pass an image of
     * the same dimensions as the source and target image.
     */
    void filter(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        Bitmap *target, Bitmap *temp = NULL,
        Float minValue = 0.0f, Float maxValue = 1.0f) const;

    /**
     * \brief Apply a separable convolution filter to the image
     *
     * Applies the provided filter while observing the specified
     * horizontal and vertical boundary conditions when looking up data
     * outside of the input domain.
     *
     * In comparison to \ref convolve(), this function operates
     * in the primal domain.
     *
     * A maximum and maximum image value can be specified to prevent to
     * prevent out-of-range values that are created by the filtering process.
     *
     * This function allocates a new output image and returns it.
     */
    ref<Bitmap> filter(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        Float minValue = 0.0f, Float maxValue = 1.0f) const;

    //! @}
    // ======================================================================

    // ======================================================================
    //! @{ \name Access to metadata associated with the bitmap
    // ======================================================================

    /// Return the bitmap's gamma identifier (-1: sRGB)
    inline Float getGamma() const { return m_gamma; }

    /// Set the bitmap's gamma identifier (-1: sRGB)
    inline void setGamma(Float gamma) { m_gamma = gamma; }

    /// Set a string-valued metadata field
    inline void setMetadataString(const std::string &key, const std::string &value) {
        m_metadata.setString(key, value);
    }

    /// Return a string-valued metadata field
    inline std::string getMetadataString(const std::string &key) const {
        return m_metadata.getAsString(key);
    }

    /// Return a \ref Properties object containing the image metadata
    inline Properties &getMetadata() { return m_metadata; }

    /// Return a \ref Properties object containing the image metadata (const version)
    inline const Properties &getMetadata() const { return m_metadata; }

    /// Set the a \ref Properties object containing the image metadata
    inline void setMetadata(const Properties &metadata) { m_metadata = metadata; }

    /**
     * \brief Explicitly set the names of the image channels
     *
     * The main use of this feature is when writing \ref EMultiChannel
     * images to OpenEXR files. In other cases, the names are ignored.
     */
    void setChannelNames(const std::vector<std::string> &names);

    /// Return the names of the image channels if explicitly specified (empty by default)
    inline std::vector<std::string> &getChannelNames() { return m_channelNames; }

    /// Return the names of the image channels if explicitly specified (empty by default)
    inline const std::vector<std::string> &getChannelNames() const { return m_channelNames; }

    //! @}
    // ======================================================================

    // ======================================================================
    //! @{ \name Bitmap raster data accessors
    // ======================================================================

    /// Access the underlying raster
    inline void *getData() { return m_data; }

    /// Access the underlying bit raster (const version)
    inline const void *getData() const { return m_data; }

    /// Access the underlying raster (for uint8 bitmaps)
    inline uint8_t *getUInt8Data() { return m_data; }

    /// Access the underlying bit raster (for uint8 bitmaps, const version)
    inline const uint8_t *getUInt8Data() const { return m_data; }

    /// Access the underlying raster data (for uint16 bitmaps)
    inline uint16_t *getUInt16Data() { return (uint16_t *) m_data; }

    /// Access the underlying raster data (for uint16 bitmaps, const version)
    inline const uint16_t *getUInt16Data() const { return (const uint16_t *) m_data; }

    /// Access the underlying raster data (for uint32 bitmaps)
    inline uint32_t *getUInt32Data() { return (uint32_t *) m_data; }

    /// Access the underlying raster data (for uint32 bitmaps, const version)
    inline const uint32_t *getUInt32Data() const { return (const uint32_t *) m_data; }

    /// Access the underlying raster data (for float16 bitmaps)
    inline half *getFloat16Data() { return (half *) m_data; }

    /// Access the underlying raster data (for float16 bitmaps, const version)
    inline const half *getFloat16Data() const { return (const half *) m_data; }

    /// Access the underlying raster data (for float32 bitmaps)
    inline float *getFloat32Data() { return (float *) m_data; }

    /// Access the underlying raster data (for float32 bitmaps, const version)
    inline const float *getFloat32Data() const { return (const float *) m_data; }

    /// Access the underlying raster data (for float64 bitmaps)
    inline double *getFloat64Data() { return (double *) m_data; }

    /// Access the underlying raster data (for float64 bitmaps, const version)
    inline const double *getFloat64Data() const { return (const double *) m_data; }

    /// Access the underlying raster data (for float32/float64 bitmaps)
    inline Float *getFloatData() { return (Float *) m_data; }

    /// Access the underlying raster data (for float/float64 bitmaps, const version)
    inline const Float *getFloatData() const { return (const Float *) m_data; }

    //! @}
    // ======================================================================

    /// Run static initialization code (sets up OpenEXR for multithreading)
    static void staticInitialization();

    /// Release any resources allocated in \ref staticInitialization
    static void staticShutdown();

    /// Return some human-readable information about this bitmap
    std::string toString() const;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Bitmap();

    /// Read a file stored using the PNG file format
    void readPNG(Stream *stream);

    /// Write a file using the PNG file format
    void writePNG(Stream *stream, int compression) const;

    /// Read a file stored using the JPEG file format
    void readJPEG(Stream *stream);

    /// Save a file using the JPEG file format
    void writeJPEG(Stream *stream, int quality = 100) const;

    /// Read a file stored using the RGBE file format
    void readRGBE(Stream *stream);

    /// Write a file using the RGBE file format
    void writeRGBE(Stream *stream) const;

    /// Read a file stored using the OpenEXR file format
    void readOpenEXR(Stream *stream, const std::string &prefix);

    /// Write a file using the OpenEXR file format
    void writeOpenEXR(Stream *stream) const;

    /// Read a file stored using the PFM file format
    void readPFM(Stream *stream);

    /// Write a file using the PFM file format
    void writePFM(Stream *stream) const;

    /// Read a file stored using the PPM file format
    void readPPM(Stream *stream);

    /// Write a file using the PPM file format
    void writePPM(Stream *stream) const;

    /// Read a file stored using the TGA file format
    void readTGA(Stream *stream);

    /// Read a file stored using the BMP file format
    void readBMP(Stream *stream);

    /// Update the internally cached channel count
    void updateChannelCount();

    /// Delegate for stream loading operations
    void readStream(EFileFormat format, Stream *stream, const std::string &prefix);
protected:
    EPixelFormat m_pixelFormat;
    EComponentFormat m_componentFormat;
    Vector2i m_size;
    uint8_t *m_data;
    Float m_gamma;
    uint8_t m_channelCount;
    bool m_ownsData;
    Properties m_metadata;
    std::vector<std::string> m_channelNames;
};

/** \brief Bitmap format conversion helper class
 *
 * This class implements efficient conversions between different
 * bitmap component formats. For instance, to transform an
 * 8-bit based image to floating point values, a suitable
 * converter can be obtained as follows:
 *
 * \code
 * FormatConverter *cvt = FormatConverter::getInstance(
 *     std::make_pair(Bitmap::EUInt8, Bitmap::EFloat)
 * );
 *
 * cvt->convert(...);
 * \endcode
 *
 * The \c convert methods in \ref Bitmap rely on this class
 * and may be more convenient to use.
 *
 * \remark This is not a color management system. Depending on the
 * target type, out-of-gamut values may be clipped component-wise.
 * If a luminance scale factor is applied, that is also done component-wise
 * (instead of scaling in a space that is based on human perception, such
 * as xyY or CIELab). If this and smarter gamut remapping are needed,
 * a library such as lcms2 will be more appropriate.
 *
 * \sa Bitmap::convert(EComponentFormat, Float);
 * \sa Bitmap::convert(Spectrum *);
 * \ingroup libcore
 */
class MTS_EXPORT_CORE FormatConverter {
public:
    typedef Bitmap::EComponentFormat                Format;
    typedef std::pair<Format, Format>               Conversion;
    typedef std::map<Conversion, FormatConverter *> ConverterMap;

    /**
     * \brief Transform pixels based on the conversion implemented
     * by this class
     *
     * Note that the alpha channel is assumed to be linear in both
     * the source and target bitmap, hence it won't be affected by
     * Gamma-related transformations.
     *
     * \param pixelFormat
     *    Pixel format of the source bitmap
     * \param sourceGamma
     *    Gamma value associated with pixels from the source bitmap.
     *    Special values: 1.0 denotes a linear space, and -1.0
     *    corresponds to sRGB.
     * \param source
     *    Pointer to the first pixel of the source bitmap
     * \param destFormat
     *    Pixel format of the destination bitmap
     * \param destGamma
     *    Gamma value associated with pixels from the destination
     *    bitmap. Special values: 1.0 denotes a linear space, and -1.0
     *    corresponds to sRGB.
     * \param dest
     *    Pointer to the first pixel of the destination bitmap
     * \param count
     *    How many pixels should be transformed?
     * \param multiplier
     *    An optional multiplicative factor that will be applied to all
     *    color/luminance/spectrum values in linear space (alpha and weight
     *    values will not be affected).
     * \param intent
     *    When converting from RGB to spectral color values, this flag
     *    specifies how ambiguities in this highly under-constrained problem
     *    should be resolved.
     * \param channelCount
     *    Number of channels (not including a weight channel, if any). Only
     *    needs to be specified for general multi-channel images.
     * \sa getConversion()
     */
    virtual void convert(
            Bitmap::EPixelFormat sourceFormat, Float sourceGamma, const void *_source,
            Bitmap::EPixelFormat destFormat, Float destGamma, void *_dest,
            size_t count, Float multiplier = 1.0f,
            Spectrum::EConversionIntent intent = Spectrum::EReflectance,
            int channelCount = -1) const= 0;

    /**
     * \brief Return the format conversion implemented by this
     * \c FormatConverter instance.
     *
     * \return A format pair, where the first element indicates the
     * source format, and the second element is the target format.
     */
    virtual Conversion getConversion() const = 0;

    /// Virtual destructor to delete instances using pointers to the base type
    virtual ~FormatConverter() {}

    /**
     * \brief Return a \ref FormatConverter instance, which can convert
     * from \c conv.first to \c conv.second.
     */
    static const FormatConverter *getInstance(Conversion con);

    /// Execute static initialization code (run once at program startup)
    static void staticInitialization();

    /// Release any resources allocated in \ref staticInitialization
    static void staticShutdown();
private:
    static ConverterMap m_converters;
};

//! \cond
namespace detail {
    template <typename T> inline Bitmap::EComponentFormat cfmt() { return Bitmap::EInvalid; }
    template <> inline Bitmap::EComponentFormat cfmt<uint8_t>() { return Bitmap::EUInt8; }
    template <> inline Bitmap::EComponentFormat cfmt<uint16_t>() { return Bitmap::EUInt16; }
    template <> inline Bitmap::EComponentFormat cfmt<uint32_t>() { return Bitmap::EUInt32; }
    template <> inline Bitmap::EComponentFormat cfmt<half>() { return Bitmap::EFloat16; }
    template <> inline Bitmap::EComponentFormat cfmt<float>() { return Bitmap::EFloat32; }
    template <> inline Bitmap::EComponentFormat cfmt<double>() { return Bitmap::EFloat64; }
};

template <typename T> inline Bitmap::EComponentFormat Bitmap::componentFormat() {
    return detail::cfmt<T>();
}

extern MTS_EXPORT_CORE std::ostream &operator<<(std::ostream &os, const Bitmap::EPixelFormat &value);
extern MTS_EXPORT_CORE std::ostream &operator<<(std::ostream &os, const Bitmap::EComponentFormat &value);
//! \endcond

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_BITMAP_H_ */
