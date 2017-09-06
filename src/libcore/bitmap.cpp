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

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/version.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <boost/algorithm/string.hpp>
#include <boost/scoped_array.hpp>
#include <boost/thread/mutex.hpp>
#include <set>

#if defined(__WINDOWS__)
#undef _CRT_SECURE_NO_WARNINGS
#define _MATH_DEFINES_DEFINED
#endif

#if defined(MTS_HAS_OPENEXR)
#if defined(_MSC_VER)
#pragma warning(disable : 4231) // nonstandard extension used : 'extern' before template explicit instantiation
#endif
#include <ImfInputFile.h>
#include <ImfStandardAttributes.h>
#include <ImfRgbaYca.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfIntAttribute.h>
#include <ImfFloatAttribute.h>
#include <ImfDoubleAttribute.h>
#include <ImfVecAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfVersion.h>
#include <ImfIO.h>
#include <ImathBox.h>
#endif

#if defined(MTS_HAS_LIBPNG)
#include <png.h>
#endif

#if defined(MTS_HAS_LIBJPEG)
extern "C" {
    #include <jpeglib.h>
    #include <jerror.h>
};
#endif

#if defined(MTS_HAS_FFTW)
#include <complex>
#include <fftw3.h>
#endif

MTS_NAMESPACE_BEGIN

namespace
{
// Safely convert between scalar types avoiding downcasting warning
template <typename T, typename S> inline T safe_cast(S a) {
    return static_cast<T>(a);
}
template <> inline half safe_cast(double a) {
    return static_cast<half>(static_cast<float>(a));
}
}

#if defined(MTS_HAS_OPENEXR)
/* ========================== *
 *     EXR helper classes     *
 * ========================== */

class EXRIStream : public Imf::IStream {
public:
    EXRIStream(Stream *stream) : IStream(stream->toString().c_str()),
        m_stream(stream) {
        m_offset = stream->getPos();
        m_size = stream->getSize();
    }

    bool read(char *c, int n) {
        m_stream->read(c, n);
        return m_stream->getPos() == m_size;
    }

    Imf::Int64 tellg() {
        return m_stream->getPos()-m_offset;
    }

    void seekg(Imf::Int64 pos) {
        m_stream->seek((size_t) pos + m_offset);
    }

    void clear() { }
private:
    ref<Stream> m_stream;
    size_t m_offset, m_size;
};

class EXROStream : public Imf::OStream {
public:
    EXROStream(Stream *stream) : OStream(stream->toString().c_str()),
        m_stream(stream) {
    }

    void write(const char *c, int n) {
        m_stream->write(c, n);
    }

    Imf::Int64 tellp() {
        return m_stream->getPos();
    }

    void seekp(Imf::Int64 pos) {
        m_stream->seek((size_t) pos);
    }

    void clear() { }
private:
    ref<Stream> m_stream;
};

inline bool chromaticitiesMatch(const Imf::Chromaticities &a, const Imf::Chromaticities &b) {
    Float diff2 = (a.red-b.red).length2()
    + (a.green-b.green).length2()
    + (a.blue-b.blue).length2()
    + (a.white-b.white).length2();

    return diff2 < 1e-8f;
}
#endif

#if defined(MTS_HAS_LIBPNG)
/* ========================== *
 *    PNG helper functions    *
 * ========================== */

static void png_flush_data(png_structp png_ptr) {
    png_voidp flush_io_ptr = png_get_io_ptr(png_ptr);
    ((Stream *) flush_io_ptr)->flush();
}

static void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length) {
    png_voidp read_io_ptr = png_get_io_ptr(png_ptr);
    ((Stream *) read_io_ptr)->read(data, length);
}

static void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length) {
    png_voidp write_io_ptr = png_get_io_ptr(png_ptr);
    ((Stream *) write_io_ptr)->write(data, length);
}

static void png_error_func(png_structp png_ptr, png_const_charp msg) {
    SLog(EError, "Fatal libpng error: %s\n", msg);
    exit(-1);
}

static void png_warn_func(png_structp png_ptr, png_const_charp msg) {
    if (strstr(msg, "iCCP: known incorrect sRGB profile") != NULL)
        return;
    SLog(EWarn, "libpng warning: %s\n", msg);
}
#endif

#if defined(MTS_HAS_LIBJPEG)
/* ========================== *
 *   JPEG helper functions    *
 * ========================== */

extern "C" {
    static const size_t jpeg_bufferSize = 0x8000;

    typedef struct {
        struct jpeg_source_mgr mgr;
        JOCTET * buffer;
        mitsuba::Stream *stream;
    } jbuf_in_t;

    typedef struct {
        struct jpeg_destination_mgr mgr;
        JOCTET * buffer;
        mitsuba::Stream *stream;
    } jbuf_out_t;

    METHODDEF(void) jpeg_init_source(j_decompress_ptr cinfo) {
        jbuf_in_t *p = (jbuf_in_t *) cinfo->src;
        p->buffer = new JOCTET[jpeg_bufferSize];
    }

    METHODDEF(boolean) jpeg_fill_input_buffer (j_decompress_ptr cinfo) {
        jbuf_in_t *p = (jbuf_in_t *) cinfo->src;
        size_t nBytes;

        try {
            p->stream->read(p->buffer, jpeg_bufferSize);
            nBytes = jpeg_bufferSize;
        } catch (const EOFException &e) {
            nBytes = e.getCompleted();
            if (nBytes == 0) {
                /* Insert a fake EOI marker */
                p->buffer[0] = (JOCTET) 0xFF;
                p->buffer[1] = (JOCTET) JPEG_EOI;
                nBytes = 2;
            }
        }

        cinfo->src->bytes_in_buffer = nBytes;
        cinfo->src->next_input_byte = p->buffer;
        return TRUE;
    }

    METHODDEF(void) jpeg_skip_input_data (j_decompress_ptr cinfo, long num_bytes) {
        if (num_bytes > 0) {
            while (num_bytes > (long) cinfo->src->bytes_in_buffer) {
                num_bytes -= (long) cinfo->src->bytes_in_buffer;
                jpeg_fill_input_buffer(cinfo);
            }
            cinfo->src->next_input_byte += (size_t) num_bytes;
            cinfo->src->bytes_in_buffer -= (size_t) num_bytes;
        }
    }

    METHODDEF(void) jpeg_term_source (j_decompress_ptr cinfo) {
        jbuf_in_t *p = (jbuf_in_t *) cinfo->src;
        delete[] p->buffer;
    }

    METHODDEF(void) jpeg_init_destination(j_compress_ptr cinfo) {
        jbuf_out_t *p = (jbuf_out_t *)cinfo->dest;

        p->buffer = new JOCTET[jpeg_bufferSize];
        p->mgr.next_output_byte = p->buffer;
        p->mgr.free_in_buffer = jpeg_bufferSize;
    }

    METHODDEF(boolean) jpeg_empty_output_buffer(j_compress_ptr cinfo) {
        jbuf_out_t *p = (jbuf_out_t *)cinfo->dest;
        p->stream->write(p->buffer, jpeg_bufferSize);
        p->mgr.next_output_byte = p->buffer;
        p->mgr.free_in_buffer = jpeg_bufferSize;
        return 1;
    }

    METHODDEF(void) jpeg_term_destination(j_compress_ptr cinfo) {
        jbuf_out_t *p = (jbuf_out_t *)cinfo->dest;
        p->stream->write(p->buffer,
            jpeg_bufferSize-p->mgr.free_in_buffer);
        delete[] p->buffer;
        p->mgr.free_in_buffer = 0;
    }

    METHODDEF(void) jpeg_error_exit (j_common_ptr cinfo) throw(std::runtime_error) {
        char msg[JMSG_LENGTH_MAX];
        (*cinfo->err->format_message) (cinfo, msg);
        SLog(EError, "Critcal libjpeg error: %s", msg);
    }
};
#endif

/* ========================== *
 *        Bitmap class        *
 * ========================== */

Bitmap::Bitmap(EPixelFormat pFormat, EComponentFormat cFormat,
        const Vector2i &size, uint8_t channelCount, uint8_t *data) : m_pixelFormat(pFormat),
        m_componentFormat(cFormat), m_size(size), m_data(data), m_channelCount(channelCount), m_ownsData(false) {
    AssertEx(size.x > 0 && size.y > 0, "Invalid bitmap size");

    if (m_componentFormat == EUInt8)
        m_gamma = -1.0f; // sRGB by default
    else
        m_gamma = 1.0f; // Linear by default

    updateChannelCount();

    if (!m_data) {
        m_data = static_cast<uint8_t *>(allocAligned(getBufferSize()));
        m_ownsData = true;
    }
}

Bitmap::Bitmap(EFileFormat format, Stream *stream, const std::string &prefix) : m_data(NULL), m_ownsData(false) {
    readStream(format, stream, prefix);
}

Bitmap::Bitmap(const fs::path &path, const std::string &prefix) : m_data(NULL), m_ownsData(false) {
    ref<FileStream> fs = new FileStream(path, FileStream::EReadOnly);
    readStream(EAuto, fs, prefix);
}

void Bitmap::readStream(EFileFormat format, Stream *stream, const std::string &prefix)  {
    if (format == EAuto) {
        /* Try to automatically detect the file format */
        size_t pos = stream->getPos();
        uint8_t start[8];
        stream->read(start, 8);

        if (start[0] == 'B' && start[1] == 'M') {
            format = EBMP;
        } else if (start[0] == '#' && start[1] == '?') {
            format = ERGBE;
        } else if (start[0] == 'P' && (start[1] == 'F' || start[1] == 'f')) {
            format = EPFM;
        } else if (start[0] == 'P' && start[1] == '6') {
            format = EPPM;
        #if defined(MTS_HAS_LIBJPEG)
        } else if (start[0] == 0xFF && start[1] == 0xD8) {
            format = EJPEG;
        #endif
        #if defined(MTS_HAS_LIBPNG)
        } else if (png_sig_cmp(start, 0, 8) == 0) {
            format = EPNG;
        #endif
        #if defined(MTS_HAS_OPENEXR)
        } else if (Imf::isImfMagic((const char *) start)) {
            format = EOpenEXR;
        #endif
        } else {
            /* Check for a TGAv2 file */
            char footer[18];
            stream->seek(stream->getSize() - 18);
            stream->read(footer, 18);
            if (footer[17] == 0 && strncmp(footer, "TRUEVISION-XFILE.", 17) == 0)
                format = ETGA;
        }
        stream->seek(pos);
    }

    switch (format) {
        case EBMP: readBMP(stream); break;
        case EJPEG: readJPEG(stream); break;
        case EOpenEXR: readOpenEXR(stream, prefix); break;
        case ERGBE: readRGBE(stream); break;
        case EPFM: readPFM(stream); break;
        case EPPM: readPPM(stream); break;
        case ETGA: readTGA(stream); break;
        case EPNG: readPNG(stream); break;
        default:
            Log(EError, "Bitmap: Invalid file format!");
    }
}

void Bitmap::write(const fs::path &path, int compression) const {
    std::string s = boost::to_lower_copy(path.string());
    EFileFormat format;
    if (boost::ends_with(s, "jpeg") || boost::ends_with(s, "jpg"))
        format = EJPEG;
    else if (boost::ends_with(s, "png"))
        format = EPNG;
    else if (boost::ends_with(s, "exr"))
        format = EOpenEXR;
    else if (boost::ends_with(s, "hdr") || boost::ends_with(s, "rgbe"))
        format = ERGBE;
    else if (boost::ends_with(s, "pfm"))
        format = EPFM;
    else if (boost::ends_with(s, "ppm"))
        format = EPPM;
    else {
        Log(EError, "No supported bitmap file extension: \"%s\"", path.string().c_str());
        return;
    }
    write(format, path, compression);
}

void Bitmap::write(EFileFormat format, const fs::path &path, int compression) const {
    ref<FileStream> fs = new FileStream(path, FileStream::ETruncReadWrite);
    write(format, fs, compression);
}

void Bitmap::write(EFileFormat format, Stream *stream, int compression) const {
    switch (format) {
        case EJPEG:
            if (compression == -1)
                compression = 100;
            writeJPEG(stream, compression);
            break;
        case EPNG:
            if (compression == -1)
                compression = 5;
            writePNG(stream, compression);
            break;
        case EOpenEXR: writeOpenEXR(stream); break;
        case ERGBE: writeRGBE(stream); break;
        case EPFM: writePFM(stream); break;
        case EPPM: writePPM(stream); break;
        default:
            Log(EError, "Bitmap::write(): Invalid file format!");
    }
}

size_t Bitmap::getBufferSize() const {
    size_t bitsPerRow = (size_t) m_size.x * m_channelCount * getBitsPerComponent();
    size_t bytesPerRow = (bitsPerRow + 7) / 8; // round up to full bytes
    return bytesPerRow * (size_t) m_size.y;
}

std::string Bitmap::getChannelName(int idx) const {
    Assert(idx < m_channelCount);
    char name = '\0';

    switch (m_pixelFormat) {
        case ELuminance: name = 'L'; break;
        case ELuminanceAlpha: name = "YA"[idx]; break;
        case ERGBA:
        case ERGB: name = "RGBA"[idx]; break;
        case EXYZA:
        case EXYZ: name = "XYZA"[idx]; break;
        case ESpectrumAlphaWeight:
        case ESpectrumAlpha:
            if (idx == m_channelCount-1)
                return m_pixelFormat == ESpectrumAlpha ? "A" : "W";
            else if (idx == m_channelCount-2 && m_pixelFormat == ESpectrumAlphaWeight)
                return "A";
        case ESpectrum: {
                std::pair<Float, Float> coverage = Spectrum::getBinCoverage(idx);
                return formatString("%.2f-%.2fnm", coverage.first, coverage.second);
            }
        default:
            Log(EError, "Unknown pixel format!");
    }

    return std::string(1, name);
}

void Bitmap::setChannelNames(const std::vector<std::string> &names) {
    if (!names.empty() && names.size() != m_channelCount)
        Log(EError, "setChannelNames(): tried to set %i channel names for an image with %i channels!",
            (int) names.size(), m_channelCount);
    m_channelNames = names;
}

void Bitmap::updateChannelCount() {
    switch (m_pixelFormat) {
        case ELuminance: m_channelCount = 1; break;
        case ELuminanceAlpha: m_channelCount = 2; break;
        case ERGB: m_channelCount = 3; break;
        case ERGBA: m_channelCount = 4; break;
        case EXYZ: m_channelCount = 3; break;
        case EXYZA: m_channelCount = 4; break;
        case ESpectrum: m_channelCount = SPECTRUM_SAMPLES; break;
        case ESpectrumAlpha: m_channelCount = SPECTRUM_SAMPLES + 1; break;
        case ESpectrumAlphaWeight: m_channelCount = SPECTRUM_SAMPLES + 2; break;
        case EMultiSpectrumAlphaWeight: break;
        case EMultiChannel: break;
        default:
            Log(EError, "Unknown pixel format!");
    }
}


int Bitmap::getBitsPerComponent() const {
    switch (m_componentFormat) {
        case EBitmask: return 1; break;
        case EUInt8: return 8; break;
        case EUInt16: return 16; break;
        case EUInt32: return 32; break;
        case EFloat16: return 16; break;
        case EFloat32: return 32; break;
        case EFloat64: return 64; break;
        default:
            Log(EError, "Unknown component format!");
            return -1;
    }
}

int Bitmap::getBytesPerComponent() const {
    switch (m_componentFormat) {
        case EUInt8: return 1; break;
        case EUInt16: return 2; break;
        case EUInt32: return 4; break;
        case EFloat16: return 2; break;
        case EFloat32: return 4; break;
        case EFloat64: return 8; break;
        case EBitmask:
            Log(EError, "Bitmask images have less than 1 byte per component!");
            return -1;
        default:
            Log(EError, "Unknown component format!");
            return -1;
    }
}

Bitmap::~Bitmap() {
    if (m_data && m_ownsData)
        freeAligned(m_data);
}

void Bitmap::clear() {
    memset(m_data, 0, getBufferSize());
}

ref<Bitmap> Bitmap::clone() const {
    ref<Bitmap> bitmap = new Bitmap(m_pixelFormat, m_componentFormat, m_size);
    memcpy(bitmap->m_data, m_data, getBufferSize());
    bitmap->m_metadata = m_metadata;
    bitmap->m_gamma = m_gamma;
    bitmap->m_channelNames = m_channelNames;
    return bitmap;
}

void Bitmap::flipVertically() {
    if (m_componentFormat == EBitmask)
        Log(EError, "Transformations involving bitmasks are currently not supported!");
    size_t rowSize = getBufferSize() / m_size.y;
    int halfHeight = m_size.y / 2;
    uint8_t *temp = (uint8_t *) alloca(rowSize);
    for (int i=0, j=m_size.y-1; i<halfHeight; ++i) {
        memcpy(temp, m_data + i * rowSize, rowSize);
        memcpy(m_data + i * rowSize, m_data + j * rowSize, rowSize);
        memcpy(m_data + j * rowSize, temp, rowSize);
        j--;
    }
}

ref<Bitmap> Bitmap::rotateFlip(ERotateFlipType type) const {
    /* Based on the GDI+ rotate/flip function in Wine */
    if (m_componentFormat == EBitmask)
        Log(EError, "Transformations involving bitmasks are currently not supported!");

    int width = m_size.x, height = m_size.y;
    bool flip_x = (type & 6) == 2 || (type & 6) == 4;
    bool flip_y = (type & 3) == 1 || (type & 3) == 2;
    bool rotate_90 = type & 1;

    if (rotate_90)
        std::swap(width, height);

    ref<Bitmap> result = new Bitmap(m_pixelFormat, m_componentFormat,
        Vector2i(width, height), m_channelCount);

    ssize_t bypp = getBytesPerPixel(),
            src_stride = m_size.x * bypp,
            dst_stride = width * bypp;

    uint8_t *dst = result->getUInt8Data();
    uint8_t *dst_row = dst, *src_row = m_data;

    if (flip_x)
        src_row += bypp * (m_size.x - 1);

    if (flip_y)
        src_row += src_stride * (m_size.y - 1);

    ssize_t src_x_step, src_y_step;
    if (rotate_90) {
        src_x_step = flip_y ? -src_stride : src_stride;
        src_y_step = flip_x ? -bypp : bypp;
    } else {
        src_x_step = flip_x ? -bypp : bypp;
        src_y_step = flip_y ? -src_stride : src_stride;
    }

    for (int y=0; y<height; y++) {
        uint8_t *src_pixel = src_row;
        uint8_t *dst_pixel = dst_row;

        for (int x=0; x<width; x++) {
            memcpy(dst_pixel, src_pixel, bypp);
            dst_pixel += bypp;
            src_pixel += src_x_step;
        }

        src_row += src_y_step;
        dst_row += dst_stride;
    }

    return result;
}

void Bitmap::copyFrom(const Bitmap *bitmap, Point2i sourceOffset,
        Point2i targetOffset, Vector2i size) {

    if (m_componentFormat == EBitmask)
        Log(EError, "Bitmap::copy(): bitmasks are not supported!");

    Assert(getPixelFormat() == bitmap->getPixelFormat() &&
           getComponentFormat() == bitmap->getComponentFormat() &&
           getChannelCount() == bitmap->getChannelCount());

    Vector2i offsetIncrease(
        std::max(0, std::max(-sourceOffset.x, -targetOffset.x)),
        std::max(0, std::max(-sourceOffset.y, -targetOffset.y))
    );

    sourceOffset += offsetIncrease;
    targetOffset += offsetIncrease;
    size -= offsetIncrease;

    Vector2i sizeDecrease(
        std::max(0, std::max(sourceOffset.x + size.x - bitmap->getWidth(), targetOffset.x + size.x - getWidth())),
        std::max(0, std::max(sourceOffset.y + size.y - bitmap->getHeight(), targetOffset.y + size.y - getHeight())));

    size -= sizeDecrease;

    if (size.x <= 0 || size.y <= 0)
        return;

    const size_t
        pixelStride  = getBytesPerPixel(),
        sourceStride = bitmap->getWidth() * pixelStride,
        targetStride = getWidth() * pixelStride;

    const uint8_t *source = bitmap->getUInt8Data() +
        (sourceOffset.x + sourceOffset.y * (size_t) bitmap->getWidth()) * pixelStride;

    uint8_t *target = m_data +
        (targetOffset.x + targetOffset.y * (size_t) m_size.x) * pixelStride;

    for (int y = 0; y < size.y; ++y) {
        memcpy(target, source, size.x * getBytesPerPixel());
        source += sourceStride;
        target += targetStride;
    }
}

void Bitmap::accumulate(const Bitmap *bitmap, Point2i sourceOffset,
        Point2i targetOffset, Vector2i size) {
    Assert(getPixelFormat() == bitmap->getPixelFormat() &&
           getComponentFormat() == bitmap->getComponentFormat() &&
           getChannelCount() == bitmap->getChannelCount());

    Vector2i offsetIncrease(
        std::max(0, std::max(-sourceOffset.x, -targetOffset.x)),
        std::max(0, std::max(-sourceOffset.y, -targetOffset.y))
    );

    sourceOffset += offsetIncrease;
    targetOffset += offsetIncrease;
    size -= offsetIncrease;

    Vector2i sizeDecrease(
        std::max(0, std::max(sourceOffset.x + size.x - bitmap->getWidth(), targetOffset.x + size.x - getWidth())),
        std::max(0, std::max(sourceOffset.y + size.y - bitmap->getHeight(), targetOffset.y + size.y - getHeight())));

    size -= sizeDecrease;

    if (size.x <= 0 || size.y <= 0)
        return;

    const size_t
        columns      = (size_t) size.x * m_channelCount,
        pixelStride  = getBytesPerPixel(),
        sourceStride = bitmap->getWidth() * pixelStride,
        targetStride = getWidth() * pixelStride;

    const uint8_t *source = bitmap->getUInt8Data() +
        (sourceOffset.x + sourceOffset.y * (size_t) bitmap->getWidth()) * pixelStride;

    uint8_t *target = m_data +
        (targetOffset.x + targetOffset.y * (size_t) m_size.x) * pixelStride;

    for (int y = 0; y < size.y; ++y) {
        switch (m_componentFormat) {
            case EUInt8:
                for (size_t i = 0; i < columns; ++i)
                    ((uint8_t *) target)[i] = (uint8_t) std::min(0xFF, ((uint8_t *) source)[i] + ((uint8_t *) target)[i]);

                break;

            case EUInt16:
                for (size_t i = 0; i < columns; ++i)
                    ((uint16_t *) target)[i] = (uint16_t) std::min(0xFFFF, ((uint16_t *) source)[i] + ((uint16_t *) target)[i]);
                break;

            case EUInt32:
                for (size_t i = 0; i < columns; ++i)
                    ((uint32_t *) target)[i] = std::min((uint32_t) 0xFFFFFFFFUL, ((uint32_t *) source)[i] + ((uint32_t *) target)[i]);
                break;

            case EFloat16:
                for (size_t i = 0; i < columns; ++i)
                    ((half *) target)[i] += ((half *) source)[i];
                break;

            case EFloat32:
                for (size_t i = 0; i < columns; ++i)
                    ((float *) target)[i] += ((float *) source)[i];
                break;

            case EFloat64:
                for (size_t i = 0; i < columns; ++i)
                    ((double *) target)[i] += ((double *) source)[i];
                break;

            default:
                Log(EError, "Unknown component format!");
        }

        source += sourceStride;
        target += targetStride;
    }
}

Spectrum Bitmap::average() const {
    if (m_gamma != 1 || (m_componentFormat != EFloat16 &&
                m_componentFormat != EFloat32 && m_componentFormat != EFloat64))
        Log(EError, "Bitmap::average() assumes a floating point image with linear gamma!");

    size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;
    Float *accum = new Float[m_channelCount];
    memset(accum, 0, sizeof(Float) * m_channelCount);

    switch (m_componentFormat) {
        case EFloat16: {
                const half *ptr = getFloat16Data();
                for (size_t i=0; i<pixelCount; ++i)
                    for (int ch=0; ch<m_channelCount; ++ch)
                        accum[ch] += (Float) *ptr++;
            }
            break;

        case EFloat32: {
                const float *ptr = getFloat32Data();
                for (size_t i=0; i<pixelCount; ++i)
                    for (int ch=0; ch<m_channelCount; ++ch)
                        accum[ch] += (Float) *ptr++;
            }
            break;

        case EFloat64: {
                const double *ptr = getFloat64Data();
                for (size_t i=0; i<pixelCount; ++i)
                    for (int ch=0; ch<m_channelCount; ++ch)
                        accum[ch] += (Float) *ptr++;
            }
            break;

        default:
            Log(EError, "average(): Unsupported component format!");
    }

    for (int ch=0; ch<m_channelCount; ++ch)
        accum[ch] /= pixelCount;

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(EFloat, EFloat)
    );

    Spectrum result;
    cvt->convert(m_pixelFormat, 1.0f, accum,
        ESpectrum, 1.0f, &result, 1);

    delete[] accum;
    return result;
}

#if defined(MTS_HAS_FFTW)
static boost::mutex __fftw_lock;
#endif

void Bitmap::convolve(const Bitmap *_kernel) {
    if (_kernel->getWidth() != _kernel->getHeight())
        Log(EError, "Bitmap::convolve(): convolution kernel must be square!");
    if (_kernel->getWidth() % 2 != 1)
        Log(EError, "Bitmap::convolve(): convolution kernel size must be odd!");
    if (_kernel->getChannelCount() != getChannelCount() && _kernel->getChannelCount() != 1)
        Log(EError, "Bitmap::convolve(): kernel and bitmap have different channel counts!");
    if (_kernel->getComponentFormat() != getComponentFormat())
        Log(EError, "Bitmap::convolve(): kernel and bitmap have different component formats!");
    if (m_componentFormat != EFloat16 && m_componentFormat != EFloat32 && m_componentFormat != EFloat64)
        Log(EError, "Bitmap::convolve(): unsupported component format! (must be float16/float32/float64)");

    int channelCountKernel = _kernel->getChannelCount();

    size_t kernelSize   = (size_t) _kernel->getWidth(),
           hKernelSize  = kernelSize / 2,
           width        = (size_t) m_size.x,
           height       = (size_t) m_size.y;

#if defined(MTS_HAS_FFTW)
    typedef std::complex<double> complex;

    size_t paddedWidth  = width + hKernelSize,
           paddedHeight = height + hKernelSize,
           paddedSize   = paddedWidth*paddedHeight;

    __fftw_lock.lock();
    complex *kernel  = (complex *) fftw_malloc(sizeof(complex) * paddedSize),
            *kernelS = (complex *) fftw_malloc(sizeof(complex) * paddedSize),
            *data    = (complex *) fftw_malloc(sizeof(complex) * paddedSize),
            *dataS   = (complex *) fftw_malloc(sizeof(complex) * paddedSize);

    if (!kernel || !kernelS || !data || !dataS) {
        __fftw_lock.unlock();
        SLog(EError, "Bitmap::convolve(): Unable to allocate temporary memory!");
    }

    /* Create a FFTW plan for a 2D DFT of this size */
    fftw_plan p = fftw_plan_dft_2d((int) paddedHeight, (int) paddedWidth,
        (fftw_complex *) kernel, (fftw_complex *) kernelS, FFTW_FORWARD, FFTW_ESTIMATE);
    __fftw_lock.unlock();

    memset(kernel, 0, sizeof(complex)*paddedSize);

    for (int ch=0; ch<m_channelCount; ++ch) {
        memset(data, 0, sizeof(complex)*paddedSize);
        switch (m_componentFormat) {
            case EFloat16:
                /* Copy and zero-pad the convolution kernel in a wraparound fashion */
                if (ch < channelCountKernel) {
                    for (size_t y=0; y<kernelSize; ++y) {
                        ssize_t wrappedY = math::modulo(hKernelSize - (ssize_t) y, (ssize_t) paddedHeight);
                        for (size_t x=0; x<kernelSize; ++x) {
                            ssize_t wrappedX = math::modulo(hKernelSize - (ssize_t) x, (ssize_t) paddedWidth);
                            kernel[wrappedX+wrappedY*paddedWidth] = _kernel->getFloat16Data()[(x+y*kernelSize)*channelCountKernel+ch];
                        }
                    }
                }

                /* Copy and zero-pad the input data */
                for (size_t y=0; y<height; ++y)
                    for (size_t x=0; x<width; ++x)
                        data[x+y*paddedWidth] = getFloat16Data()[(x+y*width)*m_channelCount + ch];
                break;

            case EFloat32:
                /* Copy and zero-pad the convolution kernel in a wraparound fashion */
                if (ch < channelCountKernel) {
                    for (size_t y=0; y<kernelSize; ++y) {
                        ssize_t wrappedY = math::modulo(hKernelSize - (ssize_t) y, (ssize_t) paddedHeight);
                        for (size_t x=0; x<kernelSize; ++x) {
                            ssize_t wrappedX = math::modulo(hKernelSize - (ssize_t) x, (ssize_t) paddedWidth);
                            kernel[wrappedX+wrappedY*paddedWidth] = _kernel->getFloat32Data()[(x+y*kernelSize)*channelCountKernel+ch];
                        }
                    }
                }

                /* Copy and zero-pad the input data */
                for (size_t y=0; y<height; ++y)
                    for (size_t x=0; x<width; ++x)
                        data[x+y*paddedWidth] = getFloat32Data()[(x+y*width)*m_channelCount + ch];
                break;

            case EFloat64:
                /* Copy and zero-pad the convolution kernel in a wraparound fashion */
                if (ch < channelCountKernel) {
                    for (size_t y=0; y<kernelSize; ++y) {
                        ssize_t wrappedY = math::modulo(hKernelSize - (ssize_t) y, (ssize_t) paddedHeight);
                        for (size_t x=0; x<kernelSize; ++x) {
                            ssize_t wrappedX = math::modulo(hKernelSize - (ssize_t) x, (ssize_t) paddedWidth);
                            kernel[wrappedX+wrappedY*paddedWidth] = _kernel->getFloat64Data()[(x+y*kernelSize)*channelCountKernel+ch];
                        }
                    }
                }

                /* Copy and zero-pad the input data */
                for (size_t y=0; y<height; ++y)
                    for (size_t x=0; x<width; ++x)
                        data[x+y*paddedWidth] = getFloat64Data()[(x+y*width)*m_channelCount + ch];
                break;

            default:
                Log(EError, "Unsupported component format!");
        }

        /* FFT the kernel */
        if (ch < channelCountKernel)
            fftw_execute(p);

        /* FFT the image */
        fftw_execute_dft(p, (fftw_complex *) data, (fftw_complex *) dataS);

        /* Multiply in frequency space -- also conjugate and scale to use the computed FFT plan backwards */
        double factor = (double) 1 / (paddedWidth*paddedHeight);
        for (size_t i=0; i<paddedSize; ++i)
            dataS[i] = factor * std::conj(dataS[i] * kernelS[i]);

        /* "IFFT" */
        fftw_execute_dft(p, (fftw_complex *) dataS, (fftw_complex *) data);

        switch (m_componentFormat) {
            case EFloat16:
                for (size_t y=0; y<height; ++y)
                    for (size_t x=0; x<width; ++x)
                        getFloat16Data()[(x+y*width)*m_channelCount+ch] = half((float) std::real(data[x+y*paddedWidth]));
                break;

            case EFloat32:
                for (size_t y=0; y<height; ++y)
                    for (size_t x=0; x<width; ++x)
                        getFloat32Data()[(x+y*width)*m_channelCount+ch] = (float) std::real(data[x+y*paddedWidth]);
                break;

            case EFloat64:
                for (size_t y=0; y<height; ++y)
                    for (size_t x=0; x<width; ++x)
                        getFloat64Data()[(x+y*width)*m_channelCount+ch] = (double) std::real(data[x+y*paddedWidth]);
                break;

            default:
                Log(EError, "Unsupported component format!");
        }
    }
    __fftw_lock.lock();
    fftw_destroy_plan(p);
    fftw_free(kernel);
    fftw_free(kernelS);
    fftw_free(data);
    fftw_free(dataS);
    __fftw_lock.unlock();
#else
    /* Brute force fallback version */
    uint8_t *output_ = static_cast<uint8_t *>(allocAligned(getBufferSize()));
    for (int ch=0; ch<m_channelCount; ++ch) {
        int chKernel = channelCountKernel > 1 ? ch : 0;

        switch (m_componentFormat) {
            case EFloat16: {
                    const half *input = getFloat16Data();
                    const half *kernel = _kernel->getFloat16Data();
                    half *output = (half *) output_;

                    for (size_t y=0; y<height; ++y) {
                        for (size_t x=0; x<width; ++x) {
                            double result = 0;
                            for (size_t ky=0; ky<kernelSize; ++ky) {
                                for (size_t kx=0; kx<kernelSize; ++kx) {
                                    ssize_t xs = x + kx - (ssize_t) hKernelSize,
                                            ys = y + ky - (ssize_t) hKernelSize;
                                    if (xs >= 0 && ys >= 0 && xs < (ssize_t) width && ys < (ssize_t) height)
                                        result += (double) kernel[(kx+ky*kernelSize)*channelCountKernel+chKernel]
                                            * (double) input[(xs+ys*width)*m_channelCount+ch];
                                }
                            }
                            output[(x+y*width)*m_channelCount+ch] = safe_cast<half>(result);
                        }
                    }
                }
                break;

            case EFloat32: {
                    const float *input = getFloat32Data();
                    const float *kernel = _kernel->getFloat32Data();
                    float *output = (float *) output_;

                    for (size_t y=0; y<height; ++y) {
                        for (size_t x=0; x<width; ++x) {
                            double result = 0;
                            for (size_t ky=0; ky<kernelSize; ++ky) {
                                for (size_t kx=0; kx<kernelSize; ++kx) {
                                    ssize_t xs = x + kx - (ssize_t) hKernelSize,
                                            ys = y + ky - (ssize_t) hKernelSize;
                                    if (xs >= 0 && ys >= 0 && xs < (ssize_t) width && ys < (ssize_t) height)
                                        result += kernel[(kx+ky*kernelSize)*channelCountKernel+chKernel]
                                            * input[(xs+ys*width)*m_channelCount+ch];
                                }
                            }
                            output[(x+y*width)*m_channelCount+ch] = (float) result;
                        }
                    }
                }
                break;

            case EFloat64: {
                    const double *input = getFloat64Data();
                    const double *kernel = _kernel->getFloat64Data();
                    double *output = (double *) output_;

                    for (size_t y=0; y<height; ++y) {
                        for (size_t x=0; x<width; ++x) {
                            double result = 0;
                            for (size_t ky=0; ky<kernelSize; ++ky) {
                                for (size_t kx=0; kx<kernelSize; ++kx) {
                                    ssize_t xs = x + kx - (ssize_t) hKernelSize,
                                            ys = y + ky - (ssize_t) hKernelSize;
                                    if (xs >= 0 && ys >= 0 && xs < (ssize_t) width && ys < (ssize_t) height)
                                        result += kernel[(kx+ky*kernelSize)*channelCountKernel+chKernel]
                                            * input[(xs+ys*width)*m_channelCount+ch];
                                }
                            }
                            output[(x+y*width)*m_channelCount+ch] = result;
                        }
                    }
                }
                break;
            default:
                Log(EError, "Unsupported component format!");
        }
    }
    if (m_ownsData)
        freeAligned(m_data);
    m_data = output_;
    m_ownsData = true;
#endif
}

void Bitmap::scale(Float value) {
    if (m_componentFormat == EBitmask)
        Log(EError, "Bitmap::scale(): bitmasks are not supported!");

    size_t nPixels = getPixelCount(), nChannels = getChannelCount();

    if (hasAlpha()) {
        switch (m_componentFormat) {
            case EUInt8: {
                    uint8_t *data = (uint8_t *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (uint8_t) std::min((Float) 0xFF,
                                std::max((Float) 0, *data * value + (Float) 0.5f));
                            ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EUInt16: {
                    uint16_t *data = (uint16_t *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (uint16_t) std::min((Float) 0xFFFF,
                                std::max((Float) 0, *data * value + (Float) 0.5f));
                            ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EUInt32: {
                    uint32_t *data = (uint32_t *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (uint32_t) std::min((Float) 0xFFFFFFFFUL,
                                std::max((Float) 0, *data * value + (Float) 0.5f));
                            ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EFloat16: {
                    half *data = (half *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = safe_cast<half> (*data * value); ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EFloat32: {
                    float *data = (float *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (float) (*data * value); ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EFloat64: {
                    double *data = (double *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (double) (*data * value); ++data;
                        }
                        ++data;
                    }
                }
                break;

            default:
                Log(EError, "Bitmap::scale(): unexpected data format!");
        }

    } else {
        size_t nEntries = nPixels * nChannels;

        switch (m_componentFormat) {
            case EUInt8: {
                    uint8_t *data = (uint8_t *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (uint8_t) std::min((Float) 0xFF,
                            std::max((Float) 0, data[i] * value + (Float) 0.5f));
                }
                break;

            case EUInt16: {
                    uint16_t *data = (uint16_t *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (uint16_t) std::min((Float) 0xFFFF,
                            std::max((Float) 0, data[i] * value + (Float) 0.5f));
                }
                break;

            case EUInt32: {
                    uint32_t *data = (uint32_t *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (uint32_t) std::min((Float) 0xFFFFFFFFUL,
                            std::max((Float) 0, data[i] * value + (Float) 0.5f));
                }
                break;

            case EFloat16: {
                    half *data = (half *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = safe_cast<half> (data[i] * value);
                }
                break;

            case EFloat32: {
                    float *data = (float *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (float) (data[i] * value);
                }
                break;

            case EFloat64: {
                    double *data = (double *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (double) (data[i] * value);
                }
                break;

            default:
                Log(EError, "Bitmap::scale(): unexpected data format!");
        }
    }
}

void Bitmap::pow(Float value) {
    if (m_componentFormat == EBitmask)
        Log(EError, "Bitmap::pow(): bitmasks are not supported!");

    size_t nPixels = getPixelCount(), nChannels = getChannelCount();

    if (hasAlpha()) {
        switch (m_componentFormat) {
            case EUInt8: {
                    uint8_t *data = (uint8_t *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (uint8_t) std::min((Float) 0xFF,
                                std::max((Float) 0, std::pow((Float) *data, value) + (Float) 0.5f));
                            ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EUInt16: {
                    uint16_t *data = (uint16_t *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (uint16_t) std::min((Float) 0xFFFF,
                                std::max((Float) 0, std::pow((Float) *data, value) + (Float) 0.5f));
                            ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EUInt32: {
                    uint32_t *data = (uint32_t *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (uint32_t) std::min((Float) 0xFFFFFFFFUL,
                                std::max((Float) 0, std::pow((Float) *data, value) + (Float) 0.5f));
                            ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EFloat16: {
                    half *data = (half *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = safe_cast<half> (std::pow((Float) *data, value)); ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EFloat32: {
                    float *data = (float *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (float) std::pow((Float) *data, value); ++data;
                        }
                        ++data;
                    }
                }
                break;

            case EFloat64: {
                    double *data = (double *) m_data;
                    for (size_t i=0; i<nPixels; ++i) {
                        for (size_t j=0; j<nChannels-1; ++j) {
                            *data = (double) std::pow((Float) *data, value); ++data;
                        }
                        ++data;
                    }
                }
                break;

            default:
                Log(EError, "Bitmap::pow(): unexpected data format!");
        }

    } else {
        size_t nEntries = nPixels * nChannels;

        switch (m_componentFormat) {
            case EUInt8: {
                    uint8_t *data = (uint8_t *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (uint8_t) std::min((Float) 0xFF,
                            std::max((Float) 0, std::pow((Float) data[i], value) + (Float) 0.5f));
                }
                break;

            case EUInt16: {
                    uint16_t *data = (uint16_t *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (uint16_t) std::min((Float) 0xFFFF,
                            std::max((Float) 0, std::pow((Float) data[i], value) + (Float) 0.5f));
                }
                break;

            case EUInt32: {
                    uint32_t *data = (uint32_t *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (uint32_t) std::min((Float) 0xFFFFFFFFUL,
                            std::max((Float) 0, std::pow((Float) data[i], value) + (Float) 0.5f));
                }
                break;

            case EFloat16: {
                    half *data = (half *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = safe_cast<half> (std::pow((Float) data[i], value));
                }
                break;

            case EFloat32: {
                    float *data = (float *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (float) std::pow((Float) data[i], value);
                }
                break;

            case EFloat64: {
                    double *data = (double *) m_data;
                    for (size_t i=0; i<nEntries; ++i)
                        data[i] = (double) std::pow((Float) data[i], value);
                }
                break;

            default:
                Log(EError, "Bitmap::pow(): unexpected data format!");
        }
    }
}

ref<Bitmap> Bitmap::arithmeticOperation(Bitmap::EArithmeticOperation operation, const Bitmap *_bitmap1, const Bitmap *_bitmap2) {
    ref<const Bitmap> bitmap1(_bitmap1), bitmap2(_bitmap2);

    /* Determine the 'fancier' pixel / component format by a maximum operation on the enum values */
    EPixelFormat pxFmt = (EPixelFormat) std::max(bitmap1->getPixelFormat(), bitmap2->getPixelFormat());
    EComponentFormat cFmt = (EComponentFormat) std::max(bitmap1->getComponentFormat(), bitmap2->getComponentFormat());

    if (cFmt == EBitmask)
        Log(EError, "Bitmap::arithmeticOperation(): bitmasks are not supported!");

    /* Make sure that the images match in size (resample if necessary) */
    Vector2i size(
        std::max(bitmap1->getWidth(), bitmap2->getWidth()),
        std::max(bitmap1->getHeight(), bitmap2->getHeight()));

    if (bitmap1->getSize() != size) {
        bitmap1 = bitmap1->resample(NULL,
                ReconstructionFilter::EClamp,
                ReconstructionFilter::EClamp, size,
                -std::numeric_limits<Float>::infinity(),
                std::numeric_limits<Float>::infinity());
    }

    if (bitmap2->getSize() != size) {
        bitmap2 = bitmap2->resample(NULL,
                ReconstructionFilter::EClamp,
                ReconstructionFilter::EClamp, size,
                -std::numeric_limits<Float>::infinity(),
                std::numeric_limits<Float>::infinity());
    }

    /* Convert the image format appropriately (no-op, if the format already matches) */
    bitmap1 = const_cast<Bitmap *>(bitmap1.get())->convert(pxFmt, cFmt);
    bitmap2 = const_cast<Bitmap *>(bitmap2.get())->convert(pxFmt, cFmt);

    ref<Bitmap> output = new Bitmap(pxFmt, cFmt, size);
    size_t nValues = output->getPixelCount() * output->getChannelCount();

    #define IMPLEMENT_OPS() \
        switch (operation) { \
            case EAddition:       for (size_t i=0; i<nValues; ++i) dst[i] = src1[i] + src2[i]; break; \
            case ESubtraction:    for (size_t i=0; i<nValues; ++i) dst[i] = src1[i] - src2[i]; break; \
            case EMultiplication: for (size_t i=0; i<nValues; ++i) dst[i] = src1[i] * src2[i]; break; \
            case EDivision:       for (size_t i=0; i<nValues; ++i) dst[i] = src1[i] / src2[i]; break; \
        }

    switch (cFmt) {
        case EUInt8: {
                const uint8_t *src1 = bitmap1->getUInt8Data();
                const uint8_t *src2 = bitmap2->getUInt8Data();
                uint8_t *dst = output->getUInt8Data();
                IMPLEMENT_OPS();
            }
            break;

        case EUInt16: {
                const uint16_t *src1 = bitmap1->getUInt16Data();
                const uint16_t *src2 = bitmap2->getUInt16Data();
                uint16_t *dst = output->getUInt16Data();
                IMPLEMENT_OPS();
            }
            break;

        case EUInt32: {
                const uint32_t *src1 = bitmap1->getUInt32Data();
                const uint32_t *src2 = bitmap2->getUInt32Data();
                uint32_t *dst = output->getUInt32Data();
                IMPLEMENT_OPS();
            }
            break;

        case EFloat16: {
                const half *src1 = bitmap1->getFloat16Data();
                const half *src2 = bitmap2->getFloat16Data();
                half *dst = output->getFloat16Data();
                IMPLEMENT_OPS();
            }
            break;

        case EFloat32: {
                const float *src1 = bitmap1->getFloat32Data();
                const float *src2 = bitmap2->getFloat32Data();
                float *dst = output->getFloat32Data();
                IMPLEMENT_OPS();
            }
            break;

        case EFloat64: {
                const double *src1 = bitmap1->getFloat64Data();
                const double *src2 = bitmap2->getFloat64Data();
                double *dst = output->getFloat64Data();
                IMPLEMENT_OPS();
            }
            break;

        default:
            Log(EError, "Bitmap::arithmeticOperation(): unexpected data format!");
    }

    #undef IMPLEMENT_OPS

    return output;
}

void Bitmap::colorBalance(Float r, Float g, Float b) {
    if (m_pixelFormat != ERGB && m_pixelFormat != ERGBA)
        Log(EError, "colorBalance(): expected a RGB or RGBA image!");
    int stride = m_pixelFormat == ERGB ? 3 : 4;
    size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;

    switch (m_componentFormat) {
        case EFloat16: {
                half *ptr = getFloat16Data();
                for (size_t i=0; i<pixelCount; ++i) {
                    ptr[0] = half((float) ptr[0] * (float) r);
                    ptr[1] = half((float) ptr[1] * (float) g);
                    ptr[2] = half((float) ptr[2] * (float) b);
                    ptr += stride;
                }
            }
            break;
        case EFloat32: {
                float *ptr = getFloat32Data();
                for (size_t i=0; i<pixelCount; ++i) {
                    ptr[0] = (float) (ptr[0] * r);
                    ptr[1] = (float) (ptr[1] * g);
                    ptr[2] = (float) (ptr[2] * b);
                    ptr += stride;
                }
            }
            break;
        case EFloat64: {
                double *ptr = getFloat64Data();
                for (size_t i=0; i<pixelCount; ++i) {
                    ptr[0] *= (double) r;
                    ptr[1] *= (double) g;
                    ptr[2] *= (double) b;
                    ptr += stride;
                }
            }
            break;
        default:
            Log(EError, "Bitmap::colorBalance(): unexpected data format!");
    }
}

void Bitmap::setPixel(const Point2i &pos, const Spectrum &value) {
    AssertEx(pos.x >= 0 && pos.x < m_size.x &&
             pos.y >= 0 && pos.y < m_size.y, "Bitmap::setPixel(): out of bounds!");

    size_t offset = ((size_t) pos.x + m_size.x * (size_t) pos.y) * getBytesPerPixel();

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(EFloat, m_componentFormat)
    );

    cvt->convert(ESpectrum, 1.0f, &value,
        m_pixelFormat, m_gamma, m_data + offset, 1);
}

void Bitmap::drawHLine(int y, int x1, int x2, const Spectrum &value) {
    if (y < 0 || y >= m_size.y)
        return;
    x1 = std::max(x1, 0); x2 = std::min(x2, m_size.x-1);

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(EFloat, m_componentFormat)
    );
    size_t pixelStride = getBytesPerPixel();
    uint8_t *source = (uint8_t *) alloca(pixelStride);
    cvt->convert(ESpectrum, 1.0f, &value,
        m_pixelFormat, m_gamma, source, 1);

    uint8_t *target = m_data + (x1 + y*m_size.x) * pixelStride;

    for (int x=x1; x<=x2; ++x) {
        memcpy(target, source, pixelStride);
        target += pixelStride;
    }
}

void Bitmap::drawVLine(int x, int y1, int y2, const Spectrum &value) {
    if (x < 0 || x >= m_size.x)
        return;
    y1 = std::max(y1, 0); y2 = std::min(y2, m_size.y-1);

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(EFloat, m_componentFormat)
    );
    size_t pixelStride = getBytesPerPixel(),
           rowStride = pixelStride * m_size.x;
    uint8_t *source = (uint8_t *) alloca(pixelStride);
    cvt->convert(ESpectrum, 1.0f, &value,
        m_pixelFormat, m_gamma, source, 1);

    uint8_t *target = m_data + (x + y1*m_size.x) * pixelStride;

    for (int y=y1; y<=y2; ++y) {
        memcpy(target, source, pixelStride);
        target += rowStride;
    }
}

void Bitmap::drawRect(const Point2i &offset, const Vector2i &size, const Spectrum &value) {
    drawHLine(offset.y, offset.x, offset.x + size.x - 1, value);
    drawHLine(offset.y + size.y - 1, offset.x, offset.x + size.x - 1, value);
    drawVLine(offset.x, offset.y, offset.y + size.y - 1, value);
    drawVLine(offset.x + size.x - 1, offset.y, offset.y + size.y - 1, value);
}

void Bitmap::fillRect(Point2i offset, Vector2i size, const Spectrum &value) {
    int sx = std::max(0, -offset.x), sy = std::max(0, -offset.y);
    size.x -= sx; size.y -= sy; offset.x += sx; offset.y += sy;

    size.x -= std::max(0, offset.x + size.x - m_size.x);
    size.y -= std::max(0, offset.y + size.y - m_size.y);

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(EFloat, m_componentFormat)
    );
    size_t pixelStride = getBytesPerPixel(),
           rowStride = pixelStride * m_size.x;
    uint8_t *source = (uint8_t *) alloca(pixelStride);

    cvt->convert(ESpectrum, 1.0f, &value,
        m_pixelFormat, m_gamma, source, 1);

    uint8_t *target = m_data + (offset.x + offset.y*m_size.x) * pixelStride;

    for (int y=0; y<size.y; ++y) {
        uint8_t *ptr = target;
        for (int x=0; x<size.x; ++x) {
            memcpy(ptr, source, pixelStride);
            ptr += pixelStride;
        }
        target += rowStride;
    }
}

Spectrum Bitmap::getPixel(const Point2i &pos) const {
    AssertEx(pos.x >= 0 && pos.x < m_size.x &&
             pos.y >= 0 && pos.y < m_size.y, "Bitmap::getPixel(): out of bounds!");

    size_t offset = ((size_t) pos.x + m_size.x * (size_t) pos.y) * getBytesPerPixel();

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(m_componentFormat, EFloat)
    );

    Spectrum result;
    cvt->convert(m_pixelFormat, m_gamma, m_data + offset,
        ESpectrum, 1.0f, &result, 1);

    return result;
}

void Bitmap::convert(Bitmap *target, Float multiplier, Spectrum::EConversionIntent intent) const {
    if (m_componentFormat == EBitmask || target->getComponentFormat() == EBitmask)
        Log(EError, "Conversions involving bitmasks are currently not supported!");
    if (m_size != target->getSize())
        Log(EError, "Bitmap::convert(): size mismatch!");

    if (m_pixelFormat == target->getPixelFormat() &&
        m_componentFormat == target->getComponentFormat() &&
        m_gamma == target->getGamma() && multiplier == 1.0f) {
        /* No conversion is necessary -- just run memcpy */
        memcpy(target->getData(), getData(), getBufferSize());
        return;
    }

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(m_componentFormat, target->getComponentFormat())
    );

    Assert(cvt != NULL);

    cvt->convert(m_pixelFormat, m_gamma, m_data,
        target->getPixelFormat(), target->getGamma(), target->getData(),
        (size_t) m_size.x * (size_t) m_size.y, multiplier, intent,
        m_channelCount);
}

ref<Bitmap> Bitmap::convert(EPixelFormat pixelFormat,
        EComponentFormat componentFormat, Float gamma, Float multiplier,
        Spectrum::EConversionIntent intent) {
    if (m_componentFormat == EBitmask || componentFormat == EBitmask)
        Log(EError, "Conversions involving bitmasks are currently not supported!");

    if (m_pixelFormat == pixelFormat &&
        m_componentFormat == componentFormat &&
        m_gamma == gamma && multiplier == 1.0f) {
        /* There is nothing to do -- return the current instance */
        return this;
    }

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(m_componentFormat, componentFormat)
    );

    Assert(cvt != NULL);

    ref<Bitmap> target = new Bitmap(pixelFormat, componentFormat, m_size,
            m_channelCount);
    target->setMetadata(m_metadata);
    if (m_channelNames.size() == (size_t) target->getChannelCount())
        target->setChannelNames(m_channelNames);
    target->setGamma(gamma);

    cvt->convert(m_pixelFormat, m_gamma, m_data,
        pixelFormat, gamma, target->getData(),
        (size_t) m_size.x * (size_t) m_size.y, multiplier, intent,
        m_channelCount);

    return target;
}

ref<Bitmap> Bitmap::convertMultiSpectrumAlphaWeight(const std::vector<EPixelFormat> &pixelFormats,
        EComponentFormat componentFormat, const std::vector<std::string> &channelNames) const {
    if (channelNames.size() > std::numeric_limits<uint8_t>::max())
        Log(EError, "convertMultiSpectrumAlphaWeight(): excessive number of channels!");
    ref<Bitmap> bitmap = new Bitmap(Bitmap::EMultiChannel, componentFormat,
            m_size, (uint8_t) channelNames.size());
    bitmap->setChannelNames(channelNames);
    convertMultiSpectrumAlphaWeight(this, getUInt8Data(), bitmap,
        bitmap->getUInt8Data(), pixelFormats, componentFormat,
        (size_t) m_size.x * (size_t) m_size.y);
    return bitmap;
}

void Bitmap::convertMultiSpectrumAlphaWeight(const Bitmap *source,
        const uint8_t *sourcePtr, const Bitmap *target, uint8_t *targetPtr,
        const std::vector<EPixelFormat> &pixelFormats,
        EComponentFormat componentFormat, size_t count) {
    if (source->getComponentFormat() != EFloat && source->getPixelFormat() != EMultiSpectrumAlphaWeight)
        Log(EError, "convertMultiSpectrumAlphaWeight(): unsupported!");

    Float *temp = new Float[count * target->getChannelCount()], *dst = temp;

    for (size_t k = 0; k<count; ++k) {
        const Float *srcData = (const Float *) sourcePtr + k * source->getChannelCount();
        Float weight = srcData[source->getChannelCount()-1],
              invWeight = weight == 0 ? 0 : (Float) 1 / weight;
        Float alpha = srcData[source->getChannelCount()-2] * invWeight;

        for (size_t i=0; i<pixelFormats.size(); ++i) {
            Spectrum value = ((Spectrum *) srcData)[i] * invWeight;
            Float tmp0, tmp1, tmp2;
            switch (pixelFormats[i]) {
                case Bitmap::ELuminance:
                    *dst++ = value.getLuminance();
                    break;
                case Bitmap::ELuminanceAlpha:
                    *dst++ = value.getLuminance();
                    *dst++ = alpha;
                    break;
                case Bitmap::EXYZ:
                    value.toXYZ(tmp0, tmp1, tmp2);
                    *dst++ = tmp0;
                    *dst++ = tmp1;
                    *dst++ = tmp2;
                    break;
                case Bitmap::EXYZA:
                    value.toXYZ(tmp0, tmp1, tmp2);
                    *dst++ = tmp0;
                    *dst++ = tmp1;
                    *dst++ = tmp2;
                    *dst++ = alpha;
                    break;
                case Bitmap::ERGB:
                    value.toLinearRGB(tmp0, tmp1, tmp2);
                    *dst++ = tmp0;
                    *dst++ = tmp1;
                    *dst++ = tmp2;
                    break;
                case Bitmap::ERGBA:
                    value.toLinearRGB(tmp0, tmp1, tmp2);
                    *dst++ = tmp0;
                    *dst++ = tmp1;
                    *dst++ = tmp2;
                    *dst++ = alpha;
                    break;
                case Bitmap::ESpectrum:
                    for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                        *dst++ = value[j];
                    break;
                case Bitmap::ESpectrumAlpha:
                    for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                        *dst++ = value[j];
                    *dst++ = alpha;
                    break;
                default:
                    Log(EError, "Unknown pixel format!");
            }
        }
    }

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(EFloat, target->getComponentFormat())
    );

    cvt->convert(Bitmap::EMultiChannel, 1.0f, temp, Bitmap::EMultiChannel, 1.0f, targetPtr,
            count, 1.0f, Spectrum::EReflectance, target->getChannelCount());

    delete[] temp;
}

void Bitmap::convert(void *target, EPixelFormat pixelFormat,
        EComponentFormat componentFormat, Float gamma, Float multiplier,
        Spectrum::EConversionIntent intent) const {
    if (m_componentFormat == EBitmask || componentFormat == EBitmask)
        Log(EError, "Conversions involving bitmasks are currently not supported!");

    if (m_pixelFormat == pixelFormat &&
        m_componentFormat == componentFormat &&
        m_gamma == gamma && multiplier == 1.0f) {
        /* No conversion is necessary -- just run memcpy */
        memcpy(target, getData(), getBufferSize());
        return;
    }

    const FormatConverter *cvt = FormatConverter::getInstance(
        std::make_pair(m_componentFormat, componentFormat)
    );

    Assert(cvt != NULL);

    cvt->convert(m_pixelFormat, m_gamma, m_data,
        pixelFormat, gamma, target,
        (size_t) m_size.x * (size_t) m_size.y, multiplier, intent,
        m_channelCount);
}

template <typename T> void tonemapReinhard(T *data, size_t pixels, Bitmap::EPixelFormat fmt,
        Float &logAvgLuminance, Float &maxLuminance, Float key, Float burn) {
    int channels = 0;

    switch (fmt) {
        case Bitmap::ERGB:
        case Bitmap::EXYZ:
            channels = 3;
            break;
        case Bitmap::ERGBA:
        case Bitmap::EXYZA:
            channels = 4;
            break;
        case Bitmap::ELuminanceAlpha:
            channels = 2;
            break;
        case Bitmap::ELuminance:
            channels = 1;
            break;
        default:
            SLog(EError, "Unsupported pixel format!");
    }

    if (logAvgLuminance <= 0 || maxLuminance <= 0) {
        /* Compute the log-average luminance if it has not already been provided */
        T *ptr = data;

        maxLuminance = 0;
        logAvgLuminance = 0;

        if (fmt == Bitmap::ERGB || fmt == Bitmap::ERGBA) {
            /* RGB[A] version */
            for (size_t i=0; i < pixels; ++i) {
                Float luminance = (Float) (ptr[0] * (Float) 0.212671 + ptr[1] * (Float) 0.715160 + ptr[2] * (Float) 0.072169);
                if (luminance == 1024) // ignore the "rendered by mitsuba banner.."
                    maxLuminance = 0.0f;
                maxLuminance = std::max(maxLuminance, luminance);
                logAvgLuminance += math::fastlog(1e-3f + luminance);
                ptr += channels;
            }
        } else if (fmt == Bitmap::EXYZ || fmt == Bitmap::EXYZA) {
            for (size_t i=0; i < pixels; ++i) {
                Float luminance = (Float) ptr[1];
                if (luminance == 1024) // ignore the "rendered by mitsuba banner.."
                    maxLuminance = 0.0f;
                maxLuminance = std::max(maxLuminance, luminance);
                logAvgLuminance += math::fastlog(1e-3f + luminance);
                ptr += channels;
            }
        } else {
            /* Monochrome version */
            for (size_t i=0; i < pixels; ++i) {
                Float luminance = (Float) *ptr;
                if (luminance == 1024) // ignore the "rendered by mitsuba banner.."
                    maxLuminance = 0.0f;
                maxLuminance = std::max(maxLuminance, luminance);
                logAvgLuminance += math::fastlog(1e-3f + luminance);
                ptr += channels;
            }
        }

        logAvgLuminance = math::fastexp(logAvgLuminance / pixels);
    }

    if (maxLuminance == 0) /* This is a black image -- stop now */
        return;

    burn = std::min((Float) 1, std::max((Float) 1e-8f, 1-burn));

    Float scale = key / logAvgLuminance,
          Lwhite = maxLuminance * scale;

    /* Having the 'burn' parameter scale as 1/b^4 provides a nicely behaved knob */
    Float invWp2 = 1 / (Lwhite * Lwhite * std::pow(burn, (Float) 4));

    if (fmt == Bitmap::ERGB || fmt == Bitmap::ERGBA) {
        /* RGB[A] version */
        for (size_t i=0; i < pixels; ++i) {
            /* Convert ITU-R Rec. BT.709 linear RGB to XYZ tristimulus values */
            Float X = static_cast<Float>(data[0] * 0.412453f + data[1] * 0.357580f + data[2] * 0.180423f);
            Float Y = static_cast<Float>(data[0] * 0.212671f + data[1] * 0.715160f + data[2] * 0.072169f);
            Float Z = static_cast<Float>(data[0] * 0.019334f + data[1] * 0.119193f + data[2] * 0.950227f);

            /* Convert to xyY */
            Float normalization = 1 / (X + Y + Z),
                x  = X * normalization,
                y  = Y * normalization,
                Lp = Y * scale;

            /* Apply the tonemapping transformation */
            Y = Lp * (1.0f + Lp*invWp2) / (1.0f + Lp);

            /* Convert back to XYZ */
            Float ratio = Y/y;
            X = ratio * x;
            Z = ratio * ((Float) 1.0f - x - y);

            /* Convert from XYZ tristimulus values to ITU-R Rec. BT.709 linear RGB */
            data[0] = safe_cast<T>(  3.240479f * X + -1.537150f * Y + -0.498535f * Z);
            data[1] = safe_cast<T>( -0.969256f * X +  1.875991f * Y +  0.041556f * Z);
            data[2] = safe_cast<T>(  0.055648f * X + -0.204043f * Y +  1.057311f * Z);

            data += channels;
        }
    } else if (fmt == Bitmap::EXYZ || fmt == Bitmap::EXYZA) {
        /* XYZ[A] version */
        for (size_t i=0; i < pixels; ++i) {
            Float X = static_cast<Float>(data[0]),
                  Y = static_cast<Float>(data[1]),
                  Z = static_cast<Float>(data[2]);

            /* Convert to xyY */
            Float normalization = 1 / (X + Y + Z),
                x  = X * normalization,
                y  = Y * normalization,
                Lp = Y * scale;

            /* Apply the tonemapping transformation */
            Y = Lp * (1.0f + Lp*invWp2) / (1.0f + Lp);

            /* Convert back to XYZ */
            Float ratio = Y/y;
            X = ratio * x;
            Z = ratio * ((Float) 1.0f - x - y);

            data[0] = safe_cast<T>(X);
            data[1] = safe_cast<T>(Y);
            data[2] = safe_cast<T>(Z);

            data += channels;
        }

    } else {
        /* Monochrome version */
        for (size_t i=0; i < pixels; ++i) {
            Float Lp = (Float) *data * scale;

            /* Apply the tonemapping transformation */
            *data = safe_cast<T> (Lp * (1.0f + Lp*invWp2) / (1.0f + Lp));

            data += channels;
        }
    }
}

void Bitmap::tonemapReinhard(Float &logAvgLuminance, Float &maxLuminance, Float key, Float burn) {
    Assert(m_pixelFormat == ERGB || m_pixelFormat == ERGBA ||
           m_pixelFormat == ELuminance || m_pixelFormat == ELuminanceAlpha);
    Assert(m_gamma == 1);

    size_t pixels = (size_t) m_size.x * (size_t) m_size.y;

    switch (m_componentFormat) {
        case EFloat16:
            mitsuba::tonemapReinhard(getFloat16Data(), pixels, m_pixelFormat, logAvgLuminance, maxLuminance, key, burn);
            break;
        case EFloat32:
            mitsuba::tonemapReinhard(getFloat32Data(), pixels, m_pixelFormat, logAvgLuminance, maxLuminance, key, burn);
            break;
        case EFloat64:
            mitsuba::tonemapReinhard(getFloat64Data(), pixels, m_pixelFormat, logAvgLuminance, maxLuminance, key, burn);
            break;
        default:
            Log(EError, "Bitmap::tonemapReinhard(): Unsupported component format!");
    }
}

std::vector<Bitmap::Layer> Bitmap::getLayers() const {
    typedef std::map<std::string, int> ChannelMap;
    if (m_channelNames.empty())
        Log(EError, "Bitmap::getLayers(): required color channel names were not available!");

    ChannelMap channels;
    for (size_t i=0; i<m_channelNames.size(); ++i)
        channels[boost::to_lower_copy(m_channelNames[i])] = (int) i;

    std::vector<Layer> layers;
    for (size_t i=0; i<m_channelNames.size(); ++i) {
        std::string name = boost::to_lower_copy(m_channelNames[i]);
        if (channels.find(name) == channels.end())
            continue;
        std::string prefix = name;
        char postfix = '\0';
        Layer layer;

        layer.name = m_channelNames[i];
        if (name.length() == 1) {
            prefix = layer.name = "";
            postfix = name[name.length() - 1];
        } if (name.length() >= 3 && name[name.length()-2] == '.') {
            prefix = name.substr(0, name.length() - 1);
            layer.name = layer.name.substr(0, layer.name.length() - 2);
            postfix = name[name.length() - 1];
        }

        ChannelMap::iterator
            itR = channels.find(prefix + "r"),
            itG = channels.find(prefix + "g"),
            itB = channels.find(prefix + "b"),
            itA = channels.find(prefix + "a"),
            itX = channels.find(prefix + "x"),
            itY = channels.find(prefix + "y"),
            itZ = channels.find(prefix + "z");

        bool maybeRGB = postfix == 'r' || postfix == 'g' || postfix == 'b' || postfix == 'a';
        bool maybeXYZ = postfix == 'x' || postfix == 'y' || postfix == 'z' || postfix == 'a';
        bool maybeY = postfix == 'y' || postfix == 'a';

        if (maybeRGB && itR != channels.end() && itG != channels.end() && itB != channels.end()) {
            layer.format = ERGB;
            layer.channels.push_back(itR->second);
            layer.channels.push_back(itG->second);
            layer.channels.push_back(itB->second);
            if (itA != channels.end()) {
                layer.format = ERGBA;
                layer.channels.push_back(itA->second);
                channels.erase(prefix + "a");
            }
            channels.erase(prefix + "r");
            channels.erase(prefix + "g");
            channels.erase(prefix + "b");
        } else if (maybeXYZ && itX != channels.end() && itY != channels.end() && itZ != channels.end()) {
            layer.format = EXYZ;
            layer.channels.push_back(itX->second);
            layer.channels.push_back(itY->second);
            layer.channels.push_back(itZ->second);
            if (itA != channels.end()) {
                layer.format = EXYZA;
                layer.channels.push_back(itA->second);
                channels.erase(prefix + "a");
            }
            channels.erase(prefix + "x");
            channels.erase(prefix + "y");
            channels.erase(prefix + "z");
        } else if (maybeY && itY != channels.end()) {
            layer.format = ELuminance;
            layer.channels.push_back(itY->second);
            if (itA != channels.end()) {
                layer.format = ELuminanceAlpha;
                layer.channels.push_back(itA->second);
                channels.erase(prefix + "a");
            }
            channels.erase(prefix + "y");
        } else {
            if (layer.name.empty())
                layer.name = m_channelNames[i];
            layer.format = ELuminance;
            layer.channels.push_back((int) i);
            channels.erase(name);
        }

        layers.push_back(layer);
    }
    return layers;
}

std::map<std::string, Bitmap *> Bitmap::split() const {
    std::map<std::string, Bitmap *> result;

    std::vector<Layer> layers = getLayers();
    for (size_t i=0; i<layers.size(); ++i) {
        const Layer &layer = layers[i];
        std::vector<std::string> channelNames;
        for (size_t j=0; j<layer.channels.size(); ++j)
            channelNames.push_back(m_channelNames[layer.channels[j]]);

        Bitmap *bitmap = NULL;
        {
            ref<Bitmap> _bitmap = this->extractChannels(layer.format, layer.channels);
            bitmap = _bitmap.get();
            bitmap->incRef();
        }
        bitmap->decRef(false);
        bitmap->setChannelNames(channelNames);

        if (result.find(layer.name) != result.end())
            Log(EError, "Internal error -- encountered two layers with the same name \"%s\"", layer.name.c_str());
        result[layer.name] = bitmap;
    }
    return result;
}

ref<Bitmap> Bitmap::extractChannels(EPixelFormat fmt, const std::vector<int> &channels) const {
    int channelCount = getChannelCount();

    for (size_t i=0; i<channels.size(); ++i)
        if (channels[i] < 0 || channels[i] >= channelCount)
            Log(EError, "Bitmap::extractChannel(%i): channel index "
                "must be between 0 and %i", channels[i], channelCount-1);

    ref<Bitmap> result = new Bitmap(fmt, m_componentFormat, m_size, (int) channels.size());
    result->setMetadata(m_metadata);
    result->setGamma(m_gamma);

    size_t componentSize = getBytesPerComponent();
    size_t stride = channelCount * componentSize;
    size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;

    const uint8_t *source = getUInt8Data();
    uint8_t *target = result->getUInt8Data();

    for (size_t px = 0; px<pixelCount; ++px) {
        for (size_t ch = 0; ch<channels.size(); ++ch)
            for (size_t c = 0; c<componentSize; ++c)
                *target++ = (source+channels[ch]*componentSize)[c];
        source += stride;
    }

    return result;
}

ref<Bitmap> Bitmap::extractChannel(int channelIndex) const {
    int channelCount = getChannelCount();

    if (channelIndex == 0 && channelCount == 1)
        return const_cast<Bitmap *>(this);

    if (channelIndex < 0 || channelIndex >= channelCount)
        Log(EError, "Bitmap::extractChannel(%i): channel index "
            "must be between 0 and %i", channelIndex, channelCount-1);

    ref<Bitmap> result = new Bitmap(ELuminance, m_componentFormat, m_size);
    result->setMetadata(m_metadata);
    result->setGamma(m_gamma);

    size_t componentSize = getBytesPerComponent();
    size_t offset = channelIndex * componentSize;
    size_t stride = channelCount * componentSize;
    size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;

    const uint8_t *source = getUInt8Data() + offset;
    uint8_t *target = result->getUInt8Data();

    for (size_t px = 0; px<pixelCount; ++px) {
        for (size_t c = 0; c<componentSize; ++c)
            *target++ = source[c];
        source += stride;
    }

    return result;
}

ref<Bitmap> Bitmap::join(EPixelFormat fmt,
            const std::vector<Bitmap *> &sourceBitmaps) {
    if (sourceBitmaps.size() == 0)
        Log(EError, "Bitmap::join(): Need at least one bitmap!");

    const Bitmap *bitmap0 = sourceBitmaps[0];
    if (bitmap0->getComponentFormat() == EBitmask)
        Log(EError, "Conversions involving bitmasks are currently not supported!");

    for (size_t i=1; i<sourceBitmaps.size(); ++i) {
        if (sourceBitmaps[i]->getSize() != bitmap0->getSize())
            Log(EError, "Bitmap::join(): Detected a size mismatch!");
        if (sourceBitmaps[i]->getGamma() != bitmap0->getGamma())
            Log(EError, "Bitmap::join(): Detected a gamma mismatch!");
        if (sourceBitmaps[i]->getComponentFormat() != bitmap0->getComponentFormat())
            Log(EError, "Bitmap::join(): Detected a component format mismatch!");
    }

    /* Determine channel names and metadata */
    int sourceChannelCount = 0, channelCount = -1;
    std::vector<std::string> channelNames;
    Properties metadata;
    for (size_t i=0; i<sourceBitmaps.size(); ++i) {
        sourceChannelCount += sourceBitmaps[i]->getChannelCount();
        const std::vector<std::string> &names = sourceBitmaps[i]->getChannelNames();
        channelNames.insert(channelNames.end(), names.begin(), names.end());
        metadata.merge(sourceBitmaps[i]->getMetadata());
    }

    if (!channelNames.empty()) {
        std::set<std::string> namesUnique;
        for (size_t i = 0; i < channelNames.size(); ++i) {
            if (!channelNames[i].empty())
                namesUnique.insert(channelNames[i]);
        }
        if (namesUnique.size() != channelNames.size())
            Log(EError, "Bitmap::join(): Error -- some of the image supplied color "
                " channel names, but it was not possible to assign a unique set of names.");
    }

    if (fmt == Bitmap::EMultiChannel)
        channelCount = sourceChannelCount;

    ref<Bitmap> result = new Bitmap(fmt, bitmap0->getComponentFormat(),
        bitmap0->getSize(), channelCount);

    channelCount = result->getChannelCount();
    if (channelCount != sourceChannelCount)
        Log(EError, "Bitmap::join(): Error -- supplied the wrong number "
            "of channels (%i instead of %i)", (int) sourceBitmaps.size(),
            result->getChannelCount());

    result->setMetadata(metadata);
    result->setGamma(bitmap0->getGamma());
    result->setChannelNames(channelNames);

    size_t pixelCount = (size_t) bitmap0->getSize().x * (size_t) bitmap0->getSize().y;
    uint8_t **pointers = (uint8_t **) alloca(sourceBitmaps.size() * sizeof(uint8_t *));
    size_t *pixelSize = (size_t *) alloca(sourceBitmaps.size() * sizeof(size_t));

    for (size_t i = 0; i<sourceBitmaps.size(); ++i) {
        pointers[i] = sourceBitmaps[i]->getUInt8Data();
        pixelSize[i] = sourceBitmaps[i]->getBytesPerPixel();
    }

    uint8_t *dest = result->getUInt8Data();

    for (size_t i = 0; i<pixelCount; ++i)
        for (size_t j = 0; j<sourceBitmaps.size(); ++j)
            for (size_t k= 0; k < pixelSize[j]; ++k)
                *dest++ = *pointers[j]++;

    return result;
}

ref<Bitmap> Bitmap::crop(const Point2i &offset, const Vector2i &size) const {
    Assert(offset.x >= 0 && offset.y >= 0 &&
           offset.x + size.x <= m_size.x &&
           offset.y + size.y <= m_size.y);

    size_t pixelStride = getBytesPerPixel();
    size_t sourceStride = pixelStride * m_size.x;
    size_t targetStride = pixelStride * size.x;

    ref<Bitmap> result = new Bitmap(m_pixelFormat, m_componentFormat,
            size, m_channelCount);

    result->setGamma(m_gamma);
    result->setChannelNames(m_channelNames);
    result->setMetadata(m_metadata);

    uint8_t *source = m_data + (offset.x + offset.y * m_size.x) * pixelStride;
    uint8_t *target = result->getUInt8Data();

    for (int y=0; y<size.y; ++y) {
        memcpy(target, source, targetStride);

        source += sourceStride;
        target += targetStride;
    }

    return result;
}

void Bitmap::applyMatrix(Float matrix_[3][3]) {
    int stride = 0;

    if (m_pixelFormat == ERGB || m_pixelFormat == EXYZ)
        stride = 3;
    else if (m_pixelFormat == ERGBA || m_pixelFormat == EXYZA)
        stride = 4;
    else
        Log(EError, "Bitmap::applyMatrix(): unsupported pixel format!");

    size_t pixels = (size_t) m_size.x * (size_t) m_size.y;

    switch (m_componentFormat) {
        case EFloat16: {
            float matrix[3][3];
            half *data = getFloat16Data();
            for (int i=0; i<3; ++i)
                for (int j=0; j<3; ++j)
                    matrix[i][j] = (float) matrix_[i][j];

            for (size_t i=0; i<pixels; ++i) {
                float result[3] = { 0.0f, 0.0f, 0.0f };
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result[i] += matrix[i][j] * (float) data[j];
                for (int i=0; i<3; ++i)
                    data[i] = (half) result[i];
                data += stride;
            }
        }
        break;

        case EFloat32: {
            float matrix[3][3], *data = getFloat32Data();
            for (int i=0; i<3; ++i)
                for (int j=0; j<3; ++j)
                    matrix[i][j] = (float) matrix_[i][j];

            for (size_t i=0; i<pixels; ++i) {
                float result[3] = { 0.0f, 0.0f, 0.0f };
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result[i] += matrix[i][j] * data[j];
                for (int i=0; i<3; ++i)
                    data[i] = result[i];
                data += stride;
            }
        }
        break;

        case EFloat64: {
            double matrix[3][3], *data = getFloat64Data();
            for (int i=0; i<3; ++i)
                for (int j=0; j<3; ++j)
                    matrix[i][j] = (double) matrix_[i][j];

            for (size_t i=0; i<pixels; ++i) {
                double result[3] = { 0.0, 0.0, 0.0 };
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        result[i] += matrix[i][j] * data[j];
                for (int i=0; i<3; ++i)
                    data[i] = result[i];
                data += stride;
            }
        }
        break;
        default:
            Log(EError, "Bitmap::applyMatrix(): unsupported component format!");
    }
}

/// Bitmap filtering / resampling utility function
template <typename Scalar> static void resample(ref<const ReconstructionFilter> rfilter,
    ReconstructionFilter::EBoundaryCondition bch,
    ReconstructionFilter::EBoundaryCondition bcv,
    const Bitmap *source, Bitmap *target, ref<Bitmap> temp, Float minValue,
    Float maxValue, bool filter) {

    if (!rfilter) {
        /* Resample using a 2-lobed Lanczos reconstruction filter */
        Properties rfilterProps("lanczos");
        rfilterProps.setInteger("lobes", 2);
        ReconstructionFilter *instance = static_cast<ReconstructionFilter *> (
            PluginManager::getInstance()->createObject(
            MTS_CLASS(ReconstructionFilter), rfilterProps));
        instance->configure();
        rfilter = instance;
    }

    if (source->getHeight() == target->getHeight() &&
        source->getWidth() == target->getWidth() && !filter) {
        memcpy(target->getData(), source->getData(), source->getBufferSize());
        return;
    }

    int channels = source->getChannelCount();
    bool clamp =
        minValue != -std::numeric_limits<Float>::infinity() ||
        maxValue !=  std::numeric_limits<Float>::infinity();

    if (source->getWidth() != target->getWidth() || filter) {
        /* Re-sample along the X direction */
        Resampler<Scalar> r(rfilter, bch, source->getWidth(), target->getWidth());

        /* Create a bitmap for intermediate storage */
        if (!temp) {
            if (source->getHeight() == target->getHeight() && !filter)
                temp = target; // write directly to the output bitmap
            else // otherwise: write to a temporary bitmap
                temp = new Bitmap(source->getPixelFormat(), source->getComponentFormat(),
                    Vector2i(target->getWidth(), source->getHeight()), channels);
        }

        if (clamp) {
            #if defined(MTS_OPENMP)
                #pragma omp parallel for
            #endif
            for (int y=0; y<source->getHeight(); ++y) {
                const Scalar *srcPtr = (Scalar *) source->getUInt8Data()
                    + y * source->getWidth() * channels;
                Scalar *trgPtr = (Scalar *) temp->getUInt8Data()
                    + y * target->getWidth() * channels;

                r.resampleAndClamp(srcPtr, 1, trgPtr, 1, channels,
                        safe_cast<Scalar>(minValue), safe_cast<Scalar>(maxValue));
            }
        } else {
            #if defined(MTS_OPENMP)
                #pragma omp parallel for
            #endif
            for (int y=0; y<source->getHeight(); ++y) {
                const Scalar *srcPtr = (Scalar *) source->getUInt8Data()
                    + y * source->getWidth() * channels;
                Scalar *trgPtr = (Scalar *) temp->getUInt8Data()
                    + y * target->getWidth() * channels;

                r.resample(srcPtr, 1, trgPtr, 1, channels);
            }
        }

        /* Now, read from the temporary bitmap */
        source = temp;
    }

    if (source->getHeight() != target->getHeight() || filter) {
        /* Re-sample along the Y direction */
        Resampler<Scalar> r(rfilter, bcv, source->getHeight(), target->getHeight());

        if (clamp) {
            #if defined(MTS_OPENMP)
                #pragma omp parallel for
            #endif
            for (int x=0; x<source->getWidth(); ++x) {
                const Scalar *srcPtr = (Scalar *) source->getUInt8Data() + x * channels;
                Scalar *trgPtr = (Scalar *) target->getUInt8Data() + x * channels;

                r.resampleAndClamp(srcPtr, source->getWidth(), trgPtr, target->getWidth(),
                    channels, safe_cast<Scalar>(minValue), safe_cast<Scalar>(maxValue));
            }
        } else {
            #if defined(MTS_OPENMP)
                #pragma omp parallel for
            #endif
            for (int x=0; x<source->getWidth(); ++x) {
                const Scalar *srcPtr = (Scalar *) source->getUInt8Data() + x * channels;
                Scalar *trgPtr = (Scalar *) target->getUInt8Data() + x * channels;

                r.resample(srcPtr, source->getWidth(), trgPtr, target->getWidth(), channels);
            }
        }
    }
}

void Bitmap::resample(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        Bitmap *target, Bitmap *temp, Float minValue, Float maxValue) const {

    Assert(getPixelFormat() == target->getPixelFormat() &&
        getComponentFormat() == target->getComponentFormat() &&
        getChannelCount() == target->getChannelCount() &&
        (!temp || temp->getSize() == Vector2i(target->getWidth(), getHeight())));


    switch (m_componentFormat) {
        case EFloat16:
            mitsuba::resample<half>(rfilter, bch, bcv, this, target, temp, minValue, maxValue, false);
            break;
        case EFloat32:
            mitsuba::resample<float>(rfilter, bch, bcv, this, target, temp, minValue, maxValue, false);
            break;
        case EFloat64:
            mitsuba::resample<double>(rfilter, bch, bcv, this, target, temp, minValue, maxValue, false);
            break;
        default:
            Log(EError, "resample(): Unsupported component type! (must be float16/32/64)");
    }
}

void Bitmap::filter(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        Bitmap *target, Bitmap *temp, Float minValue, Float maxValue) const {

    Assert(getPixelFormat() == target->getPixelFormat() &&
        getComponentFormat() == target->getComponentFormat() &&
        getChannelCount() == target->getChannelCount() &&
        getSize() == target->getSize() &&
        (!temp || temp->getSize() == getSize()));

    switch (m_componentFormat) {
        case EFloat16:
            mitsuba::resample<half>(rfilter, bch, bcv, this, target, temp, minValue, maxValue, true);
            break;
        case EFloat32:
            mitsuba::resample<float>(rfilter, bch, bcv, this, target, temp, minValue, maxValue, true);
            break;
        case EFloat64:
            mitsuba::resample<double>(rfilter, bch, bcv, this, target, temp, minValue, maxValue, true);
            break;
        default:
            Log(EError, "filter(): Unsupported component type! (must be float16/32/64)");
    }
}

ref<Bitmap> Bitmap::resample(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        const Vector2i &size, Float minValue, Float maxValue) const {
    ref<Bitmap> result = new Bitmap(m_pixelFormat, m_componentFormat, size);
    result->m_metadata = m_metadata;
    result->m_gamma = m_gamma;
    result->m_channelNames = m_channelNames;
    resample(rfilter, bch, bcv, result, NULL, minValue, maxValue);
    return result;
}

ref<Bitmap> Bitmap::filter(const ReconstructionFilter *rfilter,
        ReconstructionFilter::EBoundaryCondition bch,
        ReconstructionFilter::EBoundaryCondition bcv,
        Float minValue, Float maxValue) const {
    ref<Bitmap> result = new Bitmap(m_pixelFormat, m_componentFormat, getSize());
    result->m_metadata = m_metadata;
    result->m_gamma = m_gamma;
    result->m_channelNames = m_channelNames;
    filter(rfilter, bch, bcv, result, NULL, minValue, maxValue);
    return result;
}

bool Bitmap::operator==(const Bitmap &bitmap) const {
    return m_pixelFormat == bitmap.m_pixelFormat &&
        m_componentFormat == bitmap.m_componentFormat &&
        m_size == bitmap.m_size &&
        m_metadata == bitmap.m_metadata &&
        m_gamma == bitmap.m_gamma &&
        memcmp(bitmap.m_data, m_data, getBufferSize()) == 0;
}

std::string Bitmap::toString() const {
    std::ostringstream oss;
    oss << "Bitmap[" << endl
        << "  type = " << m_pixelFormat << "," << endl
        << "  componentFormat = " << m_componentFormat << "," << endl
        << "  size = " << m_size.toString() << "," << endl;

    if (isMultiChannel())
        oss << "  channelCount = " << (int) m_channelCount << "," << endl;

    if (!m_channelNames.empty()) {
        oss << "  channelNames = [ ";
        for (size_t i=0; i<m_channelNames.size(); ++i) {
            oss << "\"" << m_channelNames[i] << "\"";
            if (i+1<m_channelNames.size())
                oss << ", ";
        }
        oss << " ]," << endl;
    }

    std::vector<std::string> keys = m_metadata.getPropertyNames();
    if (!keys.empty()) {
        oss << "  metadata = {" << endl;
        for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ) {
            std::string value = m_metadata.getAsString(*it);
            if (value.size() > 50)
                value = value.substr(0, 50) + ".. [truncated]";

            oss << "    \"" << *it << "\" => \"" << value << "\"";
            if (++it != keys.end())
                oss << ",";
            oss << endl;
        }
        oss << "  }," << endl;
    }
    oss << "  gamma = " << m_gamma << "," << endl
        << "  data = [ " << memString(getBufferSize()) << " of image data ]" << endl
        << "]";
    return oss.str();
}

#if defined(MTS_HAS_LIBPNG)
void Bitmap::readPNG(Stream *stream) {
    png_structp png_ptr;
    png_infop info_ptr;
    volatile png_bytepp rows = NULL;

    /* Create buffers */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, &png_error_func, &png_warn_func);
    if (png_ptr == NULL) {
        Log(EError, "readPNG(): Unable to create PNG data structure");
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        png_destroy_read_struct(&png_ptr, (png_infopp) NULL, (png_infopp) NULL);
        Log(EError, "readPNG(): Unable to create PNG information structure");
    }

    /* Error handling */
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);
        if (rows)
            delete[] rows;
        Log(EError, "readPNG(): Error reading the PNG file!");
    }

    /* Set read helper function */
    png_set_read_fn(png_ptr, stream, (png_rw_ptr) png_read_data);

    int bitDepth, colorType, interlacetype, compressiontype, filtertype;
    png_read_info(png_ptr, info_ptr);
    png_uint_32 width = 0, height = 0;
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bitDepth,
            &colorType, &interlacetype, &compressiontype, &filtertype);

    /* Request various transformations from libpng as necessary */
    if (colorType == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png_ptr); // Always expand indexed files
    else if (colorType == PNG_COLOR_TYPE_GRAY && bitDepth > 1 && bitDepth < 8)
        png_set_expand_gray_1_2_4_to_8(png_ptr); // Expand 2- and 4-bit grayscale
    else if (bitDepth == 16 && Stream::getHostByteOrder() == Stream::ELittleEndian)
        png_set_swap(png_ptr); // Swap the byte order on little endian machines

    // Expand transparency
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png_ptr);

    /* Update the information based on the transformations */
    png_read_update_info(png_ptr, info_ptr);
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bitDepth,
        &colorType, &interlacetype, &compressiontype, &filtertype);
    m_size = Vector2i(width, height);

    switch (colorType) {
        case PNG_COLOR_TYPE_GRAY: m_pixelFormat = ELuminance; break;
        case PNG_COLOR_TYPE_GRAY_ALPHA: m_pixelFormat = ELuminanceAlpha; break;
        case PNG_COLOR_TYPE_RGB: m_pixelFormat = ERGB; break;
        case PNG_COLOR_TYPE_RGB_ALPHA: m_pixelFormat = ERGBA; break;
        default: Log(EError, "readPNG(): Unknown color type %i", colorType); break;
    }
    updateChannelCount();

    switch (bitDepth) {
        case 1: m_componentFormat = EBitmask; break;
        case 8: m_componentFormat = EUInt8; break;
        case 16: m_componentFormat = EUInt16; break;
        default: Log(EError, "readPNG(): Unsupported bit depth: %i", bitDepth);
    }

    /* Load any string-valued metadata */
    int textIdx = 0;
    png_textp text_ptr;
    png_get_text(png_ptr, info_ptr, &text_ptr, &textIdx);

    for (int i=0; i<textIdx; ++i, text_ptr++)
        setMetadataString(text_ptr->key, text_ptr->text);

    int intent; double gamma;
    if (png_get_sRGB(png_ptr, info_ptr, &intent)) {
        m_gamma = -1;
    } else if (png_get_gAMA(png_ptr, info_ptr, &gamma)) {
        m_gamma = (Float) 1 / (Float) gamma;
    } else {
        m_gamma = -1; // assume sRGB by default
    }

    Log(ETrace, "Loading a %ix%i PNG file", width, height);

    size_t bufferSize = getBufferSize();
    m_data = static_cast<uint8_t *>(allocAligned(bufferSize));
    m_ownsData = true;
    rows = new png_bytep[m_size.y];
    size_t rowBytes = png_get_rowbytes(png_ptr, info_ptr);
    Assert(rowBytes == getBufferSize() / m_size.y);

    for (int i=0; i<m_size.y; i++)
        rows[i] = m_data + i * rowBytes;

    png_read_image(png_ptr, rows);
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);

    delete[] rows;
}

void Bitmap::writePNG(Stream *stream, int compression) const {
    png_structp png_ptr;
    png_infop info_ptr;
    volatile png_bytepp rows = NULL;

    Log(EDebug, "Writing a %ix%i PNG file", m_size.x, m_size.y);

    int colorType, bitDepth;
    switch (m_pixelFormat) {
        case ELuminance: colorType = PNG_COLOR_TYPE_GRAY; break;
        case ELuminanceAlpha: colorType = PNG_COLOR_TYPE_GRAY_ALPHA; break;
        case ERGB: colorType = PNG_COLOR_TYPE_RGB; break;
        case ERGBA: colorType = PNG_COLOR_TYPE_RGBA; break;
        default:
            Log(EError, "writePNG(): Unsupported bitmap type!");
            return;
    }

    switch (m_componentFormat) {
        case EBitmask: bitDepth = 1; break;
        case EUInt8: bitDepth = 8; break;
        case EUInt16: bitDepth = 16; break;
        default:
            Log(EError, "writePNG(): Unsupported component type!");
            return;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, &png_error_func, &png_warn_func);
    if (png_ptr == NULL)
        Log(EError, "Error while creating PNG data structure");

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        Log(EError, "Error while creating PNG information structure");
    }

    /* Error handling */
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        Log(EError, "Error writing the PNG file");
    }

    png_set_write_fn(png_ptr, stream, (png_rw_ptr) png_write_data, (png_flush_ptr) png_flush_data);
    png_set_compression_level(png_ptr, compression);

    png_text *text = NULL;

    Properties metadata(m_metadata);
    if (!metadata.hasProperty("generatedBy"))
        metadata.setString("generatedBy", "Mitsuba version " MTS_VERSION);

    std::vector<std::string> keys = metadata.getPropertyNames();
    std::vector<std::string> values(keys.size());

    text = new png_text[keys.size()];
    memset(text, 0, sizeof(png_text) * keys.size());

    for (size_t i = 0; i<keys.size(); ++i) {
        values[i] = metadata.getAsString(keys[i]);
        text[i].key = const_cast<char *>(keys[i].c_str());
        text[i].text = const_cast<char *>(values[i].c_str());
        text[i].compression = PNG_TEXT_COMPRESSION_NONE;
    }

    png_set_text(png_ptr, info_ptr, text, (int) keys.size());

    if (m_gamma == -1)
        png_set_sRGB_gAMA_and_cHRM(png_ptr, info_ptr, PNG_sRGB_INTENT_ABSOLUTE);
    else
        png_set_gAMA(png_ptr, info_ptr, 1 / m_gamma);

    if (m_componentFormat == EUInt16 && Stream::getHostByteOrder() == Stream::ELittleEndian)
        png_set_swap(png_ptr); // Swap the byte order on little endian machines

    png_set_IHDR(png_ptr, info_ptr, m_size.x, m_size.y, bitDepth,
            colorType, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
            PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    rows = new png_bytep[m_size.y];

    size_t rowBytes = png_get_rowbytes(png_ptr, info_ptr);
    Assert(rowBytes == getBufferSize() / m_size.y);
    for (int i=0; i<m_size.y; i++)
        rows[i] = &m_data[rowBytes * i];

    png_write_image(png_ptr, rows);
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    if (text)
        delete[] text;
    delete[] rows;
}
#else
void Bitmap::readPNG(Stream *stream) {
    Log(EError, "Bitmap::readPNG(): libpng support was disabled at compile time!");
}
void Bitmap::writePNG(Stream *stream, int compression) const {
    Log(EError, "Bitmap::writePNG(): libpng support was disabled at compile time!");
}
#endif

#if defined(MTS_HAS_LIBJPEG)
void Bitmap::readJPEG(Stream *stream) {
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    jbuf_in_t jbuf;

    memset(&jbuf, 0, sizeof(jbuf_in_t));

    cinfo.err = jpeg_std_error(&jerr);
    jerr.error_exit = jpeg_error_exit;
    jpeg_create_decompress(&cinfo);
    cinfo.src = (struct jpeg_source_mgr *) &jbuf;
    jbuf.mgr.init_source = jpeg_init_source;
    jbuf.mgr.fill_input_buffer = jpeg_fill_input_buffer;
    jbuf.mgr.skip_input_data = jpeg_skip_input_data;
    jbuf.mgr.term_source = jpeg_term_source;
    jbuf.mgr.resync_to_restart = jpeg_resync_to_restart;
    jbuf.stream = stream;

    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);

    m_size = Vector2i(cinfo.output_width, cinfo.output_height);
    m_componentFormat = EUInt8;
    m_gamma = -1;

    switch (cinfo.output_components) {
        case 1: m_pixelFormat = ELuminance; break;
        case 3: m_pixelFormat = ERGB; break;
        default: Log(EError, "readJPEG(): Unsupported number of components!");
    }
    updateChannelCount();

    Log(ETrace, "Loading a %ix%i JPG file", m_size.x, m_size.y);

    size_t row_stride = (size_t) cinfo.output_width
        * (size_t) cinfo.output_components;

    m_data = static_cast<uint8_t *>(allocAligned(getBufferSize()));
    m_ownsData = true;

    boost::scoped_array<uint8_t*> scanlines(new uint8_t*[m_size.y]);
    for (int i=0; i<m_size.y; ++i)
        scanlines.get()[i] = m_data + row_stride*i;

    /* Process scanline by scanline */
    int counter = 0;
    while (cinfo.output_scanline < cinfo.output_height)
        counter += jpeg_read_scanlines(&cinfo, scanlines.get() + counter,
            m_size.y - cinfo.output_scanline);

    /* Release the libjpeg data structures */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
}

void Bitmap::writeJPEG(Stream *stream, int quality) const {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    jbuf_out_t jbuf;

    int components = 0;
    if (m_pixelFormat == ELuminance)
        components = 1;
    else if (m_pixelFormat == ERGB)
        components = 3;
    else
        Log(EError, "writeJPEG(): Invalid pixel format!");

    if (m_componentFormat != EUInt8)
        Log(EError, "writeJPEG(): Invalid component format!");

    memset(&jbuf, 0, sizeof(jbuf_out_t));
    cinfo.err = jpeg_std_error(&jerr);
    jerr.error_exit = jpeg_error_exit;
    jpeg_create_compress(&cinfo);

    cinfo.dest = (struct jpeg_destination_mgr *) &jbuf;
    jbuf.mgr.init_destination = jpeg_init_destination;
    jbuf.mgr.empty_output_buffer = jpeg_empty_output_buffer;
    jbuf.mgr.term_destination = jpeg_term_destination;
    jbuf.stream = stream;

    cinfo.image_width = m_size.x;
    cinfo.image_height = m_size.y;
    cinfo.input_components = components;
    cinfo.in_color_space = components == 1 ? JCS_GRAYSCALE : JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    Log(ETrace, "Writing a %ix%i JPEG file", m_size.x, m_size.y);

    /* Write scanline by scanline */
    for (int i=0; i<m_size.y; ++i) {
        uint8_t *source = m_data + i*m_size.x*cinfo.input_components;
        jpeg_write_scanlines(&cinfo, &source, 1);
    }

    /* Release the libjpeg data structures */
    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
}
#else
void Bitmap::readJPEG(Stream *stream) {
    Log(EError, "Bitmap::readJPEG(): libjpeg support was disabled at compile time!");
}
void Bitmap::writeJPEG(Stream *stream, int quality) const {
    Log(EError, "Bitmap::writeJPEG(): libjpeg support was disabled at compile time!");
}
#endif

#if defined(MTS_HAS_OPENEXR)
void Bitmap::readOpenEXR(Stream *stream, const std::string &_prefix) {
    EXRIStream istr(stream);
    Imf::InputFile file(istr);

    const Imf::Header &header = file.header();
    const Imf::ChannelList &channels = header.channels();
    Assert(channels.begin() != channels.end());

    const char *ch_r = NULL, *ch_g = NULL, *ch_b = NULL,
        *ch_a = NULL, *ch_x = NULL, *ch_y = NULL, *ch_z = NULL,
        *ch_ry = NULL, *ch_by = NULL, *ch_spec[SPECTRUM_SAMPLES];

    memset(ch_spec, 0, sizeof(const char *) * SPECTRUM_SAMPLES);
    std::string prefix = boost::to_lower_copy(_prefix);
    bool multichannel = false;

    /* First of all, check which layers are there */
    for (Imf::ChannelList::ConstIterator it = channels.begin(); it != channels.end(); ++it) {
        std::string name = boost::to_lower_copy(std::string(it.name()));

        /* Skip layers that have the wrong prefix */
        if (!boost::starts_with(name, prefix) && prefix != "*")
            continue;

        if (!ch_r && (name == "r" || name == "red" ||
                boost::ends_with(name, ".r") || boost::ends_with(name, ".red"))) {
            ch_r = it.name();
        } else if (!ch_g && (name == "g" || name == "green" ||
                boost::ends_with(name, ".g") || boost::ends_with(name, ".green"))) {
            ch_g = it.name();
        } else if (!ch_b && (name == "b" || name == "blue" ||
                boost::ends_with(name, ".b") || boost::ends_with(name, ".blue"))) {
            ch_b = it.name();
        } else if (!ch_a && (name == "a" || name == "alpha" ||
                boost::ends_with(name, ".a") || boost::ends_with(name, ".alpha"))) {
            ch_a = it.name();
        } else if (!ch_y && (name == "y" || name == "luminance" ||
                boost::ends_with(name, ".y") || boost::ends_with(name, ".luminance"))) {
            ch_y = it.name();
        } else if (!ch_x && (name == "x" || boost::ends_with(name, ".x"))) {
            ch_x = it.name();
        } else if (!ch_z && (name == "z" || boost::ends_with(name, ".z"))) {
            ch_z = it.name();
        } else if (!ch_ry && (name == "ry" || boost::ends_with(name, ".ry"))) {
            ch_ry = it.name();
        } else if (!ch_by && (name == "by" || boost::ends_with(name, ".by"))) {
            ch_by = it.name();
        } else {
            bool isSpectralChannel = false;
            #if SPECTRUM_SAMPLES != 3
                for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
                    std::pair<Float, Float> coverage = Spectrum::getBinCoverage(i);
                    if (!ch_spec[i] && boost::ends_with(name, formatString("%.2f-%.2fnm", coverage.first, coverage.second))) {
                        isSpectralChannel = true;
                        ch_spec[i] = it.name();
                        break;
                    }
                }
            #endif
            if (!isSpectralChannel) {
                if (_prefix == "*")
                    multichannel = true;
                else
                    Log(EWarn, "readOpenEXR(): Don't know what to do with the channel named '%s'", it.name());
            }
        }
    }

    bool spectral = true, specialColorProcessing = false,
         luminanceChromaFormat = false;
    for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
        if (!ch_spec[i])
            spectral = false;
    }

    std::string formatString;
    std::vector<const char *> sourceChannels;

    /* Now, try to categorize this image into some sort of
       generic class that we know how to deal with */
    if (multichannel) {
        for (Imf::ChannelList::ConstIterator it = channels.begin(); it != channels.end(); ++it)
            sourceChannels.push_back(it.name());
        m_channelCount = (int) sourceChannels.size();
        m_pixelFormat = EMultiChannel;
        formatString = "Multichannel";
    } else if (spectral) {
        m_pixelFormat = ESpectrum;
        formatString = "Spectrum";
        sourceChannels.insert(sourceChannels.begin(), ch_spec, ch_spec + SPECTRUM_SAMPLES);
    } else if (ch_r && ch_g && ch_b) {
        m_pixelFormat = ERGB;
        formatString = "RGB";
        sourceChannels.push_back(ch_r);
        sourceChannels.push_back(ch_g);
        sourceChannels.push_back(ch_b);
    } else if (ch_y && ch_by && ch_ry) {
        /* OpenEXR-specific luminance/chroma format */
        m_pixelFormat = ERGB;
        formatString = "YC";
        sourceChannels.push_back(ch_ry);
        sourceChannels.push_back(ch_y);
        sourceChannels.push_back(ch_by);
        luminanceChromaFormat = true;
    } else if (ch_x && ch_y && ch_z) {
        m_pixelFormat = EXYZ;
        formatString = "XYZ";
        sourceChannels.push_back(ch_x);
        sourceChannels.push_back(ch_y);
        sourceChannels.push_back(ch_z);
    } else if (ch_y) {
        m_pixelFormat = ELuminance;
        formatString = "Luminance";
        sourceChannels.push_back(ch_y);
    } else {
        Log(EError, "readOpenEXR(): Don't know how to deal with this file! There was "
            "no known pattern of color/luminance/chroma channels.");
    }

    /* Check if there is a chromaticity header entry */
    Imf::Chromaticities fileChroma;
    if (Imf::hasChromaticities(file.header()) &&
        (m_pixelFormat == ERGB || m_pixelFormat == ERGBA)) {
        fileChroma = Imf::chromaticities(file.header());

        Imf::Chromaticities ITURecBT709;
        Imf::Chromaticities XYZ(
            Imath::V2f(1.0f, 0.0f),
            Imath::V2f(0.0f, 1.0f),
            Imath::V2f(0.0f, 0.0f),
            Imath::V2f(1.0f/3.0f, 1.0f/3.0f));

        if (chromaticitiesMatch(fileChroma, ITURecBT709)) {
            /* Already in the right space -- do nothing. */
        } else if (chromaticitiesMatch(fileChroma, XYZ)) {
            /* This is an XYZ image */
            formatString = "XYZ";
            m_pixelFormat = EXYZ;
        } else {
            /* Non-standard chromaticities. Special processing is required.. */
            specialColorProcessing = true;
        }
    }

    if (ch_a) {
        m_pixelFormat = (EPixelFormat) (m_pixelFormat | 0x01);
        sourceChannels.push_back(ch_a);
        formatString += "/Alpha";
    }

    /* Load metadata if present */
    for (Imf::Header::ConstIterator it = header.begin(); it != header.end(); ++it) {
        std::string name = it.name(), typeName = it.attribute().typeName();
        const Imf::StringAttribute *sattr;
        const Imf::IntAttribute *iattr;
        const Imf::FloatAttribute *fattr;
        const Imf::DoubleAttribute *dattr;
        const Imf::V3fAttribute *vattr;
        const Imf::M44fAttribute *mattr;

        if (typeName == "string" &&
            (sattr = header.findTypedAttribute<Imf::StringAttribute>(name.c_str())))
            m_metadata.setString(name, sattr->value());
        else if (typeName == "int" &&
            (iattr = header.findTypedAttribute<Imf::IntAttribute>(name.c_str())))
            m_metadata.setInteger(name, iattr->value());
        else if (typeName == "float" &&
            (fattr = header.findTypedAttribute<Imf::FloatAttribute>(name.c_str())))
            m_metadata.setFloat(name, (Float) fattr->value());
        else if (typeName == "double" &&
            (dattr = header.findTypedAttribute<Imf::DoubleAttribute>(name.c_str())))
            m_metadata.setFloat(name, (Float) dattr->value());
        else if (typeName == "v3f" &&
            (vattr = header.findTypedAttribute<Imf::V3fAttribute>(name.c_str()))) {
            Imath::V3f vec = vattr->value();
            m_metadata.setVector(name, Vector(vec.x, vec.y, vec.z));
        } else if (typeName == "m44f" &&
            (mattr = header.findTypedAttribute<Imf::M44fAttribute>(name.c_str()))) {
            Matrix4x4 M;
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j)
                    M(i, j) = mattr->value().x[i][j];
            m_metadata.setTransform(name, Transform(M));
        }
    }

    updateChannelCount();
    m_gamma = 1.0f;
    Assert(m_channelCount == (uint8_t) sourceChannels.size());

    Imf::PixelType pxType = channels[sourceChannels[0]].type;

    size_t compSize;
    std::string encodingString;
    if (pxType == Imf::HALF) {
        m_componentFormat = EFloat16;
        compSize = sizeof(half);
        encodingString = "float16";
    } else if (pxType == Imf::FLOAT) {
        m_componentFormat = EFloat32;
        compSize = sizeof(float);
        encodingString = "float32";
    } else if (pxType == Imf::UINT) {
        m_componentFormat = EUInt32;
        compSize = sizeof(uint32_t);
        encodingString = "uint32";
    } else {
        Log(EError, "readOpenEXR(): Invalid component type (must be "
            "float16, float32, or uint32)");
        return;
    }

    /* Just how big is this image? */
    Imath::Box2i dataWindow = file.header().dataWindow();

    m_size = Vector2i(dataWindow.max.x - dataWindow.min.x + 1,
                      dataWindow.max.y - dataWindow.min.y + 1);

    /* Compute pixel / row strides */
    size_t pixelStride = compSize * m_channelCount,
           rowStride = pixelStride * m_size.x;

    /* Finally, allocate memory for it */
    m_data = static_cast<uint8_t *>(allocAligned(getBufferSize()));
    m_ownsData = true;
    char *ptr = (char *) m_data;

    ptr -= (dataWindow.min.x + dataWindow.min.y * m_size.x) * pixelStride;

    ref_vector<Bitmap> resampleBuffers((size_t) m_channelCount);
    ref<ReconstructionFilter> rfilter;

    /* Tell OpenEXR where the image data should be put */
    Imf::FrameBuffer frameBuffer;
    for (size_t i=0; i<sourceChannels.size(); ++i) {
        const char *channelName = sourceChannels[i];
        const Imf::Channel &channel = channels[channelName];
        Vector2i sampling(channel.xSampling, channel.ySampling);
        m_channelNames.push_back(channelName);

        if (channel.type != pxType)
            Log(EError, "readOpenEXR(): file has multiple channel formats, this is unsupported!");

        if (sampling == Vector2i(1)) {
            /* This is a full resolution channel. Load the ordinary way */
            frameBuffer.insert(channelName, Imf::Slice(pxType, ptr, pixelStride, rowStride));
            ptr += compSize;
        } else {
            /* Uh oh, this is a sub-sampled channel. We will need to scale it up */
            Vector2i channelSize(m_size.x / sampling.x, m_size.y / sampling.y);
            resampleBuffers[i] = new Bitmap(Bitmap::ELuminance, m_componentFormat, channelSize);
            uint8_t *resamplePtr = resampleBuffers[i]->getUInt8Data();
            resamplePtr -= (dataWindow.min.x/sampling.x + dataWindow.min.y/sampling.x * channelSize.x) * compSize;
            frameBuffer.insert(channelName, Imf::Slice(pxType, (char *) resamplePtr,
                compSize, compSize*channelSize.x, sampling.x, sampling.y));
            ptr += compSize;
        }
    }

    Log(EDebug, "Loading a %ix%i OpenEXR file (%s format, %s encoding)",
        m_size.x, m_size.y, formatString.c_str(), encodingString.c_str());

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dataWindow.min.y, dataWindow.max.y);

    for (size_t i=0; i<sourceChannels.size(); ++i) {
        if (!resampleBuffers[i])
            continue;

        if (!rfilter) {
            /* Upsample using a 2-lobed Lanczos reconstruction filter */
            Properties rfilterProps("lanczos");
            rfilterProps.setInteger("lobes", 2);
            rfilter = static_cast<ReconstructionFilter *> (
                PluginManager::getInstance()->createObject(
                MTS_CLASS(ReconstructionFilter), rfilterProps));
            rfilter->configure();
        }

        Log(EDebug, "Upsampling layer \"%s\" from %ix%i to %ix%i pixels",
            sourceChannels[i], resampleBuffers[i]->getWidth(),
            resampleBuffers[i]->getHeight(), m_size.x, m_size.y);

        resampleBuffers[i] = resampleBuffers[i]->resample(rfilter,
            ReconstructionFilter::EClamp, ReconstructionFilter::EClamp,
            m_size, -std::numeric_limits<Float>::max(), std::numeric_limits<Float>::max());

        size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;
        uint8_t *dst = m_data + compSize * i;
        uint8_t *src = resampleBuffers[i]->getUInt8Data();

        for (size_t j=0; j<pixelCount; ++j) {
            memcpy(dst, src, compSize);
            src += compSize;
            dst += pixelStride;
        }

        resampleBuffers[i] = NULL;
    }

    if (luminanceChromaFormat) {
        Imath::V3f yw = Imf::RgbaYca::computeYw(fileChroma);

        size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;
        switch (m_componentFormat) {
            case EFloat16: {
                    half *data = (half *) m_data;
                    for (size_t j=0; j<pixelCount; ++j) {
                        float ry = data[0], Y  = data[1], by = data[2],
                            r  = (ry + 1) * Y, b  = (by + 1) * Y;
                        data[0] = (half) r; data[2] = (half) b;
                        data[1] = (half) ((Y - r * yw.x - b * yw.z) / yw.y);
                        data += m_channelCount;
                    }
                }
                break;

            case EFloat32: {
                    float *data = (float *) m_data;
                    for (size_t j=0; j<pixelCount; ++j) {
                        float ry = data[0], Y  = data[1], by = data[2],
                            r  = (ry + 1) * Y, b  = (by + 1) * Y;
                        data[0] = (float) r; data[2] = (float) b;
                        data[1] = (float) ((Y - r * yw.x - b * yw.z) / yw.y);
                        data += m_channelCount;
                    }
                }
                break;

            case EUInt32: {
                    uint32_t *data = (uint32_t *) m_data;
                    double scale1 = std::numeric_limits<uint32_t>::max(),
                           scale2 = 1.0/scale1;
                    for (size_t j=0; j<pixelCount; ++j) {
                        double ry = data[0] * scale2,
                               Y  = data[1] * scale2,
                               by = data[2] * scale2,
                               r  = (ry + 1.0) * Y,
                               b  = (by + 1.0) * Y;
                        data[0] = (uint32_t) (r * scale1 + 0.5);
                        data[2] = (uint32_t) (b * scale1 + 0.5);
                        data[1] = (uint32_t) ((Y - r * yw.x - b * yw.z) / yw.y * scale1 + 0.5);
                        data += m_channelCount;
                    }
                }
                break;

            default:
                Log(EError, "Invalid component format!");
        }
    }

    if (specialColorProcessing) {
        /* Convert ITU-R Rec. BT.709 linear RGB */
        Imath::M44f M = Imf::RGBtoXYZ(fileChroma, 1) * Imf::XYZtoRGB(Imf::Chromaticities(), 1);

        size_t pixelCount = (size_t) m_size.x * (size_t) m_size.y;
        switch (m_componentFormat) {
            case EFloat16: {
                    half *data = (half *) m_data;
                    for (size_t j=0; j<pixelCount; ++j) {
                        Imath::V3f rgb = Imath::V3f(data[0], data[1], data[2]) * M;
                        data[0] = (half) rgb.x;
                        data[1] = (half) rgb.y;
                        data[2] = (half) rgb.z;
                        data += m_channelCount;
                    }
                }
                break;

            case EFloat32: {
                    float *data = (float *) m_data;
                    for (size_t j=0; j<pixelCount; ++j) {
                        Imath::V3f rgb = Imath::V3f(data[0], data[1], data[2]) * M;
                        data[0] = rgb.x; data[1] = rgb.y; data[2] = rgb.z;
                        data += m_channelCount;
                    }
                }
                break;

            case EUInt32: {
                    uint32_t *data = (uint32_t *) m_data;
                    double scale1 = std::numeric_limits<uint32_t>::max(),
                           scale2 = 1.0/scale1;
                    for (size_t j=0; j<pixelCount; ++j) {
                        Imath::V3d rgb = (Imath::V3d(data[0], data[1], data[2]) * scale2) * M;
                        data[0] = (uint32_t) (rgb.x * scale1 + 0.5);
                        data[1] = (uint32_t) (rgb.y * scale1 + 0.5);
                        data[2] = (uint32_t) (rgb.z * scale1 + 0.5);
                        data += m_channelCount;
                    }
                }
                break;

            default:
                Log(EError, "Invalid component format!");
        }
    }
}

void Bitmap::writeOpenEXR(Stream *stream) const {
    Log(EDebug, "Writing a %ix%i OpenEXR file", m_size.x, m_size.y);
    EPixelFormat pixelFormat = m_pixelFormat;

    #if SPECTRUM_SAMPLES == 3
        if (pixelFormat == ESpectrum)
            pixelFormat = ERGB;
        if (pixelFormat == ESpectrumAlpha)
            pixelFormat = ERGBA;
    #endif

    Properties metadata(m_metadata);
    if (!metadata.hasProperty("generatedBy"))
        metadata.setString("generatedBy", "Mitsuba version " MTS_VERSION);

    std::vector<std::string> keys = metadata.getPropertyNames();

    Imf::Header header(m_size.x, m_size.y);
    for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it) {
        Properties::EPropertyType type = metadata.getType(*it);

        switch (type) {
            case Properties::EString:
                header.insert(it->c_str(), Imf::StringAttribute(metadata.getString(*it)));
                break;
            case Properties::EInteger:
                header.insert(it->c_str(), Imf::IntAttribute(metadata.getInteger(*it)));
                break;
            case Properties::EFloat:
                header.insert(it->c_str(), Imf::FloatAttribute((float) metadata.getFloat(*it)));
                break;
            case Properties::EPoint: {
                    Point val = metadata.getPoint(*it);
                    header.insert(it->c_str(), Imf::V3fAttribute(
                        Imath::V3f((float) val.x, (float) val.y, (float) val.z)));
                }
                break;
            case Properties::ETransform: {
                    Matrix4x4 val = metadata.getTransform(*it).getMatrix();
                    header.insert(it->c_str(), Imf::M44fAttribute(Imath::M44f(
                        (float) val(0, 0), (float) val(0, 1), (float) val(0, 2), (float) val(0, 3),
                        (float) val(1, 0), (float) val(1, 1), (float) val(1, 2), (float) val(1, 3),
                        (float) val(2, 0), (float) val(2, 1), (float) val(2, 2), (float) val(2, 3),
                        (float) val(3, 0), (float) val(3, 1), (float) val(3, 2), (float) val(3, 3))));
                }
                break;
            default:
                header.insert(it->c_str(), Imf::StringAttribute(metadata.getAsString(*it)));
                break;
        }
    }

    if (pixelFormat == EXYZ || pixelFormat == EXYZA) {
        Imf::addChromaticities(header, Imf::Chromaticities(
            Imath::V2f(1.0f, 0.0f),
            Imath::V2f(0.0f, 1.0f),
            Imath::V2f(0.0f, 0.0f),
            Imath::V2f(1.0f/3.0f, 1.0f/3.0f)));
    } else if (pixelFormat == ERGB || pixelFormat == ERGBA) {
        Imf::addChromaticities(header, Imf::Chromaticities());
    }

    Imf::PixelType compType;
    size_t compStride;
    if (m_componentFormat == EFloat16) {
        compType = Imf::HALF;
        compStride = 2;
    } else if (m_componentFormat == EFloat32) {
        compType = Imf::FLOAT;
        compStride = 4;
    } else if (m_componentFormat == EUInt32) {
        compType = Imf::UINT;
        compStride = 4;
    } else {
        Log(EError, "writeOpenEXR(): Invalid component type (must be "
            "float16, float32, or uint32)");
        return;
    }

    if (!m_channelNames.empty() && m_channelNames.size() != (size_t) getChannelCount())
        Log(EWarn, "writeOpenEXR(): 'channelNames' has the wrong number of entries (%i, expected %i), ignoring..!",
            (int) m_channelNames.size(), (int) m_channelCount);

    bool explicitChannelNames = false;
    Imf::ChannelList &channels = header.channels();
    if (m_channelNames.size() == (size_t) getChannelCount()) {
        for (size_t i=0; i<m_channelNames.size(); ++i)
            channels.insert(m_channelNames[i].c_str(), Imf::Channel(compType));
        explicitChannelNames = true;
    } else if (pixelFormat == ELuminance || pixelFormat == ELuminanceAlpha) {
        channels.insert("Y", Imf::Channel(compType));
    } else if (pixelFormat == ERGB || pixelFormat == ERGBA ||
            pixelFormat == EXYZ || pixelFormat == EXYZA) {
        channels.insert("R", Imf::Channel(compType));
        channels.insert("G", Imf::Channel(compType));
        channels.insert("B", Imf::Channel(compType));
    } else if (pixelFormat == ESpectrum || pixelFormat == ESpectrumAlpha) {
        for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
            std::pair<Float, Float> coverage = Spectrum::getBinCoverage(i);
            std::string name = formatString("%.2f-%.2fnm", coverage.first, coverage.second);
            channels.insert(name.c_str(), Imf::Channel(compType));
        }
    } else if (pixelFormat == EMultiChannel) {
        for (int i=0; i<getChannelCount(); ++i)
            channels.insert(formatString("%i", i).c_str(), Imf::Channel(compType));
    } else {
        Log(EError, "writeOpenEXR(): Invalid pixel format!");
        return;
    }

    if ((pixelFormat == ELuminanceAlpha || pixelFormat == ERGBA ||
        pixelFormat == EXYZA || pixelFormat == ESpectrumAlpha) && !explicitChannelNames)
        channels.insert("A", Imf::Channel(compType));

    size_t pixelStride = m_channelCount * compStride,
           rowStride = pixelStride * m_size.x;
    char *ptr = (char *) m_data;

    Imf::FrameBuffer frameBuffer;

    if (explicitChannelNames) {
        for (size_t i=0; i<m_channelNames.size(); ++i) {
            frameBuffer.insert(m_channelNames[i].c_str(), Imf::Slice(compType, ptr, pixelStride, rowStride));
            ptr += compStride;
        }
    } else if (pixelFormat == ELuminance || pixelFormat == ELuminanceAlpha) {
        frameBuffer.insert("Y", Imf::Slice(compType, ptr, pixelStride, rowStride)); ptr += compStride;
    } else if (pixelFormat == ERGB || pixelFormat == ERGBA || pixelFormat == EXYZ || pixelFormat == EXYZA) {
        frameBuffer.insert("R", Imf::Slice(compType, ptr, pixelStride, rowStride)); ptr += compStride;
        frameBuffer.insert("G", Imf::Slice(compType, ptr, pixelStride, rowStride)); ptr += compStride;
        frameBuffer.insert("B", Imf::Slice(compType, ptr, pixelStride, rowStride)); ptr += compStride;
    } else if (pixelFormat == ESpectrum || pixelFormat == ESpectrumAlpha) {
        for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
            std::pair<Float, Float> coverage = Spectrum::getBinCoverage(i);
            std::string name = formatString("%.2f-%.2fnm", coverage.first, coverage.second);
            frameBuffer.insert(name.c_str(), Imf::Slice(compType, ptr, pixelStride, rowStride)); ptr += compStride;
        }
    } else if (pixelFormat == EMultiChannel) {
        for (int i=0; i<getChannelCount(); ++i) {
            frameBuffer.insert(formatString("%i", i).c_str(), Imf::Slice(compType, ptr, pixelStride, rowStride));
            ptr += compStride;
        }
    }

    if ((pixelFormat == ELuminanceAlpha || pixelFormat == ERGBA ||
         pixelFormat == EXYZA || pixelFormat == ESpectrumAlpha) && !explicitChannelNames)
        frameBuffer.insert("A", Imf::Slice(compType, ptr, pixelStride, rowStride));

    EXROStream ostr(stream);
    Imf::OutputFile file(ostr, header);
    file.setFrameBuffer(frameBuffer);
    file.writePixels(m_size.y);
}
#else
void Bitmap::readOpenEXR(Stream *stream, const std::string &_prefix) {
    Log(EError, "Bitmap::readOpenEXR(): OpenEXR support was disabled at compile time!");
}
void Bitmap::writeOpenEXR(Stream *stream) const {
    Log(EError, "Bitmap::writeOpenEXR(): OpenEXR support was disabled at compile time!");
}
#endif

void Bitmap::readTGA(Stream *stream) {
    Stream::EByteOrder byteOrder = stream->getByteOrder();
    stream->setByteOrder(Stream::ELittleEndian);

    int headerSize = stream->readUChar();
    if (stream->readUChar() != 0)
        Log(EError, "readTGA(): indexed files are not supported!");
    int colorType = stream->readUChar();
    if (colorType != 2 && colorType != 3 && colorType != 10 && colorType != 11)
        Log(EError, "readTGA(): only grayscale & RGB[A] files are supported!");

    stream->skip(9);

    int width = stream->readShort();
    int height = stream->readShort();
    uint8_t bpp = stream->readUChar();
    uint8_t descriptor = stream->readUChar();
    stream->skip(headerSize);

    m_size = Vector2i(width, height);
    m_gamma = -1;
    m_componentFormat = EUInt8;

    Log(ETrace, "Loading a %ix%i TGA file", m_size.x, m_size.y);

    bool vflip = !(descriptor & (1 << 5));
    bool greyscale = colorType == 3 || colorType == 11;
    bool rle = colorType & 8;

    if ((bpp == 8 && !greyscale) || (bpp != 8 && greyscale))
        Log(EError, "readTGA(): Invalid bit depth!");

    switch (bpp) {
        case 8: m_pixelFormat = ELuminance; break;
        case 24: m_pixelFormat = ERGB; break;
        case 32: m_pixelFormat = ERGBA; break;
        default:
            Log(EError, "readTGA(): Invalid bit depth!");
    }
    updateChannelCount();

    size_t bufferSize = getBufferSize(),
           rowSize = bufferSize / height;

    m_data = static_cast<uint8_t *>(allocAligned(bufferSize));
    m_ownsData = true;
    int channels = bpp/8;

    if (!rle) {
        for (int y=0; y<height; ++y) {
            size_t targetY = vflip ? (height - y - 1) : y;
            stream->read(m_data + targetY * rowSize, rowSize);
        }
    } else {
        /* Decode an RLE-encoded image */
        uint8_t temp[4],
                *ptr = m_data,
                *end = m_data + bufferSize;

        while (ptr != end) {
            uint8_t value = stream->readUChar();

            if (value & 0x80) {
                /* Run length packet */
                uint8_t count = (value & 0x7F) + 1;
                stream->read(temp, channels);
                for (uint32_t i=0; i<count; ++i)
                    for (int j=0; j<channels; ++j)
                    *ptr++ = temp[j];
            } else {
                /* Raw packet */
                uint32_t count = channels * (value + 1);
                for (uint32_t i=0; i<count; ++i)
                    *ptr++ = stream->readUChar();
            }
        }
        if (vflip)
            flipVertically();
    }

    if (!greyscale) {
        /* Convert BGR to RGB */
        for (size_t i=0; i<bufferSize; i += channels)
            std::swap(m_data[i], m_data[i+2]);
    }

    stream->setByteOrder(byteOrder);
}

void Bitmap::readBMP(Stream *stream) {
    Stream::EByteOrder byteOrder = stream->getByteOrder();
    stream->setByteOrder(Stream::ELittleEndian);

    uint8_t magic1 = stream->readUChar();
    uint8_t magic2 = stream->readUChar();

    if (magic1 != 'B' || magic2 != 'M')
        Log(EError, "readBMP(): Invalid header identifier!");

    stream->skip(8);

    uint32_t bmpOffset = stream->readUInt();
    uint32_t headerSize = stream->readUInt();
    int32_t width = stream->readInt();
    int32_t height = stream->readInt();
    uint16_t nplanes = stream->readUShort();
    uint16_t bpp = stream->readUShort();
    uint32_t compressionType = stream->readUInt();
    stream->skip(bmpOffset-34);

    if (headerSize != 40 || nplanes != 1 || width <= 0)
        Log(EError, "readBMP(): Unsupported BMP format encountered!");

    if (compressionType != 0)
        Log(EError, "readBMP(): Compressed files are currently not supported!");

    m_size = Vector2i(width, std::abs(height));
    m_componentFormat = EUInt8;
    m_gamma = -1.0f;

    switch (bpp) {
        case 1:
            m_pixelFormat = ELuminance;
            m_componentFormat = EBitmask;
            break;
        case 8: m_pixelFormat = ELuminance; break;
        case 16: m_pixelFormat = ELuminanceAlpha; break;
        case 24: m_pixelFormat = ERGB; break;
        case 32: m_pixelFormat = ERGBA; break;
        default:
            Log(EError, "readBMP(): Invalid bit depth (%i)!", bpp);
    }
    updateChannelCount();

    size_t bufferSize = getBufferSize();
    m_data = static_cast<uint8_t *>(allocAligned(bufferSize));
    m_ownsData = true;

    Log(ETrace, "Loading a %ix%i BMP file", m_size.x, m_size.y);

    int rowSize = (int) bufferSize / m_size.y;

    int padding = -rowSize & 3;
    bool vflip = height > 0;

    for (int y=0; y<m_size.y; ++y) {
        int targetY = vflip ? (m_size.y - y - 1) : y;
        stream->read(m_data + rowSize * targetY, rowSize);
        stream->skip(padding);
    }

    if (m_pixelFormat == ERGB || m_pixelFormat == ERGBA) {
        int channels = getChannelCount();
        for (size_t i=0; i<bufferSize; i += channels)
            std::swap(m_data[i], m_data[i+2]);
    }

    stream->setByteOrder(byteOrder);
}

/* The following is based on code by Bruce Walter */
namespace detail {
    static inline void RGBE_FromFloat(float *data, uint8_t rgbe[4]) {
        /* Find the largest contribution */
        Float max = std::max(std::max(data[0], data[1]), data[2]);
        if (max < 1e-32) {
            rgbe[0] = rgbe[1] = rgbe[2] = rgbe[3] = 0;
        } else {
            int e;
            /* Extract exponent and convert the fractional part into
               the [0..255] range. Afterwards, divide by max so that
               any color component multiplied by the result will be in [0,255] */
            max = std::frexp(max, &e) * (Float) 256 / max;
            rgbe[0] = (uint8_t) (data[0] * max);
            rgbe[1] = (uint8_t) (data[1] * max);
            rgbe[2] = (uint8_t) (data[2] * max);
            rgbe[3] = e+128; /* Exponent value in bias format */
        }
    }

    static inline void RGBE_ToFloat(uint8_t rgbe[4], float *data) {
        if (rgbe[3]) { /* nonzero pixel */
            float f = std::ldexp(1.0f, (int) rgbe[3] - (128+8));
            for (int i=0; i<3; ++i)
                data[i] = rgbe[i] * f;
        } else {
            memset(data, 0, sizeof(float)*3);
        }
    }

    /* The code below is only needed for the run-length encoded files.
       Run length encoding adds considerable complexity but does
       save some space.  For each scanline, each channel (r,g,b,e) is
       encoded separately for better compression. */
    static inline void RGBE_WriteBytes_RLE(Stream *stream, uint8_t *data, int numbytes) {
        int cur = 0;
        uint8_t buf[2];

        while (cur < numbytes) {
            int beg_run = cur;
            /* find next run of length at least 4 if one exists */
            int run_count = 0, old_run_count = 0;
            while (run_count < 4 && beg_run < numbytes) {
                beg_run += run_count;
                old_run_count = run_count;
                run_count = 1;
                while ((beg_run + run_count < numbytes) && (run_count < 127)
                    && (data[beg_run] == data[beg_run + run_count]))
                    run_count++;
            }
            /* if data before next big run is a short run then write it as such */
            if (old_run_count > 1 && old_run_count == beg_run - cur) {
                buf[0] = 128 + old_run_count;   /*write short run*/
                buf[1] = data[cur];
                stream->write(buf, 2);
                cur = beg_run;
            }
            /* write out bytes until we reach the start of the next run */
            while (cur < beg_run) {
                int nonrun_count = beg_run - cur;
                if (nonrun_count > 128)
                    nonrun_count = 128;
                buf[0] = nonrun_count;
                stream->write(buf, 1);
                stream->write(&data[cur], nonrun_count);
                cur += nonrun_count;
            }
            /* write out next run if one was found */
            if (run_count >= 4) {
                buf[0] = 128 + run_count;
                buf[1] = data[beg_run];
                stream->write(buf, 2);
                cur += run_count;
            }
        }
    }

    /* simple read routine.  will not correctly handle run length encoding */
    inline static void RGBE_ReadPixels(Stream *stream, float *data, size_t numpixels) {
        while (numpixels-- > 0) {
            uint8_t rgbe[4];
            stream->read(rgbe, 4);
            RGBE_ToFloat(rgbe, data);
            data += 3;
        }
    }
}

void Bitmap::readRGBE(Stream *stream) {
    std::string line = stream->readLine();

    if (line.length() < 2 || line[0] != '#' || line[1] != '?')
        Log(EError, "readRGBE(): Invalid header!");

    bool format_recognized = false;
    while (true) {
        line = stream->readLine();
        if (boost::starts_with(line, "FORMAT=32-bit_rle_rgbe"))
            format_recognized = true;
        if (boost::starts_with(line, "-Y ")) {
            if (sscanf(line.c_str(), "-Y %i +X %i", &m_size.y, &m_size.x) < 2)
                Log(EError, "readRGBE(): parser error!");
            break;
        }
    }

    if (!format_recognized)
        Log(EError, "readRGBE(): invalid format!");

    m_pixelFormat = ERGB;
    m_componentFormat = EFloat32;
    m_channelCount = 3;
    m_gamma = 1.0f;
    m_data = static_cast<uint8_t *>(allocAligned(getBufferSize()));
    m_ownsData = true;
    float *data = (float *) m_data;

    if (m_size.x < 8 || m_size.x > 0x7fff) {
        /* run length encoding is not allowed so read flat*/
        detail::RGBE_ReadPixels(stream, data, (size_t) m_size.x * (size_t) m_size.y);
        return;
    }

    uint8_t *buffer = new uint8_t[4*m_size.x];

    try {
        /* Read in each successive scanline */
        for (int y=0; y<m_size.y; ++y) {
            uint8_t rgbe[4];
            stream->read(rgbe, 4);

            if (rgbe[0] != 2 || rgbe[1] != 2 || rgbe[2] & 0x80) {
                /* this file is not run length encoded */
                detail::RGBE_ToFloat(rgbe, data);
                detail::RGBE_ReadPixels(stream, data + 3, (size_t) m_size.x * (size_t) m_size.y - 1);
                return;
            }

            if ((((int) rgbe[2]) << 8 | rgbe[3]) != m_size.x)
                Log(EError, "readRGBE(): wrong scanline width!");

            uint8_t *ptr = buffer;

            /* read each of the four channels for the scanline into the buffer */
            for (int i=0;i<4;i++) {
                uint8_t *ptr_end = buffer + (i+1) * m_size.x;

                while (ptr < ptr_end) {
                    uint8_t buf[2];
                    stream->read(buf, 2);

                    if (buf[0] > 128) {
                        /* a run of the same value */
                        int count = buf[0] - 128;
                        if (count == 0 || count > ptr_end - ptr)
                            Log(EError, "readRGBE(): bad scanline data!");

                        while (count-- > 0)
                            *ptr++ = buf[1];
                    } else {
                        /* a non-run */
                        int count = buf[0];
                        if (count == 0 || count > ptr_end - ptr)
                            Log(EError, "readRGBE(): bad scanline data!");
                        *ptr++ = buf[1];
                        if (--count > 0)
                            stream->read(ptr, count);
                        ptr += count;
                    }
                }
            }
            /* now convert data from buffer into floats */
            for (int i=0; i<m_size.x; i++) {
                rgbe[0] = buffer[i];
                rgbe[1] = buffer[m_size.x+i];
                rgbe[2] = buffer[2*m_size.x+i];
                rgbe[3] = buffer[3*m_size.x+i];
                detail::RGBE_ToFloat(rgbe, data);
                data += 3;
            }
        }
    } catch (...) {
        delete[] buffer;
        throw;
    }

    delete[] buffer;
}

void Bitmap::writeRGBE(Stream *stream) const {
    if (m_componentFormat != EFloat32)
        Log(EError, "writeRGBE(): component format must be EFloat32!");
    if (m_pixelFormat != ERGB && m_pixelFormat != ERGBA)
        Log(EError, "writeRGBE(): pixel format must be ERGB or ERGBA!");

    stream->writeLine("#?RGBE");

    std::vector<std::string> keys = m_metadata.getPropertyNames();
    for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it) {
        stream->writeLine(formatString("# Metadata [%s]:", it->c_str()));
        std::istringstream iss(m_metadata.getAsString(*it));
        std::string buf;
        while (std::getline(iss, buf))
            stream->writeLine(formatString("#   %s", buf.c_str()));
    }

    stream->writeLine("FORMAT=32-bit_rle_rgbe\n");
    stream->writeLine(formatString("-Y %i +X %i", m_size.y, m_size.x));

    float *data = (float *) m_data;
    if (m_size.x < 8 || m_size.x > 0x7fff) {
        /* Run length encoding is not allowed so write flat*/
        uint8_t rgbe[4];
        for (size_t i=0; i<(size_t) m_size.x * (size_t) m_size.y; ++i) {
            detail::RGBE_FromFloat(data, rgbe);
            data += (m_pixelFormat == ERGB) ? 3 : 4;
            stream->write(rgbe, 4);
        }
        return;
    }

    uint8_t *buffer = new uint8_t[4*m_size.x];
    for (int y=0; y<m_size.y; ++y) {
        uint8_t rgbe[4] = { 2, 2,
            (uint8_t) (m_size.x >> 8),
            (uint8_t) (m_size.x & 0xFF) };
        stream->write(rgbe, 4);

        for (int x=0; x<m_size.x; x++) {
            detail::RGBE_FromFloat(data, rgbe);

            buffer[x]            = rgbe[0];
            buffer[m_size.x+x]   = rgbe[1];
            buffer[2*m_size.x+x] = rgbe[2];
            buffer[3*m_size.x+x] = rgbe[3];

            data += (m_pixelFormat == ERGB) ? 3 : 4;
        }

        /* Write out each of the four channels separately run length encoded.
           First red, then green, then blue, then exponent */
        for (int i=0;i<4;i++)
            detail::RGBE_WriteBytes_RLE(stream, &buffer[i*m_size.x], m_size.x);
    }

    delete[] buffer;
}

/// Simple helper function for reading strings in PFM files
static std::string pfmReadString(Stream *stream) {
    std::string result;

    while (true) {
        char data = stream->readChar();
        if (::isspace(data))
            break;
        result += data;
    }

    return result;
}

void Bitmap::readPFM(Stream *stream) {
    char header[3];
    stream->read(header, 3);
    if (header[0] != 'P' || !(header[1] == 'F' || header[1] == 'f'))
        Log(EError, "readPFM(): Invalid header!");

    bool color = (header[1] == 'F');
    m_pixelFormat = color ? ERGB : ELuminance;
    m_componentFormat = EFloat32;
    m_channelCount = color ? 3 : 1;
    m_gamma = 1.0f;

    char *end_ptr = NULL;
    std::string widthString = pfmReadString(stream);
    m_size.x = (int) strtol(widthString.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        SLog(EError, "Could not parse image dimensions!");

    std::string heightString = pfmReadString(stream);
    m_size.y = (int) strtol(heightString.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        SLog(EError, "Could not parse image dimensions!");

    std::string scaleAndOrderString = pfmReadString(stream);
    float scaleAndOrder = (float) strtod(scaleAndOrderString.c_str(), &end_ptr);
    if (*end_ptr != '\0')
        SLog(EError, "Could not parse scale/order information!");

    m_data = static_cast<uint8_t *>(allocAligned(getBufferSize()));
    m_ownsData = true;
    float *data = (float *) m_data;

    Stream::EByteOrder backup = stream->getByteOrder();
    size_t size = getPixelCount() * m_channelCount;
    stream->setByteOrder(scaleAndOrder <= 0.0f ? Stream::ELittleEndian : Stream::EBigEndian);

    try {
        stream->readSingleArray(data, size);
    } catch (...) {
        stream->setByteOrder(backup);
        throw;
    }

    stream->setByteOrder(backup);
    const float scale = std::abs(scaleAndOrder);
    if (scale != 1) {
        for (size_t i=0; i<size; ++i)
            data[i] *= scale;
    }
    flipVertically();
}

void Bitmap::writePFM(Stream *stream) const {
    if (m_componentFormat != EFloat32)
        Log(EError, "writePFM(): component format must be EFloat32!");
    if (m_pixelFormat != ERGB && m_pixelFormat != ERGBA && m_pixelFormat != ELuminance)
        Log(EError, "writePFM(): pixel format must be ERGB, ERGBA, ELuminance, or ELuminanceAlpha!");

    /* Write the header */
    std::ostringstream oss;
    oss << 'P' << ((m_pixelFormat == ERGB || m_pixelFormat == ERGBA) ? 'F' : 'f') << '\n';
    oss << m_size.x << ' ' << m_size.y << '\n';
    oss << (Stream::getHostByteOrder() == Stream::ELittleEndian ? "-1" : "1") << '\n';
    std::string header = oss.str();
    stream->write(header.c_str(), header.length());

    float *data = (float *) m_data;
    if (m_pixelFormat == ERGB || m_pixelFormat == ELuminance) {
        size_t scanline = (size_t) m_size.x * m_channelCount;

        for (int y=0; y<m_size.y; ++y)
            stream->write(data + scanline*(m_size.y - 1 - y), scanline * sizeof(float));
    } else {
        /* For convenience: also handle images with an alpha
           channel, but strip it out while saving the data */
        size_t scanline = (size_t) m_size.x * m_channelCount;
        float *temp = (float *) alloca(scanline * sizeof(float));

        for (int y=0; y<m_size.y; ++y) {
            const float *source = data + scanline*(m_size.y - 1 - y);
            float *dest = temp;

            for (int x=0; x<m_size.x; ++x) {
                for (uint8_t j=0; j<m_channelCount-1; ++j)
                    *dest++ = *source++;
                source++;
            }

            stream->write(temp, sizeof(float) * m_size.x * (m_channelCount-1));
        }
    }
}

void Bitmap::readPPM(Stream *stream) {
    int field = 0, nChars = 0;

    std::string fields[4];

    while (field < 4) {
        char c = stream->readChar();
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            if (nChars != 0) {
                nChars = 0;
                ++field;
            }
        } else {
            fields[field] += c;
            ++nChars;
        }
    }
    if (fields[0] != "P6")
        Log(EError, "readPPM(): invalid format!");

    int intValues[3];
    for (int i=0; i<3; ++i) {
        char *end_ptr = NULL;
        intValues[i] = strtol(fields[i+1].c_str(), &end_ptr, 10);
        if (*end_ptr != '\0')
            SLog(EError, "readPPM(): unable to parse the file header!");
    }

    m_size.x = intValues[0];
    m_size.y = intValues[1];
    m_pixelFormat = ERGB;
    m_channelCount = 3;
    m_gamma = -1.0f;
    m_ownsData = true;
    m_componentFormat = intValues[2] <= 0xFF ? EUInt8 : EUInt16;
    size_t size = getBufferSize();
    m_data = static_cast<uint8_t *>(allocAligned(size));
    stream->read(m_data, size);
}

void Bitmap::writePPM(Stream *stream) const {
    if (m_pixelFormat != ERGB || (m_componentFormat != EUInt8 && m_componentFormat != EUInt16))
        Log(EError, "writePPM(): Only 8 or 16-bit RGB images are supported");
    stream->writeLine(formatString("P6\n%i\n%i\n%i\n", m_size.x, m_size.y,
        m_componentFormat == EUInt8 ? 0xFF : 0xFFFF).c_str());
    stream->write(m_data, getBufferSize());
}

void Bitmap::staticInitialization() {
#if defined(MTS_HAS_OPENEXR)
    /* Prevent races during the OpenEXR initialization */
    Imf::staticInitialize();

    /* Use multiple threads to read/write OpenEXR files */
    Imf::setGlobalThreadCount(getCoreCount());
#endif

    /* Initialize the Bitmap format conversion */
    FormatConverter::staticInitialization();

#if defined(MTS_HAS_FFTW)
    /* Initialize FFTW if enabled */
    fftw_init_threads();
    fftw_plan_with_nthreads(getCoreCount());
#endif
}

void Bitmap::staticShutdown() {
    FormatConverter::staticShutdown();

#if defined(MTS_HAS_FFTW)
    fftw_cleanup_threads();
#endif
}

ref<Bitmap> Bitmap::expand() {
    if (m_componentFormat != EBitmask)
        return this;

    ref<Bitmap> output= new Bitmap(m_pixelFormat, EUInt8, m_size);
    uint8_t *outputBuffer = output->getUInt8Data();

    size_t bytesPerRow = (m_size.x * m_channelCount + 7) / 8; // round up to full bytes
    for (int y=0; y<m_size.y; ++y) {
        uint8_t *inputBuffer = m_data + (bytesPerRow * y);

        for (int x=0; x<m_size.x; ++x) {
            int entry = x / 8, bit = x % 8;
            *outputBuffer++ = (inputBuffer[entry] & (1 << bit)) ? 255 : 0;
        }
    }

    return output;
}

void Bitmap::drawWorkUnit(const Point2i &offset, const Vector2i &size, int worker) {
    int ox = offset.x, oy = offset.y,
        ex = ox + size.x, ey = oy + size.y;
    if (size.x < 3 || size.y < 3)
        return;

    const float *color = NULL;

    /* Use desaturated colors to highlight the worker
       responsible for rendering the current image */
    const float white[]     = { 1.0f, 1.0f, 1.0f };
    const float red[]       = { 1.0f, 0.3f, 0.3f };
    const float green[]     = { 0.3f, 1.0f, 0.3f };
    const float blue[]      = { 0.3f, 0.3f, 1.0f };
    const float gray[]      = { 0.5f, 0.5f, 0.5f };
    const float yellow[]    = { 1.0f, 1.0f, 0.0f };
    const float magenta[]   = { 1.0f, 0.3f, 1.0f };
    const float turquoise[] = { 0.3f, 1.0f, 1.0f };

    switch (worker % 8) {
        case 1: color = green; break;
        case 2: color = yellow; break;
        case 3: color = blue; break;
        case 4: color = gray; break;
        case 5: color = red; break;
        case 6: color = magenta; break;
        case 7: color = turquoise; break;
        case 0:
        default:
            color = white;
            break;
    }

    float scale = .7f * (color[0] * 0.212671f + color[1] * 0.715160f + color[2] * 0.072169f);

    Spectrum spec;
    spec.fromLinearRGB(color[0]*scale,
        color[1]*scale,
        color[2]*scale);

    drawHLine(oy, ox, ox + 3, spec);
    drawHLine(oy, ex - 4, ex - 1, spec);
    drawHLine(ey - 1, ox, ox + 3, spec);
    drawHLine(ey - 1, ex - 4, ex - 1, spec);
    drawVLine(ox, oy, oy + 3, spec);
    drawVLine(ex - 1, oy, oy + 3, spec);
    drawVLine(ex - 1, ey - 4, ey - 1, spec);
    drawVLine(ox, ey - 4, ey - 1, spec);
}

std::ostream &operator<<(std::ostream &os, const Bitmap::EPixelFormat &value) {
    switch (value) {
        case Bitmap::ELuminance: os << "luminance"; break;
        case Bitmap::ELuminanceAlpha: os << "luminanceAlpha"; break;
        case Bitmap::ERGB: os << "rgb"; break;
        case Bitmap::ERGBA: os << "rgba"; break;
        case Bitmap::EXYZ: os << "xyz"; break;
        case Bitmap::EXYZA: os << "xyza"; break;
        case Bitmap::ESpectrum: os << "spectrum"; break;
        case Bitmap::ESpectrumAlpha: os << "spectrumAlpha"; break;
        case Bitmap::ESpectrumAlphaWeight: os << "spectrumAlphaWeight"; break;
        case Bitmap::EMultiSpectrumAlphaWeight: os << "multiSpectrumAlphaWeight"; break;
        case Bitmap::EMultiChannel: os << "multiChannel"; break;
        default: os << "invalid"; break;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Bitmap::EComponentFormat &value) {
    switch (value) {
        case Bitmap::EBitmask: os << "bitmask"; break;
        case Bitmap::EUInt8: os << "uint8"; break;
        case Bitmap::EUInt16: os << "uint16"; break;
        case Bitmap::EUInt32: os << "uint32"; break;
        case Bitmap::EFloat16: os << "float16"; break;
        case Bitmap::EFloat32: os << "float32"; break;
        case Bitmap::EFloat64: os << "float64"; break;
        default: os << "invalid"; break;
    }
    return os;
}

MTS_IMPLEMENT_CLASS(Bitmap, false, Object)
MTS_NAMESPACE_END
