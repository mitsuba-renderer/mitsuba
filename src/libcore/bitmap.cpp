/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

MTS_NAMESPACE_BEGIN

#if defined(MTS_HAS_OPENEXR)
/* ========================== *
 *     EXR helper classes     *
 * ========================== */

class EXRIStream : public Imf::IStream {
public:
	EXRIStream(Stream *stream) : IStream(stream->toString().c_str()),
		m_stream(stream) {
		m_offset = stream->getPos();
	}

	bool read(char *c, int n) {
		m_stream->read(c, n);
		return m_stream->isEOF();
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
	size_t m_offset;
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
		const Vector2i &size, int channelCount) : m_pixelFormat(pFormat),
		m_componentFormat(cFormat), m_size(size), m_channelCount(channelCount) {
	AssertEx(size.x > 0 && size.y > 0, "Invalid bitmap size");

	if (m_componentFormat == EUInt8)
		m_gamma = -1.0f; // sRGB by default
	else
		m_gamma = 1.0f; // Linear by default

	updateChannelCount();

	m_data = static_cast<uint8_t *>(allocAligned(getBufferSize()));
}

Bitmap::Bitmap(EFileFormat format, Stream *stream, const std::string &prefix) : m_data(NULL) {
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
		case ETGA: readTGA(stream); break;
		case EPNG: readPNG(stream); break;
		default:
			Log(EError, "Bitmap: Invalid file format!");
	}
}

void Bitmap::write(EFileFormat format, Stream *stream, int compression,
		const std::vector<std::string> *channelNames) const {
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
		case EOpenEXR: writeOpenEXR(stream, channelNames); break;
		case ERGBE: writeRGBE(stream); break;
		case EPFM: writePFM(stream); break;
		default:
			Log(EError, "Bitmap::write(): Invalid file format!");
	}
}

size_t Bitmap::getBufferSize() const {
	size_t bitsPerRow = m_size.x * m_channelCount * getBitsPerComponent();
	size_t bytesPerRow = (bitsPerRow + 7) / 8; // round up to full bytes
	return bytesPerRow * (size_t) m_size.y;
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
	if (m_data)
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
		columns      = size.x * m_channelCount,
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
		(size_t) m_size.x * (size_t) m_size.y, multiplier, intent);
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

	ref<Bitmap> target = new Bitmap(pixelFormat, componentFormat, m_size);
	target->setMetadata(m_metadata);
	target->setGamma(gamma);

	cvt->convert(m_pixelFormat, m_gamma, m_data,
		pixelFormat, gamma, target->getData(),
		(size_t) m_size.x * (size_t) m_size.y, multiplier, intent);

	return target;
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
		(size_t) m_size.x * (size_t) m_size.y, multiplier, intent);
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
			data[0] = (T)(  3.240479f * X + -1.537150f * Y + -0.498535f * Z);
			data[1] = (T)( -0.969256f * X +  1.875991f * Y +  0.041556f * Z);
			data[2] = (T)(  0.055648f * X + -0.204043f * Y +  1.057311f * Z);

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

			data[0] = (T) X;
			data[1] = (T) Y;
			data[2] = (T) Z;

			data += channels;
		}

	} else {
		/* Monochrome version */
		for (size_t i=0; i < pixels; ++i) {
			Float Lp = (Float) *data * scale;

			/* Apply the tonemapping transformation */
			*data = (T) (Lp * (1.0f + Lp*invWp2) / (1.0f + Lp));

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

ref<Bitmap> Bitmap::separateChannel(int channelIndex) {
	int channelCount = getChannelCount();

	if (channelIndex == 0 && channelCount == 1)
		return this;

	Assert(channelIndex > 0 && channelIndex < channelCount);

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
			*target++ = *source;
		source += stride;
	}

	return result;
}

ref<Bitmap> Bitmap::join(EPixelFormat fmt,
			const std::vector<Bitmap *> &sourceBitmaps) {
	const Bitmap *ch0 = sourceBitmaps.at(0);
	if (ch0->getComponentFormat() == EBitmask)
		Log(EError, "Conversions involving bitmasks are currently not supported!");

	ref<Bitmap> result = new Bitmap(fmt, ch0->getComponentFormat(),
		ch0->getSize());
	size_t channelCount = (size_t) result->getChannelCount();

	result->setMetadata(ch0->getMetadata());
	result->setGamma(ch0->getGamma());

	if (channelCount == sourceBitmaps.size())
		Log(EError, "Bitmap::join(): Error -- supplied the wrong number "
			"of channels (%i instead of %i)", (int) sourceBitmaps.size(),
			result->getChannelCount());

	for (size_t i=0; i<channelCount; ++i) {
		if (sourceBitmaps[i]->getSize() != ch0->getSize())
			Log(EError, "Bitmap::join(): Detected a size mismatch!");

		if (sourceBitmaps[i]->getComponentFormat() != ch0->getComponentFormat())
			Log(EError, "Bitmap::join(): Detected a component format mismatch!");

		if (sourceBitmaps[i]->getPixelFormat() != ELuminance)
			Log(EError, "Bitmap::join(): Detected a pixel format mismatch (expected ELuminance)!");
	}

	size_t pixelCount = (size_t) ch0->getSize().x
		* (size_t) ch0->getSize().y;
	size_t componentSize = ch0->getBytesPerComponent();

	uint8_t **pointers = (uint8_t **) alloca(channelCount * sizeof(uint8_t *));
	for (size_t i = 0; i<channelCount; ++i)
		pointers[i] = sourceBitmaps[i]->getUInt8Data();

	uint8_t *dest = result->getUInt8Data();

	for (size_t i = 0; i<pixelCount; ++i)
		for (size_t j = 0; j<channelCount; ++j)
			for (size_t k= 0; k < componentSize; ++k)
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


/// Bitmap resampling utility function
template <typename Scalar> static void resample(const ReconstructionFilter *rfilter,
	ReconstructionFilter::EBoundaryCondition bch,
	ReconstructionFilter::EBoundaryCondition bcv,
	const Bitmap *source, Bitmap *target, Float minValue, Float maxValue) {
	ref<Bitmap> temp; // Pointer to a temporary bitmap

	int channels = source->getChannelCount();

	if (source->getWidth() != target->getWidth()) {
		/* Re-sample along the X direction */
		Resampler<Scalar> r(rfilter, bch, source->getWidth(), target->getWidth());

		/* Create a bitmap for intermediate storage */
		if (source->getHeight() == target->getHeight())
			temp = target; // write directly to the output bitmap
		else // otherwise: write to a temporary bitmap
			temp = new Bitmap(source->getPixelFormat(), source->getComponentFormat(),
				Vector2i(target->getWidth(), source->getHeight()), channels);

		#if defined(MTS_OPENMP)
			#pragma omp parallel for
		#endif
		for (int y=0; y<source->getHeight(); ++y) {
			const Scalar *srcPtr = (Scalar *) source->getUInt8Data()
				+ y * source->getWidth() * channels;
			Scalar *trgPtr = (Scalar *) temp->getUInt8Data()
				+ y * target->getWidth() * channels;

			r.resampleAndClamp(srcPtr, 1, trgPtr, 1, channels,
					(Scalar) minValue, (Scalar) maxValue);
		}

		/* Now, read from the temporary bitmap */
		source = temp;
	}

	if (source->getHeight() != target->getHeight()) {
		/* Re-sample along the Y direction */
		Resampler<Scalar> r(rfilter, bcv, source->getHeight(), target->getHeight());

		#if defined(MTS_OPENMP)
			#pragma omp parallel for
		#endif
		for (int x=0; x<source->getWidth(); ++x) {
			const Scalar *srcPtr = (Scalar *) source->getUInt8Data() + x * channels;
			Scalar *trgPtr = (Scalar *) target->getUInt8Data() + x * channels;

			r.resampleAndClamp(srcPtr, source->getWidth(), trgPtr, target->getWidth(),
				channels, (Scalar) minValue, (Scalar) maxValue);
		}
	}
}

void Bitmap::resample(const ReconstructionFilter *rfilter,
		ReconstructionFilter::EBoundaryCondition bch,
		ReconstructionFilter::EBoundaryCondition bcv,
		Bitmap *target, Float minValue, Float maxValue) const {

	Assert(getPixelFormat() == target->getPixelFormat() &&
		getComponentFormat() == target->getComponentFormat() &&
		getChannelCount() == target->getChannelCount());

	switch (m_componentFormat) {
		case EFloat16:
			mitsuba::resample<half>(rfilter, bch, bcv, this, target, minValue, maxValue);
			break;
		case EFloat32:
			mitsuba::resample<float>(rfilter, bch, bcv, this, target, minValue, maxValue);
			break;
		case EFloat64:
			mitsuba::resample<double>(rfilter, bch, bcv, this, target, minValue, maxValue);
			break;
		default:
			Log(EError, "resample(): Unsupported component type! (must be float16/32/64)");
	}
}

ref<Bitmap> Bitmap::resample(const ReconstructionFilter *rfilter,
		ReconstructionFilter::EBoundaryCondition bch,
		ReconstructionFilter::EBoundaryCondition bcv,
		const Vector2i &size, Float minValue, Float maxValue) const {
	ref<Bitmap> result = new Bitmap(m_pixelFormat, m_componentFormat, size);
	result->m_metadata = m_metadata;
	result->m_gamma = m_gamma;
	resample(rfilter, bch, bcv, result, minValue, maxValue);
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
		<< "  type = " << m_pixelFormat << endl
		<< "  componentFormat = " << m_componentFormat << endl
		<< "  size = " << m_size.toString() << endl;

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
		oss << "  }" << endl;
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
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, &png_error_func, NULL);
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

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, &png_error_func, NULL);
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

	/* First of all, check which layers are there */
	for (Imf::ChannelList::ConstIterator it = channels.begin(); it != channels.end(); ++it) {
		std::string name = boost::to_lower_copy(std::string(it.name()));

		/* Skip layers that have the wrong prefix */
		if (!boost::starts_with(name, prefix))
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
			#if SPECTRAL_SAMPLES != 3
				for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
					std::pair<Float, Float> coverage = Spectrum::getBinCoverage(i);
					if (!ch_spec[i] && boost::ends_with(name, formatString("%.2f-%.2fnm", coverage.first, coverage.second))) {
						isSpectralChannel = true;
						ch_spec[i] = it.name();
						break;
					}
				}
			#endif
			if (!isSpectralChannel)
				Log(EWarn, "readOpenEXR(): Don't know how to handle the channel named '%s'", it.name());
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
	if (spectral) {
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
	Assert(m_channelCount == (int) sourceChannels.size());

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
	char *ptr = (char *) m_data;

	ptr -= (dataWindow.min.x + dataWindow.min.y * m_size.x) * pixelStride;

	ref_vector<Bitmap> resampleBuffers(m_channelCount);
	ref<ReconstructionFilter> rfilter;

	/* Tell OpenEXR where the image data should be put */
	Imf::FrameBuffer frameBuffer;
	for (size_t i=0; i<sourceChannels.size(); ++i) {
		const char *channelName = sourceChannels[i];
		const Imf::Channel &channel = channels[channelName];
		Vector2i sampling(channel.xSampling, channel.ySampling);

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

void Bitmap::writeOpenEXR(Stream *stream,
		const std::vector<std::string> *channelNames) const {
	Log(EDebug, "Writing a %ix%i OpenEXR file", m_size.x, m_size.y);
	EPixelFormat pixelFormat = m_pixelFormat;

	#if SPECTRUM_SAMPLES == 3
		if (pixelFormat == ESpectrum)
			pixelFormat = ERGB;
		if (pixelFormat == ESpectrumAlpha)
			pixelFormat = ERGBA;
	#endif

	Properties metadata(m_metadata);
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

	Imf::ChannelList &channels = header.channels();
	if (channelNames) {
		if (channelNames->size() != (size_t) m_channelCount)
			Log(EError, "writeOpenEXR(): 'channelNames' has the wrong number of entries!");
		for (size_t i=0; i<channelNames->size(); ++i)
			channels.insert((*channelNames)[i].c_str(), Imf::Channel(compType));
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

	if (pixelFormat == ELuminanceAlpha || pixelFormat == ERGBA ||
	    pixelFormat == EXYZA || pixelFormat == ESpectrumAlpha)
		channels.insert("A", Imf::Channel(compType));

	size_t pixelStride = m_channelCount * compStride,
	       rowStride = pixelStride * m_size.x;
	char *ptr = (char *) m_data;

	Imf::FrameBuffer frameBuffer;

	if (channelNames) {
		for (size_t i=0; i<channelNames->size(); ++i) {
			frameBuffer.insert((*channelNames)[i].c_str(), Imf::Slice(compType, ptr, pixelStride, rowStride));
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

	if (pixelFormat == ELuminanceAlpha || pixelFormat == ERGBA ||
		pixelFormat == EXYZA || pixelFormat == ESpectrumAlpha)
		frameBuffer.insert("A", Imf::Slice(compType, ptr, pixelStride, rowStride));

	EXROStream ostr(stream);
	Imf::OutputFile file(ostr, header);
	file.setFrameBuffer(frameBuffer);
	file.writePixels(m_size.y);
}
#else
void Bitmap::readOpenEXR(Stream *stream) {
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
				while( (beg_run + run_count < numbytes) && (run_count < 127)
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
	for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ) {
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
	Float scale = std::abs(scaleAndOrder);
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
				for (int j=0; j<m_channelCount-1; ++j)
					*dest++ = *source++;
				source++;
			}

			stream->write(temp, sizeof(float) * m_size.x * (m_channelCount-1));
		}
	}
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
}

void Bitmap::staticShutdown() {
	FormatConverter::staticShutdown();
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

std::ostream &operator<<(std::ostream &os, const Bitmap::EPixelFormat &value) {
	switch (value) {
		case Bitmap::ELuminance: os << "luminance"; break;
		case Bitmap::ELuminanceAlpha: os << "luminance-alpha"; break;
		case Bitmap::ERGB: os << "rgb"; break;
		case Bitmap::ERGBA: os << "rgba"; break;
		case Bitmap::EXYZ: os << "xyz"; break;
		case Bitmap::EXYZA: os << "xyza"; break;
		case Bitmap::ESpectrum: os << "spectrum"; break;
		case Bitmap::ESpectrumAlpha: os << "spectrum-alpha"; break;
		case Bitmap::ESpectrumAlphaWeight: os << "spectrum-alpha-weight"; break;
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
