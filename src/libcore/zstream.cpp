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

#include <mitsuba/core/zstream.h>

MTS_NAMESPACE_BEGIN

ZStream::ZStream(Stream *childStream, EStreamType streamType, int level)
		: m_childStream(childStream), m_didWrite(false) {
	m_deflateStream.zalloc = Z_NULL;
	m_deflateStream.zfree = Z_NULL;
	m_deflateStream.opaque = Z_NULL;

	int windowBits = 15 + (streamType == EGZipStream ? 16 : 0);

	int retval = deflateInit2(&m_deflateStream, level,
		Z_DEFLATED, windowBits, 8, Z_DEFAULT_STRATEGY);

	if (retval != Z_OK)
		Log(EError, "Could not initialize ZLIB: error code %i", retval);

	m_inflateStream.zalloc = Z_NULL;
	m_inflateStream.zfree = Z_NULL;
	m_inflateStream.opaque = Z_NULL;
	m_inflateStream.avail_in = 0;
	m_inflateStream.next_in = Z_NULL;

	retval = inflateInit2(&m_inflateStream, windowBits);
	if (retval != Z_OK)
		Log(EError, "Could not initialize ZLIB: error code %i", retval);
}

bool ZStream::canWrite() const {
	return m_childStream->canWrite();
}

bool ZStream::canRead() const {
	return m_childStream->canRead();
}

void ZStream::seek(size_t pos) {
	Log(EError, "seek(): unsupported in a ZLIB stream!");
}

size_t ZStream::getPos() const {
	Log(EError, "getPos(): unsupported in a ZLIB stream!");
	return 0;
}

size_t ZStream::getSize() const {
	Log(EError, "getSize(): unsupported in a ZLIB stream!");
	return 0;
}

void ZStream::truncate(size_t size) {
	Log(EError, "truncate(): unsupported in a ZLIB stream!");
}

void ZStream::flush() {
	Log(EError, "flush(): not implemented!");
}

void ZStream::write(const void *ptr, size_t size) {
	m_deflateStream.avail_in = (uInt) size;
	m_deflateStream.next_in = (uint8_t *) ptr;

	int outputSize = 0;

	do {
		m_deflateStream.avail_out = sizeof(m_deflateBuffer);
		m_deflateStream.next_out = m_deflateBuffer;

		int retval = deflate(&m_deflateStream, Z_NO_FLUSH);
		if (retval == Z_STREAM_ERROR)
			Log(EError, "deflate(): stream error!");

		outputSize = sizeof(m_deflateBuffer) - m_deflateStream.avail_out;

		m_childStream->write(m_deflateBuffer, outputSize);
	} while (outputSize != 0);

	Assert(m_deflateStream.avail_in == 0);
	m_didWrite = true;
}

void ZStream::read(void *ptr, size_t size) {
	uint8_t *targetPtr = (uint8_t *) ptr;
	while (size > 0) {
		if (m_inflateStream.avail_in == 0) {
			size_t remaining = m_childStream->getSize() - m_childStream->getPos();
			m_inflateStream.next_in = m_inflateBuffer;
			m_inflateStream.avail_in = (uInt) std::min(remaining, sizeof(m_inflateBuffer));
			if (m_inflateStream.avail_in == 0)
				Log(EError, "Read less data than expected (%i more bytes required)", size);
			m_childStream->read(m_inflateBuffer, m_inflateStream.avail_in);
		}

		m_inflateStream.avail_out = (uInt) size;
		m_inflateStream.next_out = targetPtr;

		int retval = inflate(&m_inflateStream, Z_NO_FLUSH);
		switch (retval) {
			case Z_STREAM_ERROR:
				Log(EError, "inflate(): stream error!");
			case Z_NEED_DICT:
				Log(EError, "inflate(): need dictionary!");
			case Z_DATA_ERROR:
				Log(EError, "inflate(): data error!");
			case Z_MEM_ERROR:
				Log(EError, "inflate(): memory error!");
		};

		size_t outputSize = size - (size_t) m_inflateStream.avail_out;
		targetPtr += outputSize;
		size -= outputSize;

		if (size > 0 && retval == Z_STREAM_END)
			Log(EError, "inflate(): attempting to read past the end of the stream!");
	}
}

ZStream::~ZStream() {
	if (m_didWrite) {
		m_deflateStream.avail_in = 0;
		m_deflateStream.next_in = NULL;
		int outputSize = 0;

		do {
			m_deflateStream.avail_out = sizeof(m_deflateBuffer);
			m_deflateStream.next_out = m_deflateBuffer;

			int retval = deflate(&m_deflateStream, Z_FINISH);
			if (retval == Z_STREAM_ERROR)
				Log(EError, "deflate(): stream error!");

			outputSize = sizeof(m_deflateBuffer) - m_deflateStream.avail_out;

			m_childStream->write(m_deflateBuffer, outputSize);
		} while (outputSize != 0);
	}

	deflateEnd(&m_deflateStream);
	inflateEnd(&m_inflateStream);
}

std::string ZStream::toString() const {
	std::ostringstream oss;
	oss << "ZStream[" << endl
		<< "  childStream = " << indent(m_childStream->toString()) << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(ZStream, false, Stream)
MTS_NAMESPACE_END
