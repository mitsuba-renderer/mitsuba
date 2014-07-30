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
#if !defined(__MITSUBA_CORE_ZSTREAM_H_)
#define __MITSUBA_CORE_ZSTREAM_H_

#include <mitsuba/mitsuba.h>
#include <zlib.h>

/// Buffer size used to communicate with zlib. The larger, the better.
#define ZSTREAM_BUFSIZE 32768

MTS_NAMESPACE_BEGIN

/**
 * \brief Transparent compression/decompression stream based on \c zlib.
 *
 * This class transparently decompresses and compresses reads and writes
 * to a nested stream, respectively.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE ZStream : public Stream {
public:
	// =============================================================
	//! @{ \name Constructors
	// =============================================================
	enum EStreamType {
		/// A raw deflate stream
		EDeflateStream,
		/// A gzip-compatible stream
		EGZipStream
	};

	/// Create a new compression stream
	ZStream(Stream *childStream, EStreamType streamType = EDeflateStream,
		int level = Z_DEFAULT_COMPRESSION);

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Compression stream-specific features
	// =============================================================

	/// Return the child stream of this compression stream
	inline const Stream *getChildStream() const { return m_childStream.get(); }

	/// Return the child stream of this compression stream
	inline Stream *getChildStream() { return m_childStream; }

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Implementation of the Stream interface
	// =============================================================

	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void seek(size_t pos);
	size_t getPos() const;
	size_t getSize() const;
	void truncate(size_t size);
	void flush();
	bool canWrite() const;
	bool canRead() const;

	//! @}
	// =============================================================

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	// \brief Virtual destructor
	virtual ~ZStream();
private:
	ref<Stream> m_childStream;
	z_stream m_deflateStream, m_inflateStream;
	uint8_t m_deflateBuffer[ZSTREAM_BUFSIZE];
	uint8_t m_inflateBuffer[ZSTREAM_BUFSIZE];
	bool m_didWrite;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_ZSTREAM_H_ */
