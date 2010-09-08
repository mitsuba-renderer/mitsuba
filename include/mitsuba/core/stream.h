/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__STREAM_H)
#define __STREAM_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract seekable stream class which defines
 * all functions to be implemented by stream subclasses,
 * such as <tt>FileStream</tt>, <tt>MemoryStream</tt>,
 * <tt>SocketStream</tt>, <tt>ConsoleStream</tt> or
 * <tt>SSHStream</tt>.
 *
 * Also supports platform-independent reads and writes 
 * (automatic endianness conversion).
 */
class MTS_EXPORT_CORE Stream : public Object {
public:
	/// Defines the byte order to use in this Stream
	enum EByteOrder {
		/// PowerPC, SPARC, Motorola 68K
		EBigEndian = 0,
		/// x86
		ELittleEndian = 1,
		/// Network byte order
		ENetworkByteOrder = EBigEndian
	};

	/// Create a new stream
	Stream();

	/// Read data from the stream
	virtual void read(void *ptr, size_t size) = 0;

	/// Write data into the stream
	virtual void write(const void *ptr, size_t size) = 0;

	/// Seek to a position inside the stream
	virtual void setPos(size_t pos) = 0;

	/// Truncate the stream to a given size
	virtual void truncate(size_t size) = 0;

	/// Get the current position inside the stream
	virtual size_t getPos() const = 0;

	/// Return the current size of the stream
	virtual size_t getSize() const = 0;

	/// Flush the stream's buffers
	virtual void flush() = 0;

	/// Can we write to the stream?
	virtual bool canWrite() const = 0;

	/// Can we read from the stream?
	virtual bool canRead() const = 0;

	/// Skip the given number of bytes
	void skip(size_t amount);

	/// Write a null-terminated string to the stream
	void writeString(const std::string &value);

	/// Write a string followed by a newline
	void writeLine(const std::string &value);

	/// Write a signed short (16 bit) to the stream
	void writeShort(short value);

	/// Write an unsigned short (16 bit) to the stream
	void writeUShort(unsigned short value);

	/// Write a signed int (32 bit) to the stream
	void writeInt(int value);
	
	/// Write an array of signed ints (32 bit) to the stream
	void writeIntArray(const int *values, size_t size);

	/// Write an unsigned int (32 bit) to the stream
	void writeUInt(unsigned int value);
	
	/// Write an array of unsigned ints (32 bit) to the stream
	void writeUIntArray(const unsigned int *values, size_t size);

	/// Write a signed int (32 bit) to the stream
	void writeLong(int64_t value);

	/// Write an array of signed ints (64 bit) to the stream
	void writeLongArray(const int64_t *values, size_t size);

	/// Write an unsigned int (32 bit) to the stream
	void writeULong(uint64_t value);
	
	/// Write an array of unsigned ints (64 bit) to the stream
	void writeULongArray(const uint64_t *values, size_t size);

	/// Write a signed character (8 bit) to the stream
	void writeChar(char value);

	/// Write an unsigned character (8 bit) to the stream
	void writeUChar(unsigned char value);

	/// Write a boolean (8 bit) to the stream
	inline void writeBool(bool value) { writeUChar(value); }
	
	/// Write a single-precision floating point number (32 bit) to the stream
	void writeSingle(float pFloat);
	
	/// Write a single-precision floating point array (32 bit) to the stream
	void writeSingleArray(const float *pFloat, size_t size);

	/// Write a double-precision floating point number (64 bit) to the stream
	void writeDouble(double pDouble);
	
	/// Write a double-precision floating point array (64 bit) to the stream
	void writeDoubleArray(const double *pDouble, size_t size);

	/// Write a floating point number (configured precision) to the stream
	inline void writeFloat(Float pFloat) {
#ifdef SINGLE_PRECISION
		writeSingle(pFloat);
#else
		writeDouble(pFloat);
#endif
	}

	/// Write an array of floating point values (configured precision) to the stream
	inline void writeFloatArray(const Float *pFloat, size_t size) {
#ifdef SINGLE_PRECISION
		writeSingleArray(pFloat, size);
#else
		writeDoubleArray(pFloat, size);
#endif
	}

	/// Return whether we are at the end of the stream
	bool isEOF() const;

	/// Read a line from the stream and return it as a string
	std::string readLine();

	/// Read a null-terminated string from the stream
	std::string readString();

	/// Read a signed short (16 bit) from the stream
	short readShort();

	/// Read an unsigned short (16 bit) from the stream
	unsigned short readUShort();

	/// Read a signed int (32 bit) from the stream
	int readInt();
	
	/// Read an array of signed ints (32 bit) from the stream
	void readIntArray(int *dst, size_t size);

	/// Read an unsigned int (32 bit) from the stream
	unsigned int readUInt();
	
	/// Read an array of unsigned ints (32 bit) from the stream
	void readUIntArray(unsigned int *dest, size_t size);

	/// Read a signed int (64 bit) from the stream
	int64_t readLong();

	/// Read an array of signed ints (64 bit) from the stream
	void readLongArray(int64_t *dst, size_t size);

	/// Read an unsigned int (64 bit) from the stream
	uint64_t readULong();

	/// Read an array of unsigned ints (64 bit) from the stream
	void readULongArray(uint64_t *dst, size_t size);

	/// Read a signed character (8 bit) from the stream
	char readChar();

	/// Read an unsigned character (8 bit) from the stream
	unsigned char readUChar();

	/// Read a boolean (8 bit) from the stream
	inline bool readBool() { return static_cast<bool> (readUChar()); }

	/// Read a single-precision floating point number (32 bit) from the stream
	float readSingle();

	/// Read a double-precision floating point array (64 bit) from the stream
	void readSingleArray(float *pDouble, size_t size);

	/// Read a double-precision floating point number (64 bit) from the stream
	double readDouble();

	/// Read a double-precision floating point array (64 bit) from the stream
	void readDoubleArray(double *pDouble, size_t size);

	/// Write a floating point number (configured precision) to the stream
	inline Float readFloat() {
#ifdef SINGLE_PRECISION
		return readSingle();
#else
		return readDouble();
#endif
	}

	/// Write an array of floating point values (configured precision) to the stream
	inline void readFloatArray(Float *pFloat, size_t size) {
#ifdef SINGLE_PRECISION
		readSingleArray(pFloat, size);
#else
		readDoubleArray(pFloat, size);
#endif
	}

	/**
	 * Copy content from this stream into another stream
	 * @param numBytes 
	 * 		The number of bytes to copy. When -1 is specified,
	 * 		copying proceeds until the end of the source stream.
	 */
	void copyTo(Stream *stream, int64_t numBytes = -1);

	/// Return a string representation
	virtual std::string toString() const;

	/// Set the stream byte order
	void setByteOrder(EByteOrder byteOrder);

	/// Return the stream byte order
	inline EByteOrder getByteOrder() const { return m_byteOrder; }

	/// Return the host byte order
	inline EByteOrder getHostByteOrder() { return m_hostByteOrder; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Stream() { }
private:
	static EByteOrder m_hostByteOrder;
	EByteOrder m_byteOrder;
};

MTS_NAMESPACE_END

#endif /* __STREAM_H */
