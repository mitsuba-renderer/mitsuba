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

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/* Endianness utilities */
inline short swapShort(short s) {
	unsigned char b1, b2;
	b1 = s & 255;
	b2 = (s >> 8) & 255;
	return (b1 << 8) + b2;
}

inline int swapLong(int i) {
	unsigned char b1, b2, b3, b4;

	b1 = i & 255;
	b2 = (i >> 8) & 255;
	b3 = (i >> 16) & 255;
	b4 = (i >> 24) & 255;

	return ((int) b1 << 24) + ((int) b2 << 16) + ((int) b3 << 8) + b4;
}

inline int64_t swapLongLong(int64_t i) {
	int64_t result;
	unsigned char *src = (unsigned char *) &i;
	unsigned char *dst = (unsigned char *) &result;

	dst[0] = src[7];
	dst[1] = src[6];
	dst[2] = src[5];
	dst[3] = src[4];
	dst[4] = src[3];
	dst[5] = src[2];
	dst[6] = src[1];
	dst[7] = src[0];

	return result;
}


inline float swapFloat(float f) {
	union {
		float f;
		unsigned char b[4];
	} dat1, dat2;

	dat1.f = f;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	return dat2.f;
}

inline double swapDouble(double d) {
	union {
		double d;
		unsigned char b[8];
	} dat1, dat2;

	dat1.d = d;
	dat2.b[0] = dat1.b[7];
	dat2.b[1] = dat1.b[6];
	dat2.b[2] = dat1.b[5];
	dat2.b[3] = dat1.b[4];
	dat2.b[4] = dat1.b[3];
	dat2.b[5] = dat1.b[2];
	dat2.b[6] = dat1.b[1];
	dat2.b[7] = dat1.b[0];
	return dat2.d;
}

static Stream::EByteOrder getByteOrder() {
	union {
		uint8_t  charValue[2];
		uint16_t shortValue;
	};
	charValue[0] = 1;
	charValue[1] = 0;

	if (shortValue == 1)
		return Stream::ELittleEndian;
	else
		return Stream::EBigEndian;
}

static const char *byteOrderToString(Stream::EByteOrder byteOrder) {
	if (byteOrder == Stream::ELittleEndian)
		return "little-endian";
	else if (byteOrder == Stream::EBigEndian)
		return "big-endian";
	else
		return "unknown";
}

	
Stream::EByteOrder Stream::m_hostByteOrder = mitsuba::getByteOrder();

Stream::Stream() {
	m_byteOrder = m_hostByteOrder;
}

void Stream::setByteOrder(EByteOrder pOrder) {
	m_byteOrder = pOrder;
}

void Stream::skip(size_t amount) {
	setPos(getPos() + amount);
}

void Stream::writeInt(int value) {
	if (m_byteOrder != m_hostByteOrder)
		value = swapLong(value);
	write(&value, sizeof(int));
}

void Stream::writeIntArray(const int *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		int *temp = new int[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = swapLong(data[i]);
		write(temp, sizeof(int)*size);
		delete[] temp;
	} else {
		write(data, sizeof(int)*size);
	}
}

void Stream::writeUInt(unsigned int value) {
	if (m_byteOrder != m_hostByteOrder)
		value = swapLong(value);
	write(&value, sizeof(unsigned int));
}

void Stream::writeUIntArray(const unsigned int *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		unsigned int *temp = new unsigned int[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = swapLong(data[i]);
		write(temp, sizeof(unsigned int)*size);
		delete[] temp;
	} else {
		write(data, sizeof(unsigned int)*size);
	}
}

void Stream::writeLong(int64_t value) {
	if (m_byteOrder != m_hostByteOrder)
		value = swapLongLong(value);
	write(&value, sizeof(int64_t));
}

void Stream::writeLongArray(const int64_t *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		int64_t *temp = new int64_t[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = swapLongLong(data[i]);
		write(temp, sizeof(int64_t)*size);
		delete[] temp;
	} else {
		write(data, sizeof(int64_t)*size);
	}
}

void Stream::writeULong(uint64_t value) {
	if (m_byteOrder != m_hostByteOrder)
		value = swapLongLong(value);
	write(&value, sizeof(uint64_t));
}

void Stream::writeULongArray(const uint64_t *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		uint64_t *temp = new uint64_t[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = swapLongLong(data[i]);
		write(temp, sizeof(uint64_t)*size);
		delete[] temp;
	} else {
		write(data, sizeof(uint64_t)*size);
	}
}

void Stream::writeShort(short value) {
	if (m_byteOrder != m_hostByteOrder)
		value = swapShort(value);
	write(&value, sizeof(short));
}

void Stream::writeUShort(unsigned short value) {
	if (m_byteOrder != m_hostByteOrder)
		value = swapShort(value);
	write(&value, sizeof(unsigned short));
}

void Stream::writeChar(char value) {
	write(&value, sizeof(char));
}

void Stream::writeUChar(unsigned char value) {
	write(&value, sizeof(unsigned char));
}

void Stream::writeSingle(float pFloat) {
	if (m_byteOrder != m_hostByteOrder) {
		pFloat = swapFloat(pFloat);
	}
	write(&pFloat, sizeof(float));
}

void Stream::writeSingleArray(const float *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		float *temp = new float[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = swapFloat(data[i]);
		write(temp, sizeof(float)*size);
		delete[] temp;
	} else {
		write(data, sizeof(float)*size);
	}
}

void Stream::writeDoubleArray(const double *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		double *temp = new double[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = swapDouble(data[i]);
		write(temp, sizeof(double)*size);
		delete[] temp;
	} else {
		write(data, sizeof(double)*size);
	}
}

void Stream::writeDouble(double pDouble) {
	if (m_byteOrder != m_hostByteOrder)
		pDouble = swapDouble(pDouble);
	write(&pDouble, sizeof(double));
}

void Stream::writeString(const std::string &value) {
	write(value.c_str(), value.length() + 1);
}

void Stream::writeLine(const std::string &value) {
	write(value.c_str(), value.length());
	writeChar('\n');
}

bool Stream::isEOF() const {
	return (getPos() == getSize());
}

int64_t Stream::readLong() {
	int64_t value;
	read(&value, sizeof(int64_t));
	if (m_byteOrder != m_hostByteOrder)
		value = swapLongLong(value);
	return value;
}

void Stream::readLongArray(int64_t *dest, size_t size) {
	read(dest, sizeof(int64_t)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = swapLongLong(dest[i]);
	}
}

uint64_t Stream::readULong() {
	uint64_t value;
	read(&value, sizeof(uint64_t));
	if (m_byteOrder != m_hostByteOrder)
		value = swapLongLong(value);
	return value;
}

void Stream::readULongArray(uint64_t *dest, size_t size) {
	read(dest, sizeof(uint64_t)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = swapLongLong(dest[i]);
	}
}

int Stream::readInt() {
	int value;
	read(&value, sizeof(int));
	if (m_byteOrder != m_hostByteOrder)
		value = swapLong(value);
	return value;
}

void Stream::readIntArray(int *dest, size_t size) {
	read(dest, sizeof(int)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = swapLong(dest[i]);
	}
}

unsigned int Stream::readUInt() {
	unsigned int value;
	read(&value, sizeof(unsigned int));
	if (m_byteOrder != m_hostByteOrder)
		value = swapLong(value);
	return value;
}

void Stream::readUIntArray(unsigned int *dest, size_t size) {
	read(dest, sizeof(unsigned int)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = swapLong(dest[i]);
	}
}

short Stream::readShort() {
	short value;
	read(&value, sizeof(short));
	if (m_byteOrder != m_hostByteOrder)
		value = swapShort(value);
	return value;
}

unsigned short Stream::readUShort() {
	unsigned short value;
	read(&value, sizeof(unsigned short));
	if (m_byteOrder != m_hostByteOrder)
		value = swapShort(value);
	return value;
}

char Stream::readChar() {
	char value;
	read(&value, sizeof(char));
	return value;
}

unsigned char Stream::readUChar() {
	unsigned char value;
	read(&value, sizeof(unsigned char));
	return value;
}

float Stream::readSingle() {
	float value;
	read(&value, sizeof(float));
	if (m_byteOrder != m_hostByteOrder)
		value = swapFloat(value);
	return value;
}

void Stream::readSingleArray(float *data, size_t size) {
	read(data, sizeof(float)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			data[i] = swapFloat(data[i]);
	}
}

double Stream::readDouble() {
	double value;
	read(&value, sizeof(double));
	if (m_byteOrder != m_hostByteOrder) 
		value = swapDouble(value);
	return value;
}

void Stream::readDoubleArray(double *data, size_t size) {
	read(data, sizeof(double)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			data[i] = swapDouble(data[i]);
	}
}

std::string Stream::readLine() {
	std::string retval;
	char data;
	bool nl = false;
	
	do {
		try {
			data = readChar();
		} catch (std::exception &e) {
			if (getPos() == getSize()) {
				if (retval.size() != 0) {
					return retval;
				}
			}
			throw e;
		}
		if (data != 13 && data != 10)
			retval += data;
		else if (data == 10)
			nl = true;
	} while (!nl);
	return retval;
}

std::string Stream::readString() {
	std::string retval;
	char data;
	
	do {
		data = readChar();
		if (data != 0)
			retval += data;
	} while (data != 0);
	return retval;
}

void Stream::copyTo(Stream *stream, int64_t numBytes) {
	const int block_size = 512;
	char data[block_size];
	size_t copied = 0;

	size_t amount = (numBytes == -1) ? (getSize() - getPos()) : (size_t) numBytes;
	for (size_t i=0; i<amount; i+=block_size) {
		size_t blockSize = (i + block_size) <= amount ? block_size : amount-i;
	
		read(data, blockSize);
		copied += blockSize;
		stream->write(data, blockSize);
	}
}

std::string Stream::toString() const {
	std::ostringstream oss;

	oss << "hostByteOrder="
		<< byteOrderToString(m_hostByteOrder)
		<< ", byteOrder="
		<< byteOrderToString(m_byteOrder);

	return oss.str();
}

MTS_IMPLEMENT_CLASS(Stream, true, Object)
MTS_NAMESPACE_END
