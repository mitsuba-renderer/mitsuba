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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>

MTS_NAMESPACE_BEGIN

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

Stream::EByteOrder Stream::m_hostByteOrder = mitsuba::getByteOrder();

Stream::Stream() : m_byteOrder(m_hostByteOrder) { }

void Stream::setByteOrder(EByteOrder value) {
	m_byteOrder = value;
}

void Stream::skip(size_t amount) {
	seek(getPos() + amount);
}

void Stream::writeInt(int value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(int));
}

void Stream::writeIntArray(const int *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		int *temp = new int[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(int)*size);
		delete[] temp;
	} else {
		write(data, sizeof(int)*size);
	}
}

void Stream::writeUInt(unsigned int value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(unsigned int));
}

void Stream::writeUIntArray(const unsigned int *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		unsigned int *temp = new unsigned int[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(unsigned int)*size);
		delete[] temp;
	} else {
		write(data, sizeof(unsigned int)*size);
	}
}

void Stream::writeLong(int64_t value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(int64_t));
}

void Stream::writeLongArray(const int64_t *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		int64_t *temp = new int64_t[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(int64_t)*size);
		delete[] temp;
	} else {
		write(data, sizeof(int64_t)*size);
	}
}

void Stream::writeULong(uint64_t value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(uint64_t));
}

void Stream::writeULongArray(const uint64_t *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		uint64_t *temp = new uint64_t[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(uint64_t)*size);
		delete[] temp;
	} else {
		write(data, sizeof(uint64_t)*size);
	}
}

void Stream::writeShort(short value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(short));
}

void Stream::writeShortArray(const short *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		short *temp = new short[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(short)*size);
		delete[] temp;
	} else {
		write(data, sizeof(short)*size);
	}
}

void Stream::writeUShort(unsigned short value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(unsigned short));
}

void Stream::writeUShortArray(const unsigned short *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		unsigned short *temp = new unsigned short[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(unsigned short)*size);
		delete[] temp;
	} else {
		write(data, sizeof(unsigned short)*size);
	}
}

void Stream::writeChar(char value) {
	write(&value, sizeof(char));
}

void Stream::writeUChar(unsigned char value) {
	write(&value, sizeof(unsigned char));
}

void Stream::writeHalf(half halfValue) {
	short value = halfValue.bits();
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(short));
}

void Stream::writeHalfArray(const half *data, size_t size) {
	BOOST_STATIC_ASSERT(sizeof(half) == 2);
	if (m_byteOrder != m_hostByteOrder) {
		short *temp = new short[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i].bits());
		write(temp, sizeof(short)*size);
		delete[] temp;
	} else {
		write(data, sizeof(half)*size);
	}
}

void Stream::writeSingle(float value) {
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	write(&value, sizeof(float));
}

void Stream::writeSingleArray(const float *data, size_t size) {
	if (m_byteOrder != m_hostByteOrder) {
		float *temp = new float[size];
		for (size_t i=0; i<size; ++i)
			temp[i] = endianness_swap(data[i]);
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
			temp[i] = endianness_swap(data[i]);
		write(temp, sizeof(double)*size);
		delete[] temp;
	} else {
		write(data, sizeof(double)*size);
	}
}

void Stream::writeDouble(double pDouble) {
	if (m_byteOrder != m_hostByteOrder)
		pDouble = endianness_swap(pDouble);
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
		value = endianness_swap(value);
	return value;
}

void Stream::readLongArray(int64_t *dest, size_t size) {
	read(dest, sizeof(int64_t)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = endianness_swap(dest[i]);
	}
}

uint64_t Stream::readULong() {
	uint64_t value;
	read(&value, sizeof(uint64_t));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readULongArray(uint64_t *dest, size_t size) {
	read(dest, sizeof(uint64_t)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = endianness_swap(dest[i]);
	}
}

int Stream::readInt() {
	int value;
	read(&value, sizeof(int));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readIntArray(int *dest, size_t size) {
	read(dest, sizeof(int)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = endianness_swap(dest[i]);
	}
}

unsigned int Stream::readUInt() {
	unsigned int value;
	read(&value, sizeof(unsigned int));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readUIntArray(unsigned int *dest, size_t size) {
	read(dest, sizeof(unsigned int)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = endianness_swap(dest[i]);
	}
}

short Stream::readShort() {
	short value;
	read(&value, sizeof(short));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readShortArray(short *dest, size_t size) {
	read(dest, sizeof(short)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = endianness_swap(dest[i]);
	}
}

unsigned short Stream::readUShort() {
	unsigned short value;
	read(&value, sizeof(unsigned short));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readUShortArray(unsigned short *dest, size_t size) {
	read(dest, sizeof(unsigned short)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			dest[i] = endianness_swap(dest[i]);
	}
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

half Stream::readHalf() {
	half value;
	read(&value, sizeof(half));
	if (m_byteOrder != m_hostByteOrder)
		value.setBits(endianness_swap(value.bits()));
	return value;
}

void Stream::readHalfArray(half *data, size_t size) {
	read(data, sizeof(half)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			data[i].setBits(endianness_swap(data[i].bits()));
	}
}

float Stream::readSingle() {
	float value;
	read(&value, sizeof(float));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readSingleArray(float *data, size_t size) {
	read(data, sizeof(float)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			data[i] = endianness_swap(data[i]);
	}
}


double Stream::readDouble() {
	double value;
	read(&value, sizeof(double));
	if (m_byteOrder != m_hostByteOrder)
		value = endianness_swap(value);
	return value;
}

void Stream::readDoubleArray(double *data, size_t size) {
	read(data, sizeof(double)*size);
	if (m_byteOrder != m_hostByteOrder) {
		for (size_t i=0; i<size; ++i)
			data[i] = endianness_swap(data[i]);
	}
}

std::string Stream::readLine() {
	std::string retval;
	char data;
	bool nl = false;

	do {
		try {
			read(&data, sizeof(char));
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
		read(&data, sizeof(char));
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
		<< m_hostByteOrder
		<< ", byteOrder="
		<< m_byteOrder;

	return oss.str();
}

std::ostream &operator<<(std::ostream &os, const Stream::EByteOrder &value) {
	switch (value) {
		case Stream::ELittleEndian: os << "little-endian"; break;
		case Stream::EBigEndian: os << "big-endian"; break;
		default: os << "invalid"; break;
	}
	return os;
}

MTS_IMPLEMENT_CLASS(Stream, true, Object)
MTS_NAMESPACE_END
