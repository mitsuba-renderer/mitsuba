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
#if !defined(__MITSUBA_CORE_BARRAY_H_)
#define __MITSUBA_CORE_BARRAY_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Blocked generic 2D array data type
 *
 * This class implements a blocked 2D array for cache-efficient
 * access to two-dimensional data.
 */
template <typename Value, size_t logblockSize = 2> class BlockedArray {
public:
    static const size_t blockSize = 1 << logblockSize;

    /// Create an unitialized blocked array
    BlockedArray() : m_data(NULL), m_size(-1), m_owner(false) { }

    /**
     * \brief Allocate memory for a new blocked array of
     * the specified width and height
     */
    BlockedArray(const Vector2i &size) : m_data(NULL),
            m_size(-1), m_owner(false) {
        alloc(size);
    }

    /**
     * \brief Allocate memory for a blocked array of
     * the specified width and height
     */
    void alloc(const Vector2i &size) {
        if (m_data && m_owner)
            delete[] m_data;

        m_xBlocks = (size.x + blockSize - 1) / blockSize;
        m_yBlocks = (size.y + blockSize - 1) / blockSize;
        m_data = (Value *) allocAligned(m_xBlocks * m_yBlocks
            * blockSize * blockSize * sizeof(Value));
        m_owner = true; /* We own this pointer */
        m_size = size;
    }

    /**
     * \brief Initialize the blocked array with a given pointer
     * and array size.
     *
     * This is useful in case memory has already been allocated.
     */
    void map(void *ptr, const Vector2i &size) {
        if (m_data && m_owner)
            delete[] m_data;

        m_xBlocks = (size.x + blockSize - 1) / blockSize;
        m_yBlocks = (size.y + blockSize - 1) / blockSize;
        m_data = (Value *) ptr;
        m_owner = false; /* We do not own this pointer */
        m_size = size;
    }

    /**
     * \brief Initialize the contents of the blocked array with values
     * from a non-blocked source in row-major order.
     *
     * \remark This function performs type casts when <tt>Value != AltValue</tt>
     */
    template <typename AltValue> void init(const AltValue *data) {
        for (int y=0; y<m_size.y; ++y)
            for (int x=0; x<m_size.x; ++x)
                    (*this)(x, y) = Value(*data++);
    }

    /**
     * \brief Initialize the contents of the blocked array with values
     * from a non-blocked source in row-major order and collect component-wise
     * minimum, maximum, and average information.
     *
     * Assumes that \c AltValue is some kind of \c TVector or \c TSpectrum instance.
     *
     * \remark This function performs type casts when <tt>Value != AltValue</tt>
     */
    template <typename AltValue> void init(const AltValue *data,
            AltValue &min_, AltValue &max_, AltValue &avg_) {
        typedef typename AltValue::Scalar Scalar;

        AltValue
            min(+std::numeric_limits<Scalar>::infinity()),
            max(-std::numeric_limits<Scalar>::infinity()),
            avg((Scalar) 0);

        for (int y=0; y<m_size.y; ++y) {
            for (int x=0; x<m_size.x; ++x) {
                const AltValue &value = *data++;
                for (int i=0; i<AltValue::dim; ++i) {
                    min[i]  = std::min(min[i], value[i]);
                    max[i]  = std::max(max[i], value[i]);
                    avg[i] += value[i];
                }
                (*this)(x, y) = Value(value);
            }
        }
        min_ = min;
        max_ = max;
        avg_ = avg / (Scalar) (m_size.x * m_size.y);
    }

    /**
     * \brief Copy the contents of the blocked array to a non-blocked
     * destination buffer in row-major order.
     *
     * This is effectively the opposite of \ref init().
     *
     * \remark This function performs type casts when <tt>Value != AltValue</tt>
     */
    template <typename AltValue> void copyTo(AltValue *data) const {
        for (int y=0; y<m_size.y; ++y)
            for (int x=0; x<m_size.x; ++x)
                    *data++ = AltValue((*this)(x, y));
    }


    /**
     * \brief Zero out unused memory portions
     *
     * This is useful in case we want to write the internal representation to
     * disk and avoid accessing uninitialized memory (otherwise, valgrind
     * or other similar tools will complain..)
     */
    void cleanup() {
        size_t unusedColsRight = m_xBlocks * blockSize - m_size.x,
               unusedRowsBottom = m_yBlocks * blockSize - m_size.y;

        for (int y=0; y < (int) unusedRowsBottom; ++y)
            memset(&(this->operator()(0, y + m_size.y)), 0, sizeof(Value) * m_xBlocks * blockSize);

        if (unusedColsRight > 0) {
            for (int y=0; y < (int) (m_yBlocks * blockSize); ++y)
                for (int x=0; x < (int) unusedColsRight; ++x)
                    memset(&(this->operator()(x + m_size.x, y)), 0, sizeof(Value));
        }
    }

    /// Return the size of the array
    inline const Vector2i &getSize() const { return m_size; }

    /// Return the hypothetical heap memory requirements of a blocked array for the given size
    inline static size_t bufferSize(const Vector2i &size) {
        size_t xBlocks = (size.x + blockSize - 1) / blockSize,
               yBlocks = (size.y + blockSize - 1) / blockSize,
               fullWidth = xBlocks * blockSize,
               fullHeight = yBlocks * blockSize;
        return fullWidth * fullHeight * sizeof(Value);
    }

    /// Return the size of the allocated buffer
    inline size_t getBufferSize() const {
        return m_xBlocks * m_yBlocks * blockSize * blockSize * sizeof(Value);
    }

    /// Return the width of the array
    inline int getWidth() const { return m_size.x; }

    /// Return the height of the array
    inline int getHeight() const { return m_size.y; }

    /// Release all memory
    ~BlockedArray() {
        if (m_data && m_owner)
            freeAligned(m_data);
    }

    /// Access the specified entry
    inline Value &operator()(int x, int y) {
        size_t xb = getBlock(x),  yb = getBlock(y),
               xo = getOffset(x), yo = getOffset(y);

        return m_data[
            // Offset to block
            blockSize * blockSize * (xb + yb * m_xBlocks) +
            // Offset within block
            blockSize * yo + xo
        ];
    }

    /// Access the specified entry (const version)
    inline const Value &operator()(int x, int y) const {
        size_t xb = getBlock(x),  yb = getBlock(y),
               xo = getOffset(x), yo = getOffset(y);

        return m_data[
            // Offset to block
            blockSize * blockSize * (xb + yb * m_xBlocks) +
            // Offset within block
            blockSize * yo + xo
        ];
    }

    /// Return a pointer to the internal representation
    inline Value *getData() { return m_data; }

    /// Return a pointer to the internal representation (const version)
    inline const Value *getData() const { return m_data; }
protected:
    /// Determine the index of the block which contains the given global index
    inline size_t getBlock(int a) const { return (size_t) (a >> logblockSize); }

    /// Determine the offset within the block that contains the given global index
    inline size_t getOffset(int a) const { return (size_t) (a & (blockSize - 1)); }
private:
    Value *m_data;
    Vector2i m_size;
    size_t m_xBlocks, m_yBlocks;
    bool m_owner;
};

/**
 * \brief Linear (i.e. non-blocked) generic 2D array data type
 *
 * This class implements a linearly stored 2D array. It is mainly meant
 * as a drop-in replacement for \ref BlockedArray so that the performance
 * tradeoffs can be benchmarked.
 */
template <typename Value> class LinearArray {
public:
    /// Create an unitialized linear array
    LinearArray() : m_data(NULL), m_size(-1), m_owner(false) { }

    /**
     * \brief Allocate memory for a new linear array of
     * the specified width and height
     */
    LinearArray(const Vector2i &size) : m_data(NULL),
            m_size(-1), m_owner(false) {
        alloc(size);
    }

    /**
     * \brief Allocate memory for a linear array of
     * the specified width and height
     */
    void alloc(const Vector2i &size) {
        if (m_data && m_owner)
            delete[] m_data;

        size_t arraySize = (size_t) size.x * (size_t) size.y * sizeof(Value);
        m_data = (Value *) allocAligned(arraySize);
        m_size = size;
    }

    /*
     * \brief Initialize the linear array with a given pointer
     * and array size.
     *
     * This is useful in case memory has already been allocated.
     */
    void map(void *ptr, const Vector2i &size) {
        if (m_data && m_owner)
            delete[] m_data;

        m_data = (Value *) ptr;
        m_owner = false; /* We do not own this pointer */
        m_size = size;
    }

    /**
     * \brief Initialize the contents of the linear array with values
     * from a non-blocked source in row-major order.
     *
     * \remark This function performs type casts when <tt>Value != AltValue</tt>
     */
    template <typename AltValue> void init(const AltValue *data) {
        Value *ptr = m_data;
        for (int y=0; y<m_size.y; ++y)
            for (int x=0; x<m_size.x; ++x)
                    *ptr++ = Value(*data++);
    }

    /**
     * \brief Initialize the contents of the linear array with values
     * from a non-blocked source in row-major order and collect component-wise
     * minimum, maximum, and average information.
     *
     * Assumes that \c AltValue is some kind of \c TVector or \c TSpectrum instance.
     *
     * \remark This function performs type casts when <tt>Value != AltValue</tt>
     */
    template <typename AltValue> void init(const AltValue *data,
            AltValue &min_, AltValue &max_, AltValue &avg_) {
        typedef typename AltValue::Scalar Scalar;

        AltValue
            min(+std::numeric_limits<Scalar>::infinity()),
            max(-std::numeric_limits<Scalar>::infinity()),
            avg((Scalar) 0);

        Value *ptr = m_data;

        for (int y=0; y<m_size.y; ++y) {
            for (int x=0; x<m_size.x; ++x) {
                const AltValue &value = *data++;
                for (int i=0; i<AltValue::dim; ++i) {
                    min[i]  = std::min(min[i], value[i]);
                    max[i]  = std::max(max[i], value[i]);
                    avg[i] += value[i];
                }
                *ptr++ = Value(value);
            }
        }
        min_ = min;
        max_ = max;
        avg_ = avg / (Scalar) (m_size.x * m_size.y);
    }

    /**
     * \brief Copy the contents of the linear array to a non-blocked
     * destination buffer in row-major order.
     *
     * This is effectively the opposite of \ref init().
     *
     * \remark This function performs type casts when <tt>Value != AltValue</tt>
     */
    template <typename AltValue> void copyTo(AltValue *data) const {
        Value *ptr = m_data;
        for (int y=0; y<m_size.y; ++y)
            for (int x=0; x<m_size.x; ++x)
                    *data++ = AltValue(*ptr++);
    }


    /**
     * \brief Zero out unused memory portions
     *
     * Since this is a non-blocked array, this function does nothing
     */
    void cleanup() {
    }

    /// Return the size of the array
    inline const Vector2i &getSize() const { return m_size; }

    /// Return the hypothetical heap memory requirements of a blocked array for the given size
    inline static size_t bufferSize(const Vector2i &size) {
        return (size_t) size.x * (size_t) size.y * sizeof(Value);
    }

    /// Return the size of the allocated buffer
    inline size_t getBufferSize() const {
        return (size_t) m_size.x * (size_t) m_size.y * sizeof(Value);
    }

    /// Return the width of the array
    inline int getWidth() const { return m_size.x; }

    /// Return the height of the array
    inline int getHeight() const { return m_size.y; }

    /// Release all memory
    ~LinearArray() {
        if (m_data && m_owner)
            freeAligned(m_data);
    }

    /// Access the specified entry
    inline Value &operator()(int x, int y) {
        return m_data[x + (size_t) m_size.x * y];
    }

    /// Access the specified entry (const version)
    inline const Value &operator()(int x, int y) const {
        return m_data[x + (size_t) m_size.x * y];
    }

    /// Return a pointer to the internal representation
    inline Value *getData() { return m_data; }

    /// Return a pointer to the internal representation (const version)
    inline const Value *getData() const { return m_data; }

private:
    Value *m_data;
    Vector2i m_size;
    bool m_owner;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_BARRAY_H_ */
