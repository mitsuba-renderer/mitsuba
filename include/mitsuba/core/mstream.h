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
#if !defined(__MITSUBA_CORE_MSTREAM_H_)
#define __MITSUBA_CORE_MSTREAM_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple memory buffer-based stream with automatic memory management
 *
 * The underlying memory storage of this implementation dynamically expands
 * as data is written to the stream.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE MemoryStream : public Stream {
public:
    // =============================================================
    //! @{ \name Constructors
    // =============================================================

    /// Create a new memory stream
    MemoryStream(size_t initialSize = 512);

    /**
     * \brief Create a memory stream, which operates on a
     * pre-allocated buffer.
     *
     * A memory stream created in this way will never resize the
     * underlying buffer. An exception is thrown e.g. when attempting
     * to extend its size
     *
     * \remark This constructor is not available in the python bindings
     */
    MemoryStream(void *ptr, size_t size);

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Memory stream-specific features
    // =============================================================

    /// Return the underlying data
    inline uint8_t *getData() { return m_data; }

    /// Return the underlying data (const version)
    inline const uint8_t *getData() const { return m_data; }

    /// Return the underlying data at the current position
    inline uint8_t *getCurrentData() { return m_data + m_pos; }

    /// Return the underlying data at the current position (const version)
    inline const uint8_t *getCurrentData() const { return m_data + m_pos; }

    /// Set size and position to zero without changing the underlying buffer
    void reset();

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
    void resize(size_t newSize);

    // \brief Virtual destructor
    virtual ~MemoryStream();
protected:
    size_t m_capacity;
    size_t m_size;
    size_t m_pos;
    bool m_ownsBuffer;
    uint8_t *m_data;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_MSTREAM_H_ */
