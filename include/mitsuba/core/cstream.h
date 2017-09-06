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
#if !defined(__MITSUBA_CORE_CSTREAM_H_)
#define __MITSUBA_CORE_CSTREAM_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Stream-style interface to the default stdin/stdout console streams
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE ConsoleStream : public Stream {
public:
    // =============================================================
    //! @{ \name Constructors
    // =============================================================

    /// Create a new console stream
    ConsoleStream();

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
    /// Virtual destructor
    virtual ~ConsoleStream();
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_CSTREAM_H_ */
