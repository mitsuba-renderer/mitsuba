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
#if !defined(__MITSUBA_RENDER_RANGE_H_)
#define __MITSUBA_RENDER_RANGE_H_

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief A work unit specifying a range of some quantity to be processed.
 *
 * An example usage is in \ref ParticleProcess, where this class specifies
 * sequences of particles to be traced.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER RangeWorkUnit : public WorkUnit {
public:
    inline void set(const WorkUnit *wu) {
        const RangeWorkUnit *other = static_cast<const RangeWorkUnit *>(wu);
        m_rangeStart = other->m_rangeStart;
        m_rangeEnd = other->m_rangeEnd;
    }

    inline void load(Stream *stream) {
        m_rangeStart = stream->readSize();
        m_rangeEnd = stream->readSize();
    }

    inline void save(Stream *stream) const {
        stream->writeSize(m_rangeStart);
        stream->writeSize(m_rangeEnd);
    }

    inline std::string toString() const {
        std::ostringstream oss;
        oss << "RangeWorkUnit[rangeStart=" << m_rangeStart
            << ", rangeEnd=" << m_rangeEnd << "]";
        return oss.str();
    }

    inline void setRange(size_t start, size_t end) {
        m_rangeStart = start;
        m_rangeEnd = end;
    }

    inline size_t getRangeStart() const {
        return m_rangeStart;
    }

    inline size_t getRangeEnd() const {
        return m_rangeEnd;
    }

    inline size_t getSize() const {
        return m_rangeEnd - m_rangeStart + 1;
    }

    MTS_DECLARE_CLASS()
private:
    size_t m_rangeStart, m_rangeEnd;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_RANGE_H_ */
