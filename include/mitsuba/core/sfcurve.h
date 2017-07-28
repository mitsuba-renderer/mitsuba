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
#if !defined(__MITSUBA_CORE_SFCURVE_H_)
#define __MITSUBA_CORE_SFCURVE_H_

#include <mitsuba/core/bsphere.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief 2D version of the Hilbert space-filling curve
 *
 * Based on http://voxelizator3d.wordpress.com/
 *
 * \ingroup libcore
 */
template <typename T> class HilbertCurve2D {
public:
    enum EDirection {
        ENorth,
        EEast,
        ESouth,
        EWest
    };
    typedef TVector2<T> VectorType;
    typedef TPoint2<T>  PointType;

    /// Create an empty Hilbert curve
    inline HilbertCurve2D() : m_size((T) 0), m_pos((T) 0) { }

    /// Initialize for the specified 2D size
    void initialize(const VectorType &size) {
        if (size == m_size)
            return;
        m_points.clear();
        m_points.reserve(m_size.x*m_size.y);
        m_size = size; m_pos = PointType((T) 0);
        const Float invLog2 = (Float) 1 / math::fastlog((Float) 2);
        generate(
            (int) std::ceil(invLog2 * math::fastlog((Float) std::max(m_size.x, m_size.y))),
            ENorth, EEast, ESouth, EWest);
    }

    /// Return one of the generated points
    inline const PointType &operator[](size_t idx) const {
        return m_points[idx];
    }

    /// Return a reference to the computed points
    inline const std::vector<PointType> &getPoints() const {
        return m_points;
    }

    /// Return the total number of points
    inline size_t getPointCount() const {
        return m_points.size();
    }

    /// Return the size of the underlying 2D rectangle
    inline const VectorType &getSize() const {
        return m_size;
    }
protected:
    inline void move(EDirection dir) {
        switch (dir) {
            case ENorth: m_pos.y--; break;
            case EEast:  m_pos.x++; break;
            case ESouth: m_pos.y++; break;
            case EWest:  m_pos.x--; break;
        }
    }

    void generate(int order,
            EDirection front, EDirection right,
            EDirection back, EDirection left) {
        if (order == 0) {
            if (m_pos.x < m_size.x && m_pos.y < m_size.y)
                m_points.push_back(m_pos);
        } else {
            generate(order-1, left,  back,  right, front); move(right);
            generate(order-1, front, right, back,  left);  move(back);
            generate(order-1, front, right, back,  left);  move(left);
            generate(order-1, right, front, left,  back);
        }
    }
private:
    VectorType m_size;
    PointType m_pos;
    std::vector<PointType> m_points;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SFCURVE_H_ */
