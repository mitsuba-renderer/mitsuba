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
#if !defined(__MITSUBA_BIDIR_GEODIST2_H_)
#define __MITSUBA_BIDIR_GEODIST2_H_

#include <mitsuba/bidir/path.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Clamped two-tailed geometric distribution
 *
 * This class implements a specialized clamped two-tailed geometric
 * distribution with support for sample generation and evaluation
 * of the probability mass and cumulative distribution functions.
 *
 * The precise form of the distribution is
 *
 * \f[
 * P(i) = \begin{cases}
 *     c b^{-|i-o|},&\text{if } i\ge l\text{ and }i\le r\\
 *     0,&\mathrm{otherwise}
 * \end{cases}
 * \f]
 *
 * where \f$b\in\mathbb{R}\f$ is the base parameter of the distribution,
 * \f$[l,r]\subseteq\mathbb{Z}\f$ denotes the domain of the probability
 * mass function, \f$o\in\mathbb{Z}\f$ is offset, and \f$c\f$ is a suitably
 * chosen normalization constant.
 *
 * This function is used to propose bidirectional mutations; see the
 * MLT writeup for details.
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
struct TwoTailedGeoDistr {
public:
    /// Create a new two-tailed distribution for the given base constant
    TwoTailedGeoDistr(Float base) : m_base(base) {
        m_baseNormalization = 1.0f / (Float) (base+1);
        m_invLogBase = 1.0f / std::log(base);
    }

    /// Configure the domain and center of the distribution
    void configure(int center, int start, int end) {
        m_center = center;
        m_start = start - center;
        m_end = end - center;
        m_offset = R(m_start - 1);
        m_normalization = R(m_end) - m_offset;
    }

    /// Evaluate the probability mass function at position \c i
    inline Float pmf(int i) const {
        i -= m_center;

        if (i < m_start || i > m_end)
            return 0.0f;

        return r(i) / m_normalization;
    }

    /// Evaluate the cumulative distribution function at position \c i
    inline Float cdf(int i) const {
        i -= m_center;

        if (i < m_start)
            return 0.0f;
        else if (i > m_end)
            i = m_end;

        return (R(i) - m_offset) / m_normalization;
    }

    /// Draw a position according to the probability mass function
    inline int sample(Float xi) const {
        return std::max(m_start,
            Rinv(xi * m_normalization + m_offset)) + m_center;
    }

protected:
    inline Float r(int i) const {
        return (m_base-1) * m_baseNormalization
            * std::pow(m_base, - (Float) std::abs(i));
    }

    inline Float R(int i) const {
        if (i <= 0)
            return std::pow(m_base, (Float) (i+1)) * m_baseNormalization;
        else
            return 1-std::pow(m_base, - (Float) i) * m_baseNormalization;
    }

    inline int Rinv(Float x) const {
        Float result;
        if (x < m_base * m_baseNormalization)
            result = std::log((1+m_base) * x) * m_invLogBase - 1;
        else
            result = -std::log((1+m_base) * (1-x)) * m_invLogBase;
        return (int) std::ceil(result);
    }
private:
    Float m_base, m_invLogBase, m_baseNormalization;
    Float m_normalization, m_offset;
    int m_center, m_start, m_end;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_GEODIST2_H_ */
