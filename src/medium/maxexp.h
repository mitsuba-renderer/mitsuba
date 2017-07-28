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

#if !defined(__MAXEXP_H)
#define __MAXEXP_H

#if defined(_MSC_VER)
# include <functional>
#endif

MTS_NAMESPACE_BEGIN

class MaxExpDist {
public:
    MaxExpDist(const std::vector<Float> &sigmaT)
     : m_sigmaT(sigmaT), m_cdf(sigmaT.size()+1), m_intervalStart(sigmaT.size()) {
        /* Sort the coefficients in decreasing order */
        std::sort(m_sigmaT.begin(), m_sigmaT.end(), std::greater<Float>());

        m_cdf[0] = 0;
        for (size_t i=0; i<m_sigmaT.size(); ++i) {
            if (i > 0 && m_sigmaT[i] == m_sigmaT[i-1])
                SLog(EError, "Internal error: sigmaT must vary across channels");
            /* Integrate max(f_1(t), .., f_n(t)) on [0, \infty]*/
            Float lower = (i==0) ? -1 : -std::pow((m_sigmaT[i]/m_sigmaT[i-1]),
                        -m_sigmaT[i] / (m_sigmaT[i]-m_sigmaT[i-1]));
            Float upper = (i==m_sigmaT.size()-1) ? 0 : -std::pow((m_sigmaT[i+1]/m_sigmaT[i]),
                        -m_sigmaT[i] / (m_sigmaT[i+1]-m_sigmaT[i]));
            m_cdf[i+1] = m_cdf[i] + (upper - lower);

            /* Store the interval covered by each f_i */
            m_intervalStart[i] = (i == 0) ? 0
                : math::fastlog(m_sigmaT[i]/m_sigmaT[i-1]) / (m_sigmaT[i]-m_sigmaT[i-1]);
        }

        /* Turn into a discrete CDF and keep the normalization factor */
        m_normalization = m_cdf[m_cdf.size()-1];
        m_invNormalization = 1 / m_normalization;

        for (size_t i=0; i<m_cdf.size(); ++i)
            m_cdf[i] *= m_invNormalization;
    }

    Float sample(Float u, Float &pdf) const {
        /* Find the f_i for this sample */
        const Float *lowerBound = std::lower_bound(&m_cdf[0], &m_cdf[m_cdf.size()], u);
        int index = std::max(0, (int) (lowerBound - &m_cdf[0]) - 1);
        SAssert(index >= 0 && index < (int) m_sigmaT.size());

        /* Sample according to f_i */
        Float t = -math::fastlog(math::fastexp(-m_intervalStart[index] * m_sigmaT[index])
            - m_normalization * (u - m_cdf[index])) / m_sigmaT[index];

        /* Compute the probability of this sample */
        pdf = m_sigmaT[index] * math::fastexp(-m_sigmaT[index] * t) * m_invNormalization;

        return t;
    }

    Float pdf(Float t) const {
        const Float *lowerBound = std::lower_bound(&m_intervalStart[0],
                &m_intervalStart[m_intervalStart.size()], t);
        int index = std::max(0, (int) (lowerBound - &m_intervalStart[0]) - 1);
        SAssert(index >= 0 && index < (int) m_sigmaT.size());
        return m_sigmaT[index] * math::fastexp(-m_sigmaT[index] * t) * m_invNormalization;
    }

    Float cdf(Float t) const {
        const Float *lowerBound = std::lower_bound(&m_intervalStart[0],
                &m_intervalStart[m_intervalStart.size()], t);
        int index = std::max(0, (int) (lowerBound - &m_intervalStart[0]) - 1);
        SAssert(index >= 0 && index < (int) m_sigmaT.size());

        Float lower = (index==0) ? -1 : -std::pow((m_sigmaT[index]/m_sigmaT[index-1]),
                    -m_sigmaT[index] / (m_sigmaT[index]-m_sigmaT[index-1]));
        Float upper = -math::fastexp(-m_sigmaT[index] * t);

        return m_cdf[index] + (upper - lower) * m_invNormalization;
    }
private:
    std::vector<Float> m_sigmaT;
    std::vector<Float> m_cdf;
    std::vector<Float> m_intervalStart;
    Float m_normalization, m_invNormalization;
};

MTS_NAMESPACE_END

#endif /* __MAXEXP_H */
