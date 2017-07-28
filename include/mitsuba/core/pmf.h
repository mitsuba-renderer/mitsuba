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
#if !defined(__MITSUBA_CORE_PMF_H_)
#define __MITSUBA_CORE_PMF_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Discrete probability distribution
 *
 * This data structure can be used to transform uniformly distributed
 * samples to a stored discrete probability distribution.
 *
 * \ingroup libcore
 */
struct DiscreteDistribution {
public:
    /// Allocate memory for a distribution with the given number of entries
    explicit inline DiscreteDistribution(size_t nEntries = 0) {
        reserve(nEntries);
        clear();
    }

    /// Clear all entries
    inline void clear() {
        m_cdf.clear();
        m_cdf.push_back(0.0f);
        m_normalized = false;
    }

    /// Reserve memory for a certain number of entries
    inline void reserve(size_t nEntries) {
        m_cdf.reserve(nEntries+1);
    }

    /// Append an entry with the specified discrete probability
    inline void append(Float pdfValue) {
        m_cdf.push_back(m_cdf[m_cdf.size()-1] + pdfValue);
    }

    /// Return the number of entries so far
    inline size_t size() const {
        return m_cdf.size()-1;
    }

    /// Access an entry by its index
    inline Float operator[](size_t entry) const {
        return m_cdf[entry+1] - m_cdf[entry];
    }

    /// Have the probability densities been normalized?
    inline bool isNormalized() const {
        return m_normalized;
    }

    /**
     * \brief Return the original (unnormalized) sum of all PDF entries
     *
     * This assumes that \ref normalize() has previously been called
     */
    inline Float getSum() const {
        return m_sum;
    }

    /**
     * \brief Return the normalization factor (i.e. the inverse of \ref getSum())
     *
     * This assumes that \ref normalize() has previously been called
     */
    inline Float getNormalization() const {
        return m_normalization;
    }

    /**
     * \brief Normalize the distribution
     *
     * Throws an exception when no entries were previously
     * added to the distribution.
     *
     * \return Sum of the (previously unnormalized) entries
     */
    inline Float normalize() {
        SAssert(m_cdf.size() > 1);
        m_sum = m_cdf[m_cdf.size()-1];
        if (m_sum > 0) {
            m_normalization = 1.0f / m_sum;
            for (size_t i=1; i<m_cdf.size(); ++i)
                m_cdf[i] *= m_normalization;
            m_cdf[m_cdf.size()-1] = 1.0f;
            m_normalized = true;
        } else {
            m_normalization = 0.0f;
        }
        return m_sum;
    }

    /**
     * \brief %Transform a uniformly distributed sample to the stored distribution
     *
     * \param[in] sampleValue
     *     An uniformly distributed sample on [0,1]
     * \return
     *     The discrete index associated with the sample
     */
    inline size_t sample(Float sampleValue) const {
        std::vector<Float>::const_iterator entry =
                std::lower_bound(m_cdf.begin(), m_cdf.end(), sampleValue);
        size_t index = std::min(m_cdf.size()-2,
            (size_t) std::max((ptrdiff_t) 0, entry - m_cdf.begin() - 1));

        /* Handle a rare corner-case where a entry has probability 0
           but is sampled nonetheless */
        while (operator[](index) == 0 && index < m_cdf.size()-1)
            ++index;

        return index;
    }

    /**
     * \brief %Transform a uniformly distributed sample to the stored distribution
     *
     * \param[in] sampleValue
     *     An uniformly distributed sample on [0,1]
     * \param[out] pdf
     *     Probability value of the sample
     * \return
     *     The discrete index associated with the sample
     */
    inline size_t sample(Float sampleValue, Float &pdf) const {
        size_t index = sample(sampleValue);
        pdf = operator[](index);
        return index;
    }

    /**
     * \brief %Transform a uniformly distributed sample to the stored distribution
     *
     * The original sample is value adjusted so that it can be "reused".
     *
     * \param[in, out] sampleValue
     *     An uniformly distributed sample on [0,1]
     * \return
     *     The discrete index associated with the sample
     */
    inline size_t sampleReuse(Float &sampleValue) const {
        size_t index = sample(sampleValue);
        sampleValue = (sampleValue - m_cdf[index])
            / (m_cdf[index + 1] - m_cdf[index]);
        return index;
    }

    /**
     * \brief %Transform a uniformly distributed sample.
     *
     * The original sample is value adjusted so that it can be "reused".
     *
     * \param[in,out]
     *     An uniformly distributed sample on [0,1]
     * \param[out] pdf
     *     Probability value of the sample
     * \return
     *     The discrete index associated with the sample
     */
    inline size_t sampleReuse(Float &sampleValue, Float &pdf) const {
        size_t index = sample(sampleValue, pdf);
        sampleValue = (sampleValue - m_cdf[index])
            / (m_cdf[index + 1] - m_cdf[index]);
        return index;
    }

    /**
     * \brief Turn the underlying distribution into a
     * human-readable string format
     */
    std::string toString() const {
        std::ostringstream oss;
        oss << "DiscreteDistribution[sum=" << m_sum << ", normalized="
            << (int) m_normalized << ", cdf={";
        for (size_t i=0; i<m_cdf.size(); ++i) {
            oss << m_cdf[i];
            if (i != m_cdf.size()-1)
                oss << ", ";
        }
        oss << "}]";
        return oss.str();
    }
private:
    std::vector<Float> m_cdf;
    Float m_sum, m_normalization;
    bool m_normalized;
};

namespace math {
    /// Alias sampling data structure (see \ref makeAliasTable() for details)
    template <typename QuantizedScalar, typename Index> struct AliasTableEntry {
        /// Probability of sampling the current entry
        QuantizedScalar prob;
        /// Index of the alias entry
        Index index;
    };

    /**
     * \brief Create the lookup table needed for Walker's alias sampling
     * method implemented in \ref sampleAlias(). Runs in linear time.
     *
     * The basic idea of this method is that one can "redistribute" the
     * probability mass of a distribution to make it uniform. This
     * this can be done in a way such that the probability of each entry in
     * the "flattened" PMF consists of probability mass from at most *two*
     * entries in the original PMF. That then leads to an efficient O(1)
     * sampling algorithm with a O(n) preprocessing step to set up this
     * special decomposition.
     *
     * The downside of this method is that it generally does not preserve
     * the nice stratification properties of QMC number sequences.
     *
     * \return The original (un-normalized) sum of all probabilities
     * in \c pmf.
     */
    template <typename Scalar, typename QuantizedScalar, typename Index> float makeAliasTable(
            AliasTableEntry<QuantizedScalar, Index> *tbl, Scalar *pmf, Index size) {
        /* Allocate temporary storage for classification purposes */
        Index *c = new Index[size],
              *c_short = c - 1, *c_long  = c + size;

        /* Begin by computing the normalization constant */
        Scalar sum = 0;
        for (size_t i=0; i<size; ++i)
            sum += pmf[i];

        Scalar normalization = (Scalar) 1 / sum;
        for (Index i=0; i<size; ++i) {
            /* For each entry, determine whether there is
               "too little" or "too much" probability mass */
            Scalar value = size * normalization * pmf[i];
            if (value < 1)
                *++c_short = i;
            else if (value > 1)
                *--c_long  = i;
            tbl[i].prob  = value;
            tbl[i].index = i;
        }

        /* Perform pairwise exchanges while there are entries
           with too much probability mass */
        for (Index i=0; i < size-1 && c_long - c < size; ++i) {
            Index short_index = c[i],
                  long_index  = *c_long;

            tbl[short_index].index = long_index;
            tbl[long_index].prob  -= (Scalar) 1 - tbl[short_index].prob;

            if (tbl[long_index].prob <= 1)
                ++c_long;
        }

        delete[] c;

        return sum;
    }

    /// Generate a sample in constant time using the alias method
    template <typename Scalar, typename QuantizedScalar, typename Index> Index sampleAlias(
            const AliasTableEntry<QuantizedScalar, Index> *tbl, Index size, Scalar sample) {
        Index l = std::min((Index) (sample * size), (Index) (size - 1));
        Scalar prob = (Scalar) tbl[l].prob;

        sample = sample * size - l;

        if (prob == 1 || (prob != 0 && sample < prob))
            return l;
        else
            return tbl[l].index;
    }

    /**
     * \brief Generate a sample in constant time using the alias method
     *
     * This variation shifts and scales the uniform random sample so
     * that it can be reused for another sampling operation
     */
    template <typename Scalar, typename QuantizedScalar, typename Index> Index sampleAliasReuse(
            const AliasTableEntry<QuantizedScalar, Index> *tbl, Index size, Scalar &sample) {
        Index l = std::min((Index) (sample * size), (Index) (size - 1));
        Scalar prob = (Scalar) tbl[l].prob;

        sample = sample * size - l;

        if (prob == 1 || (prob != 0 && sample < prob)) {
            sample /= prob;
            return l;
        } else {
            sample = (sample - prob) / (1 - prob);
            return tbl[l].index;
        }
    }
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_PMF_H_ */
