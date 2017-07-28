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

#include <mitsuba/core/timer.h>
#include <mitsuba/core/qmc.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Stores scrambling permutations for Van Der Corput-type
 * sequences with prime bases.
 */
class PermutationStorage : public Object {
public:
    /**
     * \brief Create new permutations
     * \param scramble
     *    Selects the desired permutation type, where <tt>-1</tt> denotes the
     *    Faure permutations; any other number causes a pseudorandom permutation
     *    to be built seeded by the value of \c scramble.
     */
    PermutationStorage(int scramble);

    /// Return the permutation corresponding to the given prime number basis
    inline uint16_t *getPermutation(uint32_t basis) const {
        return m_permutations[basis];
    }

    /// Return the inverse permutation corresponding to the given prime number basis
    inline uint16_t *getInversePermutation(uint32_t basis) const {
        return m_invPermutations[basis];
    }

    /// Return the original scramble value
    inline int getScramble() const { return m_scramble; }

private:
    /**
     * \ref Compute the Faure permutations using dynamic programming
     *
     * For reference, see "Good permutations for extreme discrepancy"
     * by Henri Faure, Journal of Number Theory, Vol. 42, 1, 1992.
     */
    void computeFaurePermutations(uint32_t maxBase, uint16_t **perm);

    /// Randomly permute an array using the TEA pseudorandom generator
    inline void shuffle(uint16_t *it1, uint16_t *it2) {
        for (uint16_t * it = it2 - 1; it > it1; --it)
            std::iter_swap(it, it1 + sampleUInt((uint32_t) (it-it1)));
    }

    /// Draw a 32-bit integer in [0, n-1] from a TEA pseudorandom generator
    inline uint64_t sampleUInt(uint32_t n) {
        /* First, construct a bitmask */
        uint64_t bitmask = n;
        bitmask |= bitmask >> 1; bitmask |= bitmask >> 2;
        bitmask |= bitmask >> 4; bitmask |= bitmask >> 8;
        bitmask |= bitmask >> 16;
        uint32_t result;

        /* Next, generate numbers until one in [0, n) is found */
        while ((result = (uint32_t) (sampleTEA(m_scramble, m_rngIndex++) & bitmask)) >= n)
            ;

        return result;
    }

    /// Invert one of the permutations
    void invertPermutation(uint32_t i);

    MTS_DECLARE_CLASS()
protected:
    virtual ~PermutationStorage();

private:
    uint16_t *m_storage, *m_invStorage;
    uint16_t **m_permutations;
    uint16_t **m_invPermutations;
    uint32_t m_rngIndex;
    int m_scramble;
};

MTS_NAMESPACE_END
