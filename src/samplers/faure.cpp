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

#include "faure.h"

MTS_NAMESPACE_BEGIN

PermutationStorage::PermutationStorage(int scramble) : m_rngIndex(0), m_scramble(scramble) {
    /* Compute the size of the final permutation table (corresponding to primes) */
    size_t finalSize = 0;
    for (size_t i=0; i<primeTableSize; ++i)
        finalSize += primeTable[i];

    /* Allocate memory for them */
    m_storage = new uint16_t[finalSize];
    m_permutations = new uint16_t*[primeTableSize];

    /* Check whether Faure or random permutations were requested */
    ref<Timer> timer = new Timer();
    if (scramble == -1) {
        /* Efficiently compute all Faure permutations using dynamic programming */
        uint16_t initialBases = primeTable[primeTableSize-1];
        size_t initialSize = ((size_t) initialBases * (size_t) (initialBases + 1))/2;
        uint16_t *initialStorage = new uint16_t[initialSize];
        uint16_t **initialPerm = new uint16_t*[initialBases+1],
                 *ptr = initialStorage;

        Log(EDebug, "Constructing Faure permutations using %s of memory",
            memString(initialSize * sizeof(uint16_t)).c_str());

        initialPerm[0] = NULL;
        for (size_t i=1; i<=initialBases; ++i) {
            initialPerm[i] = ptr;
            ptr += i;
        }
        computeFaurePermutations(initialBases, initialPerm);

        Log(EDebug, "Compactifying permutations to %s of memory",
            memString(finalSize * sizeof(uint16_t)).c_str());

        ptr = m_storage;
        for (size_t i=0; i<primeTableSize; ++i) {
            int prime = primeTable[i];
            memcpy(ptr, initialPerm[prime], prime * sizeof(uint16_t));
            m_permutations[i] = ptr;  ptr += prime;
        }

        delete[] initialStorage;
        delete[] initialPerm;
    } else {
        Log(EDebug, "Generating random permutations for the seed value = %i", scramble);

        uint16_t *ptr = m_storage;
        for (size_t i=0; i<primeTableSize; ++i) {
            int prime = primeTable[i];
            for (int j=0; j<prime; ++j)
                ptr[j] = (uint16_t) j;
            shuffle(ptr, ptr + prime);
            m_permutations[i] = ptr;  ptr += prime;
        }
    }
    Log(EDebug, "Done (took %i ms)", timer->getMilliseconds());

    /* Invert the first two permutations */
    m_invStorage = new uint16_t[5];
    m_invPermutations = new uint16_t*[2];
    m_invPermutations[0] = m_invStorage;
    m_invPermutations[1] = m_invStorage + 2;
    invertPermutation(0);
    invertPermutation(1);
}

PermutationStorage::~PermutationStorage()  {
    delete[] m_storage;
    delete[] m_invStorage;
    delete[] m_permutations;
    delete[] m_invPermutations;
}

/**
 * \ref Compute the Faure permutations using dynamic programming
 *
 * For reference, see "Good permutations for extreme discrepancy"
 * by Henri Faure, Journal of Number Theory, Vol. 42, 1, 1992.
 */
void PermutationStorage::computeFaurePermutations(uint32_t maxBase, uint16_t **perm) {
    SAssert(maxBase >= 2);

    /* Dimension 1 */
    perm[1][0] = 0;

    /* Dimension 2 */
    perm[2][0] = 0;
    perm[2][1] = 1;

    for (uint32_t b=2; b<=maxBase; ++b) {
        if (b & 1) {
            /* Odd dimension */
            uint16_t c = (b - 1) /2;

            for (uint16_t i=0; i<b; ++i) {
                if (i == c) {
                    perm[b][i] = c;
                } else {
                    uint16_t f = perm[b-1][i - (int) (i>c)];
                    perm[b][i] = f + (int) (f >= c);
                }
            }
        } else {
            /* Even dimension */
            uint16_t c = b / 2;

            for (uint16_t i=0; i<b; ++i)
                perm[b][i] = i < c ? 2*perm[c][i] : 2*perm[c][i-c]+1;;
        }
    }
}

void PermutationStorage::invertPermutation(uint32_t i) {
    uint16_t *perm = m_permutations[i],
                *invPerm = m_invPermutations[i];
    for (int j=0; j<primeTable[i]; ++j)
        invPerm[perm[j]] = j;
}

MTS_IMPLEMENT_CLASS(PermutationStorage, false, Object)
MTS_NAMESPACE_END
