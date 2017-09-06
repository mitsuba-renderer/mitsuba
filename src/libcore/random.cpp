/*
   SIMD oriented Fast Mersenne Twister (SFMT) pseudorandom number generator
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/

   Copyright (c) 2006,2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
   University. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

       * Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
       * Redistributions in binary form must reproduce the above
         copyright notice, this list of conditions and the following
         disclaimer in the documentation and/or other materials provided
         with the distribution.
       * Neither the name of the Hiroshima University nor the names of
         its contributors may be used to endorse or promote products
         derived from this software without specific prior written
         permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   M. Saito and M. Matsumoto,
     ``SIMD-oriented Fast Mersenne Twister:
       a 128-bit Pseudorandom Number Generator''
     Monte Carlo and Quasi-Monte Carlo Method 2006.
     Springer (2008) 607--622.
     DOI: 10.1007/978-3-540-74496-2_36
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and
     Computer Simulation 8. (Jan. 1998) 3--30.
*/

#include <mitsuba/core/random.h>
#include <mitsuba/core/fstream.h>

#if MTS_SSE
# define MTS_SFMT_SSE 1
#else
# define MTS_SFMT_SSE 0
#endif

#if MTS_SFMT_SSE
# include <mitsuba/core/sse.h>
#endif

#include <limits>


/************************ SFMT-19937 Parameters *******************************/

/* Mersenne Exponent. The period of the sequence
 * is a multiple of 2^MEXP-1. */
#define MEXP 19937
/* SFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define N (MEXP / 128 + 1)
/* N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define N32 (N * 4)
/* N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define N64 (N * 2)

/* MEXP dependent values */
#define POS1    122
#define SL1 18
#define SL2 1
#define SR1 11
#define SR2 1
#define MSK1    0xdfffffefU
#define MSK2    0xddfecb7fU
#define MSK3    0xbffaffffU
#define MSK4    0xbffffff6U
#define PARITY1 0x00000001U
#define PARITY2 0x00000000U
#define PARITY3 0x00000000U
#define PARITY4 0x13c9e684U

/******************************************************************************/


MTS_NAMESPACE_BEGIN

// Helper elements
namespace {

/** Compute a hash value representing the SFMT parameters */
inline uint32_t sfmtHash() {
    // Based on boost::hash_combine
    const uint32_t data[] =
        {MEXP, POS1, SL1, SL2, SR1, SR2, MSK1, MSK2, MSK3, MSK4};
    uint32_t hash = 17;
    for (int i = 0; i < 10; ++i)
        hash ^= data[i] + 0x9e3779b9U + (hash<<6) + (hash>>2);
    return hash;
}

/** 128-bit data structure */
union w128_t {
#if MTS_SFMT_SSE
    __m128i si;
#else
    uint64_t u64[2];
#endif
    uint32_t u[4];
};



#if !MTS_SFMT_SSE

/**
 * This function simulates SIMD 128-bit right shift by the standard C.
 * The 128-bit integer given in in is shifted by (shift * 8) bits.
 * This function simulates the LITTLE ENDIAN SIMD.
 * \param out the output of this function
 * \param in the 128-bit data to be shifted
 * \param shift the shift value
 */
inline static void rshift128(w128_t &out, const w128_t &in, const int shift) {
    uint64_t oh, ol;

    const uint64_t &th = in.u64[1];
    const uint64_t &tl = in.u64[0];

    oh = th >> (shift * 8);
    ol = tl >> (shift * 8);
    ol |= th << (64 - shift * 8);
    out.u64[0] = ol;
    out.u64[1] = oh;
}

/**
 * This function simulates SIMD 128-bit left shift by the standard C.
 * The 128-bit integer given in in is shifted by (shift * 8) bits.
 * This function simulates the LITTLE ENDIAN SIMD.
 * \param out the output of this function
 * \param in the 128-bit data to be shifted
 * \param shift the shift value
 */
inline static void lshift128(w128_t &out, const w128_t &in, const int shift) {
    uint64_t oh, ol;

    const uint64_t &th = in.u64[1];
    const uint64_t &tl = in.u64[0];

    oh = th << (shift * 8);
    ol = tl << (shift * 8);
    oh |= tl >> (64 - shift * 8);
    out.u64[0] = ol;
    out.u64[1] = oh;
}

#endif // !MTS_SFMT_SSE



/**
 * This function represents the recursion formula.
 * \param a a 128-bit part of the interal state array
 * \param b a 128-bit part of the interal state array
 * \param c a 128-bit part of the interal state array
 * \param d a 128-bit part of the interal state array
 * \param mask 128-bit mask
 * \return output
 */
#if MTS_SFMT_SSE
FINLINE __m128i mm_recursion(const  __m128i &a, const __m128i &b,
    const __m128i &c, const __m128i &d, const __m128i &mask) {
    __m128i v, x, y, z;

    x = _mm_load_si128(&a);
    y = _mm_srli_epi32(b, SR1);
    z = _mm_srli_si128(c, SR2);
    v = _mm_slli_epi32(d, SL1);
    z = _mm_xor_si128(z, x);
    z = _mm_xor_si128(z, v);
    x = _mm_slli_si128(x, SL2);
    y = _mm_and_si128(y, mask);
    z = _mm_xor_si128(z, x);
    z = _mm_xor_si128(z, y);
    return z;
}
#else // MTS_SFMT_SSE
inline static void do_recursion(w128_t &r, const w128_t &a, const w128_t &b,
    const w128_t &c, const w128_t &d) {
    w128_t x;
    w128_t y;

    lshift128(x, a, SL2);
    rshift128(y, c, SR2);
    r.u[0] = a.u[0] ^ x.u[0] ^ ((b.u[0] >> SR1) & MSK1) ^ y.u[0]
    ^ (d.u[0] << SL1);
    r.u[1] = a.u[1] ^ x.u[1] ^ ((b.u[1] >> SR1) & MSK2) ^ y.u[1]
    ^ (d.u[1] << SL1);
    r.u[2] = a.u[2] ^ x.u[2] ^ ((b.u[2] >> SR1) & MSK3) ^ y.u[2]
    ^ (d.u[2] << SL1);
    r.u[3] = a.u[3] ^ x.u[3] ^ ((b.u[3] >> SR1) & MSK4) ^ y.u[3]
    ^ (d.u[3] << SL1);
}
#endif // MTS_SFMT_SSE

} // namespace

// Actual state structure definition
struct Random::State {
    union {
        /** the 128-bit internal state array */
        w128_t sfmt[N];
        /** the 32bit integer version of the 128-bit internal state array */
        uint32_t psfmt32[N32];
        /** the 64bit integer version of the 128-bit internal state array */
        uint64_t psfmt64[N64];
    };

    /** index counter to the 32-bit internal state array */
    int idx;

    /** a parity check vector which certificate the period of 2^{MEXP} */
    const static uint32_t parity[4];

    /** Hash of the SFMT parameters */
    const static uint32_t s_magic;


    /* Default constructor, set the index to an invalid value */
    State() : idx(-1) {}

    /* Creates a new instance from a stream */
    void unserialize(Stream * stream) {
        stream->readULongArray(&psfmt64[0], N64);
        idx = stream->readInt();
    }

    /* Writes the current instance to a stream */
    void serialize(Stream *stream) const {
        stream->writeULongArray(&psfmt64[0], N64);
        stream->writeInt(idx);
    }

    /// Helper function to check whether the state is initialized
    inline bool isInitialized() const {
        return idx >= 0;
    }

    /**
      * \brief This function initializes the internal state array with a 64-bit
      * integer seed.
      *
      * \param seed a 64-bit integer used as the seed.
      */
    void init_gen_rand(uint64_t seed);

    /**
     * \brief This function initializes the internal state array,
     * with an array of 32-bit integers used as the seeds
     * \param init_key the array of 32-bit integers, used as a seed.
     * \param key_length the length of init_key.
     */
    void init_by_array(const uint32_t *init_key, int key_length);

    /**
    * This function generates and returns 64-bit pseudorandom number.
    * init_gen_rand or init_by_array must be called before this function.
    * The function gen_rand64 should not be called after gen_rand32,
    * unless an initialization is again executed.
    * \return 64-bit pseudorandom number
    */
    FINLINE uint64_t gen_rand64() {
        if (idx >= N32) {
            gen_rand_all();
            idx = 0;
        }

        uint64_t r = psfmt64[idx / 2];
        idx += 2;
        return r;
    }

private:

    /**
    * This function represents a function used in the initialization
    * by init_by_array
    * \param x 32-bit integer
    * \return 32-bit integer
    */
    static inline uint32_t func1(uint32_t x) {
        return (x ^ (x >> 27)) * (uint32_t)1664525UL;
    }

    /**
    * This function represents a function used in the initialization
    * by init_by_array
    * \param x 32-bit integer
    * \return 32-bit integer
    */
    static inline uint32_t func2(uint32_t x) {
        return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
    }

    /* This function certificate the period of 2^{MEXP} */
    void period_certification() {
        int inner = 0;
        int i, j;
        uint32_t work;

        for (i = 0; i < 4; ++i)
            inner ^= psfmt32[i] & parity[i];
        for (i = 16; i > 0; i >>= 1)
            inner ^= inner >> i;
        inner &= 1;
        /* check OK */
        if (inner == 1) {
            return;
        }
        /* check NG, and modification */
        for (i = 0; i < 4; ++i) {
            work = 1;
            for (j = 0; j < 32; ++j) {
                if ((work & parity[i]) != 0) {
                    psfmt32[i] ^= work;
                    return;
                }
                work = work << 1;
            }
        }
    }

    /*
     * This function fills the internal state array with pseudorandom
     * integers.
     */
    inline void gen_rand_all() {
#if MTS_SFMT_SSE
        int i;
        __m128i r, r1, r2, mask;
        mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);

        r1 = _mm_load_si128(&sfmt[N - 2].si);
        r2 = _mm_load_si128(&sfmt[N - 1].si);
        for (i = 0; i < N - POS1; ++i) {
            r = mm_recursion(sfmt[i].si, sfmt[i + POS1].si, r1, r2, mask);
            _mm_store_si128(&sfmt[i].si, r);
            r1 = r2;
            r2 = r;
        }
        for (; i < N; ++i) {
            r = mm_recursion(sfmt[i].si, sfmt[i + POS1 - N].si, r1, r2, mask);
            _mm_store_si128(&sfmt[i].si, r);
            r1 = r2;
            r2 = r;
        }
#else // MTS_SFMT_SSE
        int i;
        w128_t *r1, *r2;

        r1 = &sfmt[N - 2];
        r2 = &sfmt[N - 1];
        for (i = 0; i < N - POS1; ++i) {
            do_recursion(sfmt[i], sfmt[i], sfmt[i + POS1], *r1, *r2);
            r1 = r2;
            r2 = &sfmt[i];
        }
        for (; i < N; ++i) {
            do_recursion(sfmt[i], sfmt[i], sfmt[i + POS1 - N], *r1, *r2);
            r1 = r2;
            r2 = &sfmt[i];
        }
#endif
    } // MTS_SFMT_SSE
};

const uint32_t Random::State::parity[4] = {PARITY1, PARITY2, PARITY3, PARITY4};
const uint32_t Random::State::s_magic = sfmtHash();


void Random::State::init_gen_rand(const uint64_t seed) {
    psfmt64[0] = seed;
    for (int i = 1; i < N64; ++i) {
        psfmt64[i] = (6364136223846793005ULL * (psfmt64[i-1]
                      ^ (psfmt64[i-1] >> 62)) + i);
    }
    idx = N32;
    period_certification();
    SAssert(isInitialized());
}

void Random::State::init_by_array(const uint32_t *init_key, int key_length) {
    int i, j, count;
    uint32_t r;
    int lag;
    int mid;
    int size = N * 4;

    if (size >= 623) {
        lag = 11;
    } else if (size >= 68) {
        lag = 7;
    } else if (size >= 39) {
        lag = 5;
    } else {
        lag = 3;
    }
    mid = (size - lag) / 2;

    memset(sfmt, 0x8b, sizeof(sfmt));
    if (key_length + 1 > N32) {
        count = key_length + 1;
    } else {
        count = N32;
    }
    r = func1(psfmt32[0] ^ psfmt32[mid]
          ^ psfmt32[N32 - 1]);
    psfmt32[mid] += r;
    r += key_length;
    psfmt32[mid + lag] += r;
    psfmt32[0] = r;

    count--;
    for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
        r = func1(psfmt32[i] ^ psfmt32[(i + mid) % N32]
              ^ psfmt32[(i + N32 - 1) % N32]);
        psfmt32[(i + mid) % N32] += r;
        r += init_key[j] + i;
        psfmt32[(i + mid + lag) % N32] += r;
        psfmt32[i] = r;
        i = (i + 1) % N32;
    }
    for (; j < count; j++) {
        r = func1(psfmt32[i] ^ psfmt32[(i + mid) % N32]
              ^ psfmt32[(i + N32 - 1) % N32]);
        psfmt32[(i + mid) % N32] += r;
        r += i;
        psfmt32[(i + mid + lag) % N32] += r;
        psfmt32[i] = r;
        i = (i + 1) % N32;
    }
    for (j = 0; j < N32; j++) {
        r = func2(psfmt32[i] + psfmt32[(i + mid) % N32]
              + psfmt32[(i + N32 - 1) % N32]);
        psfmt32[(i + mid) % N32] ^= r;
        r -= i;
        psfmt32[(i + mid + lag) % N32] ^= r;
        psfmt32[i] = r;
        i = (i + 1) % N32;
    }

    idx = N32;
    period_certification();
    SAssert(isInitialized());
}

Random::Random() : mt(NULL) {
    mt = (State *) allocAligned(sizeof(State));
    Assert(mt != NULL);
#if defined(__WINDOWS__)
    seed();
#else
#if 0
    uint64_t buf[N64];
    memset(buf, 0, N64 * sizeof(uint64_t)); /* Make GCC happy */
    ref<FileStream> urandom = new FileStream("/dev/urandom", FileStream::EReadOnly);
    urandom->readULongArray(buf, N64);
    seed(buf, N64);
#else
    seed();
#endif
#endif
}

Random::Random(Random *random) : mt(NULL) {
    mt = (State *) allocAligned(sizeof(State));
    Assert(mt != NULL);
    seed(random);
}

Random::Random(uint64_t seedval) : mt(NULL) {
    mt = (State *) allocAligned(sizeof(State));
    Assert(mt != NULL);
    seed(seedval);
}

Random::~Random() {
    if (mt != NULL)
        freeAligned(mt);
}

Random::Random(Stream *stream, InstanceManager *manager)
        : SerializableObject(stream, manager), mt(NULL) {
    // Check the magic number first
    uint32_t magic = stream->readUInt();
    if (magic != State::s_magic)
        Log(EError, "Incorrected SFMT magic number: expected %08x, actual %08x",
            State::s_magic, magic);
    mt = (State *) allocAligned(sizeof(State));
    mt->unserialize(stream);
}

void Random::serialize(Stream *stream, InstanceManager *manager) const {
    stream->writeInt(State::s_magic);
    mt->serialize(stream);
}

void Random::seed(uint64_t s) {
    mt->init_gen_rand(s);
}

void Random::seed(Random *random) {
    uint64_t buf[N64];
    for (int i=0; i<N64; ++i)
        buf[i] = random->nextULong();
    seed(buf, N64);
}

void Random::set(Random *random) {
    Assert(random != NULL && random->mt != NULL && mt != NULL);
    *mt = *(random->mt);
}

void Random::seed(uint64_t *init_key, uint64_t key_length) {
    const uint64_t max_key_length = std::numeric_limits<int>::max() / 2;
    if (key_length > max_key_length) {
        Log(EWarn, "Excessive SFMT initialization data, igoring extra values.");
        key_length = max_key_length;
    }
    const uint32_t * key = reinterpret_cast<uint32_t*>(init_key);
    mt->init_by_array(key, static_cast<int>(key_length) * 2);
}


uint64_t Random::nextULong() {
    return mt->gen_rand64();
}


namespace {
    // Helper function to create bitmasks according to the size
    template <typename T> T makeBitmask(T n);

    template <> inline uint32_t makeBitmask(uint32_t n) {
        uint32_t bitmask = n;
        bitmask |= bitmask >> 1;
        bitmask |= bitmask >> 2;
        bitmask |= bitmask >> 4;
        bitmask |= bitmask >> 8;
        bitmask |= bitmask >> 16;
        return bitmask;
    }

    template <> inline uint64_t makeBitmask(uint64_t n) {
        uint64_t bitmask = n;
        bitmask |= bitmask >> 1;
        bitmask |= bitmask >> 2;
        bitmask |= bitmask >> 4;
        bitmask |= bitmask >> 8;
        bitmask |= bitmask >> 16;
        bitmask |= bitmask >> 32;
        return bitmask;
    }

    #if defined(MTS_AMBIGUOUS_SIZE_T)
    inline size_t makeBitmask(size_t n) {
        if (sizeof(size_t) == 8)
            return (size_t) makeBitmask((uint64_t) n);
        else
            return (size_t) makeBitmask((uint32_t) n);
    }
    #endif
} // namespace


uint32_t Random::nextUInt(uint32_t n) {
    /* Determine a bit mask */
    const uint32_t bitmask = makeBitmask(n);
    uint32_t result;

    /* Generate numbers until one in [0, n) is found */
    while ((result = (uint32_t) (nextULong() & bitmask)) >= n)
        ;

    return result;
}

size_t Random::nextSize(size_t n) {
    /* Determine a bit mask */
    const size_t bitmask = makeBitmask(n);
    size_t result;

    /* Generate numbers until one in [0, n) is found */
    while ((result = (size_t) (nextULong() & bitmask)) >= n)
        ;

    return result;
}

#if defined(DOUBLE_PRECISION)
Float Random::nextFloat() {
    /* Trick from MTGP: generate an uniformly distributed
       single precision number in [1,2) and subtract 1. */
    union {
        uint64_t u;
        double d;
    } x;
    x.u = (nextULong() >> 12) | 0x3ff0000000000000ULL;
    return x.d - 1.0;
}

#else

Float Random::nextFloat() {
    /* Trick from MTGP: generate an uniformly distributed
       single precision number in [1,2) and subtract 1. */
    union {
        uint32_t u;
        float f;
    } x;
    x.u = ((nextULong() & 0xFFFFFFFF) >> 9) | 0x3f800000UL;
    return x.f - 1.0f;
}
#endif

Float Random::nextStandardNormal() {
    /* Marsaglia polar method for generating two standard
       normal variates. One is subsequently thrown away */
    Float x, y, r;
    do {
        x = 2.0f * nextFloat() - 1.0f;
        y = 2.0f * nextFloat() - 1.0f;
        r = x*x + y*y;
    } while (r >= 1 || r == 0);

    return x * std::sqrt(-2.0f * math::fastlog(r) / r);
}

MTS_IMPLEMENT_CLASS_S(Random, false, SerializableObject)
MTS_NAMESPACE_END
