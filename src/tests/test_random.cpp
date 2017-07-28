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

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/mstream.h>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/common_factor.hpp>

#include <limits>
#include <vector>
#include <functional>
#include <algorithm>

MTS_NAMESPACE_BEGIN

namespace {

// Small class which represents a binary 32x32 matrix, used for the rank test.
// The matrix stores each row as an uint32_t.
class BinaryMatrix32x32
{
public:

    // Initializes the matrix with random data. This assumes the random
    // number generator generates good quality 64-bit unsigned integers
    BinaryMatrix32x32(ref<Random> & rnd) {
        for (int i = 0; i < 16; ++i) {
            m_data64[i] = rnd->nextULong();
        }
    }

    // Initializes the matrix using formatted input via "operator>>", reading
    // 32 uint32_t values consecutively
    BinaryMatrix32x32(std::istream &is) {
        for (int i = 0; i < 32; ++i) {
            is >> m_rows[i];
        }
    }

    // Copies the rows from the given data
    BinaryMatrix32x32(const uint32_t(&rows)[32]) {
        std::copy(rows, rows+32, m_rows);
    }

    // Returns a reference to the given row WITHOUT bounds checking
    uint32_t & row(size_t idx) {
        return m_rows[idx];
    }

    // Returns a reference to the given row WITHOUT bounds checking (const)
    const uint32_t & row(size_t idx) const {
        return m_rows[idx];
    }

    // Returns the rank modulus-2 of this matrix. The method generates the
    // (non-reduced) row echeleton form of the matrix using modular arithmetic
    int rank_mod2() const;

private:

    // Functor which returns true if applying "and" with the given mask
    // to the value is non-zero.
    template <typename T>
    struct nonzero_and : public std::unary_function<T, bool>
    {
        const T mask;

        // Constructs a functor with a given mask
        nonzero_and(T m) : mask(m) {}

        bool operator() (const T & value) const {
            return (value & mask) != static_cast<T>(0);
        }
    };


    // Functor which returns true is the element is different than zero
    template <typename T>
    struct is_nonzero : public std::unary_function<T, bool>
    {
        bool operator() (const T & v) const {
            return v != 0;
        }
    };


    // Actual data, providing also 64-bit access for initialization.
    union {
        uint32_t m_rows[32];
        uint64_t m_data64[16];
    };
};


int BinaryMatrix32x32::rank_mod2() const
{
    // Scratch space which will hold the row echeleton form
    std::vector<uint32_t> m(m_rows, m_rows + 32);
    typedef std::vector<uint32_t>::iterator row_iter;

    // Start assuming full rank
    int rank = 32;

    for (int col = 0, pivot_row_idx = 0; col < 32; ++col) {
        const uint32_t mask = (0x80000000U >> col);
#ifndef NDEBUG
        const uint32_t mask_valid = ~((0xffffffffU >> col) >> 1);
#endif
        // Find the index of the element to use for the current pivot
        row_iter candidate_row = std::find_if(m.begin()+pivot_row_idx, m.end(),
            nonzero_and<uint32_t>(mask));
        if (candidate_row == m.end()) {
            // This means that the element on the corresponding diagonal is zero
            // indicating a rank deficiency
            --rank;
            continue;
        }

        // Swap the pivot row to the correct position
        const uint32_t pivot_row = *candidate_row;
        std::swap(m[pivot_row_idx], *candidate_row);

        // Make the leading column in the other rows zero. Because
        // everything is mod2, it is enough to xor (addition mod2) the pivot row
        // to those which do not already have a leading zero
        for (row_iter row = m.begin()+pivot_row_idx+1; row != m.end(); ++row) {
            if ((*row & mask) != 0) {
                *row ^= pivot_row;
                assert((*row & mask) == 0);
                assert((*row & mask_valid) == 0);
            }
        }

        ++pivot_row_idx;
    }

    // The rank is the number of non-zero rows as well
    assert(rank == std::count_if(m.begin(), m.end(), is_nonzero<uint32_t>()));
    return rank;
}



// Simple array traits to get the size of a fixed array
template <typename T, size_t N>
inline size_t array_size(const T (&arr)[N]) {
    return N;
}



// Traits type so that a random number generator returns the specified type
// (e.g.) float, uint32_t, uint64_t.
template <typename T>
struct rnd_traits;

template <>
struct rnd_traits<Float> {
    static Float rnd(ref<Random> & rnd) {
        return rnd->nextFloat();
    }

    const static bool fixedRange = true;
    static Float min() { return static_cast<Float>(0); }
    static Float max() { return static_cast<Float>(1); }
};

template <>
struct rnd_traits<uint64_t> {
    static uint64_t rnd(ref<Random> & rnd) {
        return rnd->nextULong();
    }

    const static bool fixedRange = sizeof(uint64_t) != sizeof(size_t);
    static uint64_t min() { return 0L; }
    static uint64_t max() { return 0xFFFFFFFFFFFFFFFFULL; }
};

template <>
struct rnd_traits<uint32_t> {
    static uint32_t rnd(ref<Random> & rnd) {
        return rnd->nextUInt(max());
    }

    const static bool fixedRange = false;
    static uint32_t min() { return 0L; }
    static uint32_t max() { return std::numeric_limits<uint32_t>::max(); }
};



// Fake iterator which never moves and only returns contant references
// to its value.
template <typename T>
class fake_iterator : public std::iterator<std::input_iterator_tag, T>
{
public:
    const T value;

    // Constructs the fake iterator with a given value
    fake_iterator(const T & v) : value(v) {}

    // Constant reference to the value
    const T& operator*() const {
        return value;
    }

    // No-op!
    fake_iterator& operator++() {
        return *this;
    }
};



// Kolmogorov-Smirnov probability function.
// Taken from Numerical Recipes in C, 14.3
double probks(const double alam)
{
    const double EPS1=1.0e-6,EPS2=1.0e-16;
    double a2,fac=2.0,sum=0.0,term,termbf=0.0;

    a2 = -2.0*alam*alam;
    for (int j=1;j<=100;j++) {
        term=fac*exp(a2*j*j);
        sum += term;
        if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
        fac = -fac;
        termbf=fabs(term);
    }
    return 1.0;
}



// Simple functor for a uniform distribution
struct UniformCDF
{
    const boost::math::uniform dist;

    UniformCDF(double lower = 0.0, double upper = 1.0) : dist(lower, upper) {}

    template <typename RealType>
    double operator() (RealType x) const {
        return boost::math::cdf(dist, static_cast<double>(x));
    }
};



// Templated function to calculate the p-value of the chi-squared statistic
// where degrees-of-freedom = number-of-items - restrictions
template <typename ResultIter, typename ReferenceIter>
double chi_squared_pval(
    const ResultIter result_begin, const ResultIter result_end,
    const ReferenceIter reference_begin, const size_t restrictions = 1)
{
    double chi_squared = 0.0;
    size_t df = 0;
    ReferenceIter it_ref = reference_begin;
    for (ResultIter it = result_begin; it != result_end; ++it, ++it_ref)
    {
        const double expected = static_cast<double>(*it_ref);
        const double delta = static_cast<double>(*it) - expected;
        chi_squared += (delta * delta) / expected;
        ++df;
    }

    assert(restrictions < df);
    df -= restrictions;

    // Find the critical value for the given degrees of freedom
    boost::math::chi_squared chSqDist(static_cast<double>(df));
    const double pval = 1.0 - boost::math::cdf(chSqDist, chi_squared);

    SLog(EDebug, "Chi-square statistic = %g (df=%i), P-value: %.6f%%",
        chi_squared, df, pval * 100.0);
    return pval;
}



// Get the p-val of a Kolmogorov-Smirnov test. Note that the input data WILL
// be modified by sorting it! Code based on the one provided in
// Numerical Recipes in C, 14.3
template <typename RandomIter, typename CDF_Functor>
double ks_pval(RandomIter itFirst, RandomIter itLast, const CDF_Functor & cdf)
{
    const size_t n = itLast - itFirst;
    std::sort(itFirst, itLast);

    double nInv = 1.0 / static_cast<double>(n);
    double fo = 0.0;
    double d  = 0.0;
    double fn = nInv;
    for (RandomIter it = itFirst; it != itLast; ++it, fn += nInv) {
        const double ff = cdf(*it);
        const double dt = std::max(std::abs(fo - ff), std::abs(fn - ff));
        d = std::max(d, dt);
        fo = fn;
    }

    const double nSqrt = std::sqrt(static_cast<double>(n));
    const double query = (nSqrt + 0.12 + 0.11/nSqrt) * d;
    const double pVal  = probks(query);
    return pVal;
}

} // namespace



class TestRandom : public TestCase {
public:
    MTS_BEGIN_TESTCASE()
    MTS_DECLARE_TEST(test00_validate)
    MTS_DECLARE_TEST(test01_mean)
    MTS_DECLARE_TEST(test02_range)
    MTS_DECLARE_TEST(test03_circle)
    MTS_DECLARE_TEST(test04_uniform_distribution)
    MTS_DECLARE_TEST(test05_rank_validation)
    MTS_DECLARE_TEST(test05_rank)
    MTS_DECLARE_TEST(test06_relative_primes)
    MTS_DECLARE_TEST(test07_uniform_distribution_ks);
    MTS_DECLARE_TEST(test08_serialize);
    MTS_DECLARE_TEST(test09_set);
    MTS_DECLARE_TEST(benchmark);
    MTS_END_TESTCASE()

    void test00_validate();
    void test01_mean();
    void test02_range();
    void test03_circle();
    void test04_uniform_distribution();
    void test05_rank_validation();
    void test05_rank();
    void test06_relative_primes();
    void test07_uniform_distribution_ks();
    void test08_serialize();
    void test09_set();
    void benchmark();

private:

    // Sidak correction to apply to alpha when running multiple runs of a test
    template <typename T_Alpha, typename TN>
    inline static double sidak_correction(const T_Alpha alpha, TN numRuns) {
        return 1.0 - std::pow(1.0 - alpha, 1.0 / numRuns);
    }

    // Small class used to run a chi-squared tests on the uniformity of p-vals
    // using 10 bins
    class UnifromTestPVals {
    public:
        UnifromTestPVals() : m_bins(10, 0), m_count(0) {}

        inline void add(const double pVal) {
            assert(0.0 <= pVal && pVal <= 1.0);
            ++m_bins[std::min(static_cast<size_t>(pVal * 10.0),
                static_cast<size_t>(10))];
            ++m_count;
        }

        inline const size_t & count() const {
            return m_count;
        }

        inline double pval() const {
            SLog(EDebug, "P-Vals histogram: [%d %d %d %d %d %d %d %d %d %d]",
                m_bins[0], m_bins[1], m_bins[2], m_bins[3], m_bins[4],
                m_bins[5], m_bins[6], m_bins[7], m_bins[8], m_bins[9]);
            double pval = chi_squared_pval(m_bins.begin(), m_bins.end(),
                fake_iterator<double>(count()/10.0));
            return pval;
        }

        inline void reset() {
            m_count = 0;
            for (iterator it = m_bins.begin(); it != m_bins.end(); ++it) {
                *it = 0;
            }
        }

    private:
        typedef std::vector<size_t>::iterator iterator;

        std::vector<size_t> m_bins;
        size_t m_count;
    };

    // Helper template to generate random numbers based on selected types
    template <typename T>
    void test02_helper(ref<Random> & rnd);

    // Chi-squared test for uniform distributions using bins. The defaults
    // follow the recomendations in NIST SP800-22b [May 2001]
    double uniform_distribution_instance(ref<Random> & rnd, const double alpha,
        const int numBins = 10, const int numPerBin = 20);

    // Kolmogorov-Smirnov test for continuous uniform distribution.
    double uniform_distribution_ks_instance(ref<Random> & rnd, const double alpha,
        const int numSamples = 10000);

    // Chi-squared test for binary matrix rank in GF_2 (integers modulus 2)
    double binary_rank_instance(ref<Random> & rnd, const double alpha,
        const int numMatrices = 100);
};



// Basic validation of the output based on reference data. Note that this test
// depends on the implementation details of the actual random generator.
void TestRandom::test00_validate()
{
    // Reference data for the SFMT19937 generator
    static const uint64_t reference[] = {
        0xa0920029ffafd7fcULL, 0x8cb79aaa80983938ULL, 0x2931508e0eb382d3ULL,
        0xf25a7ef6e4aaaa0eULL, 0xb92c989aba00eb9fULL, 0x67708ba73386fcb1ULL,
        0xe181a1fba15a30aeULL, 0x91214d71bd952085ULL, 0xa5a1300b7ae38d06ULL,
        0xa371e58646ab9244ULL, 0xd6eea5a18cde9b23ULL, 0xb89b5cd7caef14b5ULL,
        0xdcc8deee44b378c1ULL, 0x60ff242cba453e77ULL, 0x0737a743791a92a6ULL,
        0xef6a5098247294ecULL, 0xdc2ace07791777b6ULL, 0x68b7f67f70454e6eULL,
        0xb2b0ca8770402df1ULL, 0xae02036a5d0f1c7cULL, 0x53b2ead3d002ab34ULL,
        0x2d27e108016b2109ULL, 0xfac57edfe74fa265ULL, 0xab887a5cc0d658f8ULL,
        0x2da51426f482a604ULL, 0x95ac0f5cd9c897d8ULL, 0x9a1138cad7bed05cULL,
        0xf1c63d4183a1fc15ULL, 0x7454c2074cb0f33eULL, 0xdfe6ba0cb454c658ULL,
        0xc648c98c84a02046ULL, 0x18bf0d60845c9a0aULL, 0xf3b195b59d58eb0bULL,
        0xe5913e1316ef168fULL, 0x0074ee40066a8a37ULL, 0x5a4f67458e1503d1ULL,
        0xb2c160030df27955ULL, 0x89b5824b4c622d05ULL, 0x465c079ce9d5004cULL,
        0xa35d5fc2f800c9e1ULL, 0x545bfcb0a792bff4ULL, 0xf376ba15bd381d9eULL,
        0xa5c5c8c3245296feULL, 0xef5c7c58adbe0f5aULL, 0x237cd0e052a4f3b5ULL,
        0x46791e4bafdbe873ULL, 0xeb553d5551f2e623ULL, 0x64888fc94a6273ebULL,
        0x18d1401f9122dea5ULL, 0x9a86acccff9c6886ULL, 0x341d8d19d87e7e1dULL,
        0xc392a6c561710430ULL, 0xf3c33f434622725bULL, 0x4d4eeaa9bbe61b29ULL,
        0x8a92debb09f9f647ULL, 0xff793ab7077a6db2ULL, 0xcad376f93e5122aeULL,
        0xfc4b3996eee71482ULL, 0x370744ff51033571ULL, 0x2d25f4cd4425f482ULL,
        0xfaa0fcff515731d0ULL, 0x6044efa340f7d122ULL, 0x9562f1da3878cfd1ULL,
        0x4cd76dcb48ad898fULL, 0x956c7bc007773be6ULL, 0xd8098111715d1b7aULL,
        0x7825ef73528041d7ULL, 0xec3c7358e325b9eeULL, 0xdfff80a3b811cf5eULL,
        0x15f120db46900240ULL, 0x0e91adf07c0a909aULL, 0xdb78e0341f471647ULL,
        0xbd8d94bff30f52a5ULL, 0xb252c49d83e83b8fULL, 0x6899ef6152e9fab3ULL,
        0x5e830c15eb20de44ULL, 0x0bb6534e279cd74dULL, 0x24601364bcf5c41eULL,
        0xc56bc2d6178d3c12ULL, 0xc973d66cf485e61eULL, 0xe9669ea118c7d066ULL,
        0x4e83ca3c4c61be07ULL, 0x29436c2adf39df5bULL, 0x5bcc1b89c8e60129ULL,
        0x6c12477dda60dfa9ULL, 0x10b7d143a893e5cfULL, 0x8e8788479eb9c46cULL,
        0xca1fb615177b23bfULL, 0x12bbfd94ed6c1639ULL, 0xc104cb2829008604ULL,
        0xd503c709810d19caULL, 0x680e124d816baa64ULL, 0xe3fd0b85a804f87bULL,
        0x795f496c2f0e40dcULL, 0xa279b2d3914fdf44ULL, 0x1298f51a1e686c59ULL,
        0x12da5425b8b0659aULL, 0x0883762cacf3d7deULL, 0xc2248d1624c3ab56ULL,
        0x986afca42415be51ULL, 0x77e6db9411e3acb4ULL, 0x1ec13c3a6da68977ULL,
        0xa81ae0f29d3b4f18ULL, 0xe333ddfeafad232cULL, 0xf6e225f4be604f3fULL,
        0x1822533dee81bed0ULL, 0x9bb9e99f32aa06f4ULL, 0x3f50d7b8ada62a72ULL,
        0x80991888ed430c8bULL, 0x5db79a9e513ccf90ULL, 0x80ed2f4e2b101c26ULL,
        0xf87942c8414c6924ULL, 0xa5c916b00a43057bULL, 0x3344cc2b85cc578cULL,
        0xdab43a5b8f1eba4cULL, 0x6969ae0f5d408e75ULL, 0x2ad603814ac7f2b4ULL,
        0xc48627f1158b8569ULL, 0x07506aac05b3a2e7ULL, 0x1bba4937e8fff4acULL,
        0x6a5a53b44dea2992ULL, 0x0f63fef9b670f6d4ULL, 0x850e5b4311817dfcULL,
        0x731432acbdef6cc9ULL, 0xcf6bec49d83445cbULL, 0x1c37844c43115731ULL,
        0x4b1d4d0d41359e06ULL, 0x94d5627a3d60914bULL, 0xe5e57c7c36a53731ULL,
        0x24d9f088235616b8ULL, 0x422753a71ffbd8bfULL, 0x7401083cfef42de1ULL,
        0x73fad35162e46eb3ULL, 0x6e5007979f6091baULL, 0xaf26aa1f2b13db8dULL,
        0x076b75d2a7576c7aULL, 0x5cbc4279f99c4052ULL, 0xd046b8447701b010ULL,
        0xc02a64493891db3cULL, 0xf284e47131572bdbULL, 0x94ef806e82de130aULL,
        0x91a367464bd3e33dULL, 0x98ed5c62ce5e8003ULL, 0xbe0dab7966830631ULL,
        0x4c4df25842688b93ULL, 0xdc28c5012b41de83ULL, 0xbaeed57fa456fac1ULL,
        0xc3a5c276b537a890ULL, 0x0bb622a3afff2e8fULL, 0x203c7130964ff24cULL,
        0x37cfb4a8782aaa42ULL, 0x0aa01d65bd46aa08ULL, 0x22d9af7dc9e68db1ULL,
        0xf0e31d8334e77039ULL, 0x946c91d38e8a676dULL, 0x106fb8d10abd15e9ULL,
        0xd4fd9624793346a9ULL, 0x6a2dcdc0a8483eb3ULL, 0x47a2f5afa9b1698aULL,
        0x36ac1e683ae4c990ULL, 0x0541ca4336fb5db3ULL, 0xa4aa8e9d07efefd2ULL,
        0x7f14e91b7d059861ULL, 0x7e861b1d6e6159f9ULL, 0x0d2ca4917cf25b7bULL,
        0xc7553885fc0c3c84ULL, 0xb4b47f44067cd449ULL, 0x35cd7fce20e6a5f0ULL,
        0x009b73dc6a6d133dULL, 0xa3fbb325fc4e0cdbULL, 0x33ae306162bb10faULL,
        0x07771e22071ed72cULL, 0xaa121eb3c4889037ULL, 0xf5a839494b992b60ULL,
        0x229c81e208775235ULL, 0xf2d1275cd7163061ULL, 0xe90524dbc8209329ULL,
        0x5d8462a860410a04ULL, 0xa01f44c6718cd890ULL, 0xc4b3b1b230fb3298ULL,
        0x22b8ea097d6d1e4dULL, 0x2a3383d8ed6a1a87ULL, 0x3167a536e997e940ULL,
        0xa57c58406ba86e70ULL, 0x6826abb51cc39baaULL, 0xfdc99a118439caf7ULL,
        0xbcf6bbebd3527d71ULL, 0x7a0456e2e0e2e18bULL, 0xf34a5164cc9434aeULL,
        0x2dc2730933e64c9aULL, 0x91bc4f5b718f9de8ULL, 0xdab99f2e2e71c2bcULL
    };

    ref<Random> rnd = new Random(4321);
    for (size_t i = 0; i < array_size(reference); ++i) {
        const uint64_t actual = rnd->nextULong();
        assertTrue(actual == reference[i]);
    }
}



/* Basic test: the expected value should be 0.5 */
void TestRandom::test01_mean()
{
    ref<Random> rnd = new Random;
    double sum = 0.0;
    const int N = 1000000;
    for (int i = 0; i < N; ++i) {
        sum += rnd->nextFloat();
    }
    sum /= N;
    assertEqualsEpsilon(static_cast<Float>(sum), 0.5f, 1e-3f);
}



template <typename T>
void TestRandom::test02_helper(ref<Random> & rnd)
{
    typedef rnd_traits<T> r;
    const T v = r::rnd(rnd);
    assertTrue(r::min() <= v && v < r::max());
}


/* Check that the returned values are within range */
void TestRandom::test02_range()
{
    const int N = 1000;

    ref<Random> rnd = new Random;
    for (int i = 0; i < N*N; ++i) {
        test02_helper<Float>(rnd);
        test02_helper<uint32_t>(rnd);
        test02_helper<uint64_t>(rnd);
    }

    std::srand(0x695dff6c);
    for (int k = 0; k < N; ++k) {
        const unsigned int max = static_cast<unsigned int>(std::rand()) + 1;
        for (int i = 0; i < N; ++i) {
            const uint32_t v1 = rnd->nextUInt(max);
            assertTrue(v1 < max);
            const size_t v2   = rnd->nextSize(max);
            assertTrue(v2 < max);
        }
    }
}



// Simple Monte Carlo estimation of the area of a circle of radius 0.5
void TestRandom::test03_circle()
{
    ref<Random> rnd = new Random;
    int hit = 0;
    const int N = 10000000;

    for (int i = 0; i < N; ++i) {
        const Float dx = rnd->nextFloat() - 0.5f;
        const Float dy = rnd->nextFloat() - 0.5f;
        if (dx*dx + dy*dy <= 0.25f) {
            ++hit;
        }
    }
    const Float area = static_cast<Float>(hit) / N;
    assertEqualsEpsilon(area, 0.25f * M_PI, 1e-4f);
}



/* Binned version of a Chi-squared test */
void TestRandom::test04_uniform_distribution()
{
    const int N = 1000;
    const double alpha = 0.01;
    const double alpha_pvals = 0.0001;

    const double alphaT = sidak_correction(alpha, N);
    Log(EDebug, "Original alpha: %g, adjusted: %g", alpha, alphaT);

    ref<Random> rnd = new Random;
    UnifromTestPVals pTest;

    for (int i = 0; i < N; ++i) {
        const double pval = uniform_distribution_instance(rnd, alphaT);
        pTest.add(pval);
    }

    const double pval = pTest.pval();
    Log(EDebug, "Level 2 chi-squared pval: %g", pval);
    assertFalse(pval < alpha_pvals);
}



double TestRandom::uniform_distribution_instance(ref<Random> & rnd,
    const double alpha, const int numBins, const int numPerBin)
{
    const int N = numBins * numPerBin;
    std::vector<int> bins(numBins, 0);

    for (int i = 0; i < N; ++i) {
        ++bins[rnd->nextUInt(numBins)];
    }

    // Since all the bins are expected to have the same number of elements we
    // use a fake iterator
    const double pval = chi_squared_pval(bins.begin(), bins.end(),
        fake_iterator<int>(numPerBin));

    assertFalse(pval < alpha);
    return pval;
}



void TestRandom::test05_rank()
{
    const int N = 100;
    const double alpha = 0.01;
    const double alpha_pvals = 0.0001;

    const double alphaT = sidak_correction(alpha, N);
    Log(EDebug, "Original alpha: %g, adjusted: %g", alpha, alphaT);

    ref<Random> rnd = new Random;
    UnifromTestPVals pTest;
    for (int i = 0; i < N; ++i) {
        const double pval = binary_rank_instance(rnd, alphaT);
        pTest.add(pval);
    }

    const double pval = pTest.pval();
    Log(EDebug, "Level 2 chi-squared pval: %g", pval);
    assertFalse(pval < alpha_pvals);
}



double TestRandom::binary_rank_instance(ref<Random> & rnd, const double alpha,
        const int numMatrices)
{
    // Bin together the ranks 29 and lower, since they account for less than 1%
    std::vector<int> ranks(4, 0);

    for (int i = 0; i < numMatrices; ++i) {
        BinaryMatrix32x32 matrix(rnd);
        const int rank_mod2 = matrix.rank_mod2();
        const int bin_idx = std::max(0, rank_mod2 - 29);
        ++ranks[bin_idx];
    }

    assert(ranks[0]+ranks[1]+ranks[2]+ranks[3] == numMatrices);
    double f = 100.0 / numMatrices;
    Log(EDebug, "  Count: <=29: %6.3f%%, 30: %6.3f%%, 31: %6.3f%%, 32: %6.3f%%",
        f*ranks[0], f*ranks[1], f*ranks[2], f*ranks[3]);

    // For a binary matrix of size L-by-k, the probability of having
    // rank x (modulus-2) is given by: (LaTeX code)
    //   2^{x(L+k-x)-L k}\prod _{i=0}^{x-1} \frac{\left(1-2^{i-L}\right)\left(1-2^{i-k}\right)}{1-2^{i-x}}
    //
    // Thus, for 32x32 random binary matrices:
    //   P(rank == 32) ~= 0.288788
    //   P(rank == 31) ~= 0.577576
    //   P(rank == 30) ~= 0.128350
    //   P(rank <= 29) ~= 0.005286
    const double expected[] = {
        0.005285450249787358 * numMatrices, 0.1283502644231667 * numMatrices,
        0.5775761901732048   * numMatrices, 0.2887880951538411 * numMatrices};

    // Chi-squared P-value
    const double pval = chi_squared_pval(ranks.begin(), ranks.end(), expected);

    assertFalse(pval < alpha);
    return pval;
}



// Test that the probability of two random integers m,n are relatively prime is:
//  P(gcd(m,n) == 1) = [Riemann_zeta(2)]^-1 = 6 / [pi^2] = 0.60792
// http://mathworld.wolfram.com/RelativelyPrime.html [February 2012]
//
void TestRandom::test06_relative_primes()
{
    using boost::math::gcd;
    ref<Random> rnd = new Random;
    int relative_primes = 0;
    const int N = 100000;

    for (int i = 0; i < N; ++i) {
        const uint64_t a = rnd->nextULong();
        const uint64_t b = rnd->nextULong();
        if (gcd(a, b) == 1) {
            ++relative_primes;
        }
    }
    const Float ratio = static_cast<Float>(relative_primes) / N;
    Log(EDebug, "  Relative primes ratio %g", ratio);
    assertEqualsEpsilon(ratio, static_cast<Float>(0.6079271018540266), 1e-3f);
}



void TestRandom::test07_uniform_distribution_ks()
{
    const int N = 1000;
    const double alpha = 0.01;
    const double alpha_pvals = 0.0001;

    const double alphaT = sidak_correction(alpha, N);
    Log(EDebug, "Original alpha: %g, adjusted: %g", alpha, alphaT);

    ref<Random> rnd = new Random;
    UnifromTestPVals pTest;
    for (int i = 0; i < N; ++i) {
        const double pval = uniform_distribution_ks_instance(rnd, alphaT);
        pTest.add(pval);
    }

    const double pval = pTest.pval();
    Log(EDebug, "Level 2 chi-squared pval: %g", pval);
    assertFalse(pval < alpha_pvals);
}



double TestRandom::uniform_distribution_ks_instance(ref<Random> & rnd,
    const double alpha, const int numSamples)
{
    std::vector<Float> values(numSamples);
    for (int i = 0; i != numSamples; ++i) {
        values[i] = rnd->nextFloat();
    }
    const double pval = ks_pval(values.begin(), values.end(), UniformCDF());
    Log(EDebug, "  KS pval: %g", pval);
    assertFalse(pval < alpha);
    return pval;
}



void TestRandom::test08_serialize()
{
    const int N = 2000;
    std::vector<uint64_t> vec;
    typedef std::vector<uint64_t>::const_iterator citerator;

    ref<InstanceManager> manager = new InstanceManager;
    ref<Stream> s = new MemoryStream;
    {
        // Create a new random and serialize it in its initial state
        ref<Random> rndBase = new Random(0x4daaccdcbcbe32dcULL);
        rndBase->serialize(s, manager);
        assertTrue(s->getPos() != 0);

        // Now generate some numbers which must be matched afterwards
        for (int i = 0; i < N; ++i) {
            vec.push_back(rndBase->nextULong());
        }
    }

    // Create another random from the stream and check that the numbers match
    s->seek(0);
    assertTrue(s->getPos() == 0);
    ref<Random> rnd = new Random(s, manager);
    for (citerator it = vec.begin(); it != vec.end(); ++it) {
        const uint64_t actual   = rnd->nextULong();
        const uint64_t expected = *it;
        if (actual != expected) {
            Log(EWarn, "Expected: %#llx, actual: %#llx", expected, actual);
        }
        assertTrue(actual == expected);
    }
}



// Test that the member function Random::set works as intended
void TestRandom::test09_set()
{
    const int N = 1000000;
    ref<Random> rnd1 = new Random(1234);
    ref<Random> rnd2 = new Random(5678);

    for (int i = 0; i < N; ++i) {
        const uint64_t v1 = rnd1->nextULong();
        const uint64_t v2 = rnd2->nextULong();
        assertTrue(v1 != v2);
    }

    rnd1->set(rnd2);
    for (int i = 0; i < N; ++i) {
        const uint64_t v1 = rnd1->nextULong();
        const uint64_t v2 = rnd2->nextULong();
        assertTrue(v1 == v2);
    }

    rnd1->seed(rnd2);
    for (int i = 0; i < N; ++i) {
        const uint64_t v1 = rnd1->nextULong();
        const uint64_t v2 = rnd2->nextULong();
        assertTrue(v1 != v2);
    }
}



// Simple benchmark based on the mean test
void TestRandom::benchmark()
{
    const int N1 = 1000;
    const int N2 = 1000000;
    const Float epsilon = static_cast<Float>(1e-4);
    Float estimate = 0;
    Float current  = 0;

    ref<Random> rnd = new Random;

    // Warm-up
    for (int i = 0; i < N2; ++i) {
        current += rnd->nextFloat();
    }
    estimate += current / N2;

    // Actual benchmark
    ref<Timer> timer = new Timer;
    for (int i = 0; i < N1; ++i) {
        current = 0;
        for (int j = 0; j < N2; ++j) {
            current += rnd->nextFloat();
        }
        estimate += current / N2;
    }
    const Float seconds = timer->getSeconds();

    const int N = N1 * N2;
    Log(EInfo, "Generated %.1fM random numbers in %.2f s (%.3f M-random/s)",
        1e-6 * N, seconds, 1e-6 * N / seconds);
    estimate /= (N1+1);
    assertEqualsEpsilon(estimate, static_cast<Float>(0.5), epsilon);
}



void TestRandom::test05_rank_validation()
{
    static const char * test_data = ""
" 29 "
"934570178 1818650468 408011057 2633465384 1189357396 7132812 1115149267 "
"4144625854 3394558166 2481291603 3286485936 348336557 1516718101 3513035402 "
"2518314415 1043572362 1703798231 1845098087 3389283412 924141955 1500645532 "
"894152049 1099775424 1620151351 53312758 4240093792 2840064169 651879355 "
"349410984 4127881193 1073863074 815621553"
" 28 "
"1321871512 481771948 2518452259 3803302007 2031526828 924085365 2639828258 "
"2129602798 477539808 596613155 770303758 3014185450 1755957682 858966567 "
"4076851383 2956375876 2595769487 3788320775 1997036153 3891940070 3671745174 "
"2183559570 2302115729 1266026172 293570257 940537108 4108870857 1935603129 "
"596114944 2949545818 2389748847 3864099202"
" 32 "
"3457448930 2950179555 3334461710 1829897390 3667254803 1293645867 487559943 "
"2523512728 601334015 3247436681 2455090053 3361224013 1709775768 4174078144 "
"106826005 2006073602 1249902020 587441769 1745194657 2317487761 3730927443 "
"3785224664 209380234 4131882346 25599294 3468384082 525776797 3837201811 "
"944565628 1873637896 3744067778 2771652562"
" 30 "
"95189795 261233384 143196684 2163030708 3758681608 2842486163 4247308938 "
"647338281 92460902 3719788414 4115324765 1974230608 1058518282 396959343 "
"4169649845 4132002696 3147224738 351594319 3693920125 1371608093 999466035 "
"3792322985 3218315611 191086562 1002326448 3943129189 1094604251 3856315038 "
"493681781 1780437586 2047830815 117020576"
" 31 "
"3011170709 1177324815 2228178685 311351348 1997356083 2766886690 2367570474 "
"1537934837 2018043995 3632309606 4244859096 2373774622 941330676 1504626312 "
"1175703768 2269719187 3831573389 353279281 2151473632 883714442 1171942135 "
"2103632104 957623454 4093375032 498126609 642658564 2791104065 3957581843 "
"4143394920 74227394 3274092581 2073834508"
" 31 "
"117156631 2718403680 1832840156 1374368646 3350677187 3698633219 740230718 "
"3628905823 2906357122 2363461099 250267267 1778241769 2935367473 1057527317 "
"686616291 4078677907 1282501021 711454124 1923381554 1467371560 394773832 "
"2815630891 2619127968 1054855855 1917937710 648741454 4268831639 3929865586 "
"3033057205 3500398781 2203122514 2955287792 " "";

    // Use a string stream to illustrate a way to use a text file to load
    // reference ranks and the matrices' rows as integers
    std::istringstream is(test_data);
    assertFalse(!is);
    int refRank = -1;
    size_t counter = 0;
    while (!is.eof()) {
        is >> refRank;
        if (!is) break;
        assertTrue(0 <= refRank && refRank <= 32);
        BinaryMatrix32x32 matrix(is);
        if (!is) break;
        const int rank = matrix.rank_mod2();
        assertEquals(rank, refRank);
        ++counter;
    }
    Log(EDebug, "  Successfully tested " SIZE_T_FMT " matrices for rank-mod2",
        counter);
}



MTS_EXPORT_TESTCASE(TestRandom, "Testcase for (pseudo) random number generation")

MTS_NAMESPACE_END
