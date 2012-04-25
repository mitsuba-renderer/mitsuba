/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__UTIL_H)
#define __UTIL_H

#include <boost/static_assert.hpp>

MTS_NAMESPACE_BEGIN
 
/*! \addtogroup libcore */
/*! @{ */

// -----------------------------------------------------------------------
//! @{ \name String-related utility functions
// -----------------------------------------------------------------------

/**
 * \brief Given a list of delimiters, tokenize
 * a std::string into a vector of strings
 */
extern MTS_EXPORT_CORE std::vector<std::string> tokenize(
	const std::string &string,
	const std::string &delim
);

/// Trim spaces (' ', '\\n', '\\r', '\\t') from the ends of a string
extern MTS_EXPORT_CORE std::string trim(const std::string& str);

/// Indent a string (Used for recursive toString() structure dumping)
extern MTS_EXPORT_CORE std::string indent(const std::string &string, int amount=1);

/// Wrapped snprintf
extern MTS_EXPORT_CORE std::string formatString(const char *pFmt, ...);

/**
 * \brief Convert a time difference (in ms) to a string representation
 * \param time Time value in milliseconds
 * \param precise When set to true, a higher-precision string representation
 * is generated.
 */
extern MTS_EXPORT_CORE std::string timeString(Float time, bool precise = false);

/// Turn a memory size into a human-readable string
extern MTS_EXPORT_CORE std::string memString(size_t size);

/// Return a string representation of a list of objects
template<class Iterator> std::string containerToString(const Iterator &start, const Iterator &end) {
	std::ostringstream oss;
	oss << "{" << endl;
	Iterator it = start;
	while (it != end) {
		oss << "  " << indent((*it)->toString());
		++it;
		if (it != end) 
			oss << "," << endl;
		else
			oss << endl;
	}
	oss << "}";
	return oss.str();
}

/// Simple functor for sorting string parameters by length and content
struct SimpleStringOrdering {
	bool operator()(const std::string &a, const std::string &b) const {
		if (a.length() == b.length())
			return a < b;
		return a.length() < b.length();
	}
};

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Miscellaneous
// -----------------------------------------------------------------------

/// Allocate an aligned region of memory
extern MTS_EXPORT_CORE void * __restrict allocAligned(size_t size);

#if defined(WIN32)
/// Return a string version of GetLastError()
extern std::string MTS_EXPORT_CORE lastErrorText();
#endif

/// Free an aligned region of memory
extern MTS_EXPORT_CORE void freeAligned(void *ptr);

/// Determine the number of available CPU cores
extern MTS_EXPORT_CORE int getProcessorCount();

/// Return the host name of this machine
extern MTS_EXPORT_CORE std::string getHostName();

/// Return the fully qualified domain name of this machine 
extern MTS_EXPORT_CORE std::string getFQDN();

/**
 * \brief Enable floating point exceptions (to catch NaNs, overflows, 
 * arithmetic with infinity). 
 *
 * On Intel processors, this applies to both x87 and SSE2 math
 *
 * \return \c true if floating point exceptions were active
 * before calling the function
 */
extern MTS_EXPORT_CORE bool enableFPExceptions();

/**
 * \brief Disable floating point exceptions
 *
 * \return \c true if floating point exceptions were active
 * before calling the function
 */
extern MTS_EXPORT_CORE bool disableFPExceptions();

/// Restore floating point exceptions to the specified state
extern MTS_EXPORT_CORE void restoreFPExceptions(bool state);

/// Cast between types that have an identical binary representation.
template<typename T, typename U> inline T union_cast(const U &val) {
	BOOST_STATIC_ASSERT(sizeof(T) == sizeof(U));

	union {
		U u;
		T t;
	} caster = {val};

	return caster.t;
}

/// Swaps the byte order of the underlying representation
template<typename T> inline T endianness_swap(T value) {
	union {
		T value;
		uint8_t byteValue[sizeof(T)];
	} u;

	u.value = value;
	std::reverse(&u.byteValue[0], &u.byteValue[sizeof(T)]);
	return u.value;
}

/**
 * \brief Apply an arbitrary permutation to an array in linear time
 * 
 * This algorithm is based on Donald Knuth's book
 * "The Art of Computer Programming, Volume 3: Sorting and Searching"
 * (1st edition, section 5.2, page 595)
 *
 * Given a permutation and an array of values, it applies the permutation
 * in linear time without requiring additional memory. This is based on
 * the fact that each permutation can be decomposed into a disjoint set
 * of permutations, which can then be applied individually.
 *
 * \param data
 *     Pointer to the data that should be permuted
 * \param perm
 *     Input permutation vector having the same size as \c data. After
 *     the function terminates, this vector will be set to the 
 *     identity permutation.
 */
template <typename DataType, typename IndexType> void permute_inplace(
		DataType *data, std::vector<IndexType> &perm) {
	for (size_t i=0; i<perm.size(); i++) {
		if (perm[i] != i) {
			/* The start of a new cycle has been found. Save
			   the value at this position, since it will be
			   overwritten */
			IndexType j = i;
			DataType curval = data[i];

			do {
				/* Shuffle backwards */
				IndexType k = perm[j];
				data[j] = data[k];

				/* Also fix the permutations on the way */
				perm[j] = j;
				j = k;

				/* Until the end of the cycle has been found */
			} while (perm[j] != i);

			/* Fix the final position with the saved value */
			data[j] = curval; 
			perm[j] = j;
		}
	}
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Numerical utility functions
// -----------------------------------------------------------------------

static const int primeTableSize = 1000;
/// Table of the first 1000 prime numbers
extern const int MTS_EXPORT_CORE primeTable[primeTableSize];

/// sqrt(a^2 + b^2) without underflow (like 'hypot' on compilers that support C99)
extern MTS_EXPORT_CORE Float hypot2(Float a, Float b);

/// Base-2 logarithm
extern MTS_EXPORT_CORE Float log2(Float value);

/// Friendly modulo function (always positive)
extern MTS_EXPORT_CORE int modulo(int a, int b);

/// Friendly modulo function (always positive)
extern MTS_EXPORT_CORE Float modulo(Float a, Float b);

/// Integer floor function
inline int floorToInt(Float value) {
	return (int) std::floor(value);
}
/// Base-2 logarithm (32-bit integer version)
extern MTS_EXPORT_CORE int log2i(uint32_t value);

/// Base-2 logarithm (64-bit integer version)
extern MTS_EXPORT_CORE int log2i(uint64_t value);

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline int log2i(size_t value) {
	if (sizeof(size_t) == 8)
		return log2i((uint64_t) value);
	else
		return log2i((uint32_t) value);
}
#endif

/// Check if an integer is a power of two (unsigned 32 bit version)
inline bool isPowerOfTwo(uint32_t i) { return (i & (i-1)) == 0; }

/// Check if an integer is a power of two (signed 32 bit version)
inline bool isPowerOfTwo(int32_t i) { return i > 0 && (i & (i-1)) == 0; }

/// Check if an integer is a power of two (64 bit version)
inline bool isPowerOfTwo(uint64_t i) { return (i & (i-1)) == 0; }

/// Check if an integer is a power of two (signed 64 bit version)
inline bool isPowerOfTwo(int64_t i) { return i > 0 && (i & (i-1)) == 0; }

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline bool isPowerOfTwo(size_t value) {
	if (sizeof(size_t) == 8) /// will be optimized away
		return isPowerOfTwo((uint64_t) value);
	else
		return isPowerOfTwo((uint32_t) value);
}
#endif

/// Round an integer to the next power of two
extern MTS_EXPORT_CORE uint32_t roundToPowerOfTwo(uint32_t i);

/// Round an integer to the next power of two (64 bit version)
extern MTS_EXPORT_CORE uint64_t roundToPowerOfTwo(uint64_t i);

#if defined(MTS_AMBIGUOUS_SIZE_T)
/// Round an integer to the next power of two
inline size_t roundToPowerOfTwo(size_t value) {
	if (sizeof(size_t) == 8) /// will be optimized away
		return (size_t) roundToPowerOfTwo((uint64_t) value);
	else
		return (size_t) roundToPowerOfTwo((uint32_t) value);
}
#endif

/// Windowed sinc filter (Lanczos envelope, tau=number of cycles)
extern MTS_EXPORT_CORE Float lanczosSinc(Float t, Float tau = 2);

/**
 * \brief Solve a quadratic equation of the form a*x^2 + b*x + c = 0.
 * \return \c true if a solution could be found
 */
extern MTS_EXPORT_CORE bool solveQuadratic(Float a, Float b, 
	Float c, Float &x0, Float &x1);

/**
 * \brief Solve a double-precision quadratic equation of the 
 * form a*x^2 + b*x + c = 0.
 * \return \c true if a solution could be found
 */
extern MTS_EXPORT_CORE bool solveQuadraticDouble(double a, double b, 
	double c, double &x0, double &x1);

/**
 * \brief Evaluate a cubic spline interpolant of a regularly sampled 1D function
 * 
 * This implementation uses Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param p
 *      Evaluation point of the interpolant
 * \param data
 *      Floating point array containing \c nKnots regularly spaced evaluations
 *      in the range \a [a,b] of the function to be approximated.
 * \param min 
 *      Position of the first knot
 * \param max
 *      Position of the last knot
 * \param size 
 *      Total number of knots
 * \return
 *      The interpolated value or zero when \a t lies outside of \a [a,b]
 */
extern MTS_EXPORT_CORE Float interpCubic1D(Float p, const Float *data, 
		Float min, Float max, size_t size);

/**
 * \brief Evaluate a cubic spline interpolant of a regularly sampled 2D function
 * 
 * This implementation uses a tensor product of Catmull-Rom splines, i.e. it uses 
 * finite differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param p
 *      Evaluation point of the interpolant
 * \param data
 *      Floating point array containing \c nKnots regularly spaced evaluations
 *      in the range \a [a,b] of the function to be approximated.
 * \param min
 *      Position of the first knot on each dimension
 * \param max
 *      Position of the last knot on each dimension
 * \param size
 *      Total number of knots for each dimension
 * \return
 *      The interpolated value or zero when \a t lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float interpCubic2D(const Point2 &p, const Float *data, 
		const Point2 &min, const Point2 &max, const Size2 &size);

/**
 * \brief Evaluate a cubic spline interpolant of a regularly sampled 3D function
 * 
 * This implementation uses a tensor product of Catmull-Rom splines, i.e. it uses 
 * finite differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param p
 *      Evaluation point of the interpolant
 * \param data
 *      Floating point array containing \c nKnots regularly spaced evaluations
 *      in the range \a [a,b] of the function to be approximated.
 * \param min
 *      Position of the first knot on each dimension
 * \param max
 *      Position of the last knot on each dimension
 * \param size
 *      Total number of knots for each dimension
 * \return
 *      The interpolated value or zero when \a t lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float interpCubic3D(const Point3 &p, const Float *data, 
		const Point3 &min, const Point3 &max, const Size3 &size);

/**
 * \brief Calculate the radical inverse function
 *
 * (Implementation based on "Instant Radiosity" by Alexander Keller 
 * in Computer Graphics Proceedings, Annual Conference Series, 
 * SIGGRAPH 97, pp. 49-56. 
 */
extern MTS_EXPORT_CORE Float radicalInverse(int b, size_t i);

/**
 * \brief Incrementally calculate the radical inverse function
 *
 * (Implementation based on "Instant Radiosity" by Alexander Keller 
 * in Computer Graphics Proceedings, Annual Conference Series, 
 * SIGGRAPH 97, pp. 49-56. 
 */
extern MTS_EXPORT_CORE Float radicalInverseIncremental(int b, Float x);

/** 
 * Rational approximation to the inverse normal 
 * cumulative distribution function
 * Source: http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c
 * \author Peter J. Acklam
 */
extern MTS_EXPORT_CORE double normalQuantile(double p);

//// Convert radians to degrees
inline Float radToDeg(Float value) { return value * (180.0f / M_PI); }

/// Convert degrees to radians
inline Float degToRad(Float value) { return value * (M_PI / 180.0f); }

/// Simple floating point clamping function
inline Float clamp(Float value, Float min, Float max) {
	if (value < min)
		return min;
	else if (value > max)
		return max;
	else return value;
}

/// Simple integer clamping function
inline int clamp(int value, int min, int max) {
	if (value < min)
		return min;
	else if (value > max)
		return max;
	else return value;
}

/// Linearly interpolate between two values
inline Float lerp(Float t, Float v1, Float v2) {
    return ((Float) 1 - t) * v1 + t * v2;
}

/// S-shaped smoothly varying interpolation between two values
inline Float smoothStep(Float min, Float max, Float value) {
    Float v = clamp((value - min) / (max - min), (Float) 0, (Float) 1);
    return v * v * (-2 * v  + 3);
}

/**
 * \brief Numerically well-behaved routine for computing the angle 
 * between two unit direction vectors
 *
 * This should be used wherever one is tempted to compute the
 * arc cosine of a dot product!
 *
 * Proposed by Don Hatch at
 * http://www.plunk.org/~hatch/rightway.php
 */
template <typename VectorType> inline Float unitAngle(const VectorType &u, const VectorType &v) {
	if (dot(u, v) < 0)
		return M_PI - 2 * std::asin((v+u).length()/2);
	else
		return 2 * std::asin((v-u).length()/2);
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Warping and sampling-related utility functions
// -----------------------------------------------------------------------

/**
 * \brief Solve a 2x2 linear equation system using basic linear algebra
 */
extern MTS_EXPORT_CORE bool solveLinearSystem2x2(const Float a[2][2], const Float b[2], Float x[2]);


/// Complete the set {a} to an orthonormal base
extern MTS_EXPORT_CORE void coordinateSystem(const Vector &a, Vector &b, Vector &c);

/**
 * \brief Generate (optionally jittered) stratified 1D samples
 * \param random Source of random numbers
 * \param dest A pointer to a floating point array with at least
 *             count entries
 * \param count The interval [0, 1] is split into count strata
 * \param jitter Randomly jitter the samples?
 */
extern MTS_EXPORT_CORE void stratifiedSample1D(Random *random, Float *dest, 
	int count, bool jitter);

/**
 * \brief Generate (optionally jittered) stratified 2D samples
 * \param random Source of random numbers
 * \param dest A pointer to a floating point array 
 * \param countX The X axis interval [0, 1] is split into countX strata
 * \param countY The Y axis interval [0, 1] is split into countY strata
 * \param jitter Randomly jitter the samples?
 */
extern MTS_EXPORT_CORE void stratifiedSample2D(Random *random, Point2 *dest, 
	int countX, int countY, bool jitter);

/// Generate latin hypercube samples
extern MTS_EXPORT_CORE void latinHypercube(Random *random, Float *dest, size_t nSamples, size_t nDim);

/// Convert spherical coordinates to a direction
extern MTS_EXPORT_CORE Vector sphericalDirection(Float theta, Float phi);

/// Convert a direction to spherical coordinates
extern MTS_EXPORT_CORE Point2 toSphericalCoordinates(const Vector &v);

/// Sample a vector on the unit sphere (PDF: 1/(4 * PI), wrt. solid angles)
extern MTS_EXPORT_CORE Vector squareToSphere(const Point2 &sample);

/// Sample a vector on the unit hemisphere (PDF: 1/(2 * PI), wrt. solid angles)
extern MTS_EXPORT_CORE Vector squareToHemisphere(const Point2 &sample);

/// Sample a vector on the unit hemisphere (PDF: cos(theta) / PI, wrt. solid angles)
extern MTS_EXPORT_CORE Vector squareToHemispherePSA(const Point2 &sample);

/// Sample a vector that lies in a cone of angles
extern MTS_EXPORT_CORE Vector squareToCone(Float cosCutoff, const Point2 &sample);
extern MTS_EXPORT_CORE Float squareToConePdf(Float cosCutoff);

/// Sample a vector on a 2D disk (PDF: 1/(2 * PI))
extern MTS_EXPORT_CORE Point2 squareToDisk(const Point2 &sample);

/// Low-distortion concentric square to disk mapping by Peter Shirley (PDF: 1/(2 * PI))
extern MTS_EXPORT_CORE Point2 squareToDiskConcentric(const Point2 &sample);

/// Low-distortion concentric disk to square mapping 
extern MTS_EXPORT_CORE Point2 diskToSquareConcentric(const Point2 &sample);

/// Convert an uniformly distributed square sample into barycentric coordinates
extern MTS_EXPORT_CORE Point2 squareToTriangle(const Point2 &sample);

/// Sample a point on a 2D standard normal distribution (uses the Box-Muller transformation)
extern MTS_EXPORT_CORE Point2 squareToStdNormal(const Point2 &sample);

//! @}
// -----------------------------------------------------------------------

/**
 * \brief Calculates the unpolarized fresnel reflection coefficient for a 
 * dielectric material
 *
 * \param cosThetaI
 * 		Cosine of the angle between the normal and the incident ray
 * \param cosThetaT
 * 		Cosine of the angle between the normal and the transmitted ray
 * \param etaI
 * 		Refraction coefficient at the incident direction
 * \param etaT
 * 		Refraction coefficient at the transmitted direction
 */
extern MTS_EXPORT_CORE Float fresnelDielectric(Float cosThetaI, 
		Float cosThetaT, Float etaI, Float etaT);

/**
 * \brief Calculates the unpolarized fresnel reflection coefficient on
 * an interface to a conductor.
 *
 * \param cosThetaI
 * 		Cosine of the angle between the normal and the incident ray
 * \param eta
 * 		Real refractive index (wavelength-dependent)
 * \param k
 * 		Imaginary refractive index (wavelength-dependent)
 */
extern MTS_EXPORT_CORE Spectrum fresnelConductor(Float cosThetaI, 
		const Spectrum &eta, const Spectrum &k);

/**
 * \brief Calculates the unpolarized fresnel reflection coefficient for a 
 * dielectric material. Handles incidence from either sides.
 *
 * \param cosThetaI
 * 		Cosine of the angle between the normal and the incident ray
 * \param extIOR
 * 		Refraction coefficient outside of the material
 * \param intIOR
 * 		Refraction coefficient inside the material
 */
extern MTS_EXPORT_CORE Float fresnel(Float cosThetaI, Float extIOR,
		Float intIOR);

/**
 * \brief Calculates the diffuse unpolarized fresnel reflectance of
 * a dielectric material (sometimes referred to as "Fdr"). 
 *
 * This value quantifies what fraction of completely diffuse incident 
 * illumination will be reflected by a dielectric material on average.
 *
 * \param eta
 *      Relative refraction coefficient, i.e. etaT/etaI
 * \param fast
 *      Compute an approximate value? If set to \c true, the 
 *      implementation will use a polynomial approximation with
 *      a max relative error of ~0.5% on the interval 0.5 < \c eta < 2.
 *      When \c fast=false, the code will use Gauss-Lobatto quadrature 
 *      to compute the diffuse reflectance more accurately, and for
 *      a wider range of refraction coefficients, but at a cost
 *      in terms of performance.
 */
extern MTS_EXPORT_CORE Float fresnelDiffuseReflectance(
	Float eta, bool fast = false);

/*! @} */

MTS_NAMESPACE_END

#endif /* __UTIL_H */
