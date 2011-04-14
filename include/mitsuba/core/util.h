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
 * Convert a time difference (in ms) to a string representation
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
 * Enable floating point exceptions (to catch NaNs, overflows, 
 * arithmetic with infinity). On Intel processors, this applies
 * to both x87 and SSE2 math
 *
 * \return \c true if floating point exceptions were active
 * before calling the function
 */
extern MTS_EXPORT_CORE bool enableFPExceptions();

/**
 * Disable floating point exceptions
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
 * This algorithm is based on Donald Knuth's book
 * "The Art of Computer Programming, Volume 3: Sorting and Searching"
 * (1st edition, section 5.2, page 595)
 *
 * Given a permutation and an array of values, it applies the permutation
 * in linear time without requiring additional memory. This is based on
 * the fact that each permutation can be decomposed into a disjoint set
 * of permutations, which can then be applied individually.
 */
template <typename T> void permute_inplace(T *values, std::vector<size_t> &perm) {
	for (size_t i=0; i<perm.size(); i++) {
		if (perm[i] != i) {
			/* The start of a new cycle has been found. Save
			   the value at this position, since it will be
			   overwritten */
			size_t j = i;
			T curval = values[i];

			do {
				/* Shuffle backwards */
				size_t k = perm[j];
				values[j] = values[k];

				/* Also fix the permutations on the way */
				perm[j] = j;
				j = k;

				/* Until the end of the cycle has been found */
			} while (perm[j] != i);

			/* Fix the final position with the saved value */
			values[j] = curval;
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
inline bool isPow2(uint32_t i) { return (i & (i-1)) == 0; }

/// Check if an integer is a power of two (signed 32 bit version)
inline bool isPow2(int32_t i) { return i > 0 && (i & (i-1)) == 0; }

/// Check if an integer is a power of two (64 bit version)
inline bool isPow2(uint64_t i) { return (i & (i-1)) == 0; }

/// Check if an integer is a power of two (signed 64 bit version)
inline bool isPow2(int64_t i) { return i > 0 && (i & (i-1)) == 0; }

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline bool isPow2(size_t value) {
	if (sizeof(size_t) == 8)
		return isPow2((uint64_t) value);
	else
		return isPow2((uint32_t) value);
}
#endif

/// Round an integer to the next power of two
extern MTS_EXPORT_CORE uint32_t roundToPow2(uint32_t i);

/// Round an integer to the next power of two (64 bit version)
extern MTS_EXPORT_CORE uint64_t roundToPow2(uint64_t i);

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline size_t roundToPow2(size_t value) {
	if (sizeof(size_t) == 8)
		return (size_t) roundToPow2((uint64_t) value);
	else
		return (size_t) roundToPow2((uint32_t) value);
}
#endif

//// Windowed sinc filter (Lanczos envelope, tau=number of cycles)
extern MTS_EXPORT_CORE Float lanczosSinc(Float t, Float tau = 2);

/**
 * \brief Solve a quadratic equation of the form a*x^2 + b*x + c = 0.
 * \return \c true if a solution could be found
 */
extern MTS_EXPORT_CORE bool solveQuadratic(Float a, Float b, 
	Float c, Float &x0, Float &x1);

/**
 * Calculate the radical inverse function
 * (Implementation based on "Instant Radiosity" by Alexander Keller 
 * in Computer Graphics Proceedings, Annual Conference Series, 
 * SIGGRAPH 97, pp. 49-56. 
 */
extern MTS_EXPORT_CORE Float radicalInverse(int b, int i);

/**
 * Incrementally calculate the radical inverse function
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
extern MTS_EXPORT_CORE void latinHypercube(Random *random, Float *dest, int nSamples, int nDim);

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

/// Convert an uniformly distributed square sample into barycentric coordinates
extern MTS_EXPORT_CORE Point2 squareToTriangle(const Point2 &sample);

//! @}
// -----------------------------------------------------------------------

/**
 * Calculates the unpolarized fresnel reflection coefficient for a 
 * dielectric material
 *
 * \param cosTheta1
 * 		Cosine of the angle between the normal and the incident ray
 * \param cosTheta2
 * 		Cosine of the angle between the normal and the transmitted ray
 * \param etaExt
 * 		Refraction coefficient outside of the material
 * \param etaInt
 * 		Refraction coefficient inside the material
 */
extern MTS_EXPORT_CORE Float fresnelDielectric(Float cosTheta1, 
		Float cosTheta2, Float etaExt, Float etaInt);

/**
 * Calculates the unpolarized fresnel reflection coefficient for a 
 * dielectric material. Handles incidence from either sides.
 *
 * \param cosTheta1
 * 		Cosine of the angle between the normal and the incident ray
 * \param etaExt
 * 		Refraction coefficient outside of the material
 * \param etaInt
 * 		Refraction coefficient inside the material
 */
extern MTS_EXPORT_CORE Float fresnel(Float cosTheta1, Float etaExt,
		Float etaInt);

/**
 * Calculates the unpolarized fresnel reflection coefficient on
 * an interface to a conductor.
 *
 * \param cosTheta
 * 		Cosine of the angle between the normal and the incident ray
 * \param eta
 * 		Relative refractive index (per wavelength)
 * \param k
 * 		Absorption coefficient (per wavelength)
 */
extern MTS_EXPORT_CORE Spectrum fresnelConductor(Float cosTheta, 
		const Spectrum &eta, const Spectrum &k);

/*! @} */

MTS_NAMESPACE_END

#endif /* __UTIL_H */
