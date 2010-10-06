/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__UTIL_H)
#define __UTIL_H

#include <boost/static_assert.hpp>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Miscellaneous
// -----------------------------------------------------------------------

static const int primeTableSize = 1000;
/// Table of the first 1000 prime numbers
extern const int MTS_EXPORT_CORE primeTable[primeTableSize];

/**
 * \brief Given a list of delimiters, tokenize
 * a std::string into a vector of strings
 */
extern MTS_EXPORT_CORE std::vector<std::string> tokenize(
	const std::string &string,
	const std::string &delim
);

/// Indent a string (Used for recursive toString() structure dumping)
extern MTS_EXPORT_CORE std::string indent(const std::string &string, int amount=1);

/// Convert a time difference (in ms) to a string representation
extern MTS_EXPORT_CORE std::string timeToString(Float time);

/// Trim spaces (' ', '\\n', '\\r', '\\t') from the ends of a string
extern MTS_EXPORT_CORE std::string trim(const std::string& str);

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

/// Wrapped snprintf
extern MTS_EXPORT_CORE std::string formatString(const char *pFmt, ...);

/// Base-2 logarithm
extern MTS_EXPORT_CORE Float log2(Float value);

/// Base-2 logarithm (integer version)
extern MTS_EXPORT_CORE int log2i(int value);

/// Friendly modulo function (always positive)
extern MTS_EXPORT_CORE int modulo(int a, int b);

/// Check if an integer is a power of two
extern MTS_EXPORT_CORE bool isPowerOfTwo(unsigned int i);

/// Round an integer to the next power of two
extern MTS_EXPORT_CORE unsigned int roundToPowerOfTwo(unsigned int i);

//// Windowed sinc filter (Lanczos envelope, tau=number of cycles)
extern MTS_EXPORT_CORE Float lanczosSinc(Float t, Float tau = 2);

/**
 * \brief Solve a quadratic equation of the form a*x^2 + b*x + c = 0.
 * Returns true if a solution could be found
 */
extern MTS_EXPORT_CORE bool solveQuadratic(Float a, Float b, Float c, Float &x0, Float &x1);

/**
 * \brief Similar to solveQuadratic(), but always uses double precision independent
 * of the chosen compile-time precision.
 */
extern MTS_EXPORT_CORE bool solveQuadraticDouble(double a, double b, double c, double &x0, double &x1);

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

/// Convert radians to degrees
inline Float radToDeg(Float rad) {
    return 180.0f * rad / M_PI;
}

/// Convert degrees to radians
inline Float degToRad(Float deg) {
    return deg * M_PI / 180.0f;
}

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
extern MTS_EXPORT_CORE Float fresnelDielectric(Float cosTheta1, Float cosTheta2, 
							   Float etaExt, Float etaInt);

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
extern MTS_EXPORT_CORE Float fresnel(Float cosTheta1, Float etaExt, Float etaInt);

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
extern MTS_EXPORT_CORE Spectrum fresnelConductor(Float cosTheta, const Spectrum &eta, const Spectrum &k);
 
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

/// Cast between types, which have an identical binary representation.
template<typename T, typename U> inline T union_cast(const U &val) {
	BOOST_STATIC_ASSERT(sizeof(T) == sizeof(U));

	union {
		U u;
		T t;
	} caster = {val};

	return caster.t;
}

/// Return a string representation of a list of objects
template<class T> std::string listToString(const std::vector<T> &vec) {
	std::ostringstream oss;
	oss << "{" << endl;
	for (size_t i=0; i<vec.size(); ++i) {
		oss << "  " << indent(vec[i]->toString());
		if (i != vec.size()-1)
			oss << "," << endl;
		else
			oss << endl;
	}
	oss << "}";
	return oss.str();
}

MTS_NAMESPACE_END

#endif /* __UTIL_H */
