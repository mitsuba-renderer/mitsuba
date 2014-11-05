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

#if !defined(__MITSUBA_CORE_MATH_H_)
#define __MITSUBA_CORE_MATH_H_

MTS_NAMESPACE_BEGIN

/**
 * Contains elementary 1D math functions that were either not provided by the standard,
 * or which are not consistently provided on all platforms/compilers
 */
namespace math {

/// Cross-platform implementation of the error function
extern MTS_EXPORT_CORE Float erf(Float x);

/// Cross-platform implementation of the inverse error function
extern MTS_EXPORT_CORE Float erfinv(Float x);

/// sqrt(a^2 + b^2) without range issues (like 'hypot' on compilers that support C99, single precision)
extern MTS_EXPORT_CORE float hypot2(float a, float b);

/// sqrt(a^2 + b^2) without range issues (like 'hypot' on compilers that support C99, double precision)
extern MTS_EXPORT_CORE double hypot2(double a, double b);

/// Base-2 logarithm (single precision)
extern MTS_EXPORT_CORE float log2(float value);

/// Base-2 logarithm (double precision)
extern MTS_EXPORT_CORE double log2(double value);

/// Generic clamping function
template <typename Scalar> inline Scalar clamp(Scalar value, Scalar min, Scalar max) {
	return std::min(max, std::max(min, value));
}

/// Linearly interpolate between two values
template <typename Scalar> inline Scalar lerp(Scalar t, Scalar v1, Scalar v2) {
    return ((Scalar) 1 - t) * v1 + t * v2;
}

/// S-shaped smoothly varying interpolation between two values
template <typename Scalar> inline Scalar smoothStep(Scalar min, Scalar max, Scalar value) {
    Scalar v = clamp((value - min) / (max - min), (Scalar) 0, (Scalar) 1);
    return v * v * (-2 * v  + 3);
}

/// Always-positive modulo function (assumes b > 0)
inline int32_t modulo(int32_t a, int32_t b) {
	int32_t r = a % b;
	return (r < 0) ? r+b : r;
}

/// Always-positive modulo function (assumes b > 0)
inline int64_t modulo(int64_t a, int64_t b) {
	int64_t r = a % b;
	return (r < 0) ? r+b : r;
}

#if defined(MTS_AMBIGUOUS_SIZE_T)
inline ssize_t modulo(ssize_t a, ssize_t b) {
	if (sizeof(ssize_t) == 8)
		return modulo((int64_t) a, (int64_t) b);
	else
		return modulo((int32_t) a, (int32_t) b);
}
#endif

/// Always-positive modulo function, single precision version (assumes b > 0)
inline float modulo(float a, float b) {
	float r = std::fmod(a, b);
	return (r < 0.0f) ? r+b : r;
}

/// Always-positive modulo function, double precision version (assumes b > 0)
inline double modulo(double a, double b) {
	double r = std::fmod(a, b);
	return (r < 0.0) ? r+b : r;
}

/// Integer floor function (single precision)
template <typename Scalar> inline int floorToInt(Scalar value) { return (int) std::floor(value); }

/// Integer ceil function (single precision)
template <typename Scalar> inline int ceilToInt(Scalar value) { return (int) std::ceil(value); }

/// Integer round function (single precision)
inline int roundToInt(float value)  {
	#if defined(__MSVC__)
		return (int) (value < 0.0f ? std::ceil(value - 0.5f) : std::floor(value + 0.5f));
	#else
		return (int) ::roundf(value);
	#endif
}

/// Integer round function (double precision)
inline int roundToInt(double value) {
	#if defined(__MSVC__)
		return (int) (value < 0.0 ? std::ceil(value - 0.5) : std::floor(value + 0.5));
	#else
		return (int) ::round(value);
	#endif
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

#if defined(__LINUX__) && defined(__x86_64__)
	/*
	   The Linux/x86_64 single precision implementations of 'exp'
	   and 'log' suffer from a serious performance regression.
	   It is about 5x faster to use the double-precision versions
	   with the extra overhead of the involved FP conversion.

	   Until this is fixed, the following aliases make sure that
	   the fastest implementation is used in every case.
	 */
	inline float fastexp(float value) {
		return (float) ::exp((double) value);
	}

	inline double fastexp(double value) {
		return ::exp(value);
	}

	inline float fastlog(float value) {
		return (float) ::log((double) value);
	}

	inline double fastlog(double value) {
		return ::log(value);
	}
#else
	inline float fastexp(float value) {
		return ::expf(value);
	}

	inline double fastexp(double value) {
		return ::exp(value);
	}

	inline float fastlog(float value) {
		return ::logf(value);
	}

	inline double fastlog(double value) {
		return ::log(value);
	}
#endif

#if defined(_GNU_SOURCE)
	inline void sincos(float theta, float *sin, float *cos) {
		::sincosf(theta, sin, cos);
	}

	inline void sincos(double theta, double *sin, double *cos) {
		::sincos(theta, sin, cos);
	}

#else
	inline void sincos(float theta, float *_sin, float *_cos) {
		*_sin = sinf(theta);
		*_cos = cosf(theta);
	}

	inline void sincos(double theta, double *_sin, double *_cos) {
		*_sin = sin(theta);
		*_cos = cos(theta);
	}
#endif

	/// Arcsine variant that gracefully handles arguments > 1 that are due to roundoff errors
	inline float safe_asin(float value) {
		return std::asin(std::min(1.0f, std::max(-1.0f, value)));
	}

	/// Arcsine variant that gracefully handles arguments > 1 that are due to roundoff errors
	inline double safe_asin(double value) {
		return std::asin(std::min(1.0, std::max(-1.0, value)));
	}

	/// Arccosine variant that gracefully handles arguments > 1 that are due to roundoff errors
	inline float safe_acos(float value) {
		return std::acos(std::min(1.0f, std::max(-1.0f, value)));
	}

	/// Arccosine variant that gracefully handles arguments > 1 that are due to roundoff errors
	inline double safe_acos(double value) {
		return std::acos(std::min(1.0, std::max(-1.0, value)));
	}

	/// Square root variant that gracefully handles arguments < 0 that are due to roundoff errors
	inline float safe_sqrt(float value) {
		return std::sqrt(std::max(0.0f, value));
	}

	/// Square root variant that gracefully handles arguments < 0 that are due to roundoff errors
	inline double safe_sqrt(double value) {
		return std::sqrt(std::max(0.0, value));
	}

	/// Simple signum function -- note that it returns the FP sign of the input (and never zero)
	inline Float signum(Float value) {
		#if defined(__WINDOWS__)
			return (Float) _copysign(1.0, value);
		#elif defined(SINGLE_PRECISION)
			return copysignf((float) 1.0, value);
		#elif defined(DOUBLE_PRECISION)
			return copysign((double) 1.0, value);
		#endif
	}

	/// Cast to single precision and round up if not exactly representable (passthrough)
	inline float castflt_up(float val) { return val; }

	/// Cast to single precision and round up if not exactly representable
	inline float castflt_up(double val) {
		union {
			float a;
			int b;
		};

		a = (float) val;
		if ((double) a < val)
			b += a < 0 ? -1 : 1;
		return a;
	}

	/// Cast to single precision and round down if not exactly representable (passthrough)
	inline float castflt_down(float val) { return val; }

	/// Cast to single precision and round down if not exactly representable
	inline float castflt_down(double val) {
		union {
			float a;
			int b;
		};

		a = (float) val;
		if ((double) a > val)
			b += a > 0 ? -1 : 1;
		return a;
	}
}; /* namespace math */

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_MATH_H_ */
