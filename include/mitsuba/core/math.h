/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

namespace math {
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
}; /* namespace math */

MTS_NAMESPACE_END

#if defined(_MSC_VER) && _MSC_VER < 1800
extern "C" {
	extern MTS_EXPORT_CORE float nextafterf(float x, float y);
	extern MTS_EXPORT_CORE double nextafter(double x, double y);
};
#endif

#endif /* __MITSUBA_CORE_MATH_H_ */
