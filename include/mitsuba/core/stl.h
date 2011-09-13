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

#ifndef __MITSUBA_STL_H
#define __MITSUBA_STL_H
 
/* Include some SGI STL extensions, which might be missing */
#ifdef __GNUC__
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::compose1;
#else
#include <functional>

/// \cond
// (Don't include in the documentation)
namespace std {
	template <class _Pair> struct _Select1st : public unary_function<_Pair, typename _Pair::first_type> {
		const typename _Pair::first_type& operator()(const _Pair& __x) const {
			return __x.first;
		}
	};
	
	template <class _Pair> struct _Select2nd : public unary_function<_Pair, typename _Pair::second_type> {
		const typename _Pair::second_type& operator()(const _Pair& __x) const {
			return __x.second;
		}
	};

	template <class _Pair> struct select1st : public _Select1st<_Pair> {};
	template <class _Pair> struct select2nd : public _Select2nd<_Pair> {};

	template <class _Operation1, class _Operation2> class unary_compose : public unary_function<typename _Operation2::argument_type, typename _Operation1::result_type> {
	protected:
		_Operation1 _M_fn1;
		_Operation2 _M_fn2;
	public:
		unary_compose(const _Operation1& __x, const _Operation2& __y) : _M_fn1(__x), _M_fn2(__y) {}
		typename _Operation1::result_type operator()(const typename _Operation2::argument_type& __x) const {
			return _M_fn1(_M_fn2(__x));
		}
	};

	template <class _Operation1, class _Operation2> inline unary_compose<_Operation1,_Operation2> compose1(const _Operation1& __fn1, const _Operation2& __fn2) {
		return unary_compose<_Operation1,_Operation2>(__fn1, __fn2);
	}

#if defined(_MSC_VER)
	#include <float.h>

	#define snprintf _snprintf
	#define vsnprintf _vsnprintf
	
	inline int isinf(float value) {
		int type = ::_fpclass(value);
		if (type == _FPCLASS_PINF || type == _FPCLASS_NINF)
			return 1;
		return 0;
	}
	
	inline char tolower(char c) {
		return ::tolower(c);
	}

	inline char toupper(char c) {
		return ::toupper(c);
	}

	inline int isinf(double value) {
		int type = ::_fpclass(value);
		if (type == _FPCLASS_PINF || type == _FPCLASS_NINF)
			return 1;
		return 0;
	}
#endif
};
using std::select2nd;
using std::compose1;
#endif

namespace std {
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
};

#if defined(WIN32)
inline bool mts_isnan(float f) {
	int classification = ::_fpclass(f);
	return classification == _FPCLASS_QNAN 
		|| classification == _FPCLASS_SNAN;
}

inline bool mts_isnan(double f) {
	int classification = ::_fpclass(f);
	return classification == _FPCLASS_QNAN
		|| classification == _FPCLASS_SNAN;
}
extern "C" {
	extern MTS_EXPORT_CORE float nextafterf(float x, float y);
	extern MTS_EXPORT_CORE double nextafter(double x, double y);
};
#elif defined(__clang__)
inline bool mts_isnan(float f) {
	return std::isnan(f);
}

inline bool mts_isnan(double f) {
	return std::isnan(f);
}
#else
inline bool mts_isnan(float f) {
	return std::fpclassify(f) == FP_NAN;
}

inline bool mts_isnan(double f) {
	return std::fpclassify(f) == FP_NAN;
}
#endif
/// @endcond
#endif

