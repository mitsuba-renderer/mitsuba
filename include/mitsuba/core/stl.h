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
#if !defined(__MITSUBA_CORE_STL_H_)
#define __MITSUBA_CORE_STL_H_

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

	inline char tolower(char c) {
		return ::tolower(c);
	}

	inline char toupper(char c) {
		return ::toupper(c);
	}

	inline bool isnan(float f) {
		return _isnan(f);
	}

	inline bool isnan(double f) {
		return _isnan(f);
	}

	inline bool isfinite(float f) {
		return _finite(f);
	}

	inline bool isfinite(double f) {
		return _finite(f);
	}

	inline bool isinf(float f) {
		return !_finite(f);
	}

	inline bool isinf(double f) {
		return !_finite(f);
	}
#endif
};
using std::select2nd;
using std::compose1;
#endif
/// \endcond

#endif /* __MITSUBA_CORE_STL_H_ */
