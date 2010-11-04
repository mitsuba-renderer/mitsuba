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

#ifndef __MITSUBA_STL_H
#define __MITSUBA_STL_H
 
/* Include some SGI STL extensions, which might be missing */
#ifdef __GNUC__
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::compose1;
#else
#include <functional>

namespace std {
	/// \cond
	// (Don't include in the documentation)
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
	/// @endcond

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

/* Forward declarations */
MTS_NAMESPACE_BEGIN
extern MTS_EXPORT_CORE void * __restrict allocAligned(size_t size);
extern MTS_EXPORT_CORE void freeAligned(void *ptr);

/**
 * \brief Aligned memory allocator for use with SSE2-based code
 * \headerfile mitsuba/core/stl.h mitsuba/mitsuba.h
 * 
 * Basic implementaiton, which forwards all calls to \ref allocAligned.
 */
template <typename T> class aligned_allocator {
public:
	typedef size_t    size_type;
	typedef ptrdiff_t difference_type;
	typedef T*        pointer;
	typedef const T*  const_pointer;
	typedef T&        reference;
	typedef const T&  const_reference;
	typedef T         value_type;

	/// \cond
	template <class U> struct rebind {
		typedef aligned_allocator<U> other;
	};
	/// \endcond

	pointer address (reference value) const {
		return &value;
	}

	const_pointer address (const_reference value) const {
		return &value;
	}
	
	aligned_allocator() throw() { }	

	aligned_allocator(const aligned_allocator&) throw() { }

	template <class U> aligned_allocator (const aligned_allocator<U>&) throw()  { }

	~aligned_allocator() throw() { }

	size_type max_size () const throw() {
		return INT_MAX;
	}


	T* __restrict allocate (size_type num, const_pointer *hint = 0) {
		return (T *) mitsuba::allocAligned(num*sizeof(T));
	}

	void construct (pointer p, const T& value) {
		*p=value;
	}

	void destroy (pointer p) {
		p->~T();
	};
	
	void deallocate (pointer p, size_type num) {
		freeAligned(p);
	}
};

MTS_NAMESPACE_END

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
#endif

