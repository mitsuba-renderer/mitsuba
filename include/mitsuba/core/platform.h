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
#if !defined(__MITSUBA_CORE_PLATFORM_H_)
#define __MITSUBA_CORE_PLATFORM_H_

/// Disable BOOST's autolinking feature
#define BOOST_ALL_NO_LIB 1

#if !defined(_OPENMP) && !defined(MTS_NO_OPENMP)
    #define MTS_NO_OPENMP
#endif

#if defined(_MSC_VER)
    #define __MSVC__
    #define __WINDOWS__

    // Don't complain about perfectly fine ISO C++
    #define _CRT_SECURE_NO_WARNINGS
    #define _CRT_NONSTDC_NO_DEPRECATE
    #define _CRT_SECURE_NO_DEPRECATE

    #define _WIN32_WINNT 0x0501 // Windows XP
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN

    #pragma warning(disable : 4251) // 'field' : class 'A' needs to have dll-interface to be used by clients of class 'B'
    #pragma warning(disable : 4800) // 'type' : forcing value to bool 'true' or 'false' (performance warning)
    #pragma warning(disable : 4996) // Secure SCL warnings

    #include <stdint.h>

    #if _MSC_VER >= 1400
        #include <memory.h>
        #include <string.h>
        #include <math.h>
        #pragma intrinsic(memset, memcmp, memcpy, strlen, strcmp, strcpy, _strset, strcat, fabs, abs)
    #endif

    #if _MSC_VER >= 1600
        #ifdef SINGLE_PRECISION
                #pragma detect_mismatch( "MTS_FLOAT_PRECISION", "SINGLE")
        #elif  DOUBLE_PRECISION
                #pragma detect_mismatch( "MTS_FLOAT_PRECISION", "DOUBLE")
        #endif
        #define MTS_STRINGIFY(s) #s
        #define MTS_XSTRINGIFY(s) MTS_STRINGIFY(s)
        #pragma detect_mismatch("MTS_SPECTRUM_SAMPLES", MTS_XSTRINGIFY(SPECTRUM_SAMPLES))
    #endif
#elif defined(__APPLE__)
    #define __OSX__
#elif defined(__linux)
    #define __LINUX__
    #if !defined(_GNU_SOURCE)
        #define _GNU_SOURCE
    #endif
#else
    #error Unknown OS
#endif

#ifdef __MSVC__
    #define MTS_DONT_EXPORT // not supported on MSVC
    #define SIZE_T_FMT "%Iu"
    #define BOOST_FILESYSTEM_NO_LIB
    #define BOOST_SYSTEM_NO_LIB
    #define MTS_EXPORT __declspec(dllexport)
    #define MTS_IMPORT __declspec(dllimport)
    #define MTS_MAY_ALIAS // not supported on Windows
#else
    #define MTS_EXPORT __attribute__ ((visibility("default")))
    #define MTS_IMPORT
    #define MTS_MAY_ALIAS __attribute__ ((__may_alias__))

    #include <stdint.h>

    #define SIZE_T_FMT "%zd"
#endif

#define MTS_MODULE_CORE 1
#define MTS_MODULE_RENDER 2
#define MTS_MODULE_HW 3
#define MTS_MODULE_BIDIR 4
#define MTS_MODULE_PYTHON 5

#if MTS_BUILD_MODULE == MTS_MODULE_CORE
    #define MTS_EXPORT_CORE MTS_EXPORT
#else
    #define MTS_EXPORT_CORE MTS_IMPORT
#endif
#if MTS_BUILD_MODULE == MTS_MODULE_RENDER
    #define MTS_EXPORT_RENDER MTS_EXPORT
#else
    #define MTS_EXPORT_RENDER MTS_IMPORT
#endif
#if MTS_BUILD_MODULE == MTS_MODULE_HW
    #define MTS_EXPORT_HW MTS_EXPORT
#else
    #define MTS_EXPORT_HW MTS_IMPORT
#endif
#if MTS_BUILD_MODULE == MTS_MODULE_BIDIR
    #define MTS_EXPORT_BIDIR MTS_EXPORT
#else
    #define MTS_EXPORT_BIDIR MTS_IMPORT
#endif
#if MTS_BUILD_MODULE == MTS_MODULE_PYTHON
    #define MTS_EXPORT_PYTHON MTS_EXPORT
#else
    #define MTS_EXPORT_PYTHON MTS_IMPORT
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__LP64__) || defined(_LP64) || defined(WIN64)
#define __64BIT__
#else
#define __32BIT__
#endif

#if !defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
#define __LITTLE_ENDIAN__ 1 // Little endian by default
#endif

#define MTS_NAMESPACE_BEGIN namespace mitsuba {
#define MTS_NAMESPACE_END }

#if defined(__GNUC__)
#define FINLINE                inline __attribute__((always_inline))
#define NOINLINE               __attribute__((noinline))
#define EXPECT_TAKEN(a)        __builtin_expect(!!(a), true)
#define EXPECT_NOT_TAKEN(a)    __builtin_expect(!!(a), false)
#elif defined(__MSVC__)
#define FINLINE                __forceinline
#define NOINLINE               __declspec(noinline)
#define MM_ALIGN16             __declspec(align(16))
#define EXPECT_TAKEN(a)        (a)
#define EXPECT_NOT_TAKEN(a)    (a)
#else
#error Unsupported compiler!
#endif

#ifdef MTS_SSE
#define SSE_STR "SSE2 enabled"
#else
#define SSE_STR "SSE2 disabled"
#endif

/* The default OpenMP implementation on OSX is seriously broken,
   for instance it segfaults when launching OpenMP threads
   from context other than the main application thread */
#if defined(__OSX__) && !defined(__INTEL_COMPILER) && !defined(MTS_NO_OPENMP)
#define MTS_NO_OPENMP
#endif

#if !defined(MTS_NO_OPENMP)
#define MTS_OPENMP
#endif

/* Compile with Boost::Filesystem v3 */
#define BOOST_FILESYSTEM_VERSION 3

#include <string>

MTS_NAMESPACE_BEGIN

#if defined(DOUBLE_PRECISION)
typedef double Float;
#elif defined(SINGLE_PRECISION)
typedef float Float;
#else
#error No precision flag was defined!
#endif

#if defined(__OSX__)
extern MTS_EXPORT_CORE void __mts_autorelease_init();
extern MTS_EXPORT_CORE void __mts_autorelease_shutdown();
extern MTS_EXPORT_CORE void __mts_autorelease_begin();
extern MTS_EXPORT_CORE void __mts_autorelease_end();
extern MTS_EXPORT_CORE std::string __mts_bundlepath();
extern MTS_EXPORT_CORE void __mts_chdir_to_bundlepath();
extern MTS_EXPORT_CORE void __mts_init_cocoa();
extern MTS_EXPORT_CORE void __mts_set_appdefaults();
#define MTS_AUTORELEASE_BEGIN() __mts_autorelease_begin();
#define MTS_AUTORELEASE_END() __mts_autorelease_end();
#define MTS_AMBIGUOUS_SIZE_T 1
#else
#define MTS_AUTORELEASE_BEGIN()
#define MTS_AUTORELEASE_END()
#endif

MTS_NAMESPACE_END


/// \cond
// Try to make MSVC++ behave a bit more like C++
// with an underlying C99 implementation
// (and don't include this in the documentation)
#if defined(_MSC_VER)

#include <float.h>

#define snprintf _snprintf
#define vsnprintf _vsnprintf

namespace mitsuba {
#if defined(__64BIT__)
typedef long long ssize_t;
#else
typedef long ssize_t;
#endif
};

namespace std {
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
};
#endif
/// \endcond

#endif /* __MITSUBA_CORE_PLATFORM_H_ */
