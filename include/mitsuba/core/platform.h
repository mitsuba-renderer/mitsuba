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

#if !defined(__PLATFORM_H)
#define __PLATFORM_H

/// Disable BOOST's autolinking feature
#define BOOST_ALL_NO_LIB 1

#if defined(_MSC_VER)
	/* Disable MSVC STL debug + security features (slow..!) */
	#ifdef _SECURE_SCL_THROWS
		#undef _SECURE_SCL_THROWS
	#endif
	#ifdef _SCL_SECURE_NO_WARNINGS
		#undef _SCL_SECURE_NO_WARNINGS
	#endif
	#ifdef _HAS_ITERATOR_DEBUGGING
		#undef _HAS_ITERATOR_DEBUGGING
	#endif
	#undef _STLP_DEBUG

	#define _SECURE_SCL_THROWS 0
	#define _HAS_ITERATOR_DEBUGGING 0
	#define _SCL_SECURE_NO_WARNINGS 0
	#define __MSVC__
	#define __WINDOWS__

	// Don't complain about perfectly fine ISO C++
	#define _CRT_SECURE_NO_WARNINGS
	#define _CRT_NONSTDC_NO_DEPRECATE
	#define _CRT_SECURE_NO_DEPRECATE

	#define _WIN32_WINNT 0x0501
	#define NOMINMAX
	#define WIN32_LEAN_AND_MEAN
	#include <winsock2.h> // IPv6 support
	#include <windows.h>

	#pragma warning(disable : 4251) // 'field' : class 'A' needs to have dll-interface to be used by clients of class 'B'
	#pragma warning(disable : 4800) // 'type' : forcing value to bool 'true' or 'false' (performance warning)
	#pragma warning(disable : 4996) // Secure SCL warnings
	#if _MSC_VER < 1600
		#include <stdint_msvc.h>    // Does not exist in MSVC. Use a replacement
	#else
		#include <stdint.h>
	#endif
	#if _MSC_VER >= 1400
		#include <memory.h>
		#include <string.h>
		#include <math.h>
		#pragma intrinsic(memset, memcmp, memcpy, strlen, strcmp, strcpy, _strset, strcat, fabs, abs)
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
	#define MTS_MODULE_CORE 1
	#define MTS_MODULE_RENDER 2
	#define MTS_MODULE_HW 3
	#define MTS_MODULE_BIDIR 4
	#define MTS_MODULE_PYTHON 5

	#define MTS_EXPORT __declspec(dllexport)
	#define MTS_IMPORT __declspec(dllimport)

	#if MTS_BUILD_MODULE == MTS_MODULE_CORE
		#define MTS_EXPORT_CORE __declspec(dllexport)
	#else
		#define MTS_EXPORT_CORE __declspec(dllimport)
	#endif
	#if MTS_BUILD_MODULE == MTS_MODULE_RENDER
		#define MTS_EXPORT_RENDER __declspec(dllexport)
	#else
		#define MTS_EXPORT_RENDER __declspec(dllimport)
	#endif
	#if MTS_BUILD_MODULE == MTS_MODULE_HW
		#define MTS_EXPORT_HW __declspec(dllexport)
	#else
		#define MTS_EXPORT_HW __declspec(dllimport)
	#endif
	#if MTS_BUILD_MODULE == MTS_MODULE_BIDIR
		#define MTS_EXPORT_BIDIR __declspec(dllexport)
	#else
		#define MTS_EXPORT_BIDIR __declspec(dllimport)
	#endif
	#if MTS_BUILD_MODULE == MTS_MODULE_PYTHON
		#define MTS_EXPORT_PYTHON __declspec(dllexport)
	#else
		#define MTS_EXPORT_PYTHON __declspec(dllimport)
	#endif

	#define SIZE_T_FMT "%Iu"
	#define BOOST_FILESYSTEM_NO_LIB 
	#define BOOST_SYSTEM_NO_LIB 
#else
	#define MTS_EXPORT
	#define MTS_EXPORT_CORE
	#define MTS_EXPORT_RENDER
	#define MTS_EXPORT_HW
	#define MTS_EXPORT_BIDIR
	#define MTS_EXPORT_PYTHON
	#include <stdint.h>

	#define SIZE_T_FMT "%zd"
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

/* The default OpenMP implementation on OSX is seriously broken,
 * for instance it segfaults when launching OpenMP threads
 * from any other context than the main application thread
 */
#if defined(__OSX__) && !defined(__INTEL_COMPILER)
#define MTS_BROKEN_OPENMP 1
#else
#define MTS_BROKEN_OPENMP 0
#endif

/* Compile with Boost filesystem v2. At some point,
 * the transition to v3 should be made, but as of now
 * many Linux distributions still ship with Boost 1.42,
 * which does not support version 3.
 */
#define BOOST_FILESYSTEM_VERSION 2

/* Use ELF support for thread-local storage on Linux? This
 * is potentially faster but causes problems when dynamically
 * loading Mitsuba from Python, so let's keep it disabled for now
 */
#define MTS_USE_ELF_TLS 0

#include <string>

MTS_NAMESPACE_BEGIN
#if defined(__OSX__)
extern void __mts_autorelease_init();
extern void __mts_autorelease_shutdown();
extern void __mts_autorelease_begin();
extern void __mts_autorelease_end();
extern std::string __mts_bundlepath();
extern void __mts_chdir_to_bundlepath();
extern void __mts_init_cocoa();
#define MTS_AUTORELEASE_BEGIN() __mts_autorelease_begin();
#define MTS_AUTORELEASE_END() __mts_autorelease_end();
#define MTS_AMBIGUOUS_SIZE_T 1
#else
#define MTS_AUTORELEASE_BEGIN() 
#define MTS_AUTORELEASE_END() 
#endif
MTS_NAMESPACE_END

#endif /* __PLATFORM_H */
