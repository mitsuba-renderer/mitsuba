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

#if !defined(__PLATFORM_H)
#define __PLATFORM_H

#if defined(_MSC_VER)
	/* Disable MSVC STL debug + security features (slow..!) */
	#ifdef _SECURE_SCL
		#undef _SECURE_SCL
	#endif
	#ifdef _SECURE_SCL_THROWS
		#undef _SECURE_SCL_THROWS
	#endif
	#ifdef _HAS_ITERATOR_DEBUGGING
		#undef _HAS_ITERATOR_DEBUGGING
	#endif
	#undef _STLP_DEBUG

	#define _SECURE_SCL 0
	#define _SECURE_SCL_THROWS 0
	#define _HAS_ITERATOR_DEBUGGING 0
	#define __MSVC__
	#define __WINDOWS__

	// Don't complain about perfectly fine ISO C++
	#define _CRT_SECURE_NO_WARNINGS
	#define _CRT_NONSTDC_NO_DEPRECATE
	#define _CRT_SECURE_NO_DEPRECATE

	#define _WIN32_WINNT 0x0501
	#define NOMINMAX
	#include <winsock2.h> // IPv6 support
	#include <windows.h>

	#pragma warning(disable : 4251) // 'field' : class 'A' needs to have dll-interface to be used by clients of class 'B'
	#pragma warning(disable : 4800) // 'type' : forcing value to bool 'true' or 'false' (performance warning)
	#include <stdint_msvc.h>        // Does not exist in MSVC. Use a replacement
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
#else
	#error Unknown OS
#endif

#ifdef __MSVC__
	#define MTS_MODULE_CORE 1
	#define MTS_MODULE_RENDER 2
	#define MTS_MODULE_HW 3

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

	#define SIZE_T_FMT "%Iu"
#else
	#define MTS_EXPORT
	#define MTS_EXPORT_CORE
	#define MTS_EXPORT_RENDER
	#define MTS_EXPORT_HW
	#include <stdint.h>

	#define SIZE_T_FMT "%zd"
#endif

#define MTS_NAMESPACE_BEGIN namespace mitsuba {
#define MTS_NAMESPACE_END }

#include <string>

MTS_NAMESPACE_BEGIN
#if defined(__OSX__)
extern void __ubi_autorelease_init();
extern void __ubi_autorelease_shutdown();
extern void __ubi_autorelease_begin();
extern void __ubi_autorelease_end();
extern std::string __ubi_bundlepath();
extern void __ubi_chdir_to_bundlepath();
#define MTS_AUTORELEASE_BEGIN() __ubi_autorelease_begin();
#define MTS_AUTORELEASE_END() __ubi_autorelease_end();
#else
#define MTS_AUTORELEASE_BEGIN() 
#define MTS_AUTORELEASE_END() 
#endif
MTS_NAMESPACE_END

#endif /* __PLATFORM_H */
