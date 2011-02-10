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

#if !defined(__ATOMIC_H)
#define __ATOMIC_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

#if defined(__OSX__)
#include <libkern/OSAtomic.h>
#endif

/* The following implementations are based on PBRT */

template <typename T> inline bool atomicCompareAndExchangePtr(T **v, T *newValue, T *oldValue) {
#if defined(WIN32)
    return InterlockedCompareExchangePointer(
		reinterpret_cast<volatile PVOID *>(v), newValue, oldValue) == oldValue;
#else
  #if !defined(__clang__)
	return __sync_bool_compare_and_swap(v, oldValue, newValue);
  #else
  #if __SIZEOF_POINTER__ == 8 
	return __sync_bool_compare_and_swap(
		reinterpret_cast<long long volatile *>(v), reinterpret_cast<long long>(oldValue), 
		reinterpret_cast<long long>(newValue));
  #else
	return __sync_bool_compare_and_swap(
		reinterpret_cast<long volatile *>(v), reinterpret_cast<long>(oldValue), 
		reinterpret_cast<long>(newValue));
  #endif
  #endif
#endif
}

inline bool atomicCompareAndExchange(volatile int32_t *v, int32_t newValue, int32_t oldValue) {
#if defined(WIN32)
    return InterlockedCompareExchange(
		reint32_terpret_cast<volatile int32_t *>(v), newValue, oldValue) == oldValue;
#else
	return __sync_bool_compare_and_swap(v, oldValue, newValue);
#endif
}

inline float atomicAdd(volatile float *val, float delta) {
	/* Atomic FP addition from PBRT */
	union bits { float f; int32_t i; };
	bits oldVal, newVal;
	do {
        // On IA32/x64, adding a PAUSE instruction in compare/exchange loops
        // is recommended to improve performance.  (And it does!)
#if (defined(__i386__) || defined(__amd64__))
        __asm__ __volatile__ ("pause\n");
#endif
		oldVal.f = *val;
		newVal.f = oldVal.f + delta;
	} while (atomicCompareAndExchange((volatile int32_t *) val, newVal.i, oldVal.i) != oldVal.i);
	return newVal.f;
}

MTS_NAMESPACE_END

#endif /* __ATOMIC_H */
