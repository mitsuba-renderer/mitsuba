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

template <typename T> inline bool atomicCompareAndExchangePtr(T **v, T *newValue, T *oldValue) {
#if defined(WIN32)
    return InterlockedCompareExchangePointer(
		reinterpret_cast<volatile PVOID *>(v), newValue, oldValue) == oldValue;
#else
	return __sync_bool_compare_and_swap(v, oldValue, newValue);
#endif
}

MTS_NAMESPACE_END

#endif /* __ATOMIC_H */
