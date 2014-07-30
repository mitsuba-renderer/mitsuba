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
#if !defined(__MITSUBA_CORE_NORMAL_H_)
#define __MITSUBA_CORE_NORMAL_H_

#include <mitsuba/core/vector.h>

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/normal.h mitsuba/mitsuba.h
 * \brief Three-dimensional normal data structure
 *
 * Internally represented using floating point numbers of the chosen
 * compile-time precision. The main difference of this data structure
 * when compared to \ref TVector3<Float> is in how instances of
 * \ref Normal are treated by linear transformations.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct Normal : public TVector3<Float> {
	/** \brief Construct a new normal without initializing it.
	 *
	 * This construtor is useful when the normal will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
	Normal() { }

	/// Initialize the vector with the specified X and Z components
	Normal(Float x, Float y, Float z) : TVector3<Float>(x, y, z) { }

	/// Initialize all components of the the normal with the specified value
	explicit Normal(Float val) : TVector3<Float>(val) { }

	/// Unserialize a normal from a binary data stream
	Normal(Stream *stream) {
		x = stream->readElement<Float>();
		y = stream->readElement<Float>();
		z = stream->readElement<Float>();
	}

	/// Construct a normal from a vector data structure
	Normal(const TVector3<Float> &v) : TVector3<Float>(v.x, v.y, v.z) { }

	/// Assign a vector to this normal
	void operator=(const TVector3<Float> &v) {
		x = v.x; y = v.y; z = v.z;
	}
};

inline Normal normalize(const Normal &n) {
	return n / n.length();
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_NORMAL_H_ */
