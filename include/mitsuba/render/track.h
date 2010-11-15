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

#if !defined(__ANIMATION_TRACK_H)
#define __ANIMATION_TRACK_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

template <typename T> class AnimationTrack {
public:
	enum EType {
		EInvalid = 0,
		ELocationX, ELocationY, ELocationZ, ELocationXYZ,
		EScaleX, EScaleY, EScaleZ, EScaleXYZ,
		ERotationX, ERotationY, ERotationZ, ERotationXYZ
	};

	AnimationTrack(EType type, size_t nKeyframes) 
		: m_type(type), m_times(nKeyframes), m_values(nKeyframes) { }

	/// Return the type of this track
	inline EType getType() const { return m_type; }

	/// Return the number of keyframes
	inline size_t getSize() const { return m_times.size(); }

	/// Set the time value of a certain keyframe
	inline void setTime(size_t idx, Float time) { m_times[idx] = time; }

	/// Return the time value of a certain keyframe
	inline Float getTime(size_t idx) const { return m_times[idx]; }

	/// Set the value of a certain keyframe
	inline void setValue(size_t idx, const T &value) { m_values[idx] = value; }
	
	/// Return the value of a certain keyframe
	inline const T &getValue(size_t idx) const { return m_values[idx]; }

	/// Evaluate the animation track at an arbitrary time value
	inline T lookup(Float time) const {
	}
private:
	EType m_type;
	std::vector<Float> m_times;
	std::vector<T> m_values;
};

MTS_NAMESPACE_END

#endif /* __ANIMATION_TRACK_H */
