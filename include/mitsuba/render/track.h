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

#include <mitsuba/core/quat.h>

MTS_NAMESPACE_BEGIN
	
template <typename T> class AnimationTrack;

/// Base class of animation tracks
class MTS_EXPORT_RENDER AbstractAnimationTrack : public Object {
	template<typename T> friend class AnimationTrack;
public:
	enum EType {
		EInvalid = 0,
		ELocationX, ELocationY, ELocationZ, ELocationXYZ, 
		EScaleX, EScaleY, EScaleZ, EScaleXYZ,
		ERotationX, ERotationY, ERotationZ,
		ERotationQuat
	};

	/// Return the type of this track
	inline EType getType() const { return m_type; }

	/// Set the time value of a certain keyframe
	inline void setTime(size_t idx, Float time) { m_times[idx] = time; }

	/// Return the time value of a certain keyframe
	inline Float getTime(size_t idx) const { return m_times[idx]; }
	
	/// Return the number of keyframes
	inline size_t getSize() const { return m_times.size(); }

	MTS_DECLARE_CLASS()
protected:
	AbstractAnimationTrack(EType type, size_t nKeyframes) 
		: m_type(type), m_times(nKeyframes) { }

protected:
	EType m_type;
	std::vector<Float> m_times;
};

/// Parameterizable animation track
template <typename T> class AnimationTrack : public AbstractAnimationTrack {
public:
	typedef T value_type;

	AnimationTrack(EType type, size_t nKeyframes) 
		: AbstractAnimationTrack(type, nKeyframes), m_values(nKeyframes) { }

	/// Set the value of a certain keyframe
	inline void setValue(size_t idx, const value_type &value) { m_values[idx] = value; }
	
	/// Return the value of a certain keyframe
	inline const value_type &getValue(size_t idx) const { return m_values[idx]; }

	/// Evaluate the animation track using linear interpolation
	value_type lerp(int idx0, int idx1, Float t) const;

	/// Evaluate the animation track at an arbitrary time value
	inline value_type lookup(Float time) const {
		std::vector<Float>::const_iterator entry = 
				std::lower_bound(m_times.begin(), m_times.end(), time);
		int idx0 = (int) (entry - m_times.begin()) - 1;
		int idx1 = idx0 + 1;
		idx0 = std::max(idx0, 0); idx1 = std::min(idx1, (int) m_times.size() - 1);
		Float t = 0.5f;
		if (m_times[idx0] != m_times[idx1])
			t = (time-m_times[idx0]) / (m_times[idx1]-m_times[idx0]);
		return lerp(idx0, idx1, t);
	}
private:
	std::vector<value_type> m_values;
};
	
template<typename T> T AnimationTrack<T>::lerp(int idx0, int idx1, Float t) const {
	return m_values[idx0] * (1-t) + m_values[idx1] * t;
}

/// Partial specialization for quaternions (uses \ref slerp())
template<> Quaternion AnimationTrack<Quaternion>::lerp(int idx0, int idx1, Float t) const {
	return slerp(m_values[idx0], m_values[idx1], t);
}

MTS_NAMESPACE_END

#endif /* __ANIMATION_TRACK_H */
