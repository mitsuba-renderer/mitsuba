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

#if !defined(__ANIMATION_TRACK_H)
#define __ANIMATION_TRACK_H

#include <mitsuba/core/quat.h>

MTS_NAMESPACE_BEGIN
	
template <typename T> class AnimationTrack;

/**
 * \brief Base class of animation tracks
 * \ingroup librender
 */
class MTS_EXPORT_RENDER AbstractAnimationTrack : public Object {
	template<typename T> friend class AnimationTrack;
public:
	enum EType {
		EInvalid = 0,
		ETranslationX, ETranslationY, ETranslationZ, ETranslationXYZ, 
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

	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream) const = 0;

	MTS_DECLARE_CLASS()
protected:
	AbstractAnimationTrack(EType type, size_t nKeyframes) 
		: m_type(type), m_times(nKeyframes) { }

	virtual ~AbstractAnimationTrack() { }
protected:
	EType m_type;
	std::vector<Float> m_times;
};

/**
 * \brief Parameterizable animation track
 * \ingroup librender
 */
template <typename T> class AnimationTrack : public AbstractAnimationTrack {
public:
	typedef T value_type;

	AnimationTrack(EType type, size_t nKeyframes) 
		: AbstractAnimationTrack(type, nKeyframes), m_values(nKeyframes) { }

	AnimationTrack(EType type, Stream *stream) 
		: AbstractAnimationTrack(type, (size_t) stream->readUInt()) {
		m_values.resize(m_times.size());
		stream->readFloatArray(&m_times[0], m_times.size());
		for (size_t i=0; i<m_values.size(); ++i)
			unserialize(stream, m_values[i]);
	}

	/// Set the value of a certain keyframe
	inline void setValue(size_t idx, const value_type &value) { m_values[idx] = value; }
	
	/// Return the value of a certain keyframe
	inline const value_type &getValue(size_t idx) const { return m_values[idx]; }

	/// Serialize to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeUInt(m_type);
		stream->writeUInt((uint32_t) m_times.size());
		stream->writeFloatArray(&m_times[0], m_times.size());
		for (size_t i=0; i<m_values.size(); ++i)
			serialize(stream, m_values[i]);
	}
			
	/// Evaluate the animation track at an arbitrary time value
	inline value_type eval(Float time) const {
		SAssert(m_times.size() > 0);
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
protected:
	/// Evaluate the animation track using linear interpolation
	value_type lerp(int idx0, int idx1, Float t) const;

	void unserialize(Stream *stream, value_type &value) {
		value = stream->readElement<value_type>();
	}

	void serialize(Stream *stream, const value_type &value) const {
		stream->writeElement<value_type>(value);
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

template<> void AnimationTrack<Point>::unserialize(Stream *stream, Point &value) {
	value = Point(stream);
}

template<> void AnimationTrack<Point>::serialize(Stream *stream, const Point &value) const {
	value.serialize(stream);
}

template<> void AnimationTrack<Vector>::unserialize(Stream *stream, Vector &value) {
	value = Vector(stream);
}

template<> void AnimationTrack<Vector>::serialize(Stream *stream, const Vector &value) const {
	value.serialize(stream);
}

template<> void AnimationTrack<Quaternion>::unserialize(Stream *stream, Quaternion &value) {
	value = Quaternion(stream);
}

template<> void AnimationTrack<Quaternion>::serialize(Stream *stream, const Quaternion &value) const {
	value.serialize(stream);
}

/**
 * \brief Animated transformation with an underlying keyframe representation
 * \ingroup librender
 */
class MTS_EXPORT_RENDER AnimatedTransform : public Object {
public:
	/// Create a new animated transform
	AnimatedTransform() { }
	
	/// Unseraizlie a animated transform
	AnimatedTransform(Stream *stream);

	/// Return the number of associated animation tracks
	inline size_t getTrackCount() const { return m_tracks.size(); }

	/// Look up one of the tracks by index
	inline const AbstractAnimationTrack *getTrack(size_t idx) const { return m_tracks[idx]; }

	/// Append an animation track
	void addTrack(AbstractAnimationTrack *track);

	/// Compute the transformation at the specified time value
	void eval(Float t, Transform &trafo) const;

	/// Serialize to a binary data stream
	void serialize(Stream *stream) const;

	/// Return the extents along the time axis
	void computeTimeBounds(Float &min, Float &max) const {
		min = std::numeric_limits<Float>::infinity();
		max = -std::numeric_limits<Float>::infinity();

		for (size_t i=0; i<m_tracks.size(); ++i) {
			AbstractAnimationTrack *track = m_tracks[i];
			size_t size = track->getSize();
			SAssert(size > 0);
			min = std::min(min, track->getTime(0));
			max = std::max(max, track->getTime(size-1));
		}
	}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~AnimatedTransform();
private:
	std::vector<AbstractAnimationTrack *> m_tracks;
};

MTS_NAMESPACE_END

#endif /* __ANIMATION_TRACK_H */
