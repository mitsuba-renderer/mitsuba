/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
#if !defined(__MITSUBA_RENDER_TRACK_H_)
#define __MITSUBA_RENDER_TRACK_H_

#include <mitsuba/core/quat.h>
#include <mitsuba/core/simplecache.h>

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
		EInvalid        = 0,
		ETranslationX   = 1,
		ETranslationY   = 2,
		ETranslationZ   = 3,
		ETranslationXYZ = 4,
		EScaleX         = 5,
		EScaleY         = 6,
		EScaleZ         = 7,
		EScaleXYZ       = 8,
		ERotationX      = 9,
		ERotationY      = 10,
		ERotationZ      = 11,
		ERotationQuat   = 12
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
		: AbstractAnimationTrack(type, stream->readSize()) {
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
	inline void serialize(Stream *stream) const {
		stream->writeUInt(m_type);
		stream->writeSize(m_times.size());
		stream->writeFloatArray(&m_times[0], m_times.size());
		for (size_t i=0; i<m_values.size(); ++i)
			serialize(stream, m_values[i]);
	}

	/// Evaluate the animation track at an arbitrary time value
	inline value_type eval(Float time) const {
		SAssert(m_times.size() > 0);
		std::vector<Float>::const_iterator entry =
				std::lower_bound(m_times.begin(), m_times.end(), time);
		size_t idx0 = (size_t) std::max(
				(ptrdiff_t) (entry - m_times.begin()) - 1,
				(ptrdiff_t) 0);
		size_t idx1 = std::min(idx0+1, m_times.size()-1);
		Float t = 0.5f;
		if (m_times[idx0] != m_times[idx1]) {
			time = std::max(m_times[idx0], std::min(m_times[idx1], time));
			t = (time-m_times[idx0]) / (m_times[idx1]-m_times[idx0]);
		}
		return lerp(idx0, idx1, t);
	}
protected:
	/// Evaluate the animation track using linear interpolation
	inline value_type lerp(size_t idx0, size_t idx1, Float t) const;

	inline void unserialize(Stream *stream, value_type &value) {
		value = stream->readElement<value_type>();
	}

	inline void serialize(Stream *stream, const value_type &value) const {
		stream->writeElement<value_type>(value);
	}
private:
	std::vector<value_type> m_values;
};

template<typename T> inline T AnimationTrack<T>::lerp(size_t idx0, size_t idx1, Float t) const {
	return m_values[idx0] * (1-t) + m_values[idx1] * t;
}

/// Partial specialization for quaternions (uses \ref slerp())
template<> inline Quaternion AnimationTrack<Quaternion>::lerp(size_t idx0, size_t idx1, Float t) const {
	return slerp(m_values[idx0], m_values[idx1], t);
}

template<> inline void AnimationTrack<Point>::unserialize(Stream *stream, Point &value) {
	value = Point(stream);
}

template<> inline void AnimationTrack<Point>::serialize(Stream *stream, const Point &value) const {
	value.serialize(stream);
}

template<> inline void AnimationTrack<Vector>::unserialize(Stream *stream, Vector &value) {
	value = Vector(stream);
}

template<> inline void AnimationTrack<Vector>::serialize(Stream *stream, const Vector &value) const {
	value.serialize(stream);
}

template<> inline void AnimationTrack<Quaternion>::unserialize(Stream *stream, Quaternion &value) {
	value = Quaternion(stream);
}

template<> inline void AnimationTrack<Quaternion>::serialize(Stream *stream, const Quaternion &value) const {
	value.serialize(stream);
}

/**
 * \brief Animated transformation with an underlying keyframe representation
 * \ingroup librender
 */
class MTS_EXPORT_RENDER AnimatedTransform : public Object {
protected:
	/// Internal functor used by \ref eval() and \ref SimpleCache
	struct MTS_EXPORT_RENDER TransformFunctor {
	public:
		inline TransformFunctor(const std::vector<AbstractAnimationTrack *> &tracks)
			: m_tracks(tracks) {}

		void operator()(const Float &time, Transform &trafo) const;
	private:
		const std::vector<AbstractAnimationTrack *> &m_tracks;
	};
public:
	/**
	 * \brief Create a new animated transformation
	 *
	 * When the transformation is constant (i.e. there are no
	 * animation tracks), the supplied parameter specifies the
	 * target value.
	 */
	AnimatedTransform(const Transform &trafo = Transform())
		: m_transform(trafo) { }

	/// Unserialized an animated transformation from a binary data stream
	AnimatedTransform(Stream *stream);

	/// Return the number of associated animation tracks
	inline size_t getTrackCount() const { return m_tracks.size(); }

	/// Look up one of the tracks by index
	inline const AbstractAnimationTrack *getTrack(size_t idx) const { return m_tracks[idx]; }

	/// Append an animation track
	void addTrack(AbstractAnimationTrack *track);

	/**
	 * \brief Compute the transformation for the specified time value
	 *
	 * Note that the returned reference leads to a thread-local cache.
	 * This means that it will become invalidated at the next call
	 * to this function.
	 */
	inline const Transform &eval(Float t) const {
		if (m_tracks.size() == 0)
			return m_transform;
		else
			return m_cache.get(TransformFunctor(m_tracks), t);
	}

	/// Is the animation static?
	inline bool isStatic() const { return m_tracks.size() == 0; }

	/// Transform a point by an affine / non-projective matrix
	inline Point transformAffine(Float t, const Point &p) const {
		return eval(t).transformAffine(p);
	}

	/// Transform a point by an affine / non-projective matrix (no temporaries)
	inline void transformAffine(Float t, const Point &p, Point &dest) const {
		eval(t).transformAffine(p, dest);
	}

	/// Transform a ray by an affine / non-projective matrix
	inline Ray transformAffine(Float t, const Ray &r) const {
		return eval(t).transformAffine(r);
	}

	/// Transform a ray by an affine / non-projective matrix (no temporaries)
	inline void transformAffine(Float t, const Ray &r, Ray &dest) const {
		eval(t).transformAffine(r, dest);
	}

	/// Matrix-vector multiplication for points in 3d space
	inline Point operator()(Float t, const Point &p) const {
		return eval(t).transformAffine(p);
	}

	/// Matrix-vector multiplication for points in 3d space (no temporaries)
	inline void operator()(Float t, const Point &p, Point &dest) const {
		eval(t).operator()(p, dest);
	}

	/// Matrix-vector multiplication for vectors in 3d space
	inline Vector operator()(Float t, const Vector &v) const {
		return eval(t).operator()(v);
	}

	/// Matrix-vector multiplication for vectors in 3d space (no temporaries)
	inline void operator()(Float t, const Vector &v, Vector &dest) const {
		eval(t).operator()(v, dest);
	}

	/// Matrix-vector multiplication for normals in 3d space
	inline Normal operator()(Float t, const Normal &n) const {
		return eval(t).operator()(n);
	}

	/// Matrix-vector multiplication for normals in 3d space (no temporaries)
	inline void operator()(Float t, const Normal &n, Normal &dest) const {
		eval(t).operator()(n, dest);
	}

	/// \brief Transform a ray
	inline Ray operator()(Float t, const Ray &r) const {
		return eval(t).operator()(r);
	}

	/// Transform a ray (no temporaries)
	inline void operator()(Float t, const Ray &r, Ray &dest) const {
		eval(t).operator()(r, dest);
	}

	/// Serialize to a binary data stream
	void serialize(Stream *stream) const;

	/// Return the extents along the time axis
	AABB1 getTimeBounds() const;

	/// Return an axis-aligned box bounding the amount of translation
	AABB getTranslationBounds() const;

	/// Compute the spatial bounds of a transformed (static) AABB
	AABB getSpatialBounds(const AABB &aabb) const;

	/// Return a human-readable string description
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~AnimatedTransform();
private:
	std::vector<AbstractAnimationTrack *> m_tracks;
	mutable SimpleCache<Float, Transform> m_cache;
	Transform m_transform;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_TRACK_H_ */
