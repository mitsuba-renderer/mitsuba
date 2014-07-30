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
#if !defined(__MITSUBA_CORE_TRACK_H_)
#define __MITSUBA_CORE_TRACK_H_

#include <mitsuba/core/quat.h>
#include <mitsuba/core/simplecache.h>
#include <set>

MTS_NAMESPACE_BEGIN

template <typename T> class AnimationTrack;

/**
 * \brief Base class of animation tracks
 * \ingroup librender
 */
class MTS_EXPORT_CORE AbstractAnimationTrack : public Object {
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

	/// Clone this track
	virtual AbstractAnimationTrack *clone() const = 0;

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
	typedef T ValueType;

	AnimationTrack(EType type, size_t nKeyframes = 0)
		: AbstractAnimationTrack(type, nKeyframes), m_values(nKeyframes) { }

	AnimationTrack(EType type, Stream *stream)
		: AbstractAnimationTrack(type, stream->readSize()) {
		m_values.resize(m_times.size());
		stream->readFloatArray(&m_times[0], m_times.size());
		for (size_t i=0; i<m_values.size(); ++i)
			unserialize(stream, m_values[i]);
	}

	/// Copy constructor
	AnimationTrack(const AnimationTrack *track)
		 : AbstractAnimationTrack(track->getType(), track->getSize()) {
		m_times = track->m_times;
		m_values = track->m_values;
	}

	/// Set the value of a certain keyframe
	inline void setValue(size_t idx, const ValueType &value) { m_values[idx] = value; }

	/// Return the value of a certain keyframe
	inline const ValueType &getValue(size_t idx) const { return m_values[idx]; }

	/// Reserve space for a certain number of entries
	inline void reserve(size_t count) { m_times.reserve(count); m_values.reserve(count); }

	/// Append a value
	inline void append(Float time, const ValueType &value) {
		m_times.push_back(time);
		m_values.push_back(value);
	}

	/// Clone this instance
	AbstractAnimationTrack *clone() const {
		return new AnimationTrack(this);
	}

	/// Prepend a transformation to every entry of this track
	void prependTransformation(const ValueType &value) {
		for (size_t i=0; i<m_values.size(); ++i)
			m_values[i] = concatenateTransformations(m_values[i], value);
	}

	/// Append a transformation to every entry of this track
	void appendTransformation(const ValueType &value) {
		for (size_t i=0; i<m_values.size(); ++i)
			m_values[i] = concatenateTransformations(value, m_values[i]);
	}

	/// Serialize to a binary data stream
	inline void serialize(Stream *stream) const {
		stream->writeUInt(m_type);
		stream->writeSize(m_times.size());
		stream->writeFloatArray(&m_times[0], m_times.size());
		for (size_t i=0; i<m_values.size(); ++i)
			serialize(stream, m_values[i]);
	}

	/// Evaluate the animation track at an arbitrary time value
	inline ValueType eval(Float time) const {
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

private:
	struct SortPredicate {
		inline bool operator()(const std::pair<Float, ValueType> &p1,
		                       const std::pair<Float, ValueType> &p2) const {
			return p1.first < p2.first;
		}
	};

	struct UniqueTimePredicate {
		inline bool operator()(const std::pair<Float, ValueType> &p1,
		                       const std::pair<Float, ValueType> &p2) const {
			return p1.first == p2.first;
		}
	};

public:
	/**
	 * \brief Sort all animation tracks and remove
	 * unnecessary data (for user-provided input)
	 *
	 * \return \c false if this animation track was deemed to be "trivial"
	 * after the cleanup (for instance, it only contains (0,0,0) translation operations)
	 */
	bool sortAndSimplify() {
		SAssert(m_values.size() == m_times.size());
		if (m_values.size() == 0)
			return false;

		std::vector< std::pair<Float, ValueType> > temp(m_values.size());
		for (size_t i=0; i<m_values.size(); ++i)
			temp[i] = std::make_pair(m_times[i], m_values[i]);
		std::sort(temp.begin(), temp.end(), SortPredicate());

		m_times.clear(); m_values.clear();
		m_times.push_back(temp[0].first);
		m_values.push_back(temp[0].second);

		for (size_t i=1; i<temp.size(); ++i) {
			Float time = temp[i].first;
			const ValueType &value = temp[i].second;

			if (m_times.back() == time)
				SLog(EError, "Duplicate time value in animated transformation!");

			/* Ignore irrelevant keys */
			if (i+1 < temp.size() && value == temp[i+1].second &&
					value == m_values.back())
				continue;
			else if (i+1 == temp.size() && value == m_values.back())
				continue;

			m_times.push_back(time);
			m_values.push_back(value);
		}

		return !(m_values.size() == 0 || (m_values.size() == 1 && isNoOp(m_values[0])));
	}
protected:
	/// Evaluate the animation track using linear interpolation
	inline ValueType lerp(size_t idx0, size_t idx1, Float t) const;

	/// Is this a "no-op" transformation?
	inline bool isNoOp(const ValueType &value) const;

	/// Concatenate two transformations
	inline ValueType concatenateTransformations(
			const ValueType &value1, const ValueType &value2) const;

	inline void unserialize(Stream *stream, ValueType &value) {
		value = stream->readElement<ValueType>();
	}

	inline void serialize(Stream *stream, const ValueType &value) const {
		stream->writeElement<ValueType>(value);
	}
private:
	std::vector<ValueType> m_values;
};

template<typename T> inline T AnimationTrack<T>::lerp(size_t idx0, size_t idx1, Float t) const {
	return m_values[idx0] * (1-t) + m_values[idx1] * t;
}

/// Partial specialization for quaternions (uses \ref slerp())
template<> inline Quaternion AnimationTrack<Quaternion>::lerp(size_t idx0, size_t idx1, Float t) const {
	return slerp(m_values[idx0], m_values[idx1], t);
}

template<typename T> inline T AnimationTrack<T>::concatenateTransformations(
		const T &value1, const T &value2) const {
	return value1 * value2;
}

template<> inline Vector AnimationTrack<Vector>::concatenateTransformations(
		const Vector &value1, const Vector &value2) const {
	if (m_type == ETranslationXYZ)
		return value1 + value2;
	else
		return Vector(value1.x * value2.x, value1.y * value2.y, value1.z * value2.z);
}

template<> inline Point AnimationTrack<Point>::concatenateTransformations(
		const Point &value1, const Point &value2) const {
	return value1 + value2;
}

template<> inline Float AnimationTrack<Float>::concatenateTransformations(
		const Float &value1, const Float &value2) const {
	if (m_type == ETranslationX || m_type == ETranslationY || m_type == ETranslationZ)
		return value1 + value2;
	else
		return value1 * value2;
}

template<typename T> inline bool AnimationTrack<T>::isNoOp(const ValueType &value) const {
	return false;
}

template<> inline bool AnimationTrack<Float>::isNoOp(const Float &value) const {
	if ((m_type == ETranslationX || m_type == ETranslationY || m_type == ETranslationZ) && value == 0)
		return true;
	else if ((m_type == ERotationX || m_type == ERotationY || m_type == ERotationZ) && value == 0)
		return true;
	else if ((m_type == EScaleX || m_type == EScaleY || m_type == EScaleZ) && value == 1)
		return true;
	return false;
}

template<> inline bool AnimationTrack<Vector>::isNoOp(const Vector &value) const {
	if (m_type == ETranslationXYZ && value.isZero())
		return true;
	else if (m_type == EScaleXYZ && (value.x == 1 && value.y == 1 && value.z == 1))
		return true;
	return false;
}

template<> inline bool AnimationTrack<Quaternion>::isNoOp(const Quaternion &value) const {
	return value.isIdentity();
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
class MTS_EXPORT_CORE AnimatedTransform : public Object {
private:
	/// Internal functor used by \ref eval() and \ref SimpleCache
	struct MTS_EXPORT_CORE TransformFunctor {
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

	/// Copy constructor
	AnimatedTransform(const AnimatedTransform *trafo);

	/// Return the number of associated animation tracks
	inline size_t getTrackCount() const { return m_tracks.size(); }

	/// Find a track of the given type
	AbstractAnimationTrack *findTrack(AbstractAnimationTrack::EType type);

	/// Find a track of the given type
	const AbstractAnimationTrack *findTrack(AbstractAnimationTrack::EType type) const;

	/// Look up one of the tracks by index
	inline AbstractAnimationTrack *getTrack(size_t idx) { return m_tracks[idx]; }

	/// Look up one of the tracks by index (const version)
	inline const AbstractAnimationTrack *getTrack(size_t idx) const { return m_tracks[idx]; }

	/// Return the used keyframes as a set
	void collectKeyframes(std::set<Float> &result) const;

	/// Append an animation track
	void addTrack(AbstractAnimationTrack *track);

	/**
	 * \brief Convenience function, which appends a linear transformation to the track
	 *
	 * Internally, a polar decomposition is used to split the transformation into scale,
	 * translation, and rotation, which are all separately interpolated.
	 *
	 * \remark Remember to run \ref sortAndSimplify() after adding all transformations.
	 */
	void appendTransform(Float time, const Transform &trafo);

	/**
	 * \brief Compute the transformation for the specified time value
	 *
	 * Note that the returned reference leads to a thread-local cache.
	 * This means that it will become invalidated at the next call
	 * to this function.
	 */
	inline const Transform &eval(Float t) const {
		if (EXPECT_TAKEN(m_tracks.size() == 0))
			return m_transform;
		else
			return m_cache.get(TransformFunctor(m_tracks), t);
	}

	/// Is the animation static?
	inline bool isStatic() const { return m_tracks.size() == 0; }

	/**
	 * \brief Sort all animation tracks and remove unnecessary
	 * data (for user-provided input)
	 */
	void sortAndSimplify();

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

	/// Prepend a scale transformation to the transform (this is often useful)
	void prependScale(const Vector &scale);

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

#endif /* __MITSUBA_CORE_TRACK_H_ */
