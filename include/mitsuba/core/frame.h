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
#if !defined(__MITSUBA_CORE_FRAME_H_)
#define __MITSUBA_CORE_FRAME_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Stores a three-dimensional orthonormal coordinate frame
 *
 * This class is mostly used to quickly convert between different
 * cartesian coordinate systems and to efficiently compute certain
 * quantities (e.g. \ref cosTheta(), \ref tanTheta, ..).
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct Frame {
    Vector s, t;
    Normal n;

    /// Default constructor -- performs no initialization!
    inline Frame() { }

    /// Given a normal and tangent vectors, construct a new coordinate frame
    inline Frame(const Vector &s, const Vector &t, const Normal &n)
     : s(s), t(t), n(n) {
    }

    /// Construct a frame from the given orthonormal vectors
    inline Frame(const Vector &x, const Vector &y, const Vector &z)
     : s(x), t(y), n(z) {
    }

    /// Construct a new coordinate frame from a single vector
    inline Frame(const Vector &n) : n(n) {
        coordinateSystem(n, s, t);
    }

    /// Unserialize from a binary data stream
    inline Frame(Stream *stream) {
        s = Vector(stream);
        t = Vector(stream);
        n = Normal(stream);
    }

    /// Serialize to a binary data stream
    inline void serialize(Stream *stream) const {
        s.serialize(stream);
        t.serialize(stream);
        n.serialize(stream);
    }

    /// Convert from world coordinates to local coordinates
    inline Vector toLocal(const Vector &v) const {
        return Vector(
            dot(v, s),
            dot(v, t),
            dot(v, n)
        );
    }

    /// Convert from local coordinates to world coordinates
    inline Vector toWorld(const Vector &v) const {
        return s * v.x + t * v.y + n * v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared cosine of the angle between the normal and v */
    inline static Float cosTheta2(const Vector &v) {
        return v.z * v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the cosine of the angle between the normal and v */
    inline static Float cosTheta(const Vector &v) {
        return v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the u and v coordinates of the vector 'v' */
    inline static Vector2 uv(const Vector &v) {
        return Vector2(v.x, v.y);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared sine of the angle between the normal and v */
    inline static Float sinTheta2(const Vector &v) {
        return 1.0f - v.z * v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the angle between the normal and v */
    inline static Float sinTheta(const Vector &v) {
        Float temp = sinTheta2(v);
        if (temp <= 0.0f)
            return 0.0f;
        return std::sqrt(temp);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the tangent of the angle between the normal and v */
    inline static Float tanTheta(const Vector &v) {
        Float temp = 1 - v.z*v.z;
        if (temp <= 0.0f)
            return 0.0f;
        return std::sqrt(temp) / v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared tangent of the angle between the normal and v */
    inline static Float tanTheta2(const Vector &v) {
        Float temp = 1 - v.z*v.z;
        if (temp <= 0.0f)
            return 0.0f;
        return temp / (v.z * v.z);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the phi parameter in spherical coordinates */
    inline static Float sinPhi(const Vector &v) {
        Float sinTheta = Frame::sinTheta(v);
        if (sinTheta == 0.0f)
            return 1.0f;
        return math::clamp(v.y / sinTheta, (Float) -1.0f, (Float) 1.0f);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the cosine of the phi parameter in spherical coordinates */
    inline static Float cosPhi(const Vector &v) {
        Float sinTheta = Frame::sinTheta(v);
        if (sinTheta == 0.0f)
            return 1.0f;
        return math::clamp(v.x / sinTheta, (Float) -1.0f, (Float) 1.0f);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared sine of the phi parameter in  spherical
     * coordinates */
    inline static Float sinPhi2(const Vector &v) {
        return math::clamp(v.y * v.y / sinTheta2(v), (Float) 0.0f, (Float) 1.0f);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared cosine of the phi parameter in  spherical
     * coordinates */
    inline static Float cosPhi2(const Vector &v) {
        return math::clamp(v.x * v.x / sinTheta2(v), (Float) 0.0f, (Float) 1.0f);
    }

    /// Equality test
    inline bool operator==(const Frame &frame) const {
        return frame.s == s && frame.t == t && frame.n == n;
    }

    /// Inequality test
    inline bool operator!=(const Frame &frame) const {
        return !operator==(frame);
    }

    /// Return a string representation of this frame
    inline std::string toString() const {
        std::ostringstream oss;
        oss << "Frame[" << std::endl
            << "  s = " << s.toString() << "," << std::endl
            << "  t = " << t.toString() << "," << std::endl
            << "  n = " << n.toString() << std::endl
            << "]";
        return oss.str();
    }
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_FRAME_H_ */
