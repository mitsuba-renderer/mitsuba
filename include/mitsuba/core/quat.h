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
#if !defined(__MITSUBA_CORE_QUAT_H_)
#define __MITSUBA_CORE_QUAT_H_

#include <mitsuba/core/transform.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Parameterizable quaternion data structure
 * \ingroup libcore
 */
template <typename T> struct TQuaternion {
    typedef T          Scalar;

    /// Used by \ref TQuaternion::fromEulerAngles
    enum EEulerAngleConvention {
        EEulerXYZ = 0,
        EEulerXZY,
        EEulerYXZ,
        EEulerYZX,
        EEulerZXY,
        EEulerZYX
    };

    /// Imaginary component
    TVector3<T> v;

    /// Real component
    T w;

    /// Create a unit quaternion
    TQuaternion() : v(0.0f), w(1) { }

    /**
     * Initialize the quaternion with the specified
     * real and imaginary components
     */
    TQuaternion(const TVector3<T> &v, T w) : v(v), w(w) {  }

    /// Unserialize a quaternion from a binary data stream
    explicit TQuaternion(Stream *stream) {
        v = TVector3<T>(stream);
        w = stream->readElement<T>();
    }

    /// Add two quaternions and return the result
    TQuaternion operator+(const TQuaternion &q) const {
        return TQuaternion(v + q.v, w + q.w);
    }

    /// Subtract two quaternions and return the result
    TQuaternion operator-(const TQuaternion &q) const {
        return TQuaternion(v - q.v, w - q.w);
    }

    /// Add another quaternions to the current one
    TQuaternion& operator+=(const TQuaternion &q) {
        v += q.v; w += q.w;
        return *this;
    }

    /// Subtract a quaternion
    TQuaternion& operator-=(const TQuaternion &q) {
        v -= q.v; w -= q.w;
        return *this;
    }

    /// Unary negation operator
    TQuaternion operator-() const {
        return TQuaternion(-v, -w);
    }

    /// Multiply the quaternion by the given scalar and return the result
    TQuaternion operator*(T f) const {
        return TQuaternion(v*f, w*f);
    }

    /// Multiply the quaternion by the given scalar
    TQuaternion &operator*=(T f) {
        v *= f; w *= f;
        return *this;
    }

    /// Divide the quaternion by the given scalar and return the result
    TQuaternion operator/(T f) const {
#ifdef MTS_DEBUG
        if (f == 0)
            SLog(EWarn, "Quaternion: Division by zero!");
#endif
        T recip = (T) 1 / f;
        return TQuaternion(v * recip, w * recip);
    }

    /// Divide the quaternion by the given scalar
    TQuaternion &operator/=(T f) {
#ifdef MTS_DEBUG
        if (f == 0)
            SLog(EWarn, "Quaternion: Division by zero!");
#endif
        T recip = (T) 1 / f;
        v *= recip; w *= recip;
        return *this;
    }

    /// Quaternion multiplication
    TQuaternion &operator*=(const TQuaternion &q) {
        Scalar tmp = w * q.w - dot(v, q.v);
        v = cross(v, q.v) + q.w * v + w * q.v;
        w = tmp;
        return *this;
    }

    /// Quaternion multiplication (creates a temporary)
    TQuaternion operator*(const TQuaternion &q) const {
        return TQuaternion(cross(v, q.v) + q.w * v + w * q.v,
            w * q.w - dot(v, q.v));
    }

    /// Equality test
    bool operator==(const TQuaternion &q) const {
        return v == q.v && w == q.w;
    }

    /// Inequality test
    bool operator!=(const TQuaternion &q) const {
        return v != q.v || w != q.w;
    }

    /// Identity test
    bool isIdentity() const {
        return v.isZero() && w == 1;
    }

    /// Return the rotation axis of this quaternion
    inline TVector3<T> axis() const { return normalize(v); }

    /// Return the rotation angle of this quaternion (in radians)
    inline Scalar angle() const { return 2 * math::safe_acos(w); }

    /**
     * \brief Compute the exponential of a quaternion with
     * scalar part w = 0.
     *
     * Based on code the appendix of
     * "Quaternion Calculus for Computer Graphics" by Ken Shoemake
     */
    TQuaternion exp() const {
        T theta = v.length();
        T c = std::cos(theta);

        if (theta > Epsilon)
            return TQuaternion(v * (std::sin(theta) / theta), c);
        else
            return TQuaternion(v, c);
    }

    /**
     * \brief Compute the natural logarithm of a unit quaternion
     *
     * Based on code the appendix of
     * "Quaternion Calculus for Computer Graphics" by Ken Shoemake
     */
    TQuaternion log() const {
        T scale = v.length();
        T theta = std::atan2(scale, w);

        if (scale > 0)
            scale = theta/scale;

        return TQuaternion<T>(v * scale, 0.0f);
    }

    /**
     * \brief Construct an unit quaternion, which represents a rotation
     * around \a axis by \a angle radians.
     */
    static TQuaternion fromAxisAngle(const Vector &axis, Float angle) {
        T sinValue = std::sin(angle/2.0f), cosValue = std::cos(angle/2.0f);
        return TQuaternion(normalize(axis) * sinValue, cosValue);
    }

    /**
     * \brief Construct an unit quaternion, which rotates unit direction
     * \a from onto \a to.
     */
    static TQuaternion fromDirectionPair(const Vector &from, const Vector &to) {
        Float dp = dot(from, to);
        if (dp > 1-Epsilon) {
            // there is nothing to do
            return TQuaternion();
        } else if (dp < -(1-Epsilon)) {
            // Use a better-conditioned method for opposite directions
            Vector rotAxis = cross(from, Vector(1, 0, 0));
            Float length = rotAxis.length();
            if (length < Epsilon) {
                rotAxis = cross(from, Vector(0, 1, 0));
                length = rotAxis.length();
            }
            rotAxis /= length;
            return TQuaternion(rotAxis, 0);
        } else {
            // Find cos(theta) and sin(theta) using half-angle formulae
            Float cosTheta = std::sqrt(0.5f * (1 + dp));
            Float sinTheta = std::sqrt(0.5f * (1 - dp));
            Vector rotAxis = normalize(cross(from, to));
            return TQuaternion(rotAxis * sinTheta, cosTheta);
        }
    }

    inline static TQuaternion fromTransform(const Transform &trafo) {
        return fromMatrix(trafo.getMatrix());
    }

    /**
     * \brief Construct an unit quaternion matching the supplied
     * rotation matrix.
     */
    static TQuaternion fromMatrix(const Matrix4x4 &m) {
        // Implementation from PBRT, originally based on the matrix
        // and quaternion FAQ (http://www.j3d.org/matrix_faq/matrfaq_latest.html)
        T trace = m(0, 0) + m(1, 1) + m(2, 2);
        TVector3<T> v; T w;

        if (trace > Epsilon) {
            T s = std::sqrt(trace + 1.0f);
            w = s * 0.5f;
            s = 0.5f / s;
            v.x = (m(2, 1) - m(1, 2)) * s;
            v.y = (m(0, 2) - m(2, 0)) * s;
            v.z = (m(1, 0) - m(0, 1)) * s;
        } else {
            const int nxt[3] = {1, 2, 0};
            T q[3];
            int i = 0;
            if (m(1, 1) > m(0, 0)) i = 1;
            if (m(2, 2) > m(i, i)) i = 2;
            int j = nxt[i];
            int k = nxt[j];
            T s = std::sqrt((m(i, i) - (m(j, j) + m(k, k))) + 1.0f);
            q[i] = s * 0.5f;
            if (s != 0.f) s = 0.5f / s;
            w = (m(k, j) - m(j, k)) * s;
            q[j] = (m(j, i) + m(i, j)) * s;
            q[k] = (m(k, i) + m(i, k)) * s;
            v.x = q[0];
            v.y = q[1];
            v.z = q[2];
        }
        return TQuaternion(v, w);
    }

    /**
     * \brief Construct an unit quaternion matching the supplied
     * rotation expressed in Euler angles (in radians)
     */
    static TQuaternion fromEulerAngles(EEulerAngleConvention conv,
            Float x, Float y, Float z) {
        Quaternion qx = fromAxisAngle(Vector(1.0, 0.0, 0.0), x);
        Quaternion qy = fromAxisAngle(Vector(0.0, 1.0, 0.0), y);
        Quaternion qz = fromAxisAngle(Vector(0.0, 0.0, 1.0), z);

        switch (conv) {
            case EEulerXYZ:
                return qz * qy * qx;
            case EEulerXZY:
                return qy * qz * qx;
            case EEulerYXZ:
                return qz * qx * qy;
            case EEulerYZX:
                return qx * qz * qy;
            case EEulerZXY:
                return qy * qx * qz;
            case EEulerZYX:
                return qx * qy * qz;
            default:
                SLog(EError, "Internal error!");
                return TQuaternion();
        }
    }

    /// Compute the rotation matrix for the given quaternion
    Transform toTransform() const {
        /// Implementation from PBRT
        Float xx = v.x * v.x, yy = v.y * v.y, zz = v.z * v.z;
        Float xy = v.x * v.y, xz = v.x * v.z, yz = v.y * v.z;
        Float wx = v.x * w,   wy = v.y * w,   wz = v.z * w;

        Matrix4x4 m;
        m.m[0][0] = 1.f - 2.f * (yy + zz);
        m.m[0][1] =       2.f * (xy + wz);
        m.m[0][2] =       2.f * (xz - wy);
        m.m[0][3] = 0.0f;
        m.m[1][0] =       2.f * (xy - wz);
        m.m[1][1] = 1.f - 2.f * (xx + zz);
        m.m[1][2] =       2.f * (yz + wx);
        m.m[1][3] = 0.0f;
        m.m[2][0] =       2.f * (xz + wy);
        m.m[2][1] =       2.f * (yz - wx);
        m.m[2][2] = 1.f - 2.f * (xx + yy);
        m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f;
        m.m[3][1] = 0.0f;
        m.m[3][2] = 0.0f;
        m.m[3][3] = 1.0f;

        Matrix4x4 transp;
        m.transpose(transp);
        return Transform(transp, m);
    }

    /// Serialize this quaternion to a binary data stream
    void serialize(Stream *stream) const {
        v.serialize(stream);
        stream->writeElement<T>(w);
    }

    /// Return a readable string representation of this quaternion
    std::string toString() const {
        std::ostringstream oss;
        oss << "Quaternion[v=" << v.toString() << ", w=" << w << "]";
        return oss.str();
    }
};

template <typename T> inline TQuaternion<T> operator*(T f, const TQuaternion<T> &v) {
    return v*f;
}

template <typename T> inline T dot(const TQuaternion<T> &q1, const TQuaternion<T> &q2) {
    return dot(q1.v, q2.v) + q1.w * q2.w;
}

template <typename T> inline TQuaternion<T> normalize(const TQuaternion<T> &q) {
    return q / std::sqrt(dot(q, q));
}

template <typename T> inline TQuaternion<T> slerp(const TQuaternion<T> &q1,
    const TQuaternion<T> &_q2, Float t) {
    TQuaternion<T> q2(_q2);

    T cosTheta = dot(q1, q2);
    if (cosTheta < 0) {
        /* Take the short way! */
        q2 = -q2;
        cosTheta = -cosTheta;
    }
    if (cosTheta > .9995f) {
        // Revert to plain linear interpolation
        return normalize(q1 * (1.0f - t) +  q2 * t);
    } else {
        Float theta = math::safe_acos(math::clamp(cosTheta, (Float) -1.0f, (Float) 1.0f));
        Float thetap = theta * t;
        TQuaternion<T> qperp = normalize(q2 - q1 * cosTheta);
        return q1 * std::cos(thetap) + qperp * std::sin(thetap);
    }
}


MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_QUAT_H_ */
