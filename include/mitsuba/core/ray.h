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
#if !defined(__MITSUBA_CORE_RAY_H_)
#define __MITSUBA_CORE_RAY_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple n-dimensional ray data structure with
 * minimum / maximum extent information.
 *
 * The somewhat peculiar ordering of the attributes is due
 * to alignment purposes and should not be changed.
 *
 * \ingroup libcore
 * \ingroup libpython
*/
template <typename _PointType, typename _VectorType> struct TRay {
    typedef _PointType                  PointType;
    typedef _VectorType                 VectorType;
    typedef typename PointType::Scalar  Scalar;

    /* The somewhat peculiar ordering of the attributes is for
       alignment purposes in the 3D case and should not be changed. */

    PointType o;     ///< Ray origin
    Scalar mint;     ///< Minimum range for intersection tests
    VectorType d;    ///< Ray direction
    Scalar maxt;     ///< Maximum range for intersection tests
    VectorType dRcp; ///< Componentwise reciprocals of the ray direction
    Float time;  ///< Time value associated with this ray

    /// Construct a new ray
    inline TRay() : mint(Epsilon),
        maxt(std::numeric_limits<Scalar>::infinity()), time(0) {
    }

    /// Copy constructor (1)
    inline TRay(const TRay &ray)
     : o(ray.o), mint(ray.mint), d(ray.d), maxt(ray.maxt),
       dRcp(ray.dRcp), time(ray.time) {
    }

    /// Copy constructor (2)
    inline TRay(const TRay &ray, Scalar mint, Scalar maxt)
     : o(ray.o), mint(mint), d(ray.d), maxt(maxt),
       dRcp(ray.dRcp), time(ray.time) { }

    /// Construct a new ray, while not specifying a direction yet
    inline TRay(const PointType &o, Scalar time) : o(o), mint(Epsilon),
      maxt(std::numeric_limits<Scalar>::infinity()), time(time) { }

    /// Construct a new ray
    inline TRay(const PointType &o, const VectorType &d, Scalar time)
        : o(o), mint(Epsilon),  d(d),
          maxt(std::numeric_limits<Scalar>::infinity()), time(time) {
#ifdef MTS_DEBUG_FP
        bool state = disableFPExceptions();
#endif
        for (int i=0; i<3; ++i)
            dRcp[i] = (Scalar) 1 / d[i];
#ifdef MTS_DEBUG_FP
        restoreFPExceptions(state);
#endif
    }

    /// Construct a new ray
    inline TRay(const PointType &o, const VectorType &d, Scalar mint, Scalar maxt,
        Scalar time) : o(o), mint(mint),  d(d), maxt(maxt), time(time) {
#ifdef MTS_DEBUG_FP
        bool state = disableFPExceptions();
#endif
        for (int i=0; i<3; ++i)
            dRcp[i] = (Scalar) 1 / d[i];
#ifdef MTS_DEBUG_FP
        restoreFPExceptions(state);
#endif
    }

    /// Set the origin
    inline void setOrigin(const PointType &pos) { o = pos; }

    /// Set the time
    inline void setTime(Scalar tval) { time = tval; }

    /// Set the direction and update the reciprocal
    inline void setDirection(const VectorType &dir) {
        d = dir;
#ifdef MTS_DEBUG_FP
        bool state = disableFPExceptions();
#endif
        for (int i=0; i<3; ++i)
            dRcp[i] = (Scalar) 1 / dir[i];
#ifdef MTS_DEBUG_FP
        restoreFPExceptions(state);
#endif
    }

    /**
     * \brief Return the position of a point along the ray
     *
     * \remark In the Python bindings, this operator is
     * exposed as a function named \c eval -- i.e.
     * position lookups should be written as \c ray.eval(t)
     */
    inline PointType operator() (Scalar t) const { return o + t * d; }

    /// Return a string representation of this ray
    inline std::string toString() const {
        std::ostringstream oss;
        oss << "Ray[origin=" << o.toString() << ", direction="
            << d.toString() << ", mint=" << mint
            << ", maxt=" << maxt << ", time=" << time << "]";
        return oss.str();
    }
};

/** \brief %Ray differential -- enhances the basic ray class with
   information about the rays of adjacent pixels on the view plane
   \ingroup libcore
*/
struct RayDifferential : public Ray {
    Point rxOrigin, ryOrigin;
    Vector rxDirection, ryDirection;
    bool hasDifferentials;

    inline RayDifferential()
        : hasDifferentials(false) {
    }

    inline RayDifferential(const Point &p, const Vector &d, Float time)
        : Ray(p, d, time), hasDifferentials(false) {
    }

    inline explicit RayDifferential(const Ray &ray)
        : Ray(ray), hasDifferentials(false) {
    }

    inline RayDifferential(const RayDifferential &ray)
        : Ray(ray), rxOrigin(ray.rxOrigin), ryOrigin(ray.ryOrigin),
          rxDirection(ray.rxDirection), ryDirection(ray.ryDirection),
          hasDifferentials(ray.hasDifferentials) {
    }

    void scaleDifferential(Float amount) {
        rxOrigin = o + (rxOrigin - o) * amount;
        ryOrigin = o + (ryOrigin - o) * amount;
        rxDirection = d + (rxDirection - d) * amount;
        ryDirection = d + (ryDirection - d) * amount;
    }

    void scaleDifferentialUV(const Vector2 amountUV) {
        rxOrigin = o + (rxOrigin - o) * amountUV.x;
        ryOrigin = o + (ryOrigin - o) * amountUV.y;
        rxDirection = d + (rxDirection - d) * amountUV.x;
        ryDirection = d + (ryDirection - d) * amountUV.y;
    }

    inline void operator=(const RayDifferential &ray) {
        o = ray.o;
        mint = ray.mint;
        d = ray.d;
        maxt = ray.maxt;
#ifdef MTS_DEBUG_FP
        bool state = disableFPExceptions();
#endif
        dRcp = ray.dRcp;
#ifdef MTS_DEBUG_FP
        restoreFPExceptions(state);
#endif
        hasDifferentials = ray.hasDifferentials;
        rxOrigin = ray.rxOrigin;
        ryOrigin = ray.ryOrigin;
        rxDirection = ray.rxDirection;
        ryDirection = ray.ryDirection;
    }

    inline void operator=(const Ray &ray) {
        o = ray.o;
        mint = ray.mint;
        d = ray.d;
        maxt = ray.maxt;
#ifdef MTS_DEBUG_FP
        bool state = disableFPExceptions();
#endif
        dRcp = ray.dRcp;
#ifdef MTS_DEBUG_FP
        restoreFPExceptions(state);
#endif
        hasDifferentials = false;
    }

    /// Return a string representation of this ray
    inline std::string toString() const {
        std::ostringstream oss;
        oss << "RayDifferential[" << endl
            << "  origin = " << o.toString() << "," << endl
            << "  direction  = " << d.toString() << "," << endl
            << "  mint = " << mint << "," << endl
            << "  maxt = " << maxt << "," << endl
            << "  time = " << time << "," << endl
            << "  hasDifferentials = " << hasDifferentials << "," << endl
            << "  rxOrigin = " << rxOrigin.toString() << "," << endl
            << "  ryOrigin = " << ryOrigin.toString() << "," << endl
            << "  rxDirection = " << rxDirection.toString() << "," << endl
            << "  ryDirection = " << ryDirection.toString() << endl
            << "]" << endl;
        return oss.str();
    }
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_RAY_H_ */
