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

#if !defined(__RAY_H)
#define __RAY_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple three-dimensional ray data structure with
 * minimum / maximum extent information.
 *
 * The somewhat peculiar ordering of the attributes is due
 * to alignment purposes and should not be changed.
 *
 * \ingroup libcore
 * \ingroup libpython
*/
struct Ray {
	Point o;     ///< Ray origin
	Float mint;  ///< Minimum range for intersection tests
	Vector d;    ///< Ray direction
	Float maxt;  ///< Maximum range for intersection tests
	Vector dRcp; ///< Componentwise reciprocals of the ray direction
	Float time;  ///< Time value associated with this ray

	/// Construct a new ray
	inline Ray() : mint(Epsilon), 
		maxt(std::numeric_limits<Float>::infinity()), time(0.0f) {
	}

	/// Copy constructor (1)
	inline Ray(const Ray &ray) 
	 : o(ray.o), mint(ray.mint), d(ray.d), maxt(ray.maxt), 
	   dRcp(ray.dRcp), time(ray.time) {
	}

	/// Copy constructor (2)
	inline Ray(const Ray &ray, Float mint, Float maxt) 
	 : o(ray.o), mint(mint), d(ray.d), maxt(maxt), dRcp(ray.dRcp), time(ray.time) {
	}

	/// Construct a new ray, while not specifying a direction yet
	inline Ray(Point o, Float time) : o(o), mint(Epsilon), maxt(std::numeric_limits<Float>::infinity()), time(time) {
	}

	/// Construct a new ray
	inline Ray(Point o, Vector _d, Float time)
		: o(o), mint(Epsilon),  d(_d), maxt(std::numeric_limits<Float>::infinity()), time(time) {
#ifdef MTS_DEBUG_FP
		bool state = disableFPExceptions();
#endif
		dRcp.x = (Float) 1.0f / _d.x;
		dRcp.y = (Float) 1.0f / _d.y;
		dRcp.z = (Float) 1.0f / _d.z;
#ifdef MTS_DEBUG_FP
		restoreFPExceptions(state);
#endif
	}

	/// Construct a new ray
	inline Ray(Point o, Vector _d, Float mint, Float maxt, Float time)
		: o(o), mint(mint),  d(_d), maxt(maxt), time(time) {
#ifdef MTS_DEBUG_FP
		bool state = disableFPExceptions();
#endif
		dRcp.x = (Float) 1.0f / _d.x;
		dRcp.y = (Float) 1.0f / _d.y;
		dRcp.z = (Float) 1.0f / _d.z;
#ifdef MTS_DEBUG_FP
		restoreFPExceptions(state);
#endif
	}

	/// Set the origin
	inline void setOrigin(const Point &oVal) { o = oVal; }

	/// Set the origin
	inline void setTime(const Float &tval) { time = tval; }
	
	/// Set the direction and update the reciprocal
	inline void setDirection(const Vector &dVal) {
		d = dVal;
#ifdef MTS_DEBUG_FP
		bool state = disableFPExceptions();
#endif
		dRcp.x = (Float) 1.0f / dVal.x;
		dRcp.y = (Float) 1.0f / dVal.y;
		dRcp.z = (Float) 1.0f / dVal.z;
#ifdef MTS_DEBUG_FP
		restoreFPExceptions(state);
#endif
	}

	/**
	 * \brief Return 3D coordinates of a point along the ray
	 *
	 * \remark In the Python bindings, this operator is 
	 * exposed as a function named \c eval -- i.e. 
	 * position lookups should be written as \c ray.eval(t)
	 */
	inline Point operator() (Float t) const { return o + t * d; }

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
			<< "  rxOrigin = " << rxOrigin.toString() << "," << endl
			<< "  ryOrigin = " << ryOrigin.toString() << "," << endl
			<< "  rxDirection = " << rxDirection.toString() << "," << endl
			<< "  ryDirection = " << ryDirection.toString() << endl
			<< "]" << endl;
		return oss.str();
	}
};

#if defined(MTS_SSE)
/** \brief SIMD quad-packed ray for coherent ray tracing */
struct RayPacket4 {
	QuadVector o, d;
	QuadVector dRcp;
	uint8_t signs[4][4];

	inline RayPacket4() {
	}

	inline bool load(const Ray *rays) {
		for (int i=0; i<4; i++) {
			for (int axis=0; axis<3; axis++) {
				o[axis].f[i] = rays[i].o[axis];
				d[axis].f[i] = rays[i].d[axis];
				dRcp[axis].f[i] = rays[i].dRcp[axis];
				signs[axis][i] = rays[i].d[axis] < 0 ? 1 : 0;
				if (signs[axis][i] != signs[axis][0])
					return false;
			}
		}
		return true;
	}
};

struct RayInterval4 {
	SSEVector mint;
	SSEVector maxt;

	inline RayInterval4() {
		mint = SSEConstants::eps;
		maxt = SSEConstants::p_inf;
	}

	inline RayInterval4(const Ray *rays) {
		for (int i=0; i<4; i++) {
			mint.f[i] = rays[i].mint;
			maxt.f[i] = rays[i].maxt;
		}
	}
};

struct Intersection4 {
	SSEVector t;
	SSEVector u;
	SSEVector v;
	SSEVector primIndex;
	SSEVector shapeIndex;

	inline Intersection4() {
		t          = SSEConstants::p_inf;
		u          = SSEConstants::zero;
		v          = SSEConstants::zero;
		primIndex  = SSEConstants::ffffffff;
		shapeIndex = SSEConstants::ffffffff;
	}
};

#endif

MTS_NAMESPACE_END

#endif /* __RAY_H */
