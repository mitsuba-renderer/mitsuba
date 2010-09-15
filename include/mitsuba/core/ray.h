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

#if !defined(__RAY_H)
#define __RAY_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple three-dimensional ray class with
   minimum / maximum extent information */
struct Ray {
	Point o;     ///< Ray origin
	Float mint;  ///< Minimum range for intersection tests
	Vector d;    ///< Ray direction
	Float maxt;  ///< Maximum range for intersection tests
	Vector dRcp; ///< Componentwise reciprocals of the ray direction

	/// Construct a new ray
	inline Ray() : mint(Epsilon), maxt(std::numeric_limits<Float>::infinity()) {
	}

	/// Copy constructor (1)
	inline Ray(const Ray &ray) 
	 : o(ray.o), mint(ray.mint), d(ray.d), maxt(ray.maxt), dRcp(ray.dRcp) {
	}

	/// Copy constructor (2)
	inline Ray(const Ray &ray, Float mint, Float maxt) 
	 : o(ray.o), mint(mint), d(ray.d), maxt(maxt), dRcp(ray.dRcp) {
	}
	
	/// Construct a new ray
	inline Ray(Point o, Vector _d)
		: o(o), mint(Epsilon),  d(_d), maxt(std::numeric_limits<Float>::infinity()) {
#ifdef MTS_DEBUG_FP
		disable_fpexcept();
#endif
		dRcp.x = (Float) 1.0f / _d.x;
		dRcp.y = (Float) 1.0f / _d.y;
		dRcp.z = (Float) 1.0f / _d.z;
#ifdef MTS_DEBUG_FP
		enable_fpexcept();
#endif
	}

	/// Construct a new ray
	inline Ray(Point o, Vector _d, Float mint, Float maxt)
		: o(o), mint(mint),  d(_d), maxt(maxt) {
#ifdef MTS_DEBUG_FP
		disable_fpexcept();
#endif
		dRcp.x = (Float) 1.0f / _d.x;
		dRcp.y = (Float) 1.0f / _d.y;
		dRcp.z = (Float) 1.0f / _d.z;
#ifdef MTS_DEBUG_FP
		enable_fpexcept();
#endif
	}

	/// Set the origin
	inline void setOrigin(const Point &oVal) { o = oVal; }
	
	/// Set the direction and update the reciprocal
	inline void setDirection(const Vector &dVal) {
		d = dVal;
#ifdef MTS_DEBUG_FP
		disable_fpexcept();
#endif
		dRcp.x = (Float) 1.0f / dVal.x;
		dRcp.y = (Float) 1.0f / dVal.y;
		dRcp.z = (Float) 1.0f / dVal.z;
#ifdef MTS_DEBUG_FP
		enable_fpexcept();
#endif
	}

	/// Return 3d coordinates of a point on the ray
	inline Point operator() (Float t) const { return o + t * d; }

	/// Return a string representation of this ray
	inline std::string toString() const {
		std::ostringstream oss;
		oss << "Ray[orig=" << o.toString() << ", dest=" << d.toString() << "]";
		return oss.str();
	}
};

/** \brief %Ray differential -- enhances the basic ray class with 
   information about the rays of adjacent pixels on the view plane */
struct RayDifferential : public Ray {
	bool hasDifferentials;
	Ray rx, ry;

	inline RayDifferential() 
		: hasDifferentials(false) {
	}

	inline RayDifferential(const Point &p, const Vector &d) 
		: Ray(p, d), hasDifferentials(false) {
	}

	inline explicit RayDifferential(const Ray &ray) 
		: Ray(ray), hasDifferentials(false) {
	}

	inline RayDifferential(const RayDifferential &ray) 
		: Ray(ray), hasDifferentials(ray.hasDifferentials), rx(ray.rx), ry(ray.ry) {
	}
	
	inline void operator=(const RayDifferential &ray) {
		o = ray.o;
		mint = ray.mint;
		d = ray.d;
		maxt = ray.maxt;
		dRcp = ray.dRcp;
		hasDifferentials = ray.hasDifferentials;
		rx = ray.rx;
		ry = ray.ry;
	}

	inline void operator=(const Ray &ray) {
		o = ray.o;
		mint = ray.mint;
		d = ray.d;
		maxt = ray.maxt;
		dRcp = ray.dRcp;
		hasDifferentials = false;
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
