#include <mitsuba/render/shape.h>
#include "miterseg.h"

MTS_NAMESPACE_BEGIN

MiterHairSegment::MiterHairSegment(ref<Hair> parent, int iv) : Shape(Properties()) {
	m_hair = parent;
	m_iv = iv;
	configure();
}


void MiterHairSegment::configure() {
	Shape::configure();

	m_aabb.reset(); m_bsphere.radius = 0;

	const Point start = firstVertex();
	const Point end = secondVertex();
	const Vector segment = end - start;

	// Caution this is false for collapsing segments (where the intersection of the miter planes intersects the cylinder)
	m_surfaceArea = 2*M_PI * m_hair->radius() * segment.length();
	m_invSurfaceArea = 1.0f / m_surfaceArea;

	m_aabb = getWorldAABB(0, 1);
	m_bsphere.center = m_aabb.getCenter();
	for (int i=0; i<8; ++i)
		m_bsphere.expandBy(m_aabb.getCorner(i));
}


AABB MiterHairSegment::getWorldAABB(Float t0, Float t1) const {
	AABB result;

	// The bounding box is conservatively the bbox of two spheres at the
	// endpoints of the segment.  Each sphere's radius is the hair
	// radius divided by the cosine of the steepest miter angle, making
	// it a bounding sphere for the ellipsoidal boundary for that end
	// of the segment.
	// Side note: There is a possible problem here, that the miter angle
	// can get arbitrarily steep if two segments form a very sharp angle.
	// This may need to be addressed at some point.

	const Point start = firstVertex();
	const Point end = secondVertex();
	const Vector segment = end - start;

	const Point a = start + t0 * segment;
	const Point b = start + t1 * segment;

	// cosine of steepest miter angle
	const Float cos0 = dot(firstMiterNormal(), tangent());
	const Float cos1 = dot(secondMiterNormal(), tangent());
	const Float maxInvCos = 1.0 / std::min(cos0, cos1);

	const Float expandRadius = m_hair->radius() * maxInvCos;
	const Vector expandVec(expandRadius, expandRadius, expandRadius);

	result.expandBy(a - expandVec);
	result.expandBy(a + expandVec);
	result.expandBy(b - expandVec);
	result.expandBy(b + expandVec);

	return result;
}

AABB MiterHairSegment::getClippedAABB(const AABB &aabb) const {
    AABB result(m_aabb);
	result.clip(aabb);
	return result;
	/* The following is broken, I believe...
	Float nearT, farT;
	AABB result(m_aabb);
	result.clip(aabb);

	Point a = firstVertex();
	Point b = secondVertex();

	if (!result.rayIntersect(Ray(a, normalize(b-a)), nearT, farT))
		return result; // that could be improved

	nearT = std::max(nearT, (Float) 0);
	farT = std::min(farT, m_length);
	result = getWorldAABB(nearT, farT);
	result.clip(aabb);

	return result;*/
}

bool MiterHairSegment::rayIntersect(const Ray &ray, Float start, Float end, Float &t) const {
	Float nearT, farT;

	/* First compute the intersection with the infinite cylinder */

	// Projection of ray onto subspace normal to axis
	Vector axis = tangent();
	Vector relOrigin = ray.o - firstVertex();
	Vector projOrigin = relOrigin - dot(axis, relOrigin) * axis;
	Vector projDirection = ray.d - dot(axis, ray.d) * axis;

	// Quadratic to intersect circle in projection
	const Float A = projDirection.lengthSquared();
	const Float B = 2 * dot(projOrigin, projDirection);
	const Float radius = m_hair->radius();
	const Float C = projOrigin.lengthSquared() - radius*radius;

	if (!solveQuadratic(A, B, C, nearT, farT))
		return false;

	if (nearT > end || farT < start)
		return false;

	/* Next check the intersection points against the miter planes */

	Point pointNear = ray(nearT);
	Point pointFar = ray(farT);
	if (dot(pointNear - firstVertex(), firstMiterNormal()) >= 0 &&
	    dot(pointNear - secondVertex(), secondMiterNormal()) <= 0 &&
	    nearT >= start) {
		t = nearT;
	} else if (dot(pointFar - firstVertex(), firstMiterNormal()) >= 0 &&
		       dot(pointFar - secondVertex(), secondMiterNormal()) <= 0) {
		if (farT > end)
			return false;
		t = farT;
	} else {
		return false;
	}

	return true;
}

bool MiterHairSegment::rayIntersect(const Ray &ray, Intersection &its) const {
	if (!rayIntersect(ray, ray.mint, ray.maxt, its.t))
		return false;
	its.p = ray(its.t);


	/* For now I don't compute texture coordinates at all. */
	its.uv = Point2(0,0);
	its.dpdu = Vector(0,0,0);
	its.dpdv = Vector(0,0,0);

	its.geoFrame.s = tangent();
	const Vector relHitPoint = its.p - firstVertex();
	const Vector axis = tangent();
	its.geoFrame.n = Normal(relHitPoint - dot(axis, relHitPoint) * axis);
	its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
	its.shFrame = its.geoFrame;
	its.wi = its.toLocal(-ray.d);
	its.hasUVPartials = false;
	its.shape = this;

	/* Intersection refinement step */
	// Do I need this?
	/*
	Vector2 localDir(normalize(Vector2(local.x, local.y)));
	Vector rel = its.p - m_objectToWorld(Point(m_radius  * localDir.x,
			m_radius * localDir.y, local.z));
	Float correction = -dot(rel, its.geoFrame.n)/dot(ray.d, its.geoFrame.n);

	its.t += correction;
	if (its.t < ray.mint || its.t > ray.maxt) {
		its.t = std::numeric_limits<Float>::infinity();
		return false;
	}

	its.p += ray.d * correction;
	*/

	return true;
}

#if defined(MTS_SSE)
/* SSE-accelerated packet tracing is not supported for cylinders at the moment */
__m128 MiterHairSegment::rayIntersectPacket(const RayPacket4 &packet, const
		__m128 start, __m128 end, __m128 inactive, Intersection4 &its) const {
	SSEVector result(_mm_setzero_ps()), mint(start), maxt(end), mask(inactive);
	Float t;

	for (int i=0; i<4; i++) {
		Ray ray;
		for (int axis=0; axis<3; axis++) {
			ray.o[axis] = packet.o[axis].f[i];
			ray.d[axis] = packet.d[axis].f[i];
		}
		if (mask.i[i] != 0)
			continue;
		if (rayIntersect(ray, mint.f[i], maxt.f[i], t)) {
			result.i[i] = 0xFFFFFFFF;
			its.t.f[i] = t;
		}
	}
	return result.ps;
}
#endif

Float MiterHairSegment::sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
	/* Luminaire sampling not supported */
	Log(EError, "Area sampling not supported by MiterHairSegment");
	return 0;
	/*
	Point p = Point(m_radius * std::cos(sample.y),
			m_radius * std::sin(sample.y), 
			sample.x * m_length);
	sRec.p = m_objectToWorld(p);
	sRec.n = normalize(m_objectToWorld(Normal(p.x, p.y, 0.0f)));
	return m_invSurfaceArea;
	*/
}

std::string MiterHairSegment::toString() const {
	std::ostringstream oss;
	oss << "MiterHairSegment [" << endl
			<< "  index = " << m_iv << endl
			<< "]";
	return oss.str();
}


MTS_IMPLEMENT_CLASS_S(MiterHairSegment, false, Shape)
MTS_NAMESPACE_END
