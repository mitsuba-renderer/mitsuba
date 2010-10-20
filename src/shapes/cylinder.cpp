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

#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/random.h>

MTS_NAMESPACE_BEGIN

class Cylinder : public Shape {
private:
	Transform m_objectToWorld;
	Transform m_worldToObject;
	Float m_radius, m_length, m_invSurfaceArea;
public:
	Cylinder(const Properties &props) : Shape(props) {
		/**
		 * There are two ways of instantiating cylinders: either,
		 * one can specify a linear transformation to from the
		 * unit sphere using the 'toWorld' parameter, or one
		 * can explicitly specify two points and a radius.
		 */
		if (props.hasProperty("p1") && props.hasProperty("p2")
				&& props.hasProperty("radius")) {
			Point p1 = props.getPoint("p1"), p2 = props.getPoint("p2");
			Vector rel = p2 - p1;
			Float radius = props.getFloat("radius");
			Float length = rel.length();

			m_objectToWorld = 
				Transform::translate(Vector(p1)) *
				Transform::fromFrame(Frame(rel/length));
			m_radius = radius;
			m_length = length;
		} else {
			Transform objectToWorld = props.getTransform("toWorld", Transform());
			m_radius = objectToWorld(Vector(1,0,0)).length();
			m_length = objectToWorld(Vector(0,0,1)).length();
			// Remove the scale from the object-to-world trasnsform
			m_objectToWorld = objectToWorld * Transform::scale(
					Vector(1/m_radius, 1/m_radius, 1/m_length));
		}
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
		Assert(m_length > 0 && m_radius > 0);
	}

	Cylinder(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		m_objectToWorld = Transform(stream);
		m_radius = stream->readFloat();
		m_length = stream->readFloat();
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
		stream->writeFloat(m_radius);
		stream->writeFloat(m_length);
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		Ray ray;

		/* Transform into the local coordinate system and normalize */
		m_worldToObject(_ray, ray);

		const Float
			ox = ray.o.x,
			oy = ray.o.y,
			dx = ray.d.x, 
			dy = ray.d.y;

		const Float A = dx*dx + dy*dy;
		const Float B = 2 * (dx*ox + dy*oy);
		const Float C = ox*ox + oy*oy - m_radius*m_radius;

		Float nearT, farT;
		if (!solveQuadratic(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;

		const Float zPosNear = ray.o.z + ray.d.z * nearT;
		const Float zPosFar = ray.o.z + ray.d.z * farT;
		if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
			t = nearT;
		} else if (zPosFar >= 0 && zPosFar <= m_length) {
			if (farT > maxt)
				return false;
			t = farT;
		} else {
			return false;
		}

		return true;
	}


	void fillIntersectionRecord(const Ray &ray, Float t, 
			const void *temp, Intersection &its) const {
		its.p = ray(t);

		Point local = m_worldToObject(its.p);
		Float phi = std::atan2(local.y, local.x);
		if (phi < 0)
			phi += 2*M_PI;
		its.uv.x = local.z / m_length;
		its.uv.y = phi / (2*M_PI);

		Vector dpdu = Vector(-local.y, local.x, 0) * (2*M_PI);
		Vector dpdv = Vector(0, 0, m_length);
		its.dpdu = m_objectToWorld(dpdu);
		its.dpdv = m_objectToWorld(dpdv);
		its.geoFrame.n = Normal(normalize(m_objectToWorld(cross(dpdu, dpdv))));
		its.geoFrame.s = normalize(its.dpdu);
		its.geoFrame.t = normalize(its.dpdv);
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.hasUVPartials = false;
		its.shape = this;
	}

	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
		Point p = Point(m_radius * std::cos(sample.y), 
			m_radius * std::sin(sample.y), 
			sample.x * m_length);
		sRec.p = m_objectToWorld(p);
		sRec.n = normalize(m_objectToWorld(Normal(p.x, p.y, 0.0f)));
		return m_invSurfaceArea;
	}

	inline AABB getAABB() const {
		Vector x1 = m_objectToWorld(Vector(m_radius, 0, 0));
		Vector x2 = m_objectToWorld(Vector(0, m_radius, 0));
		Point p0 = m_objectToWorld(Point(0, 0, 0));
		Point p1 = m_objectToWorld(Point(0, 0, m_length));
		AABB result;

		/* To bound the cylinder, it is sufficient to find the
		   smallest box containing the two circles at the endpoints.
		   This can be done component-wise as follows */

		for (int i=0; i<3; ++i) {
			Float range = std::sqrt(x1[i]*x1[i] + x2[i]*x2[i]);

			result.min[i] = std::min(std::min(result.min[i], 
						p0[i]-range), p1[i]-range);
			result.max[i] = std::max(std::max(result.max[i], 
						p0[i]+range), p1[i]+range);
		}

		return result;
	}

	inline Float sqr(Float f) {
		return f*f;
	}


	/**
	 * Compute the ellipse created by the intersection of an infinite
	 * cylinder and a plane. Returns false in the degenerate case.
	 * Based on:
	 * www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
	 */
	inline bool intersectCylPlane(Point planePt, Normal planeNrml,
			Point cylPt, Vector cylD, Float radius) {
		if (absDot(planeNrml, cylD) < Epsilon)
			return false;

		Vector A = normalize(cylD - dot(planeNrml, cylD)*planeNrml),
			   B = cross(planeNrml, A),
			   delta = planePt - cylPt,
			   deltaProj = delta - cylD*dot(delta, cylD);

		Float c0 = 1-sqr(dot(A, cylD));
		Float c1 = 1-sqr(dot(B, cylD));
		Float c2 = 2*dot(A, deltaProj);
		Float c3 = 2*dot(B, deltaProj);
		Float c4 = dot(delta, deltaProj) - radius*radius;

		Float lambda = (c2*c2/(4*c0) + c3*c3/(4*c1) - c4)/(c0*c1);

		Float alpha0 = -c2/(2*c0),
			  beta0 = -c3/(2*c1),
			  L_alpha = std::sqrt(c1*lambda),
			  L_beta = std::sqrt(c0*lambda);

		Point center = planePt + alpha0 * A + beta0 * B;
		Vector axis1 = L_alpha * A;
		Vector axis2 = L_beta * B;
		return true;
	}

	AABB getClippedAABB(const AABB &box) const {
		/* Compute a base bounding box */
		AABB base(getAABB());
		base.clip(box);

		/* Now forget about the cylinder ends and 
		   intersect an infinite cylinder with each AABB face */
		return box;
	}

#if 0
	inline AABB getAABB(Float start, Float end) const {
		AABB result;
		const Float r = m_radius;
		const Point a = m_objectToWorld(Point(0, 0, start));
		const Point b = m_objectToWorld(Point(0, 0, end));

		result.expandBy(a - Vector(r, r, r));
		result.expandBy(a + Vector(r, r, r));
		result.expandBy(b - Vector(r, r, r));
		result.expandBy(b + Vector(r, r, r));
		return result;
	}

	AABB getAABB() const {
		/* Very approximate .. */
		return getAABB(0, m_length);
	}

	AABB getClippedAABB(const AABB &box) const {
		/* This is incorrect! */
		Float nearT, farT;
		AABB result(getAABB(0, m_length));
		result.clip(box);

		Point a = m_objectToWorld(Point(0, 0, 0));
		Point b = m_objectToWorld(Point(0, 0, m_length));

		if (!result.rayIntersect(Ray(a, normalize(b-a)), nearT, farT))
			return result; // that could be improved

		nearT = std::max(nearT, (Float) 0);
		farT = std::min(farT, m_length);
		result = getAABB(nearT, farT);
		result.clip(box);

		return result;
	}
#endif

	Float getSurfaceArea() const {
		return 2*M_PI*m_radius*m_length;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Cylinder[" << endl
			<< "  radius = " << m_radius << ", " << endl
			<< "  length = " << m_length << ", " << endl
			<< "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
			<< "  luminaire = " << indent(m_luminaire.toString()) << "," << endl
			<< "  subsurface = " << indent(m_subsurface.toString())
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(Cylinder, false, Shape)
MTS_EXPORT_PLUGIN(Cylinder, "Cylinder intersection primitive");
MTS_NAMESPACE_END
