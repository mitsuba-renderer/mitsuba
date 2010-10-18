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
		return getAABB(0, m_length);
	}

	AABB getClippedAABB(const AABB &box) const {
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
