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

MTS_NAMESPACE_BEGIN

class Cylinder : public Shape {
private:
	Float m_radius, m_length;
public:
	Cylinder(const Properties &props) : Shape(props) {
		m_radius = props.getFloat("radius", 1.0f);
		m_length = props.getFloat("length", 1.0f);
		if (m_objectToWorld.hasScale()) 
			Log(EError, "The scale needs to be specified using the 'radius' parameter!");
	}

	Cylinder(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		m_radius = stream->readFloat();
		m_length = stream->readFloat();
	}

	void configure() {
		Shape::configure();

		m_aabb.reset(); m_bsphere.radius = 0;
		m_surfaceArea = 2*M_PI*m_radius*m_length;
		m_invSurfaceArea = 1.0f / m_surfaceArea;

		m_aabb = getWorldAABB(0, m_length);
		m_bsphere.center = m_aabb.getCenter();
		for (int i=0; i<8; ++i) 
			m_bsphere.expandBy(m_aabb.getCorner(i));
	}

	bool isClippable() const {
		return true;
	}
	
	AABB getWorldAABB(Float start, Float end) const {
		AABB result;
		const Float r = m_radius;
		const Point a = m_objectToWorld(Point(0,0,start));
		const Point b = m_objectToWorld(Point(0,0,end));

		result.expandBy(a - Vector(r, r, r));
		result.expandBy(a + Vector(r, r, r));
		result.expandBy(b - Vector(r, r, r));
		result.expandBy(b + Vector(r, r, r));
		return result;
	}

	AABB getClippedAABB(const AABB &aabb) const {
		Float nearT, farT;
		AABB result(m_aabb);
		result.clip(aabb);
		
		Point a = m_objectToWorld(Point(0,0,0));
		Point b = m_objectToWorld(Point(0,0,m_length));

		if (!result.rayIntersect(Ray(a, normalize(b-a)), nearT, farT))
			return result; // that could be improved

		nearT = std::max(nearT, (Float) 0);
		farT = std::min(farT, m_length);
		result = getWorldAABB(nearT, farT);
		result.clip(aabb);

		return result;
	}

	bool rayIntersect(const Ray &_ray, Float start, Float end, Float &t) const {
		double nearT, farT; Ray ray;

		/* Transform into the local coordinate system and normalize */
		m_worldToObject(_ray, ray);
		
		const double
			ox = ray.o.x,
			oy = ray.o.y,
			dx = ray.d.x, 
			dy = ray.d.y;

		const double A = dx*dx + dy*dy;
		const double B = 2 * (dx*ox + dy*oy);
		const double C = ox*ox + oy*oy - m_radius*m_radius;

		if (!solveQuadraticDouble(A, B, C, nearT, farT))
			return false;

		if (nearT > end || farT < start)
			return false;

		const double zPosNear = ray.o.z + ray.d.z * nearT;
		const double zPosFar = ray.o.z + ray.d.z * farT;
		if (zPosNear >= 0 && zPosNear <= m_length && nearT >= start) {
			t = (Float) nearT;
		} else if (zPosFar >= 0 && zPosFar <= m_length) {
			if (farT > end)
				return false;
			t = (Float) farT;
		} else {
			return false;
		}

		return true;
	}

	bool rayIntersect(const Ray &ray, Intersection &its) const {
		if (!rayIntersect(ray, ray.mint, ray.maxt, its.t)) 
			return false;
		its.p = ray(its.t);

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
		its.geoFrame.s = normalize(its.dpdu - its.geoFrame.n
			* dot(its.geoFrame.n, its.dpdu));
		its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.hasUVPartials = false;
		its.shape = this;

		return true;
	}

#if defined(MTS_SSE)
	/* SSE-accelerated packet tracing is not supported for cylinders at the moment */
	__m128 rayIntersectPacket(const RayPacket4 &packet, const
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

	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
		Point p = Point(m_radius * std::cos(sample.y), 
			m_radius * std::sin(sample.y), 
			sample.x * m_length);
		sRec.p = m_objectToWorld(p);
		sRec.n = normalize(m_objectToWorld(Normal(p.x, p.y, 0.0f)));
		return m_invSurfaceArea;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		stream->writeFloat(m_radius);
		stream->writeFloat(m_length);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Cylinder[" << endl
			<< "  radius = " << m_radius << ", " << endl
			<< "  length = " << m_length << ", " << endl
			<< "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
			<< "  aabb = " << m_aabb.toString() << "," << endl
			<< "  bsphere = " << m_bsphere.toString() << "," << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
			<< "  luminaire = " << indent(m_luminaire.toString()) << "," << endl
			<< "  subsurface = " << indent(m_subsurface.toString()) << "," << endl
			<< "  surfaceArea = " << m_surfaceArea << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(Cylinder, false, Shape)
MTS_EXPORT_PLUGIN(Cylinder, "Cylinder intersection primitive");
MTS_NAMESPACE_END
