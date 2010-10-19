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
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/**
 * Sphere primitive.
 */

class Sphere : public Shape {
public:
	Sphere(const Properties &props) : Shape(props) {
		/**
		 * There are two ways of instantiating spheres: either,
		 * one can specify a linear transformation to from the
		 * unit sphere using the 'toWorld' parameter, or one
		 * can explicitly specify a radius and center.
		 */
		if (props.hasProperty("center") && props.hasProperty("radius")) {
			m_objectToWorld = 
				Transform::translate(Vector(props.getPoint("center")));
			m_radius = props.getFloat("radius");
		} else {
			Transform objectToWorld = props.getTransform("toWorld", Transform());
			m_radius = objectToWorld(Vector(1,0,0)).length();
			// Remove the scale from the object-to-world trasnsform
			m_objectToWorld = objectToWorld * Transform::scale(Vector(1/m_radius));
		}

		/// Are the sphere normals pointing inwards? default: no
		m_inverted = props.getBoolean("inverted", false);
		m_center = m_objectToWorld(Point(0,0,0));
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(4*M_PI*m_radius*m_radius);
	}

	Sphere(Stream *stream, InstanceManager *manager) 
			: Shape(stream, manager) {
		m_objectToWorld = Transform(stream);
		m_radius = stream->readFloat();
		m_center = Point(stream);
		m_inverted = stream->readBool();
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(4*M_PI*m_radius*m_radius);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
		stream->writeFloat(m_radius);
		m_center.serialize(stream);
		stream->writeBool(m_inverted);
	}

	AABB getAABB() const {
		AABB aabb;
		Float absRadius = std::abs(m_radius);
		aabb.min = m_center - Vector(absRadius);
		aabb.max = m_center + Vector(absRadius);
		return aabb;
	}

	Float getSurfaceArea() const {
		return 4*M_PI*m_radius*m_radius;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *tmp) const {
		Vector ro = ray.o - m_center;

		/* Transform into the local coordinate system and normalize */
		Float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
		Float B = 2 * (ray.d.x*ro.x + ray.d.y*ro.y + ray.d.z*ro.z);
		Float C = ro.x*ro.x + ro.y*ro.y +
				ro.z*ro.z - m_radius*m_radius;

		Float nearT, farT;
		if (!solveQuadratic(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;
		if (nearT < mint) {
			if (farT > maxt)
				return false;
			t = farT;		
		} else {
			t = nearT;
		}

		return true;
	}

	void fillIntersectionRecord(const Ray &ray, Float t, 
			const void *temp, Intersection &its) const {
		its.t = t;
		its.p = ray(t);
		Vector local = m_worldToObject(its.p - m_center);
		Float theta = std::acos(std::min(std::max(local.z/m_radius, 
				-(Float) 1), (Float) 1));
		Float phi = std::atan2(local.y, local.x);

		if (phi < 0)
			phi += 2*M_PI;

		its.uv.x = phi * (0.5 * INV_PI);
		its.uv.y = theta * INV_PI;
		its.dpdu = m_objectToWorld(Vector(-local.y, local.x, 0) * (2*M_PI));
    	its.geoFrame.n = normalize(its.p - m_center);
		Float zrad = std::sqrt(local.x*local.x + local.y*local.y);

		if (zrad > 0) {
			Float invZRad = 1.0f / zrad,
				  cosPhi = local.x * invZRad,
				  sinPhi = local.y * invZRad;
			its.dpdv = m_objectToWorld(Vector(local.z * cosPhi, local.z * sinPhi,
					-std::sin(theta)*m_radius) * M_PI);
    		its.geoFrame.s = normalize(its.dpdu);
			its.geoFrame.t = normalize(its.dpdv);
		} else {
			// avoid a singularity
			const Float cosPhi = 0, sinPhi = 1;
			its.dpdv = m_objectToWorld(Vector(local.z * cosPhi, local.z * sinPhi,
					-std::sin(theta)*m_radius) * M_PI);
			coordinateSystem(its.geoFrame.n, its.geoFrame.s, its.geoFrame.t);
		}

		if (m_inverted)
			its.geoFrame.n *= -1;

 		its.shFrame = its.geoFrame;
 		its.wi = its.toLocal(-ray.d);
		its.shape = this;
 		its.hasUVPartials = false;
	}

	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
		Vector v = squareToSphere(sample);
		sRec.n = Normal(v);
		sRec.p = Point(v * m_radius) + m_center;
		return 1.0f / (4*M_PI*m_radius*m_radius);
	}

	/**
	 * Improved sampling strategy given in
	 * "Monte Carlo techniques for direct lighting calculations" by
	 * Shirley, P. and Wang, C. and Zimmerman, K. (TOG 1996)
	 */
	Float sampleSolidAngle(ShapeSamplingRecord &sRec, const Point &p, const Point2 &sample) const {
		Vector w = m_center - p; Float invDistW = 1 / w.length();
		Float squareTerm = std::abs(m_radius * invDistW); // Support negative radii

		if (squareTerm >= 1-Epsilon) {
			/* We're inside the sphere - switch to uniform sampling */
			Vector d(squareToSphere(sample));

			sRec.p = m_center + d * m_radius;
			sRec.n = Normal(d);

			Vector lumToPoint = p - sRec.p;
			Float distSquared = lumToPoint.lengthSquared(), dp = dot(lumToPoint, sRec.n);

			if (dp > 0)
				return m_invSurfaceArea * distSquared * std::sqrt(distSquared) / dp;
			else
				return 0;
		}

		Float cosThetaMax = std::sqrt(std::max((Float) 0, 1 - squareTerm*squareTerm));

		Vector d = Frame(w*invDistW).toWorld(
			squareToCone(cosThetaMax, sample));

		Ray ray(p, d);
		Float t;
		if (!rayIntersect(ray, 0, std::numeric_limits<Float>::infinity(), t, NULL)) {
			// This can happen sometimes due to roundoff errors - just fail to 
			// generate a sample in this case.
			return 0;
		}

		sRec.p = ray(t);
		sRec.n = Normal(normalize(sRec.p-m_center));

		return 1 / ((2*M_PI) * (1-cosThetaMax));
	}

	Float pdfSolidAngle(const ShapeSamplingRecord &sRec, const Point &p) const {
		Vector w = p - m_center; Float invDistW = 1 / w.length();
		Float squareTerm = std::abs(m_radius * invDistW);

		if (squareTerm >= 1-Epsilon) {
			/* We're inside the sphere - switch to uniform sampling */
			Vector lumToPoint = p - sRec.p;
			Float distSquared = lumToPoint.lengthSquared(), dp = dot(lumToPoint, sRec.n);
			if (dp > 0)
				return m_invSurfaceArea * distSquared * std::sqrt(distSquared) / dp;
			else
				return 0;
		}

		Float cosThetaMax = std::sqrt(std::max((Float) 0, 1 - squareTerm*squareTerm));
		return squareToConePdf(cosThetaMax);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Sphere[" << endl
			<< "  radius = " << m_radius << ", " << endl
			<< "  center = " << m_center.toString() << ", " << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
			<< "  luminaire = " << indent(m_luminaire.toString()) << "," << endl
			<< "  subsurface = " << indent(m_subsurface.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_objectToWorld;
	Transform m_worldToObject;
	Point m_center;
	Float m_radius;
	Float m_invSurfaceArea;
	bool m_inverted;
};

MTS_IMPLEMENT_CLASS_S(Sphere, false, Shape)
MTS_EXPORT_PLUGIN(Sphere, "Sphere intersection primitive");
MTS_NAMESPACE_END
