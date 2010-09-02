#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

/**
 * Sphere primitive.
 */

class Sphere : public Shape {
private:
	Point m_center;
	Float m_radius;
public:
	Sphere(const Properties &props) : Shape(props) {
		m_radius = props.getFloat("radius", 1.0f); // Negative radius -> inside-out sphere
		if (m_objectToWorld.hasScale()) 
			Log(EError, "The scale needs to be specified using the 'radius' parameter!");
	}
	
	Sphere(Stream *stream, InstanceManager *manager) 
			: Shape(stream, manager) {
		m_radius = stream->readFloat();
		configure();
	}

	void configure() {
		Shape::configure();

		m_surfaceArea = 4*M_PI*m_radius*m_radius;
		m_invSurfaceArea = 1.0f / m_surfaceArea;
		m_center = m_objectToWorld(Point(0,0,0));
		Float absRadius = std::abs(m_radius);
		m_aabb.min = m_center - Vector(absRadius, absRadius, absRadius);
		m_aabb.max = m_center + Vector(absRadius, absRadius, absRadius);
		m_bsphere.center = m_center;
		m_bsphere.radius = absRadius;
	}

	bool rayIntersect(const Ray &ray, Float start, Float end, Float &t) const {
		/* Transform into the local coordinate system and normalize */
		double nearT, farT;
		const double ox = (double) ray.o.x - (double) m_center.x, 
			oy = (double) ray.o.y - (double) m_center.y, 
			oz = (double) ray.o.z - (double) m_center.z;
		const double dx = ray.d.x, dy = ray.d.y, dz = ray.d.z;
		const double A = dx*dx + dy*dy + dz*dz;
		const double B = 2 * (dx*ox + dy*oy + dz*oz);
		const double C = ox*ox + oy*oy + oz*oz - m_radius * m_radius;

		if (!solveQuadraticDouble(A, B, C, nearT, farT))
			return false;

		if (nearT > end || farT < start)
			return false;
		if (nearT < start) {
			if (farT > end)
				return false;
			t = (Float) farT;		
		} else {
			t = (Float) nearT;
		}
		return true;
	}

	bool rayIntersect(const Ray &ray, Intersection &its) const {
		if (!rayIntersect(ray, ray.mint, ray.maxt, its.t)) 
			return false;
		its.p = ray(its.t);

		Point local = m_worldToObject(its.p);
		Float absRadius = std::abs(m_radius);
		Float theta = std::acos(local.z / absRadius);
		Float phi = std::atan2(local.y, local.x);
		if (phi < 0)
			phi += 2*M_PI;
		its.uv.x = theta / M_PI;
		its.uv.y = phi / (2*M_PI);
		its.shape = this;
		Float zrad = std::sqrt(local.x*local.x + local.y*local.y);

		its.geoFrame.n = Normal(normalize(its.p - m_center));
		if (m_radius < 0)
			its.geoFrame.n *= -1;

		if (zrad > 0) {
			Float invZRad = 1.0f / zrad;
			Float cosPhi = local.x * invZRad;
			Float sinPhi = local.y * invZRad;
			its.dpdu = m_objectToWorld(Vector(-local.y, local.x, 0) * (2*M_PI));
			its.dpdv = m_objectToWorld(Vector(local.z * cosPhi, local.z * sinPhi,
					- absRadius * std::sin(theta)) * M_PI);
			its.geoFrame.s = normalize(its.dpdu - its.geoFrame.n
				* dot(its.geoFrame.n, its.dpdu));
			its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
		} else {
			// avoid a singularity
			Float cosPhi = 0, sinPhi = 1;
			its.dpdv = m_objectToWorld(Vector(local.z * cosPhi, local.z * sinPhi,
				-absRadius*std::sin(theta))* M_PI);
			its.dpdu = cross(its.dpdv, its.geoFrame.n);
			coordinateSystem(its.geoFrame.n, its.geoFrame.s, its.geoFrame.t);
		}

 		its.shFrame = its.geoFrame;
 		its.wi = its.toLocal(-ray.d);
 		its.hasUVPartials = false;

		return true;
	}

#if defined(MTS_SSE)
	/* SSE-accelerated packet tracing is not supported for spheres at the moment */
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
		Vector v = squareToSphere(sample);
		sRec.n = m_objectToWorld(Normal(v));
		sRec.p = m_objectToWorld(Point(v * m_radius));
		return m_invSurfaceArea;
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
		if (!rayIntersect(ray, 0, std::numeric_limits<Float>::infinity(), t)) {
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

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		stream->writeFloat(m_radius);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Sphere[" << endl
			<< "  radius = " << m_radius << ", " << endl
			<< "  center = " << m_center.toString() << ", " << endl
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

MTS_IMPLEMENT_CLASS_S(Sphere, false, Shape)
MTS_EXPORT_PLUGIN(Sphere, "Sphere intersection primitive");
MTS_NAMESPACE_END
