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

#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{sphere}{Sphere intersection primitive}
 * \order{1}
 * \parameters{
 *     \parameter{center}{\Point}{
 *	     Center of the sphere in object-space \default{(0, 0, 0)}
 *	   }
 *     \parameter{radius}{\Float}{
 *	     Radius of the sphere in object-space units \default{1}
 *	   }
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *	      Specifies an optional linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 *     \parameter{flipNormals}{\Boolean}{
 *	      Is the sphere inverted, i.e. should the normal vectors
 *		  be flipped? \default{\code{false}, i.e. the normals point outside}
 *	   }
 * }
 *
 * \renderings{
 *     \rendering{Basic example, see \lstref{sphere-basic}}
 *         {shape_sphere_basic}
 *     \rendering{A textured sphere with the default parameterization}
 *         {shape_sphere_parameterization}
 * }
 *
 * This shape plugin describes a simple sphere intersection primitive. It should
 * always be preferred over sphere approximations modeled using triangles.
 *
 * \begin{xml}[caption={A sphere can either be configured using a linear
 * \code{toWorld} transformation or the \code{center} and \code{radius} parameters (or both).
 *    The above two declarations are equivalent.}, label=lst:sphere-basic]
 * <shape type="sphere">
 *     <transform name="toWorld">
 *         <scale value="2"/>
 *         <translate x="1" y="0" z="0"/>
 *     </transform>
 *     <bsdf type="diffuse"/>
 * </shape>
 *
 * <shape type="sphere">
 *     <point name="center" x="1" y="0" z="0"/>
 *     <float name="radius" value="2"/>
 *     <bsdf type="diffuse"/>
 * </shape>
 * \end{xml}
 * When a \pluginref{sphere} shape is turned into an \pluginref{area} light source,
 * Mitsuba switches to an efficient sampling strategy \cite{Shirley91Direct} that
 * has particularly low variance. This makes it a good default choice for lighting
 * new scenes (\figref{spherelight}).
 * \renderings{
 *     \rendering{Spherical area light modeled using triangles}
 *         {shape_sphere_arealum_tri}
 *     \rendering{Spherical area light modeled using the \pluginref{sphere} plugin}
 *         {shape_sphere_arealum_analytic}
 *
 *     \caption{
 *         \label{fig:spherelight}
 *         Area lights built from the combination of the \pluginref{area}
 *         and \pluginref{sphere} plugins produce renderings that have an
 *         overall lower variance.
 *     }
 * }
 * \begin{xml}[caption=Instantiation of a sphere emitter]
 * <shape type="sphere">
 *     <point name="center" x="0" y="1" z="0"/>
 *     <float name="radius" value="1"/>

 *     <emitter type="area">
 *         <blackbody name="intensity" temperature="7000K"/>
 *     </emitter>
 * </shape>
 * \end{xml}
 */
class Sphere : public Shape {
public:
	Sphere(const Properties &props) : Shape(props) {
		m_objectToWorld =
			Transform::translate(Vector(props.getPoint("center", Point(0.0f))));
		m_radius = props.getFloat("radius", 1.0f);

		if (props.hasProperty("toWorld")) {
			Transform objectToWorld = props.getTransform("toWorld");
			Float radius = objectToWorld(Vector(1,0,0)).length();
			// Remove the scale from the object-to-world transform
			m_objectToWorld =
				  objectToWorld
				* Transform::scale(Vector(1/radius))
				* m_objectToWorld;
			m_radius *= radius;
		}

		/// Are the sphere normals pointing inwards? default: no
		m_flipNormals = props.getBoolean("flipNormals", false);
		m_center = m_objectToWorld(Point(0,0,0));
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(4*M_PI*m_radius*m_radius);

		if (m_radius <= 0)
			Log(EError, "Cannot create spheres of radius <= 0");
	}

	Sphere(Stream *stream, InstanceManager *manager)
			: Shape(stream, manager) {
		m_objectToWorld = Transform(stream);
		m_radius = stream->readFloat();
		m_center = Point(stream);
		m_flipNormals = stream->readBool();
		m_worldToObject = m_objectToWorld.inverse();
		m_invSurfaceArea = 1/(4*M_PI*m_radius*m_radius);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
		stream->writeFloat(m_radius);
		m_center.serialize(stream);
		stream->writeBool(m_flipNormals);
	}

	AABB getAABB() const {
		AABB aabb;
		aabb.min = m_center - Vector(m_radius);
		aabb.max = m_center + Vector(m_radius);
		return aabb;
	}

	Float getSurfaceArea() const {
		return 4*M_PI*m_radius*m_radius;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *tmp) const {
		Vector3d o = Vector3d(ray.o) - Vector3d(m_center);
		Vector3d d(ray.d);

		double A = d.lengthSquared();
		double B = 2 * dot(o, d);
		double C = o.lengthSquared() - m_radius*m_radius;

		double nearT, farT;
		if (!solveQuadraticDouble(A, B, C, nearT, farT))
			return false;

		if (!(nearT <= maxt && farT >= mint)) /* NaN-aware conditionals */
			return false;

		if (nearT < mint) {
			if (farT > maxt)
				return false;
			t = (Float) farT;
		} else {
			t = (Float) nearT;
		}

		return true;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Vector3d o = Vector3d(ray.o) - Vector3d(m_center);
		Vector3d d(ray.d);

		double A = d.lengthSquared();
		double B = 2 * dot(o, d);
		double C = o.lengthSquared() - m_radius*m_radius;

		double nearT, farT;
		if (!solveQuadraticDouble(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;
		if (nearT < mint && farT > maxt)
			return false;

		return true;
	}

	void fillIntersectionRecord(const Ray &ray,
			const void *temp, Intersection &its) const {
		its.p = ray(its.t);

		#if defined(SINGLE_PRECISION)
			/* Re-project onto the sphere to limit cancellation effects */
			its.p = m_center + normalize(its.p - m_center) * m_radius;
		#endif

		Vector local = m_worldToObject(its.p - m_center);
		Float theta = math::safe_acos(local.z/m_radius);
		Float phi = std::atan2(local.y, local.x);

		if (phi < 0)
			phi += 2*M_PI;

		its.uv.x = phi * (0.5f * INV_PI);
		its.uv.y = theta * INV_PI;
		its.dpdu = m_objectToWorld(Vector(-local.y, local.x, 0) * (2*M_PI));
		its.geoFrame.n = normalize(its.p - m_center);
		Float zrad = std::sqrt(local.x*local.x + local.y*local.y);
		its.shape = this;

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

		if (m_flipNormals)
			its.geoFrame.n *= -1;

 		its.shFrame.n = its.geoFrame.n;
 		its.hasUVPartials = false;
		its.instance = NULL;
		its.time = ray.time;
	}

	void samplePosition(PositionSamplingRecord &pRec, const Point2 &sample) const {
		Vector v = warp::squareToUniformSphere(sample);

		pRec.p = Point(v * m_radius) + m_center;
		pRec.n = Normal(v);

		if (m_flipNormals)
			pRec.n *= -1;

		pRec.pdf = m_invSurfaceArea;
		pRec.measure = EArea;
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_invSurfaceArea;
	}

	void getNormalDerivative(const Intersection &its,
			Vector &dndu, Vector &dndv, bool shadingFrame) const {
		Float invRadius = (m_flipNormals ? -1.0f : 1.0f) / m_radius;
		dndu = its.dpdu * invRadius;
		dndv = its.dpdv * invRadius;
	}

	/**
	 * Improved sampling strategy given in
	 * "Monte Carlo techniques for direct lighting calculations" by
	 * Shirley, P. and Wang, C. and Zimmerman, K. (TOG 1996)
	 */
	void sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
		const Vector refToCenter = m_center - dRec.ref;
		const Float refDist2 = refToCenter.lengthSquared();
		const Float invRefDist = static_cast<Float>(1) / std::sqrt(refDist2);

		/* Sine of the angle of the cone containing the
		   sphere as seen from 'dRec.ref' */
		const Float sinAlpha = m_radius * invRefDist;

		if (sinAlpha < 1-Epsilon) {
			/* The reference point lies outside of the sphere.
			   => sample based on the projected cone. */

			Float cosAlpha = math::safe_sqrt(1.0f - sinAlpha * sinAlpha);

			dRec.d = Frame(refToCenter * invRefDist).toWorld(
				warp::squareToUniformCone(cosAlpha, sample));
			dRec.pdf = warp::squareToUniformConePdf(cosAlpha);

			/* Distance to the projection of the sphere center
			   onto the ray (dRec.ref, dRec.d) */
			const Float projDist = dot(refToCenter, dRec.d);

			/* To avoid numerical problems move the query point to the
			   intersection of the of the original direction ray and a plane
			   with normal refToCenter which goes through the sphere's center */
			const Float baseT = refDist2 / projDist;
			const Point query = dRec.ref + dRec.d * baseT;

			const Vector queryToCenter = m_center - query;
			const Float queryDist2     = queryToCenter.lengthSquared();
			const Float queryProjDist  = dot(queryToCenter, dRec.d);

			/* Try to find the intersection point between the
			   sampled ray and the sphere. */
			Float A = 1.0f, B = -2*queryProjDist,
				  C = queryDist2 - m_radius*m_radius;

			Float nearT, farT;
			if (!solveQuadratic(A, B, C, nearT, farT)) {
				/* The intersection couldn't be found due to roundoff errors..
				   Don't give up -- one workaround is to project the closest
				   ray position onto the sphere */
				nearT = queryProjDist;
			}

			dRec.dist = baseT + nearT;
			dRec.n = normalize(dRec.d*nearT - queryToCenter);
			dRec.p = m_center + dRec.n * m_radius;
		} else {
			/* The reference point lies inside the sphere
			   => use uniform sphere sampling. */
			Vector d = warp::squareToUniformSphere(sample);

			dRec.p = m_center + d * m_radius;
			dRec.n = Normal(d);
			dRec.d = dRec.p - dRec.ref;

			Float dist2 = dRec.d.lengthSquared();
			dRec.dist = std::sqrt(dist2);
			dRec.d /= dRec.dist;
			dRec.pdf = m_invSurfaceArea * dist2
				/ absDot(dRec.d, dRec.n);
		}

		if (m_flipNormals)
			dRec.n *= -1;

		dRec.measure = ESolidAngle;
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		const Vector refToCenter = m_center - dRec.ref;
		const Float invRefDist = (Float) 1.0f / refToCenter.length();

		/* Sine of the angle of the cone containing the
		   sphere as seen from 'dRec.ref' */
		const Float sinAlpha = m_radius * invRefDist;

		if (sinAlpha < 1-Epsilon) {
			/* The reference point lies outside the sphere */
			Float cosAlpha = math::safe_sqrt(1 - sinAlpha*sinAlpha);
			Float pdfSA = warp::squareToUniformConePdf(cosAlpha);

			if (dRec.measure == ESolidAngle)
				return pdfSA;
			else if (dRec.measure == EArea)
				return pdfSA * absDot(dRec.d, dRec.n)
					/ (dRec.dist*dRec.dist);
			else
				return 0.0f;
		} else {
			/* The reference point lies inside the sphere */
			if (dRec.measure == ESolidAngle)
				return m_invSurfaceArea * dRec.dist * dRec.dist
					/ absDot(dRec.d, dRec.n);
			else if (dRec.measure == EArea)
				return m_invSurfaceArea;
			else
				return 0.0f;
		}
	}

	ref<TriMesh> createTriMesh() {
		/// Choice of discretization
		const uint32_t thetaSteps = 20;
		const uint32_t phiSteps = thetaSteps * 2;
		const Float dTheta = M_PI / (thetaSteps-1);
		const Float dPhi   = (2*M_PI) / (phiSteps-1);

		/// Precompute cosine and sine tables
		Float *cosPhi = new Float[phiSteps];
		Float *sinPhi = new Float[phiSteps];
		for (uint32_t i=0; i<phiSteps; ++i) {
			sinPhi[i] = std::sin(i*dPhi);
			cosPhi[i] = std::cos(i*dPhi);
		}

		size_t numTris = 2 * (phiSteps-1) * (thetaSteps-1);
		size_t numVertices = thetaSteps * phiSteps;

		ref<TriMesh> mesh = new TriMesh("Sphere approximation",
			numTris, numVertices, true, true, false);

		Point *vertices = mesh->getVertexPositions();
		Normal *normals = mesh->getVertexNormals();
		Point2 *texcoords = mesh->getVertexTexcoords();
		Triangle *triangles = mesh->getTriangles();
		uint32_t vertexIdx = 0;
		for (uint32_t theta=0; theta<thetaSteps; ++theta) {
			Float sinTheta = std::sin(theta * dTheta);
			Float cosTheta = std::cos(theta * dTheta);

			for (uint32_t phi=0; phi<phiSteps; ++phi) {
				Vector v(
					sinTheta * cosPhi[phi],
					sinTheta * sinPhi[phi],
					cosTheta
				);
				texcoords[vertexIdx] = Point2(phi * dPhi * INV_TWOPI, theta * dTheta * INV_PI);
				vertices[vertexIdx] = m_objectToWorld(Point(v*m_radius));
				normals[vertexIdx++] = m_objectToWorld(Normal(v));
			}
		}
		Assert(vertexIdx == numVertices);

		uint32_t triangleIdx = 0;
		for (uint32_t theta=1; theta<thetaSteps; ++theta) {
			for (uint32_t phi=0; phi<phiSteps-1; ++phi) {
				uint32_t nextPhi = phi + 1;
				uint32_t idx0 = phiSteps*theta + phi;
				uint32_t idx1 = phiSteps*theta + nextPhi;
				uint32_t idx2 = phiSteps*(theta-1) + phi;
				uint32_t idx3 = phiSteps*(theta-1) + nextPhi;

				triangles[triangleIdx].idx[0] = idx0;
				triangles[triangleIdx].idx[1] = idx2;
				triangles[triangleIdx].idx[2] = idx1;
				triangleIdx++;
				triangles[triangleIdx].idx[0] = idx1;
				triangles[triangleIdx].idx[1] = idx2;
				triangles[triangleIdx].idx[2] = idx3;
				triangleIdx++;
			}
		}
		Assert(triangleIdx == numTris);
		delete[] cosPhi;
		delete[] sinPhi;
		mesh->copyAttachments(this);
		mesh->configure();

		return mesh.get();
	}

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return 1;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Sphere[" << endl
			<< "  radius = " << m_radius << "," << endl
			<< "  center = " << m_center.toString() << "," << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl;
		if (isMediumTransition())
			oss << "  interiorMedium = " << indent(m_interiorMedium.toString()) << "," << endl
				<< "  exteriorMedium = " << indent(m_exteriorMedium.toString()) << "," << endl;
		oss << "  emitter = " << indent(m_emitter.toString()) << "," << endl
			<< "  sensor = " << indent(m_sensor.toString()) << "," << endl
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
	bool m_flipNormals;
};

MTS_IMPLEMENT_CLASS_S(Sphere, false, Shape)
MTS_EXPORT_PLUGIN(Sphere, "Sphere intersection primitive");
MTS_NAMESPACE_END
