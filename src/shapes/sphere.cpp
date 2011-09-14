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

#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>

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
 *     \parameter{toWorld}{\Transform}{
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
 * When using a sphere as the base object of an \pluginref{area} luminaire,
 * Mitsuba will switch to a special sphere luminaire sampling strategy 
 * \cite{Shirley91Direct} that works much better than the default approach. 
 * The resulting variance reduction makes it preferable to model most light
 * sources as sphere luminaires (\figref{spherelight}).
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
 *
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
 * \begin{xml}[caption=Instantiation of a sphere luminaire]
 * <shape type="sphere">
 *     <point name="center" x="0" y="1" z="0"/>
 *     <float name="radius" value="1"/>

 *     <luminaire type="area">
 *         <blackbody name="intensity" temperature="7000K"/>
 *     </luminaire>
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
		Float absRadius = std::abs(m_radius);
		aabb.min = m_center - Vector(absRadius);
		aabb.max = m_center + Vector(absRadius);
		return aabb;
	}

	Float getSurfaceArea() const {
		return 4*M_PI*m_radius*m_radius;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *tmp) const {
		Vector o = ray.o - m_center;
		Float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
		Float B = 2 * (ray.d.x*o.x + ray.d.y*o.y + ray.d.z*o.z);
		Float C = o.x*o.x + o.y*o.y + o.z*o.z - m_radius*m_radius;

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

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Vector o = ray.o - m_center;
		Float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
		Float B = 2 * (ray.d.x*o.x + ray.d.y*o.y + ray.d.z*o.z);
		Float C = o.x*o.x + o.y*o.y + o.z*o.z - m_radius*m_radius;

		Float nearT, farT;
		if (!solveQuadratic(A, B, C, nearT, farT))
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
		Vector local = m_worldToObject(its.p - m_center);
		Float theta = std::acos(std::min(std::max(local.z/m_radius, 
				-(Float) 1), (Float) 1));
		Float phi = std::atan2(local.y, local.x);

		if (phi < 0)
			phi += 2*M_PI;

		its.uv.x = phi * (0.5f * INV_PI);
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

		if (m_flipNormals)
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

	Float pdfArea(const ShapeSamplingRecord &sRec) const {
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

		Ray ray(p, d, 0.0f);
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

				if (idx0 != idx1) {
					triangles[triangleIdx].idx[0] = idx0;
					triangles[triangleIdx].idx[1] = idx2;
					triangles[triangleIdx].idx[2] = idx1;
					triangleIdx++;
				}
				if (idx2 != idx3) {
					triangles[triangleIdx].idx[0] = idx1;
					triangles[triangleIdx].idx[1] = idx2;
					triangles[triangleIdx].idx[2] = idx3;
					triangleIdx++;
				}
			}
		}
		Assert(triangleIdx == numTris);
		delete[] cosPhi;
		delete[] sinPhi;
		mesh->setBSDF(m_bsdf);
		mesh->setLuminaire(m_luminaire);
		mesh->configure();

		return mesh.get();
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
	bool m_flipNormals;
};

MTS_IMPLEMENT_CLASS_S(Sphere, false, Shape)
MTS_EXPORT_PLUGIN(Sphere, "Sphere intersection primitive");
MTS_NAMESPACE_END
