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

/*!\plugin{disk}{Disk intersection primitive}
 * \order{4}
 * \parameters{
 *     \parameter{toWorld}{\Transform}{
 *	      Specifies a linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 *     \parameter{flipNormals}{\Boolean}{
 *	      Is the disk inverted, i.e. should the normal vectors
 *		  be flipped? \default{\code{false}}
 *	   }
 *     \vspace{-8mm}
 * }
 * \renderings{
 *     \rendering{Rendering with an disk emitter and a textured disk, showing
 *     the default parameterization. (\lstref{disk})}{shape_disk}
 * }
 * 
 * \vspace{-1mm}
 * This shape plugin describes a simple disk intersection primitive. It is
 * usually preferable over discrete approximations made from triangles.
 *
 * By default, the disk has unit radius and is located at the origin. Its
 * surface normal points into the positive $Z$ direction.
 * To change the disk scale, rotation, or translation, use the 
 * \code{toWorld} parameter.
 *
 * \begin{xml}[caption={A simple example involving two disk instances}, label=lst:disk]
 * <scene version=$\MtsVer$>
 *     <shape type="disk">
 *         <bsdf type="diffuse">
 *             <texture name="reflectance" type="checkerboard">
 *                 <float name="uvscale" value="5"/>
 *             </texture>
 *         </bsdf>
 *     </shape>
 *     <shape type="disk">
 *         <transform name="toWorld">
 *             <rotate x="1" angle="90"/>
 *             <scale value="0.3"/>
 *             <translate y="1" z="0.3"/>
 *         </transform>
 *         <luminaire type="area">
 *             <spectrum name="intensity" value="4"/>
 *         </luminaire>
 *     </shape>
 * </scene>
 * \end{xml}
 */
class Disk : public Shape {
public:
	Disk(const Properties &props) : Shape(props) {
		m_objectToWorld = props.getTransform("toWorld", Transform());
		if (props.getBoolean("flipNormals", false))
			m_objectToWorld = m_objectToWorld * Transform::scale(Vector(1, 1, -1));
		m_worldToObject = m_objectToWorld.inverse();
	}

	Disk(Stream *stream, InstanceManager *manager) 
			: Shape(stream, manager) {
		m_objectToWorld = Transform(stream);
		m_worldToObject = m_objectToWorld.inverse();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
	}

	void configure() {
		Shape::configure();

		m_normal = normalize(m_objectToWorld(Normal(0, 0, 1)));

		Vector dpdu = m_objectToWorld(Vector(1, 0, 0));
		Vector dpdv = m_objectToWorld(Vector(0, 1, 0));

		if (std::abs(dot(dpdu, dpdv)) > Epsilon)
			Log(EError, "Error: 'toWorld' transformation contains shear!");

		if (std::abs(dpdu.length() / dpdv.length() - 1) > Epsilon)
			Log(EError, "Error: 'toWorld' transformation contains a non-uniform scale!");

		m_surfaceArea = M_PI * dpdu.length() * dpdu.length();
	}

	AABB getAABB() const {
		AABB aabb;
		aabb.expandBy(m_objectToWorld(Point( 1,  0, 0)));
		aabb.expandBy(m_objectToWorld(Point(-1,  0, 0)));
		aabb.expandBy(m_objectToWorld(Point( 0,  1, 0)));
		aabb.expandBy(m_objectToWorld(Point( 0, -1, 0)));
		return aabb;
	}

	Float getSurfaceArea() const {
		return m_surfaceArea;
	}

	inline bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		Ray ray;
		m_worldToObject.transformAffine(_ray, ray);
		Float hit = -ray.o.z/ray.d.z;

		if (hit < mint || hit > maxt)
			return false;

		Point local = ray(hit);

		if (local.x * local.x + local.y * local.y > 1)
			return false;

		t = hit;

		if (temp) {
			Float *data = static_cast<Float *>(temp);
			data[0] = local.x;
			data[1] = local.y;
		}

		return true;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Float t;
		return Disk::rayIntersect(ray, mint, maxt, t, NULL);
	}

	void fillIntersectionRecord(const Ray &ray, 
			const void *temp, Intersection &its) const {
		const Float *data = static_cast<const Float *>(temp);

		Float r = std::sqrt(data[0] * data[0] + data[1] * data[1]),
			  invR = (r == 0) ? 0.0f : (1.0f / r);

		Float phi = std::atan2(data[1], data[0]);
		if (phi < 0)
			phi += 2*M_PI;

		Float cosPhi = data[0] * invR, sinPhi = data[1] * invR;

		its.dpdu = m_objectToWorld(Vector(cosPhi, sinPhi, 0));
		its.dpdv = m_objectToWorld(Vector(-sinPhi, cosPhi, 0));

		its.shFrame = its.geoFrame = Frame(
			normalize(its.dpdu), normalize(its.dpdv), m_normal);
		its.uv = Point2(r, phi * INV_TWOPI);
		its.p = ray(its.t);
		its.wi = its.toLocal(-ray.d);
		its.shape = this;
 		its.hasUVPartials = false;
	}

	ref<TriMesh> createTriMesh() {
		const unsigned int phiSteps = 40;

		ref<TriMesh> mesh = new TriMesh(getName(),
			phiSteps-1, 2*phiSteps, true, true, false);

		Point *vertices = mesh->getVertexPositions();
		Normal *normals = mesh->getVertexNormals();
		Point2 *texcoords = mesh->getVertexTexcoords();
		Triangle *triangles = mesh->getTriangles();

		Float dphi = (2 * M_PI) / (Float) (phiSteps-1);

		Point center = m_objectToWorld(Point(0.0f));
		for (size_t i=0; i<phiSteps; ++i) {
			Float phi = i*dphi;
			vertices[i] = center;
			vertices[phiSteps+i] = m_objectToWorld(
				Point(std::cos(phi), std::sin(phi), 0)
			);

			normals[i] = m_normal;
			normals[phiSteps+i] = m_normal;
			texcoords[i] = Point2(0.0f, phi * INV_TWOPI);
			texcoords[phiSteps+i] = Point2(1.0f, phi * INV_TWOPI);
		}

		for (size_t i=0; i<phiSteps-1; ++i) {
			triangles[i].idx[0] = i;
			triangles[i].idx[1] = i+phiSteps;
			triangles[i].idx[2] = i+phiSteps+1;
		}

		mesh->setBSDF(m_bsdf);
		mesh->setLuminaire(m_luminaire);
		mesh->configure();

		return mesh.get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Disk[" << endl
			<< "  objectToWorld = " << indent(m_objectToWorld.toString()) << ", " << endl
			<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
			<< "  luminaire = " << indent(m_luminaire.toString()) << "," << endl
			<< "  subsurface = " << indent(m_subsurface.toString()) << endl
			<< "]";
		return oss.str();
	}

	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
		sRec.n = m_normal;
		Point2 p = squareToDiskConcentric(sample);
		sRec.p = m_objectToWorld(Point3(p.x, p.y, 0));
		return 1.0f / m_surfaceArea;
	}

	Float pdfArea(const ShapeSamplingRecord &sRec) const {
		return 1.0f / m_surfaceArea;
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_objectToWorld;
	Transform m_worldToObject;
	Normal m_normal;
	Float m_surfaceArea;
};

MTS_IMPLEMENT_CLASS_S(Disk, false, Shape)
MTS_EXPORT_PLUGIN(Disk, "Disk intersection primitive");
MTS_NAMESPACE_END
