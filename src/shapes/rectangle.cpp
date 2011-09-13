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

/*!\plugin{rectangle}{Rectangle intersection primitive}
 * \parameters{
 *     \parameter{toWorld}{\Transform}{
 *	      Specifies a linear object-to-world transformation.
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 * }
 * 
 * This shape plugin describes a simple rectangular intersection primitive.
 *
 * It is mainly provided as a convenience for those cases when loading an 
 * external mesh with two triangles is simply too tedious.
 */
class Rectangle : public Shape {
public:
	Rectangle(const Properties &props) : Shape(props) {
		m_objectToWorld = props.getTransform("toWorld", Transform());
		m_worldToObject = m_objectToWorld.inverse();
	}

	Rectangle(Stream *stream, InstanceManager *manager) 
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
		m_dpdu = m_objectToWorld(Vector(1, 0, 0));
		m_dpdv = m_objectToWorld(Vector(0, 1, 0));
		Normal normal = normalize(m_objectToWorld(Normal(0, 0, 1)));
		m_frame = Frame(normalize(m_dpdu), normalize(m_dpdv), normal);

		m_surfaceArea = 4 * m_dpdu.length() * m_dpdv.length();
		if (std::abs(dot(m_dpdu, m_dpdv)) > Epsilon)
			Log(EError, "Error: 'toWorld' transformation contains shear!");
	}

	AABB getAABB() const {
		AABB aabb;
		aabb.expandBy(m_objectToWorld(Point(-1, -1, 0)));
		aabb.expandBy(m_objectToWorld(Point( 1, -1, 0)));
		aabb.expandBy(m_objectToWorld(Point( 1,  1, 0)));
		aabb.expandBy(m_objectToWorld(Point(-1,  1, 0)));
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

		if (std::abs(local.x) > 1 || std::abs(local.y) > 1)
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
		return Rectangle::rayIntersect(ray, mint, maxt, t, NULL);
	}

	void fillIntersectionRecord(const Ray &ray, 
			const void *temp, Intersection &its) const {
		const Float *data = static_cast<const Float *>(temp);
		its.shFrame = its.geoFrame = m_frame;
		its.uv = Point2(data[0], data[1]);
		its.p = ray(its.t);
		its.wi = its.toLocal(-ray.d);
		its.shape = this;
 		its.hasUVPartials = false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Rectangle[" << endl
			<< "  objectToWorld = " << indent(m_objectToWorld.toString()) << ", " << endl
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
	Frame m_frame;
	Vector m_dpdu, m_dpdv;
	Float m_surfaceArea;
};

MTS_IMPLEMENT_CLASS_S(Rectangle, false, Shape)
MTS_EXPORT_PLUGIN(Rectangle, "Rectangle intersection primitive");
MTS_NAMESPACE_END
