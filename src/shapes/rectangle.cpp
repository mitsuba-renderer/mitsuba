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
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld.serialize(stream);
	}

	void configure() {
		Vector du = m_objectToWorld(Vector(1, 0, 0));
		Vector dv = m_objectToWorld(Vector(0, 1, 0));

		m_surfaceArea = 4 * du.length() * dv.length();
		if (std::abs(dot(du, dv)) > Epsilon)
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

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *) const {
		return true;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		return true;
	}

	void fillIntersectionRecord(const Ray &ray, 
			const void *temp, Intersection &its) const {
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
	Float m_surfaceArea;
};

MTS_IMPLEMENT_CLASS_S(Rectangle, false, Shape)
MTS_EXPORT_PLUGIN(Rectangle, "Rectangle intersection primitive");
MTS_NAMESPACE_END
