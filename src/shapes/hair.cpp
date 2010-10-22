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
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

class Hair : public Shape {
public:
	Hair(const Properties &props) : Shape(props) {
		fs::path path = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		m_radius = (Float) props.getFloat("radius", 0.05f);

		Log(EInfo, "Loading hair geometry from \"%s\" ..", path.leaf().c_str());

		fs::ifstream is(path);
		if (is.fail())
			Log(EError, "Could not open \"%s\"!", path.file_string().c_str());

		std::string line;
		bool newFiber = true;
		Point p;
		size_t segmentCount = 0;

		while (is.good()) {
			std::getline(is, line);
			if (line.length() > 0 && line[0] == '#')
				continue;
			if (line.length() == 0) {
				newFiber = true;
			} else {
				std::istringstream iss(line);
				iss >> p.x >> p.y >> p.z;
				if (!iss.fail()) {
					if (newFiber)
						segmentCount++;
					m_vertices.push_back(p);
					m_startFiber.push_back(newFiber);
					newFiber = false;
				}
			}
		}
		m_startFiber.push_back(true);

		for (size_t i=0; i<m_vertices.size()-1; ++i) {
			const Point a = firstVertex();
			const Point b = secondVertex();

			// cosine of steepest miter angle
			const Float cos0 = dot(firstMiterNormal(), tangent());
			const Float cos1 = dot(secondMiterNormal(), tangent());
			const Float maxInvCos = 1.0 / std::min(cos0, cos1);
			const Vector expandVec(m_radius * maxIvCos);

			result.expandBy(a - expandVec);
			result.expandBy(a + expandVec);
			result.expandBy(b - expandVec);
			result.expandBy(b + expandVec);

		}

		Log(EDebug, "Read " SIZE_T_FMT " hair vertices, " SIZE_T_FMT " segments,", 
			m_vertices.size(), segmentCount);
	}

	Hair(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
	}
	
	inline Vector tangent(int iv) const {
		return normalize(m_vertices[iv+1] - m_vertices[iv]);
	}

	inline Vector firstMiterNormal(int iv) const {
		if (!m_startFiber[iv])
			return normalize(tangent(iv - 1) + tangent(iv));
		else
			return tangent(iv);
	}

	inline Vector secondMiterNormal(int iv) const {
		if (!m_startFiber[iv+2])
			return normalize(tangent(iv) + tangent(iv+1));
		else
			return tangent();
	}


	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		return false;
	}


	void fillIntersectionRecord(const Ray &ray, Float t, 
			const void *temp, Intersection &its) const {
	}

	inline AABB getAABB() const {
		return m_aabb;
	}

	Float getSurfaceArea() const {
		Log(EError, "Hair::getSurfaceArea(): Not implemented.");
		return -1;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Hair[" << endl
			<< "  radius = " << m_radius << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<bool> m_startFiber;
	std::vector<int> m_segIndex;
	std::vector<Point> m_vertices;
	AABB m_aabb;
	Float m_radius;
};

MTS_IMPLEMENT_CLASS_S(Hair, false, Shape)
MTS_EXPORT_PLUGIN(Hair, "Hair intersection primitive");
MTS_NAMESPACE_END
