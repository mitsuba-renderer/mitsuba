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
#include <mitsuba/render/gkdtree.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN


/**
 * \brief Acceleration structure for cylindrical hair
 * segments with miter joins.
 */
class HairKDTree : public GenericKDTree<HairKDTree> {
	friend class GenericKDTree<HairKDTree>;
public:
	HairKDTree(std::vector<Point> &vertices, 
			std::vector<bool> &startFiber, Float radius)
			: m_radius(radius) {
		// Take the vertex & start fiber arrays (without copying)
		m_vertices.swap(vertices);
		m_startFiber.swap(startFiber);

		// Compute the index of the first vertex in each segment.
		m_segIndex.reserve(m_vertices.size());
		for (size_t i=0; i<m_vertices.size()-1; i++)
			if (!m_startFiber[i+1])
				m_segIndex.push_back(i);

		Log(EDebug, "Building a kd-tree for " SIZE_T_FMT " hair vertices, "
			SIZE_T_FMT " segments,", m_vertices.size(), m_segIndex.size());

		setStopPrims(0);
		buildInternal();
	}

	inline const AABB &getAABB() const {
		return m_aabb;
	}

	inline const std::vector<Point> &getVertices() const {
		return m_vertices;
	}

	inline const std::vector<bool> &getStartFiber() const {
		return m_startFiber;
	}

	inline Float getRadius() const {
		return m_radius;
	}
protected:
	AABB getAABB(int index) const {
		uint32_t iv = m_segIndex[index];

		// cosine of steepest miter angle
		const Float cos0 = dot(firstMiterNormal(iv), tangent(iv));
		const Float cos1 = dot(secondMiterNormal(iv), tangent(iv));
		const Float maxInvCos = 1.0 / std::min(cos0, cos1);
		const Vector expandVec(m_radius * maxInvCos);

		const Point a = m_vertices[iv];
		const Point b = m_vertices[iv+1];
		AABB aabb;
		aabb.expandBy(a - expandVec);
		aabb.expandBy(a + expandVec);
		aabb.expandBy(b - expandVec);
		aabb.expandBy(b + expandVec);
		return aabb;
	}

	AABB getClippedAABB(int index, const AABB &box) const {
		AABB aabb(getAABB(index));
		aabb.clip(box);
		return aabb;
	}

	inline int getPrimitiveCount() const {
		return m_segIndex.size();
	}

protected:
	/* Some utility functions */
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
			return tangent(iv);
	}

	EIntersectionResult intersect(const Ray &ray, index_type idx, 
		Float mint, Float maxt, Float &t, void *tmp) {
		return ENo;
	}

	MTS_DECLARE_CLASS()
protected:
	std::vector<Point> m_vertices;
	std::vector<bool> m_startFiber;
	std::vector<uint32_t> m_segIndex;
	Float m_radius;
};

class Hair : public Shape {
public:
	Hair(const Properties &props) : Shape(props) {
		fs::path path = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		Float radius = (Float) props.getFloat("radius", 0.05f);

		Log(EInfo, "Loading hair geometry from \"%s\" ..", path.leaf().c_str());

		fs::ifstream is(path);
		if (is.fail())
			Log(EError, "Could not open \"%s\"!", path.file_string().c_str());

		std::string line;
		bool newFiber = true;
		Point p;
		std::vector<Point> vertices;
		std::vector<bool> startFiber;

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
					vertices.push_back(p);
					startFiber.push_back(newFiber);
					newFiber = false;
				}
			}
		}

		startFiber.push_back(true);
		m_kdtree = new HairKDTree(vertices, startFiber, radius);
	}


	Hair(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		Float radius = stream->readFloat();
		size_t vertexCount = (size_t) stream->readUInt();

		std::vector<Point> vertices(vertexCount);
		std::vector<bool> startFiber(vertexCount+1);
		stream->readFloatArray((Float *) &vertices[0], vertexCount * 3);

		for (size_t i=0; i<vertexCount; ++i) 
			startFiber[i] = stream->readBool();
		startFiber[vertexCount] = true;

		m_kdtree = new HairKDTree(vertices, startFiber, radius);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		const std::vector<Point> &vertices = m_kdtree->getVertices();
		const std::vector<bool> &startFiber = m_kdtree->getStartFiber();

		stream->writeFloat(m_kdtree->getRadius());
		stream->writeUInt((uint32_t) vertices.size());
		stream->writeFloatArray((Float *) &vertices[0], vertices.size() * 3);
		for (size_t i=0; i<vertices.size(); ++i)
			stream->writeBool(startFiber[i]);
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		return false;
	}


	void fillIntersectionRecord(const Ray &ray, Float t, 
		const void *temp, Intersection &its) const {
	}

	AABB getAABB() const {
		return m_kdtree->getAABB();
	}

	Float getSurfaceArea() const {
		Log(EError, "Hair::getSurfaceArea(): Not implemented.");
		return -1;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Hair[" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<HairKDTree> m_kdtree;
};

MTS_IMPLEMENT_CLASS(HairKDTree, false, GenericKDTree)
MTS_IMPLEMENT_CLASS_S(Hair, false, Shape)
MTS_EXPORT_PLUGIN(Hair, "Hair intersection primitive");
MTS_NAMESPACE_END
