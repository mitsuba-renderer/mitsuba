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
 * \brief Space-efficient acceleration structure for cylindrical hair
 * segments with miter joints.
 */
class HairKDTree : public GenericKDTree<HairKDTree> {
	friend class GenericKDTree<HairKDTree>;
public:
	HairKDTree(std::vector<Point> &vertices, 
			std::vector<bool> &vertexStartsFiber, Float radius)
			: m_radius(radius) {
		/* Take the supplied vertex & start fiber arrays (without copying) */
		m_vertices.swap(vertices);
		m_vertexStartsFiber.swap(vertexStartsFiber);

		/* Compute the index of the first vertex in each segment. */
		m_segIndex.reserve(m_vertices.size());
		for (size_t i=0; i<m_vertices.size()-1; i++)
			if (!m_vertexStartsFiber[i+1])
				m_segIndex.push_back(i);

		Log(EDebug, "Building a kd-tree for " SIZE_T_FMT " hair vertices, "
			SIZE_T_FMT " segments,", m_vertices.size(), m_segIndex.size());

		/* Ray-cylinder intersections are expensive. Use only the
		   SAH cost as the tree subdivision stopping criterion, 
		   not the number of primitives */
		setStopPrims(0);
		setTraversalCost(10);
		setIntersectionCost(30);
		buildInternal();

		/* Optimization: replace all primitive indices by the
		   associated vertex indices (this avoids an extra 
		   indirection during traversal later on) */
		for (size_type i=0; i<m_indexCount; ++i)
			m_indices[i] = m_segIndex[m_indices[i]];

		/* Free the segIndex array, it is not needed anymore */
		std::vector<index_type>().swap(m_segIndex);
	}

	/// Return the AABB of the hair kd-tree
	inline const AABB &getAABB() const {
		return m_aabb;
	}

	/// Return the list of vertices underlying the hair kd-tree
	inline const std::vector<Point> &getVertices() const {
		return m_vertices;
	}

	/**
	 * Return a boolean list specifying whether a vertex 
	 * marks the beginning of a new fiber
	 */
	inline const std::vector<bool> &getStartFiber() const {
		return m_vertexStartsFiber;
	}

	/// Return the radius of the hairs stored in the kd-tree
	inline Float getRadius() const {
		return m_radius;
	}

	/// Intersect a ray with all segments stored in the kd-tree
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt, 
			Float &t, void *temp) const {
		Float tempT = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<false>(ray, mint, maxt, tempT, temp)) {
					t = tempT;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * \brief Intersect a ray with all segments stored in the kd-tree
	 * (Visiblity query version)
	 */
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt) const {
		Float tempT = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL)) 
					return true;
			}
		}
		return false;
	}

	/// Compute the AABB of a segment (only used during tree construction)
	AABB getAABB(int index) const {
		index_type iv = m_segIndex.at(index);

		// cosine of steepest miter angle
		const Float cos0 = dot(firstMiterNormal(iv), tangent(iv));
		const Float cos1 = dot(secondMiterNormal(iv), tangent(iv));
		const Float maxInvCos = 1.0 / std::min(cos0, cos1);
		const Vector expandVec(m_radius * maxInvCos);

		const Point a = firstVertex(iv);
		const Point b = secondVertex(iv);

		AABB aabb;
		aabb.expandBy(a - expandVec);
		aabb.expandBy(a + expandVec);
		aabb.expandBy(b - expandVec);
		aabb.expandBy(b + expandVec);
		return aabb;
	}

	/// Compute the clipped AABB of a segment (only used during tree construction)
	AABB getClippedAABB(int index, const AABB &box) const {
		AABB aabb(getAABB(index));
		aabb.clip(box);
		return aabb;
	}

	/// Return the total number of segments
	inline int getPrimitiveCount() const {
		return m_segIndex.size();
	}

	inline EIntersectionResult intersect(const Ray &ray, index_type iv, 
		Float mint, Float maxt, Float &t, void *tmp) const {
		/* First compute the intersection with the infinite cylinder */
		Float nearT, farT;
		Vector axis = tangent(iv);

		// Projection of ray onto subspace normal to axis
		Vector relOrigin = ray.o - firstVertex(iv);
		Vector projOrigin = relOrigin - dot(axis, relOrigin) * axis;
		Vector projDirection = ray.d - dot(axis, ray.d) * axis;

		// Quadratic to intersect circle in projection
		const Float A = projDirection.lengthSquared();
		const Float B = 2 * dot(projOrigin, projDirection);
		const Float C = projOrigin.lengthSquared() - m_radius*m_radius;

		if (!solveQuadratic(A, B, C, nearT, farT))
			return ENever;

		if (nearT > maxt || farT < mint)
			return ENo;

		/* Next check the intersection points against the miter planes */
		Point pointNear = ray(nearT);
		Point pointFar = ray(farT);
		if (dot(pointNear - firstVertex(iv), firstMiterNormal(iv)) >= 0 &&
			dot(pointNear - secondVertex(iv), secondMiterNormal(iv)) <= 0 &&
			nearT >= mint) {
			t = nearT;
		} else if (dot(pointFar - firstVertex(iv), firstMiterNormal(iv)) >= 0 &&
				dot(pointFar - secondVertex(iv), secondMiterNormal(iv)) <= 0) {
			if (farT > maxt)
				return ENo;
			t = farT;
		} else {
			return ENo;
		}

		index_type *storage = static_cast<index_type *>(tmp);
		if (storage)
			*storage = iv;

		return EYes;
	}
	
	inline EIntersectionResult intersect(const Ray &ray, index_type iv, 
		Float mint, Float maxt) const {
		Float tempT;
		return intersect(ray, iv, mint, maxt, tempT, NULL);
	}

	/* Some utility functions */
	inline Point firstVertex(index_type iv) const { return m_vertices[iv]; }
	inline Point secondVertex(index_type iv) const { return m_vertices[iv+1]; }
	inline Point prevVertex(index_type iv) const { return m_vertices[iv-1]; }
	inline Point nextVertex(index_type iv) const { return m_vertices[iv+2]; }

	inline bool prevSegmentExists(index_type iv) const { return !m_vertexStartsFiber[iv]; }
	inline bool nextSegmentExists(index_type iv) const { return !m_vertexStartsFiber[iv+2]; }

	inline Vector tangent(index_type iv) const { return normalize(secondVertex(iv) - firstVertex(iv)); }
	inline Vector prevTangent(index_type iv) const { return normalize(firstVertex(iv) - prevVertex(iv)); }
	inline Vector nextTangent(index_type iv) const { return normalize(nextVertex(iv) - secondVertex(iv)); }

	inline Vector firstMiterNormal(index_type iv) const {
		if (prevSegmentExists(iv))
			return normalize(prevTangent(iv) + tangent(iv));
		else
			return tangent(iv);
	}

	inline Vector secondMiterNormal(index_type iv) const {
		if (nextSegmentExists(iv))
			return normalize(tangent(iv) + nextTangent(iv));
		else
			return tangent(iv);
	}

	MTS_DECLARE_CLASS()
protected:
	std::vector<Point> m_vertices;
	std::vector<bool> m_vertexStartsFiber;
	std::vector<index_type> m_segIndex;
	Float m_radius;
};

class Hair : public Shape {
public:
	Hair(const Properties &props) : Shape(props) {
		fs::path path = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		Float radius = (Float) props.getFloat("radius", 0.05f);

		/* Object-space -> World-space transformation */
		Transform objectToWorld = props.getTransform("toWorld", Transform());

		Log(EInfo, "Loading hair geometry from \"%s\" ..", path.leaf().c_str());

		fs::ifstream is(path);
		if (is.fail())
			Log(EError, "Could not open \"%s\"!", path.file_string().c_str());

		std::string line;
		bool newFiber = true;
		Point p, lastP(0.0f);
		std::vector<Point> vertices;
		std::vector<bool> vertexStartsFiber;
		size_t nDegenerate = 0;

		while (is.good()) {
			std::getline(is, line);
			if (line.length() > 0 && line[0] == '#') {
				newFiber = true;
				continue;
			}
			std::istringstream iss(line);
			iss >> p.x >> p.y >> p.z;
			if (!iss.fail()) {
				if (newFiber || p != lastP) {
					vertices.push_back(objectToWorld(p));
					vertexStartsFiber.push_back(newFiber);
					lastP = p;
				} else {
					nDegenerate++;
				}
				newFiber = false;
			} else {
				newFiber = true;
			}
		}

		if (nDegenerate > 0)
			Log(EInfo, "Encountered " SIZE_T_FMT 
				" degenerate segments!", nDegenerate);

		vertexStartsFiber.push_back(true);

		m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
	}

	Hair(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		Float radius = stream->readFloat();
		size_t vertexCount = (size_t) stream->readUInt();

		std::vector<Point> vertices(vertexCount);
		std::vector<bool> vertexStartsFiber(vertexCount+1);
		stream->readFloatArray((Float *) &vertices[0], vertexCount * 3);

		for (size_t i=0; i<vertexCount; ++i) 
			vertexStartsFiber[i] = stream->readBool();
		vertexStartsFiber[vertexCount] = true;

		m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		const std::vector<Point> &vertices = m_kdtree->getVertices();
		const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();

		stream->writeFloat(m_kdtree->getRadius());
		stream->writeUInt((uint32_t) vertices.size());
		stream->writeFloatArray((Float *) &vertices[0], vertices.size() * 3);
		for (size_t i=0; i<vertices.size(); ++i)
			stream->writeBool(vertexStartsFiber[i]);
	}

	bool rayIntersect(const Ray &ray, Float mint, 
			Float maxt, Float &t, void *temp) const {
		return m_kdtree->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		return m_kdtree->rayIntersect(ray, mint, maxt);
	}

	void fillIntersectionRecord(const Ray &ray, Float t, 
		const void *temp, Intersection &its) const {
		its.p = ray(t);

		/* No UV coordinates for now */
		its.uv = Point2(0,0);
		its.dpdu = Vector(0,0,0);
		its.dpdv = Vector(0,0,0);

		const HairKDTree::index_type *storage = 
			static_cast<const HairKDTree::index_type *>(temp);
		HairKDTree::index_type iv = *storage;

		its.geoFrame.s = m_kdtree->tangent(iv);
		const Vector relHitPoint = its.p - m_kdtree->firstVertex(iv);
		const Vector axis = m_kdtree->tangent(iv);
		its.geoFrame.n = Normal(relHitPoint - dot(axis, relHitPoint) * axis);
		its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.hasUVPartials = false;
		its.shape = this;
	}

	const AbstractKDTree *getKDTree() const {
		return m_kdtree.get();
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
