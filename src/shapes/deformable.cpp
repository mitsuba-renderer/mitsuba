/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
#include <mitsuba/render/sahkdtree4.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/mmap.h>

#define SHAPE_PER_SEGMENT 1
#define NO_CLIPPING_SUPPORT 1

MTS_NAMESPACE_BEGIN

class SpaceTimeKDTree : public SAHKDTree4D<SpaceTimeKDTree> {
	friend class GenericKDTree<AABB4, SurfaceAreaHeuristic4, SpaceTimeKDTree>;
	friend class SAHKDTree4D<SpaceTimeKDTree>;
public:
	/// Temporarily holds some intersection information
	struct IntersectionCache {
		Point p[3];
		Float u, v;
	};

	SpaceTimeKDTree(const std::vector<Float> &frameTimes, std::vector<float *> &positions,
			Triangle *triangles, size_t vertexCount,  size_t triangleCount)
		: m_frameTimes(frameTimes), m_positions(positions), m_triangles(triangles),
		  m_vertexCount(vertexCount), m_triangleCount(triangleCount) {

		Log(EInfo, "Total amount of vertex data: %s",
				memString(vertexCount*frameTimes.size()*sizeof(float)*3).c_str());

		//setClip(false);
		//setExactPrimitiveThreshold(10);
		buildInternal();

		/* Collect some statistics */
		std::stack<const KDNode *> stack;

		stack.push(m_nodes);
		size_t spatialSplits = 0, timeSplits = 0;
		while (!stack.empty()) {
			const KDNode *node = stack.top();
			stack.pop();
			if (!node->isLeaf()) {
				if (node->getAxis() == 3) {
					timeSplits++;
				} else {
					spatialSplits++;
				}
				stack.push((const KDNode *) node->getLeft());
				stack.push((const KDNode *) node->getRight());
			}
		}

		KDLog(EInfo, "Spacetime kd-tree statistics");
		KDLog(EInfo, "  Time interval  = [%f, %f]" , m_tightAABB.min.w, m_tightAABB.max.w);
		KDLog(EInfo, "  Spatial splits = " SIZE_T_FMT, spatialSplits);
		KDLog(EInfo, "  Time splits    = " SIZE_T_FMT, timeSplits);
		KDLog(EInfo, "");

		m_spatialAABB = AABB(
			Point(m_aabb.min.x, m_aabb.min.y, m_aabb.min.z),
			Point(m_aabb.max.x, m_aabb.max.y, m_aabb.max.z)
		);
	}

	/// Return one of the points stored in the point cache
	inline Point getPoint(uint32_t frame, uint32_t index) const {
		float *ptr = m_positions[frame] + index*3;
#if defined(__LITTLE_ENDIAN__)
		return Point(
			(Float) endianness_swap(ptr[0]),
			(Float) endianness_swap(ptr[1]),
			(Float) endianness_swap(ptr[2]));
#else
		return Point((Float) ptr[0], (Float) ptr[1], (Float) ptr[2]);
#endif
	}

	// ========================================================================
	//    Implementation of functions required by the parent class
	// ========================================================================

	/// Return the total number of primitives that are organized in the tree
	inline SizeType getPrimitiveCount() const {
#ifdef SHAPE_PER_SEGMENT
		return m_triangleCount * (m_frameTimes.size() - 1);
#else
		return m_triangleCount;
#endif
	}

	/// Return the 4D extents for one of the primitives contained in the tree
	AABB4 getAABB(IndexType index) const {
#ifdef SHAPE_PER_SEGMENT
		int frameIdx = index / m_triangleCount;
		int triangleIdx  = index % m_triangleCount;
		const Triangle &tri = m_triangles[triangleIdx];

		AABB aabb;
		for (int i=0; i<3; ++i) {
			aabb.expandBy(getPoint(frameIdx, tri.idx[i]));
			aabb.expandBy(getPoint(frameIdx+1, tri.idx[i]));
		}

		return AABB4(
			Point4(aabb.min.x, aabb.min.y, aabb.min.z, m_frameTimes[frameIdx]),
			Point4(aabb.max.x, aabb.max.y, aabb.max.z, m_frameTimes[frameIdx+1])
		);
#else
		AABB aabb;
		const Triangle &tri = m_triangles[index];
		for (size_t i=0; i<m_frameTimes.size(); ++i)
			for (int j=0; j<3; ++j)
				aabb.expandBy(getPoint(i, tri.idx[j]));
		return AABB4(
			Point4(aabb.min.x, aabb.min.y, aabb.min.z, m_frameTimes[0]),
			Point4(aabb.max.x, aabb.max.y, aabb.max.z, m_frameTimes[m_frameTimes.size()-1])
		);
#endif
	}

	/// Return a clipped 4D AABB for one of the primitives contained in the tree
	AABB4 getClippedAABB(int index, const AABB4 &box) const {
		AABB clip(
			Point(box.min.x, box.min.y, box.min.z),
			Point(box.max.x, box.max.y, box.max.z)
		);
#ifdef NO_CLIPPING_SUPPORT
		AABB4 aabb = getAABB(index);
		aabb.clip(box);
		return aabb;
#elif SHAPE_PER_SEGMENT
		int frameIdx = index / m_triangleCount;
		int triangleIdx  = index % m_triangleCount;

		AABB aabb(m_triangles[triangleIdx].getClippedAABB(m_positions[frameIdx], clip)); /// XXX broken
		aabb.expandBy(m_triangles[triangleIdx].getClippedAABB(m_positions[frameIdx+1], clip));
		if (aabb.isValid())
			return AABB4(
				Point4(aabb.min.x, aabb.min.y, aabb.min.z, box.min.w),
				Point4(aabb.max.x, aabb.max.y, aabb.max.z, box.max.w));
		else
			return AABB4();
#else
		int startIndex = std::max((int) (std::lower_bound(m_frameTimes.begin(), m_frameTimes.end(),
				box.min.w) - m_frameTimes.begin()) - 1, 0);
		int endIndex = (int) (std::lower_bound(m_frameTimes.begin(), m_frameTimes.end(),
				box.max.w) - m_frameTimes.begin());
		AABB4 result;
		const Triangle &tri = m_triangles[index];

		for (int i=startIndex; i<=endIndex; ++i) {
			Point p0 = getPoint(i, tri.idx[0]);
			Point p1 = getPoint(i, tri.idx[1]);
			Point p2 = getPoint(i, tri.idx[2]);
			AABB aabb(Triangle::getClippedAABB(p0, p1, p2, clip));
			if (aabb.isValid()) {
				result.expandBy(Point4(aabb.min.x, aabb.min.y, aabb.min.z, m_frameTimes[i]));
				result.expandBy(Point4(aabb.max.x, aabb.max.y, aabb.max.z, m_frameTimes[i]));
			}
		}
		result.clip(box);
		return result;
#endif
	}

	/// Cast a normal (i.e. non-shadow) ray against a specific animated triangle
	inline bool intersect(const Ray &ray, IndexType idx,
			Float mint, Float maxt, Float &t, void *tmp) const {
#if SHAPE_PER_SEGMENT
		IndexType frameIdx = idx / m_triangleCount;
		IndexType triangleIdx = idx % m_triangleCount;
#else
		IndexType triangleIdx = idx;
		IndexType frameIdx = (IndexType) std::max((int) (std::lower_bound(
			m_frameTimes.begin(), m_frameTimes.end(), ray.time) -
			m_frameTimes.begin()) - 1, 0);
#endif
		const Triangle &tri = m_triangles[triangleIdx];

		Float alpha = (ray.time - m_frameTimes[frameIdx])
			/ (m_frameTimes[frameIdx + 1] - m_frameTimes[frameIdx]);

		if (alpha < 0 || alpha > 1)
			return false;

		/* Compute interpolated positions */
		Point p[3];
		for (int i=0; i<3; ++i)
			p[i] = (1 - alpha) * getPoint(frameIdx, tri.idx[i])
				+ alpha * getPoint(frameIdx+1, tri.idx[i]);

		Float tempU, tempV, tempT;
		if (!Triangle::rayIntersect(p[0], p[1], p[2], ray, tempU, tempV, tempT))
			return false;
		if (tempT < mint || tempT > maxt)
			return false;

		if (tmp != NULL) {
			IntersectionCache *cache =
				static_cast<IntersectionCache *>(tmp);
			t = tempT;
			memcpy(cache->p, p, sizeof(Point)*3);
			cache->u = tempU;
			cache->v = tempV;
		}
		return true;
	}

	/// Cast a shadow ray against a specific triangle
	inline bool intersect(const Ray &ray, IndexType idx,
			Float mint, Float maxt) const {
		Float tempT;
		/* No optimized version for shadow rays yet */
		return intersect(ray, idx, mint, maxt, tempT, NULL);
	}

	// ========================================================================
	//   Miscellaneous
	// ========================================================================

	/// Intersect a ray with all primitives stored in the kd-tree
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt,
			Float &t, void *temp) const {
		Float tempT = std::numeric_limits<Float>::infinity();
		Float mint, maxt;

		if (m_spatialAABB.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint && ray.time >= m_aabb.min.w && ray.time <= m_aabb.max.w)) {
				if (rayIntersectHavran<false>(ray, mint, maxt, tempT, temp)) {
					t = tempT;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * \brief Intersect a ray with all primitives stored in the kd-tree
	 * (Visiblity query version)
	 */
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt) const {
		Float tempT = std::numeric_limits<Float>::infinity();
		Float mint, maxt;

		if (m_spatialAABB.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint && ray.time >= m_aabb.min.w && ray.time <= m_aabb.max.w))
				if (rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL))
					return true;
		}
		return false;
	}

	inline const Triangle *getTriangles() const {
		return m_triangles;
	}

	/// Return an AABB with the spatial extents
	inline const AABB &getSpatialAABB() const {
		return m_spatialAABB;
	}

	MTS_DECLARE_CLASS()
protected:
	std::vector<Float> m_frameTimes;
	std::vector<float *> m_positions;
	Triangle *m_triangles;
	size_t m_vertexCount;
	size_t m_triangleCount;
	AABB m_spatialAABB;
};

class Deformable : public Shape {
public:
	Deformable(const Properties &props) : Shape(props) {
		FileResolver *fResolver = Thread::getThread()->getFileResolver();
		fs::path path = fResolver->resolve(props.getString("filename"));
		if (path.extension() != ".mdd")
			Log(EError, "Point cache files must have the extension \".mdd\"");

		m_mmap = new MemoryMappedFile(path);

		ref<MemoryStream> mStream = new MemoryStream((uint8_t *) m_mmap->getData(),
				m_mmap->getSize());
		mStream->setByteOrder(Stream::EBigEndian);

		uint32_t frameCount = mStream->readUInt();
		m_vertexCount = mStream->readUInt();

		Log(EInfo, "Point cache has %i frames and %i vertices", frameCount, m_vertexCount);

		Float clipStart = props.getFloat("clipStart", 0),
			  clipEnd   = props.getFloat("clipEnd", 0);

		std::vector<Float> frameTimes;
		std::vector<float *> positions;

		for (uint32_t i=0; i<frameCount; ++i)
			frameTimes.push_back((Float) mStream->readSingle());

		for (uint32_t i=0; i<frameCount; ++i) {
			positions.push_back(reinterpret_cast<float *>(mStream->getCurrentData()));
			mStream->skip(m_vertexCount * 3 * sizeof(float));
		}

		if (clipStart != clipEnd) {
			m_positions.reserve(positions.size());
			m_frameTimes.reserve(frameTimes.size());
			for (uint32_t i=0; i<frameCount; ++i) {
				if (frameTimes[i] >= clipStart && frameTimes[i] <= clipEnd) {
					m_frameTimes.push_back(frameTimes[i]);
					m_positions.push_back(positions[i]);
				}
			}
			if (m_frameTimes.empty())
				Log(EError, "After clipping to the time range [%f, %f] no frames were left!",
					clipStart, clipEnd);
			Log(EInfo, "Clipped away %u/%u frames", frameCount - (uint32_t) m_frameTimes.size(), frameCount);
		} else {
			m_positions = positions;
			m_frameTimes = frameTimes;
		}
	}

	Deformable(Stream *stream, InstanceManager *manager)
		: Shape(stream, manager) {
		/// TBD
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		/// TBD
	}

	void configure() {
		Shape::configure();

		if (m_mesh == NULL)
			Log(EError, "A nested triangle mesh is required so that "
				"connectivity information can be extracted!");
		if (m_mesh->getVertexCount() != m_vertexCount)
			Log(EError, "Point cache and nested geometry have mismatched "
				"numbers of vertices!");

		m_kdtree = new SpaceTimeKDTree(m_frameTimes, m_positions, m_mesh->getTriangles(),
				m_vertexCount, m_mesh->getTriangleCount());
		m_aabb = m_kdtree->getSpatialAABB();
	}

	bool rayIntersect(const Ray &ray, Float mint,
			Float maxt, Float &t, void *temp) const {
		return m_kdtree->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		return m_kdtree->rayIntersect(ray, mint, maxt);
	}

	void fillIntersectionRecord(const Ray &ray,
			const void *temp, Intersection &its) const {
		const SpaceTimeKDTree::IntersectionCache *cache
			= reinterpret_cast<const SpaceTimeKDTree::IntersectionCache *>(temp);

		const Vector b(1 - cache->u - cache->v, cache->u, cache->v);
		const Point p0 = cache->p[0];
		const Point p1 = cache->p[1];
		const Point p2 = cache->p[2];

		Normal faceNormal(cross(p1-p0, p2-p0));
		Float length = faceNormal.length();
		if (!faceNormal.isZero())
			faceNormal /= length;

		/* Just the basic attributes for now and geometric normals */
		its.p = ray(its.t);
		its.geoFrame = Frame(faceNormal);
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.shape = this;
		its.instance = this;
		its.hasUVPartials = false;
		its.time = ray.time;
	}

	AABB getAABB() const {
		return m_kdtree->getSpatialAABB();
	}

	size_t getPrimitiveCount() const {
		return m_mesh->getTriangleCount();
	}

	size_t getEffectivePrimitiveCount() const {
		return m_mesh->getTriangleCount();
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(TriMesh::m_theClass)) {
			Assert(m_mesh == NULL);
			m_mesh = static_cast<TriMesh *>(child);
			if (m_mesh->getVertexCount() != m_vertexCount)
				Log(EError, "Geometry mismatch! MDD file contains %u vertices. "
					"The attached shape uses %u!", m_vertexCount, m_mesh->getVertexCount());
		} else if (cClass->derivesFrom(Shape::m_theClass) && static_cast<Shape *>(child)->isCompound()) {
			size_t index = 0;
			Shape *shape = static_cast<Shape *>(child);
			do {
				ref<Shape> element = shape->getElement(index++);
				if (element == NULL)
					break;
				addChild(name, element);
			} while (true);
		} else {
			Shape::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Deformable[" << endl
			<< "  mesh = " << indent(m_mesh.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<MemoryMappedFile> m_mmap;
	ref<SpaceTimeKDTree> m_kdtree;
	std::vector<Float> m_frameTimes;
	std::vector<float *> m_positions;
	ref<TriMesh> m_mesh;
	uint32_t m_vertexCount;
	AABB m_aabb;
};

MTS_IMPLEMENT_CLASS(SpaceTimeKDTree, false, KDTreeBase)
MTS_IMPLEMENT_CLASS_S(Deformable, false, Shape)
MTS_EXPORT_PLUGIN(Deformable, "Deformable shape");
MTS_NAMESPACE_END
