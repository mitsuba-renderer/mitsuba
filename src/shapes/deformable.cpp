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
        IndexType frameIndex;
        Float alpha;
        IndexType shapeIndex, primIndex;
        Float u, v;
    };

    SpaceTimeKDTree(const std::vector<Float> &times) : m_times(times) { }

    SpaceTimeKDTree(Stream *stream, InstanceManager *manager) {
        size_t times = (size_t) stream->readUInt();
        m_times.resize(times);
        m_meshes.resize(times);
        for (size_t i=0; i<times; ++i) {
            m_times[i] = stream->readFloat();
            size_t count = (size_t) stream->readUInt();
            std::vector<const TriMesh *> &meshes = m_meshes.at(i);
            meshes.resize(count);
            for (size_t j=0; j<count; ++j)
                meshes[j] = static_cast<TriMesh *>(manager->getInstance(stream));
        }
    }

    ~SpaceTimeKDTree() {
        for (size_t i=0; i<m_meshes.size(); ++i)
            for (size_t j=0; j<m_meshes[i].size(); ++j)
                m_meshes[i][j]->decRef();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        stream->writeUInt((uint32_t) m_times.size());

        for (size_t i=0; i<m_times.size(); ++i) {
            const std::vector<const TriMesh *> &meshes = m_meshes.at(i);
            stream->writeFloat(m_times[i]);
            stream->writeUInt((uint32_t) meshes.size());
            for (size_t j=0; j<meshes.size(); ++j)
                manager->serialize(stream, meshes[j]);
        }
    }

    void addShape(Shape *shape) {
        std::vector<const TriMesh *> vec;
        size_t index = 0;
        if (shape->getClass()->derivesFrom(MTS_CLASS(TriMesh))) {
            shape->incRef();
            vec.push_back(static_cast<TriMesh *>(shape));
        } else if (shape->isCompound()) {
            do {
                ref<Shape> element = shape->getElement(index++);
                if (element == NULL)
                    break;
                if (!element->getClass()->derivesFrom(MTS_CLASS(TriMesh)))
                    Log(EError, "Can only add triangle meshes to the 'deformable' plugin");
                element->incRef();
                vec.push_back(static_cast<TriMesh *>(element.get()));
            } while (true);
        }

        if (vec.empty())
            Log(EError, "Can only add triangle meshes to the 'deformable' plugin");
        else
            m_meshes.push_back(vec);
    }


    void build() {
        if (m_meshes.size() < 2)
            Log(EError, "The deformable shape requires at least two sub-shapes!");

        if (m_meshes.size() != m_times.size()) {
            Log(EError, "The number of arguments to the 'times' parameter (%u) must "
                "match the number of sub-shapes (%u).", m_times.size(), m_meshes.size());
        }

        for (size_t i=1; i<m_meshes.size(); ++i) {
            const std::vector<const TriMesh *> &meshes = m_meshes[i];

            if (m_times[i] <= m_times[i-1])
                Log(EError, "Frame times must be increasing!");

            if (meshes.size() != m_meshes[0].size())
                Log(EError, "The number of compound shapes for each time value must be identical!");

            for (size_t j=0;j<m_meshes[0].size(); ++j) {
                const TriMesh *mesh0 = m_meshes[0][j];
                const TriMesh *mesh1 = m_meshes[i][j];

                if (mesh0->getTriangleCount() != mesh1->getTriangleCount())
                    Log(EError, "All sub-meshes must have the exact same number of triangles");
                if (mesh0->getVertexCount() != mesh1->getVertexCount())
                    Log(EError, "All sub-meshes must have the exact same number of triangles");
                if (memcmp(mesh0->getTriangles(), mesh1->getTriangles(), sizeof(Triangle) * mesh0->getTriangleCount()) != 0)
                    Log(EError, "All sub-meshes must have the exact same face topology");
            }
        }


        m_shapeMap.resize(m_meshes[0].size()+1);
        m_shapeMap[0] = 0;
        for (size_t i=0; i<m_meshes[0].size(); ++i)
            m_shapeMap[i+1] = m_shapeMap[i] + (SizeType) m_meshes[0][i]->getTriangleCount();

        this->setClip(false);
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

    inline IndexType findShape(IndexType &index) const {
        std::vector<IndexType>::const_iterator it = std::lower_bound(
                m_shapeMap.begin(), m_shapeMap.end(), index + 1) - 1;
        index -= *it;
        return (IndexType) (it - m_shapeMap.begin());
    }

    inline IndexType findFrame(Float time) const {
        return (IndexType) std::min(std::max((int) (std::lower_bound(
            m_times.begin(), m_times.end(), time) - m_times.begin()) - 1, 0), (int) m_times.size()-1);
    }

    // ========================================================================
    //    Implementation of functions required by the parent class
    // ========================================================================

    /// Return the total number of primitives that are organized in the tree
    inline SizeType getPrimitiveCount() const {
        return m_shapeMap[m_shapeMap.size()-1];
    }

    /// Return the 4D extents for one of the primitives contained in the tree
    AABB4 getAABB(IndexType index) const {
        IndexType shapeIndex = findShape(index);
        const Triangle &tri = m_meshes[0][shapeIndex]->getTriangles()[index];

        AABB aabb;
        for (size_t frame=0; frame<m_times.size(); ++frame) {
            const Point *pos = m_meshes[frame][shapeIndex]->getVertexPositions();
            for (int j=0; j<3; ++j)
                aabb.expandBy(pos[tri.idx[j]]);
        }

        return AABB4(
            Point4(aabb.min.x, aabb.min.y, aabb.min.z, m_times[0]),
            Point4(aabb.max.x, aabb.max.y, aabb.max.z, m_times[m_times.size()-1])
        );
    }

    /// Return a clipped 4D AABB for one of the primitives contained in the tree
    AABB4 getClippedAABB(IndexType index, const AABB4 &box) const {
        AABB clip(
            Point(box.min.x, box.min.y, box.min.z),
            Point(box.max.x, box.max.y, box.max.z)
        );
#if 0
        int startIndex = findFrame(box.min.w),
            endIndex   = findFrameUpperBound(box.max.w);

        AABB4 result;
        IndexType shapeIndex = findShape(index);
        const Triangle &tri = m_meshes[0][shapeIndex]->getTriangles()[index];

        for (int frame=startIndex; frame<=endIndex; ++frame) {
            const Point *pos = m_meshes[frame][shapeIndex]->getVertexPositions();
            AABB aabb = tri.getClippedAABB(pos, clip);
            if (aabb.isValid()) {
                result.expandBy(Point4(aabb.min.x, aabb.min.y, aabb.min.z, m_times[frame]));
                result.expandBy(Point4(aabb.max.x, aabb.max.y, aabb.max.z, m_times[frame]));
            }
        }
#endif

        AABB4 result = getAABB(index);
        result.clip(box);
        return result;
    }

    /// Cast a normal (i.e. non-shadow) ray against a specific animated triangle
    inline bool intersect(const Ray &ray, IndexType index,
            Float mint, Float maxt, Float &t, void *tmp) const {
        IntersectionCache *cache = static_cast<IntersectionCache *>(tmp);
        IndexType shapeIndex = findShape(index);
        IndexType frameIndex = findFrame(ray.time);
        Float alpha = std::max((Float) 0.0f, std::min((Float) 1.0f,
            (ray.time - m_times[frameIndex])
            / (m_times[frameIndex + 1] - m_times[frameIndex])));

        const Triangle &tri = m_meshes[0][shapeIndex]->getTriangles()[index];

        const Point *pos0 = m_meshes[frameIndex  ][shapeIndex]->getVertexPositions();
        const Point *pos1 = m_meshes[frameIndex+1][shapeIndex]->getVertexPositions();

        /* Compute interpolated positions */
        Point p[3];
        for (int i=0; i<3; ++i)
            p[i] = (1 - alpha) * pos0[tri.idx[i]] + alpha * pos1[tri.idx[i]];

        Float tempU, tempV, tempT;
        if (!Triangle::rayIntersect(p[0], p[1], p[2], ray, tempU, tempV, tempT))
            return false;

        if (tempT < mint || tempT > maxt)
            return false;

        t = tempT;
        cache->frameIndex = frameIndex;
        cache->alpha = alpha;
        cache->shapeIndex = shapeIndex;
        cache->primIndex = index;
        cache->u = tempU;
        cache->v = tempV;
        return true;
    }

    /// Cast a shadow ray against a specific triangle
    inline bool intersect(const Ray &ray, IndexType index, Float mint, Float maxt) const {
        IndexType shapeIndex = findShape(index);
        const Triangle &tri = m_meshes[0][shapeIndex]->getTriangles()[index];

        IndexType frameIndex = findFrame(ray.time);
        Float alpha = std::max((Float) 0.0f, std::min((Float) 1.0f,
            (ray.time - m_times[frameIndex])
            / (m_times[frameIndex + 1] - m_times[frameIndex])));

        const Point *pos0 = m_meshes[frameIndex  ][shapeIndex]->getVertexPositions();
        const Point *pos1 = m_meshes[frameIndex+1][shapeIndex]->getVertexPositions();

        /* Compute interpolated positions */
        Point p[3];
        for (int i=0; i<3; ++i)
            p[i] = (1 - alpha) * pos0[tri.idx[i]] + alpha * pos1[tri.idx[i]];

        Float tempU, tempV, tempT;
        if (!Triangle::rayIntersect(p[0], p[1], p[2], ray, tempU, tempV, tempT))
            return false;

        if (tempT < mint || tempT > maxt)
            return false;

        return true;
    }

    // ========================================================================
    //   Miscellaneous
    // ========================================================================

    /// Intersect a ray with all primitives stored in the kd-tree
    inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt,
            Float &t, void *temp) const {
        IntersectionCache *cache = static_cast<IntersectionCache *>(temp);

        Float tempT = std::numeric_limits<Float>::infinity();
        Float mint, maxt;

        if (m_spatialAABB.rayIntersect(ray, mint, maxt)) {
            if (_mint > mint) mint = _mint;
            if (_maxt < maxt) maxt = _maxt;

            if (EXPECT_TAKEN(maxt > mint)) {
                if (rayIntersectHavran<false>(ray, mint, maxt, tempT, cache)) {
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

            if (EXPECT_TAKEN(maxt > mint)) {
                if (rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL))
                    return true;
            }
        }
        return false;
    }

    /// Return an AABB with the spatial extents
    inline const AABB &getSpatialAABB() const {
        return m_spatialAABB;
    }

    /// Return the number of key-framed time values
    inline size_t getTimeCount() const {
        return m_times.size();
    }

    inline const std::vector<Float> getTimes() const {
        return m_times;
    }

    inline const TriMesh *getMesh(IndexType frameIndex, IndexType shapeIndex) const {
        return m_meshes[frameIndex][shapeIndex];
    }

    inline Triangle getTriangle(IndexType shapeIndex, IndexType primIndex) const {
        return m_meshes[0][shapeIndex]->getTriangles()[primIndex];
    }

    inline const std::vector<std::vector<const TriMesh *> > &getMeshes() const {
        return m_meshes;
    }

    MTS_DECLARE_CLASS()
protected:
    std::vector<Float> m_times;
    std::vector<std::vector<const TriMesh *> > m_meshes;
    std::vector<IndexType> m_shapeMap;
    AABB m_spatialAABB;
    Float m_traceTime;
};

class Deformable : public Shape {
public:
    Deformable(const Properties &props) : Shape(props) {
        std::vector<std::string> times_str =
            tokenize(props.getString("times", ""), " ,;");
        std::vector<Float> times(times_str.size());

        char *end_ptr = NULL;
        for (size_t i=0; i<times_str.size(); ++i) {
            Float value = (Float) strtod(times_str[i].c_str(), &end_ptr);
            if (*end_ptr != '\0')
                SLog(EError, "Could not parse the times parameter!");
            times[i] = value;
        }
        m_kdtree = new SpaceTimeKDTree(times);
    }

    Deformable(Stream *stream, InstanceManager *manager)
        : Shape(stream, manager) {
        m_kdtree = new SpaceTimeKDTree(stream, manager);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);
        m_kdtree->serialize(stream, manager);
    }

    void configure() {
        m_kdtree->build();
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
            = static_cast<const SpaceTimeKDTree::IntersectionCache *>(temp);
        const TriMesh *trimesh0 = m_kdtree->getMesh(cache->frameIndex,   cache->shapeIndex);
        const TriMesh *trimesh1 = m_kdtree->getMesh(cache->frameIndex+1, cache->shapeIndex);
        const Vector b(1 - cache->u - cache->v, cache->u, cache->v);
        const Triangle tri = m_kdtree->getTriangle(cache->shapeIndex, cache->primIndex);
        const uint32_t idx0 = tri.idx[0], idx1 = tri.idx[1], idx2 = tri.idx[2];
        const Float alpha = cache->alpha;

        const Point *vertexPositions0 = trimesh0->getVertexPositions();
        const Point *vertexPositions1 = trimesh1->getVertexPositions();
        const Normal *vertexNormals0 = trimesh0->getVertexNormals();
        const Normal *vertexNormals1 = trimesh1->getVertexNormals();
        const Point2 *vertexTexcoords0 = trimesh0->getVertexTexcoords();
        const Point2 *vertexTexcoords1 = trimesh1->getVertexTexcoords();
        const Color3 *vertexColors0 = trimesh0->getVertexColors();
        const Color3 *vertexColors1 = trimesh1->getVertexColors();
        const TangentSpace *vertexTangents0 = trimesh0->getUVTangents();
        const TangentSpace *vertexTangents1 = trimesh1->getUVTangents();

        const Point p0 = vertexPositions0[idx0] * (1-alpha) + vertexPositions1[idx0] * alpha;
        const Point p1 = vertexPositions0[idx1] * (1-alpha) + vertexPositions1[idx1] * alpha;
        const Point p2 = vertexPositions0[idx2] * (1-alpha) + vertexPositions1[idx2] * alpha;

        its.p = p0 * b.x + p1 * b.y + p2 * b.z;

        Vector side1(p1-p0), side2(p2-p0);
        Normal faceNormal(cross(side1, side2));
        Float length = faceNormal.length();
        if (!faceNormal.isZero())
            faceNormal /= length;

        if (EXPECT_NOT_TAKEN(vertexTangents0 && vertexTangents1)) {
            const TangentSpace &ts0 = vertexTangents0[cache->primIndex];
            const TangentSpace &ts1 = vertexTangents1[cache->primIndex];
            its.dpdu = (1-alpha) * ts0.dpdu + alpha * ts1.dpdu;
            its.dpdv = (1-alpha) * ts0.dpdv + alpha * ts1.dpdv;
        } else {
            its.dpdu = side1;
            its.dpdv = side2;
        }

        if (EXPECT_TAKEN(vertexNormals0)) {
            Normal
                n0 = (1-alpha) * vertexNormals0[idx0] + alpha * vertexNormals1[idx0],
                n1 = (1-alpha) * vertexNormals0[idx1] + alpha * vertexNormals1[idx1],
                n2 = (1-alpha) * vertexNormals0[idx2] + alpha * vertexNormals1[idx2];

            its.shFrame.n = normalize(n0 * b.x + n1 * b.y + n2 * b.z);

            /* Ensure that the geometric & shading normals face the same direction */
            if (dot(faceNormal, its.shFrame.n) < 0)
                faceNormal = -faceNormal;
        } else {
            its.shFrame.n = faceNormal;
        }
        its.geoFrame = Frame(faceNormal);

        if (EXPECT_TAKEN(vertexTexcoords0)) {
            Point2
                t0 = (1-alpha) * vertexTexcoords0[idx0] + alpha * vertexTexcoords1[idx0],
                t1 = (1-alpha) * vertexTexcoords0[idx1] + alpha * vertexTexcoords1[idx1],
                t2 = (1-alpha) * vertexTexcoords0[idx2] + alpha * vertexTexcoords1[idx2];
            its.uv = t0 * b.x + t1 * b.y + t2 * b.z;
        } else {
            its.uv = Point2(b.y, b.z);
        }

        if (EXPECT_NOT_TAKEN(vertexColors0)) {
            Color3
                c0 = (1-alpha) * vertexColors0[idx0] + alpha * vertexColors1[idx0],
                c1 = (1-alpha) * vertexColors0[idx1] + alpha * vertexColors1[idx1],
                c2 = (1-alpha) * vertexColors0[idx2] + alpha * vertexColors1[idx2];
            Color3 result(c0 * b.x + c1 * b.y + c2 * b.z);
            its.color.fromLinearRGB(result[0], result[1],
                result[2], Spectrum::EReflectance);
        }

        its.shape = m_kdtree->getMesh(0, cache->shapeIndex);
        its.hasUVPartials = false;
        its.primIndex = cache->primIndex;
        its.other = cache->shapeIndex;
        its.instance = this;
        its.time = ray.time;
    }

    void getNormalDerivative(const Intersection &its,
            Vector &dndu, Vector &dndv, bool shadingFrame) const {

        const std::vector<Float> &times = m_kdtree->getTimes();
        int frameIndex = m_kdtree->findFrame(its.time);
        Float alpha = std::max((Float) 0.0f, std::min((Float) 1.0f,
            (its.time - times[frameIndex])
            / (times[frameIndex + 1] - times[frameIndex])));

        uint32_t primIndex = its.primIndex, shapeIndex = its.other;
        const TriMesh *trimesh0 = m_kdtree->getMesh(frameIndex,   shapeIndex);
        const TriMesh *trimesh1 = m_kdtree->getMesh(frameIndex+1, shapeIndex);
        const Point *vertexPositions0 = trimesh0->getVertexPositions();
        const Point *vertexPositions1 = trimesh1->getVertexPositions();
        const Point2 *vertexTexcoords0 = trimesh0->getVertexTexcoords();
        const Point2 *vertexTexcoords1 = trimesh1->getVertexTexcoords();
        const Normal *vertexNormals0 = trimesh0->getVertexNormals();
        const Normal *vertexNormals1 = trimesh1->getVertexNormals();

        if (!vertexNormals0 || !vertexNormals1) {
            dndu = dndv = Vector(0.0f);
        } else {
            const Triangle &tri = trimesh0->getTriangles()[primIndex];
            uint32_t idx0 = tri.idx[0],
                     idx1 = tri.idx[1],
                     idx2 = tri.idx[2];

            const Point
                p0 = (1-alpha)*vertexPositions0[idx0] + alpha*vertexPositions1[idx0],
                p1 = (1-alpha)*vertexPositions0[idx1] + alpha*vertexPositions1[idx1],
                p2 = (1-alpha)*vertexPositions0[idx2] + alpha*vertexPositions1[idx2];

            /* Recompute the barycentric coordinates, since 'its.uv' may have been
               overwritten with coordinates of the texture "parameterization". */
            Vector rel = its.p - p0, du = p1 - p0, dv = p2 - p0;

            Float b1  = dot(du, rel), b2 = dot(dv, rel), /* Normal equations */
                  a11 = dot(du, du), a12 = dot(du, dv),
                  a22 = dot(dv, dv),
                  det = a11 * a22 - a12 * a12;

            if (det == 0) {
                dndu = dndv = Vector(0.0f);
                return;
            }

            Float invDet = 1.0f / det,
                  u = ( a22 * b1 - a12 * b2) * invDet,
                  v = (-a12 * b1 + a11 * b2) * invDet,
                  w = 1 - u - v;

            const Normal
                n0 = normalize((1-alpha)*vertexNormals0[idx0] + alpha*vertexNormals1[idx0]),
                n1 = normalize((1-alpha)*vertexNormals0[idx1] + alpha*vertexNormals1[idx1]),
                n2 = normalize((1-alpha)*vertexNormals0[idx2] + alpha*vertexNormals1[idx2]);

            /* Now compute the derivative of "normalize(u*n1 + v*n2 + (1-u-v)*n0)"
               with respect to [u, v] in the local triangle parameterization.

               Since d/du [f(u)/|f(u)|] = [d/du f(u)]/|f(u)|
                 - f(u)/|f(u)|^3 <f(u), d/du f(u)>, this results in
            */

            Normal N(u * n1 + v * n2 + w * n0);
            Float il = 1.0f / N.length(); N *= il;

            dndu = (n1 - n0) * il; dndu -= N * dot(N, dndu);
            dndv = (n2 - n0) * il; dndv -= N * dot(N, dndv);

            if (vertexTexcoords0 && vertexTexcoords1) {
                /* Compute derivatives with respect to a specified texture
                   UV parameterization.  */
                const Point2
                    uv0 = (1-alpha)*vertexTexcoords0[idx0] + alpha*vertexTexcoords1[idx0],
                    uv1 = (1-alpha)*vertexTexcoords0[idx1] + alpha*vertexTexcoords1[idx1],
                    uv2 = (1-alpha)*vertexTexcoords0[idx2] + alpha*vertexTexcoords1[idx2];

                Vector2 duv1 = uv1 - uv0, duv2 = uv2 - uv0;

                det = duv1.x * duv2.y - duv1.y * duv2.x;

                if (det == 0) {
                    dndu = dndv = Vector(0.0f);
                    return;
                }

                invDet = 1.0f / det;
                Vector dndu_ = ( duv2.y * dndu - duv1.y * dndv) * invDet;
                Vector dndv_ = (-duv2.x * dndu + duv1.x * dndv) * invDet;
                dndu = dndu_; dndv = dndv_;
            }
        }
    }


    void adjustTime(Intersection &its, Float time) const {
        SpaceTimeKDTree::IntersectionCache cache;

        const std::vector<Float> &times = m_kdtree->getTimes();

        cache.primIndex = its.primIndex;
        cache.shapeIndex = its.other;
        cache.frameIndex = m_kdtree->findFrame(its.time);
        cache.alpha = std::max((Float) 0.0f, std::min((Float) 1.0f,
            (its.time - times[cache.frameIndex])
            / (times[cache.frameIndex + 1] - times[cache.frameIndex])));

        const TriMesh *trimesh0 = m_kdtree->getMesh(cache.frameIndex,   cache.shapeIndex);
        const TriMesh *trimesh1 = m_kdtree->getMesh(cache.frameIndex+1, cache.shapeIndex);
        const Point *vertexPositions0 = trimesh0->getVertexPositions();
        const Point *vertexPositions1 = trimesh1->getVertexPositions();

        const Triangle tri = m_kdtree->getTriangle(cache.shapeIndex, cache.primIndex);
        const uint32_t idx0 = tri.idx[0], idx1 = tri.idx[1], idx2 = tri.idx[2];
        const Point p0 = vertexPositions0[idx0] * (1-cache.alpha) + vertexPositions1[idx0] * cache.alpha;
        const Point p1 = vertexPositions0[idx1] * (1-cache.alpha) + vertexPositions1[idx1] * cache.alpha;
        const Point p2 = vertexPositions0[idx2] * (1-cache.alpha) + vertexPositions1[idx2] * cache.alpha;

        Vector rel = its.p - p0, du = p1 - p0, dv = p2 - p0;

        Float b1  = dot(du, rel), b2 = dot(dv, rel),
              a11 = dot(du, du), a12 = dot(du, dv),
              a22 = dot(dv, dv),
              invDet = 1.0f / (a11 * a22 - a12 * a12);

        cache.u = ( a22 * b1 - a12 * b2) * invDet,
        cache.v = (-a12 * b1 + a11 * b2) * invDet;

        cache.frameIndex = m_kdtree->findFrame(time);
        cache.alpha = std::max((Float) 0.0f, std::min((Float) 1.0f,
            (time - times[cache.frameIndex])
            / (times[cache.frameIndex + 1] - times[cache.frameIndex])));

        fillIntersectionRecord(Ray(Point(0.0f), its.toWorld(-its.wi), time), &cache, its);
    }


    AABB getAABB() const {
        return m_kdtree->getSpatialAABB();
    }

    size_t getPrimitiveCount() const {
        return m_kdtree->getPrimitiveCount();
    }

    size_t getEffectivePrimitiveCount() const {
        return m_kdtree->getPrimitiveCount();
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Shape)))
            m_kdtree->addShape(static_cast<Shape *>(child));
        else
            Shape::addChild(name, child);
    }

    ref<TriMesh> createTriMesh() {
        return const_cast<TriMesh *>(m_kdtree->getMesh(0, 0));
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Deformable[" << endl
            << "   primitiveCount = " << m_kdtree->getPrimitiveCount() << "," << endl
            << "   timeCount = " << m_kdtree->getTimeCount() << "," << endl
            << "   aabb = " << indent(m_kdtree->getSpatialAABB().toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<SpaceTimeKDTree> m_kdtree;
};

MTS_IMPLEMENT_CLASS_S(SpaceTimeKDTree, false, KDTreeBase)
MTS_IMPLEMENT_CLASS_S(Deformable, false, Shape)
MTS_EXPORT_PLUGIN(Deformable, "Deformable shape");
MTS_NAMESPACE_END
