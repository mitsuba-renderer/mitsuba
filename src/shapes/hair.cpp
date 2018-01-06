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

#include "hair.h"
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/sahkdtree3.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

#define MTS_HAIR_USE_FANCY_CLIPPING 1

MTS_NAMESPACE_BEGIN

/*!\plugin{hair}{Hair intersection shape}
 * \order{11}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the hair data file that should be loaded
 *     }
 *     \parameter{radius}{\Float}{
 *       Radius of the hair segments in world-space units
 *       \default{0.025, which assumes that the scene
 *       is modeled in millimeters.}.
 *     }
 *     \parameter{angleThreshold}{\Float}{
 *       For performance reasons, the plugin will merge adjacent hair
 *       segments when the angle of their tangent directions is below
 *       than this value (in degrees). \default{1}.
 *     }
 *     \parameter{reduction}{\Float}{
 *       When the reduction ratio is set to a value between zero and one, the hair
 *       plugin stochastically culls this portion of the input data (where
 *       1 corresponds to removing all hairs). To approximately preserve the
 *       appearance in renderings, the hair radius is enlarged (see Cook et al.
 *       \cite{Cook2007Stochastic}). This parameter is convenient for fast
 *       previews. \default{0, i.e. all geometry is rendered}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *        Specifies an optional linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none, i.e. object space $=$ world space}
 *     }
 * }
 * \renderings{
 *     \centering
 *     \fbox{\includegraphics[width=6cm]{images/shape_hair}}\hspace{4.5cm}
 *     \caption{A close-up of the hair shape rendered with a diffuse
 *     scattering model (an actual hair scattering model will
 *     be needed for realistic appearance)}
 * }
 * The plugin implements a space-efficient acceleration structure for
 * hairs made from many straight cylindrical hair segments with miter
 * joints. The underlying idea is that intersections with straight cylindrical
 * hairs can be found quite efficiently, and curved hairs are easily
 * approximated using a series of such segments.
 *
 * The plugin supports two different input formats: a simple (but not
 * particularly efficient) ASCII format containing the coordinates of a
 * hair vertex on every line. An empty line marks the beginning of a
 * new hair. The following snippet is an example of this format:\newpage
 * \begin{xml}
 * .....
 * -18.5498 -21.7669 22.8138
 * -18.6358 -21.3581 22.9262
 * -18.7359 -20.9494 23.0256
 *
 * -30.6367 -21.8369 6.78397
 * -30.7289 -21.4145 6.76688
 * -30.8226 -20.9933 6.73948
 * .....
 * \end{xml}
 *
 * There is also a binary format, which starts with the identifier
 * ``\texttt{BINARY\_HAIR}'' (11 bytes), followed by the number of
 * vertices as a 32-bit little endian integer.
 * The remainder of the file consists of the vertex positions stored as
 * single-precision XYZ coordinates (again in little-endian byte ordering).
 * To mark the beginning of a new hair strand, a single $+\infty$ floating
 * point value can be inserted between the vertex data.
 */

class HairKDTree : public SAHKDTree3D<HairKDTree> {
    friend class GenericKDTree<AABB, SurfaceAreaHeuristic3, HairKDTree>;
    friend class SAHKDTree3D<HairKDTree>;
public:
    using SAHKDTree3D<HairKDTree>::IndexType;
    using SAHKDTree3D<HairKDTree>::SizeType;

    HairKDTree(std::vector<Point> &vertices,
            std::vector<bool> &vertexStartsFiber, Float radius)
            : m_radius(radius) {
        /* Take the supplied vertex & start fiber arrays (without copying) */
        m_vertices.swap(vertices);
        m_vertexStartsFiber.swap(vertexStartsFiber);
        m_hairCount = 0;

        /* Compute the index of the first vertex in each segment. */
        m_segIndex.reserve(m_vertices.size());
        for (size_t i=0; i<m_vertices.size()-1; i++) {
            if (m_vertexStartsFiber[i])
                m_hairCount++;
            if (!m_vertexStartsFiber[i+1])
                m_segIndex.push_back((IndexType) i);
        }
        m_segmentCount = m_segIndex.size();

        Log(EDebug, "Building a kd-tree for " SIZE_T_FMT " hair vertices, "
            SIZE_T_FMT " segments, " SIZE_T_FMT " hairs",
            m_vertices.size(), m_segmentCount, m_hairCount);

        /* Ray-cylinder intersections are expensive. Use only the
           SAH cost as the tree subdivision stopping criterion,
           not the number of primitives */
        setStopPrims(1);

        /* Some other defaults that work well in practice */
        setTraversalCost(10);
        setQueryCost(15);
        setExactPrimitiveThreshold(16384);
        setClip(true);
        setRetract(true);
        setEmptySpaceBonus(0.9f);

        buildInternal();

        Log(EDebug, "Total amount of storage (kd-tree & vertex data): %s",
            memString(m_nodeCount * sizeof(KDNode)
            + m_indexCount * sizeof(IndexType)
            + vertices.size() * sizeof(Point)
            + vertexStartsFiber.size() / 8).c_str());

        /* Optimization: replace all primitive indices by the
           associated vertex indices (this avoids an extra
           indirection during traversal later on) */
        for (SizeType i=0; i<m_indexCount; ++i)
            m_indices[i] = m_segIndex[m_indices[i]];

        /* Free the segIndex array, it is not needed anymore */
        std::vector<IndexType>().swap(m_segIndex);
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

    /// Return the total number of segments
    inline size_t getSegmentCount() const {
        return m_segmentCount;
    }

    /// Return the total number of hairs
    inline size_t getHairCount() const {
        return m_hairCount;
    }

    /// Return the total number of vertices
    inline size_t getVertexCount() const {
        return m_vertices.size();
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

#if MTS_HAIR_USE_FANCY_CLIPPING == 1
    /**
     * Compute the ellipse created by the intersection of an infinite
     * cylinder and a plane. Returns false in the degenerate case.
     * Based on:
     * www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
     */
    bool intersectCylPlane(Point planePt, Normal planeNrml,
            Point cylPt, Vector cylD, Float radius, Point &center,
            Vector *axes, Float *lengths) const {
        if (absDot(planeNrml, cylD) < Epsilon)
            return false;

        Assert(std::abs(planeNrml.length()-1) <Epsilon);
        Vector B, A = cylD - dot(cylD, planeNrml)*planeNrml;

        Float length = A.length();
        if (length > Epsilon && planeNrml != cylD) {
            A /= length;
            B = cross(planeNrml, A);
        } else {
            coordinateSystem(planeNrml, A, B);
        }

        Vector delta = planePt - cylPt,
               deltaProj = delta - cylD*dot(delta, cylD);

        Float aDotD = dot(A, cylD);
        Float bDotD = dot(B, cylD);
        Float c0 = 1-aDotD*aDotD;
        Float c1 = 1-bDotD*bDotD;
        Float c2 = 2*dot(A, deltaProj);
        Float c3 = 2*dot(B, deltaProj);
        Float c4 = dot(delta, deltaProj) - radius*radius;

        Float lambda = (c2*c2/(4*c0) + c3*c3/(4*c1) - c4)/(c0*c1);

        Float alpha0 = -c2/(2*c0),
              beta0 = -c3/(2*c1);

        lengths[0] = std::sqrt(c1*lambda),
        lengths[1] = std::sqrt(c0*lambda);

        center = planePt + alpha0 * A + beta0 * B;
        axes[0] = A;
        axes[1] = B;
        return true;
    }

    /**
     * \brief Intersect an infinite cylinder with an
     * AABB face and bound the resulting clipped ellipse
     */
    AABB intersectCylFace(int axis,
            const Point &min, const Point &max,
            const Point &cylPt, const Vector &cylD) const {
        int axis1 = (axis + 1) % 3;
        int axis2 = (axis + 2) % 3;

        Normal planeNrml(0.0f);
        planeNrml[axis] = 1;

        Point ellipseCenter;
        Vector ellipseAxes[2];
        Float ellipseLengths[2];

        AABB aabb;
        if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius * (1 + Epsilon),
            ellipseCenter, ellipseAxes, ellipseLengths)) {
            /* Degenerate case -- return an invalid AABB. This is
               not a problem, since one of the other faces will provide
               enough information to arrive at a correct clipped AABB */
            return aabb;
        }

        /* Intersect the ellipse against the sides of the AABB face */
        for (int i=0; i<4; ++i) {
            Point p1, p2;
            p1[axis] = p2[axis] = min[axis];
            p1[axis1] = ((i+1) & 2) ? min[axis1] : max[axis1];
            p1[axis2] = ((i+0) & 2) ? min[axis2] : max[axis2];
            p2[axis1] = ((i+2) & 2) ? min[axis1] : max[axis1];
            p2[axis2] = ((i+1) & 2) ? min[axis2] : max[axis2];

            Point2 p1l(
                dot(p1 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
                dot(p1 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);
            Point2 p2l(
                dot(p2 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
                dot(p2 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);

            Vector2 rel = p2l-p1l;
            Float A = dot(rel, rel);
            Float B = 2*dot(Vector2(p1l), rel);
            Float C = dot(Vector2(p1l), Vector2(p1l))-1;

            Float x0, x1;
            if (solveQuadratic(A, B, C, x0, x1)) {
                if (x0 >= 0 && x0 <= 1)
                    aabb.expandBy(p1+(p2-p1)*x0);
                if (x1 >= 0 && x1 <= 1)
                    aabb.expandBy(p1+(p2-p1)*x1);
            }
        }

        ellipseAxes[0] *= ellipseLengths[0];
        ellipseAxes[1] *= ellipseLengths[1];
        AABB faceBounds(min, max);

        /* Find the componentwise maxima of the ellipse */
        for (int i=0; i<2; ++i) {
            int j = (i==0) ? axis1 : axis2;
            Float alpha = ellipseAxes[0][j];
            Float beta = ellipseAxes[1][j];
            Float tmp = 1 / std::sqrt(alpha*alpha + beta*beta);
            Float cosTheta = alpha * tmp, sinTheta = beta*tmp;

            Point p1 = ellipseCenter + cosTheta*ellipseAxes[0] + sinTheta*ellipseAxes[1];
            Point p2 = ellipseCenter - cosTheta*ellipseAxes[0] - sinTheta*ellipseAxes[1];

            if (faceBounds.contains(p1))
                aabb.expandBy(p1);
            if (faceBounds.contains(p2))
                aabb.expandBy(p2);
        }

        return aabb;
    }

    AABB getAABB(IndexType index) const {
        IndexType iv = m_segIndex[index];
        Point center;
        Vector axes[2];
        Float lengths[2];

        bool success = intersectCylPlane(firstVertex(iv), firstMiterNormal(iv),
            firstVertex(iv), tangent(iv), m_radius * (1-Epsilon), center, axes, lengths);
        Assert(success);

        AABB result;
        axes[0] *= lengths[0]; axes[1] *= lengths[1];
        for (int i=0; i<3; ++i) {
            Float range = std::sqrt(axes[0][i]*axes[0][i] + axes[1][i]*axes[1][i]);
            result.min[i] = std::min(result.min[i], center[i]-range);
            result.max[i] = std::max(result.max[i], center[i]+range);
        }

        success = intersectCylPlane(secondVertex(iv), secondMiterNormal(iv),
            secondVertex(iv), tangent(iv), m_radius * (1-Epsilon), center, axes, lengths);
        Assert(success);

        axes[0] *= lengths[0]; axes[1] *= lengths[1];
        for (int i=0; i<3; ++i) {
            Float range = std::sqrt(axes[0][i]*axes[0][i] + axes[1][i]*axes[1][i]);
            result.min[i] = std::min(result.min[i], center[i]-range);
            result.max[i] = std::max(result.max[i], center[i]+range);
        }
        return result;
    }

    AABB getClippedAABB(IndexType index, const AABB &box) const {
        /* Compute a base bounding box */
        AABB base(getAABB(index));
        base.clip(box);

        IndexType iv = m_segIndex[index];
        Point cylPt = firstVertex(iv);
        Vector cylD = tangent(iv);

        /* Now forget about the cylinder ends and
           intersect an infinite cylinder with each AABB face */
        AABB clippedAABB;
        clippedAABB.expandBy(intersectCylFace(0,
                Point(base.min.x, base.min.y, base.min.z),
                Point(base.min.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(0,
                Point(base.max.x, base.min.y, base.min.z),
                Point(base.max.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(1,
                Point(base.min.x, base.min.y, base.min.z),
                Point(base.max.x, base.min.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(1,
                Point(base.min.x, base.max.y, base.min.z),
                Point(base.max.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(2,
                Point(base.min.x, base.min.y, base.min.z),
                Point(base.max.x, base.max.y, base.min.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(2,
                Point(base.min.x, base.min.y, base.max.z),
                Point(base.max.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.clip(base);

        return clippedAABB;
    }
#else
    /// Compute the AABB of a segment (only used during tree construction)
    AABB getAABB(IndexType index) const {
        IndexType iv = m_segIndex[index];

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
    AABB getClippedAABB(IndexType index, const AABB &box) const {
        AABB aabb(getAABB(index));
        aabb.clip(box);
        return aabb;
    }
#endif

    /// Return the total number of segments
    inline SizeType getPrimitiveCount() const {
        return (SizeType) m_segIndex.size();
    }

    struct IntersectionStorage {
        IndexType iv;
        Point p;
    };

    inline bool intersect(const Ray &ray, IndexType iv,
        Float mint, Float maxt, Float &t, void *tmp) const {
        /* First compute the intersection with the infinite cylinder */
        Vector3d axis = tangentDouble(iv);

        // Projection of ray onto subspace normal to axis
        Point3d rayO(ray.o);
        Vector3d rayD(ray.d);
        Point3d v1 = firstVertexDouble(iv);

        Vector3d relOrigin = rayO - v1;
        Vector3d projOrigin = relOrigin - dot(axis, relOrigin) * axis;
        Vector3d projDirection = rayD - dot(axis, rayD) * axis;

        // Quadratic to intersect circle in projection
        const double A = projDirection.lengthSquared();
        const double B = 2 * dot(projOrigin, projDirection);
        const double C = projOrigin.lengthSquared() - m_radius*m_radius;

        double nearT, farT;
        if (!solveQuadraticDouble(A, B, C, nearT, farT))
            return false;

        if (!(nearT <= maxt && farT >= mint)) /* NaN-aware conditionals */
            return false;

        /* Next check the intersection points against the miter planes */
        Point3d pointNear = rayO + rayD * nearT;
        Point3d pointFar = rayO + rayD * farT;

        Vector3d n1 = firstMiterNormalDouble(iv);
        Vector3d n2 = secondMiterNormalDouble(iv);
        Point3d v2 = secondVertexDouble(iv);
        IntersectionStorage *storage = static_cast<IntersectionStorage *>(tmp);
        Point p;

        if (dot(pointNear - v1, n1) >= 0 &&
            dot(pointNear - v2, n2) <= 0 &&
            nearT >= mint) {
            p = Point(rayO + rayD * nearT);
            t = (Float) nearT;
        } else if (dot(pointFar - v1, n1) >= 0 &&
                   dot(pointFar - v2, n2) <= 0) {
            if (farT > maxt)
                return false;
            p = Point(rayO + rayD * farT);
            t = (Float) farT;
        } else {
            return false;
        }

        if (storage) {
            storage->iv = iv;
            storage-> p = p;
        }

        return true;
    }

    inline bool intersect(const Ray &ray, IndexType iv,
        Float mint, Float maxt) const {
        Float tempT;
        return intersect(ray, iv, mint, maxt, tempT, NULL);
    }

    /* Some utility functions */
    inline Point firstVertex(IndexType iv) const { return m_vertices[iv]; }
    inline Point3d firstVertexDouble(IndexType iv) const { return Point3d(m_vertices[iv]); }
    inline Point secondVertex(IndexType iv) const { return m_vertices[iv+1]; }
    inline Point3d secondVertexDouble(IndexType iv) const { return Point3d(m_vertices[iv+1]); }
    inline Point prevVertex(IndexType iv) const { return m_vertices[iv-1]; }
    inline Point3d prevVertexDouble(IndexType iv) const { return Point3d(m_vertices[iv-1]); }
    inline Point nextVertex(IndexType iv) const { return m_vertices[iv+2]; }
    inline Point3d nextVertexDouble(IndexType iv) const { return Point3d(m_vertices[iv+2]); }

    inline bool prevSegmentExists(IndexType iv) const { return !m_vertexStartsFiber[iv]; }
    inline bool nextSegmentExists(IndexType iv) const { return !m_vertexStartsFiber[iv+2]; }

    inline Vector tangent(IndexType iv) const { return normalize(secondVertex(iv) - firstVertex(iv)); }
    inline Vector3d tangentDouble(IndexType iv) const { return normalize(Vector3d(secondVertex(iv)) - Vector3d(firstVertex(iv))); }
    inline Vector prevTangent(IndexType iv) const { return normalize(firstVertex(iv) - prevVertex(iv)); }
    inline Vector3d prevTangentDouble(IndexType iv) const { return normalize(firstVertexDouble(iv) - prevVertexDouble(iv)); }
    inline Vector nextTangent(IndexType iv) const { return normalize(nextVertex(iv) - secondVertex(iv)); }
    inline Vector3d nextTangentDouble(IndexType iv) const { return normalize(nextVertexDouble(iv) - secondVertexDouble(iv)); }

    inline Vector firstMiterNormal(IndexType iv) const {
        if (prevSegmentExists(iv))
            return normalize(prevTangent(iv) + tangent(iv));
        else
            return tangent(iv);
    }

    inline Vector secondMiterNormal(IndexType iv) const {
        if (nextSegmentExists(iv))
            return normalize(tangent(iv) + nextTangent(iv));
        else
            return tangent(iv);
    }

    inline Vector3d firstMiterNormalDouble(IndexType iv) const {
        if (prevSegmentExists(iv))
            return normalize(prevTangentDouble(iv) + tangentDouble(iv));
        else
            return tangentDouble(iv);
    }

    inline Vector3d secondMiterNormalDouble(IndexType iv) const {
        if (nextSegmentExists(iv))
            return normalize(tangentDouble(iv) + nextTangentDouble(iv));
        else
            return tangentDouble(iv);
    }


    MTS_DECLARE_CLASS()
protected:
    std::vector<Point> m_vertices;
    std::vector<bool> m_vertexStartsFiber;
    std::vector<IndexType> m_segIndex;
    size_t m_segmentCount;
    size_t m_hairCount;
    Float m_radius;
};

HairShape::HairShape(const Properties &props) : Shape(props) {
    fs::path path = Thread::getThread()->getFileResolver()->resolve(
        props.getString("filename"));
    Float radius = props.getFloat("radius", 0.025f);
    /* Skip segments, whose tangent differs by less than one degree
       compared to the previous one */
    Float angleThreshold = degToRad(props.getFloat("angleThreshold", 1.0f));
    Float dpThresh = std::cos(angleThreshold);

    /* When set to a value n>1, the hair shape object will reduce
       the input by only loading every n-th hair */
    Float reduction = props.getFloat("reduction", 0);
    if (reduction < 0 || reduction >= 1) {
        Log(EError, "The 'reduction' parameter must have a value in [0, 1)!");
    } else if (reduction > 0) {
        Float correction = 1.0f / (1-reduction);
        Log(EDebug, "Reducing the amount of geometry by %.2f%%, scaling radii by %f.",
            reduction * 100, correction);
        radius *= correction;
    }
    ref<Random> random = new Random();

    /* Object-space -> World-space transformation */
    Transform objectToWorld = props.getTransform("toWorld", Transform());
    radius *= objectToWorld(Vector(0, 0, 1)).length();

    Log(EInfo, "Loading hair geometry from \"%s\" ..", path.filename().string().c_str());
    ref<Timer> timer = new Timer();

    ref<FileStream> binaryStream = new FileStream(path, FileStream::EReadOnly);
    binaryStream->setByteOrder(Stream::ELittleEndian);

    const char *binaryHeader = "BINARY_HAIR";
    char temp[11];

    bool binaryFormat = true;
    binaryStream->read(temp, 11);
    if (memcmp(temp, binaryHeader, 11) != 0)
        binaryFormat = false;

    std::vector<Point> vertices;
    std::vector<bool> vertexStartsFiber;
    Vector tangent(0.0f);
    size_t nDegenerate = 0, nSkipped = 0;
    Point p, lastP(0.0f);
    bool ignore = false;

    if (binaryFormat) {
        size_t vertexCount = binaryStream->readUInt();
        Log(EInfo, "Loading " SIZE_T_FMT " hair vertices ..", vertexCount);
        vertices.reserve(vertexCount);
        vertexStartsFiber.reserve(vertexCount);

        bool newFiber = true;
        size_t verticesRead = 0;

        while (verticesRead != vertexCount) {
            Float value = binaryStream->readSingle();
            if (std::isinf(value)) {
                p.x = binaryStream->readSingle();
                p.y = binaryStream->readSingle();
                p.z = binaryStream->readSingle();
                newFiber = true;
                if (reduction > 0)
                    ignore = random->nextFloat() < reduction;
            } else {
                p.x = value;
                p.y = binaryStream->readSingle();
                p.z = binaryStream->readSingle();
            }
            //cout << "Read " << verticesRead << " vertices (vs goal " << vertexCount << ") .." << endl;
            p = objectToWorld(p);
            verticesRead++;

            if (ignore) {
                // Do nothing
                ++nSkipped;
            } else if (newFiber) {
                vertices.push_back(p);
                vertexStartsFiber.push_back(newFiber);
                lastP = p;
                tangent = Vector(0.0f);
            } else if (p != lastP) {
                if (tangent.isZero()) {
                    vertices.push_back(p);
                    vertexStartsFiber.push_back(newFiber);
                    tangent = normalize(p - lastP);
                    lastP = p;
                } else {
                    Vector nextTangent = normalize(p - lastP);
                    if (dot(nextTangent, tangent) > dpThresh) {
                        /* Too small of a difference in the tangent value,
                           just overwrite the previous vertex by the current one */
                        tangent = normalize(p - vertices[vertices.size()-2]);
                        vertices[vertices.size()-1] = p;
                        ++nSkipped;
                    } else {
                        vertices.push_back(p);
                        vertexStartsFiber.push_back(newFiber);
                        tangent = nextTangent;
                    }
                    lastP = p;
                }
            } else {
                nDegenerate++;
            }
            newFiber = false;
        }
    } else {
        std::string line;
        bool newFiber = true;

        fs::ifstream is(path);
        if (is.fail())
            Log(EError, "Could not open \"%s\"!", path.string().c_str());
        while (is.good()) {
            std::getline(is, line);
            if (line.length() > 0 && line[0] == '#') {
                newFiber = true;
                continue;
            }
            std::istringstream iss(line);
            iss >> p.x >> p.y >> p.z;
            if (!iss.fail()) {
                p = objectToWorld(p);
                if (ignore) {
                    // Do nothing
                    ++nSkipped;
                } else if (newFiber) {
                    vertices.push_back(p);
                    vertexStartsFiber.push_back(newFiber);
                    lastP = p;
                    tangent = Vector(0.0f);
                } else if (p != lastP) {
                    if (tangent.isZero()) {
                        vertices.push_back(p);
                        vertexStartsFiber.push_back(newFiber);
                        tangent = normalize(p - lastP);
                        lastP = p;
                    } else {
                        Vector nextTangent = normalize(p - lastP);
                        if (dot(nextTangent, tangent) > dpThresh) {
                            /* Too small of a difference in the tangent value,
                               just overwrite the previous vertex by the current one */
                            tangent = normalize(p - vertices[vertices.size()-2]);
                            vertices[vertices.size()-1] = p;
                            ++nSkipped;
                        } else {
                            vertices.push_back(p);
                            vertexStartsFiber.push_back(newFiber);
                            tangent = nextTangent;
                        }
                        lastP = p;
                    }
                } else {
                    nDegenerate++;
                }
                newFiber = false;
            } else {
                newFiber = true;
                if (reduction > 0)
                    ignore = random->nextFloat() < reduction;
            }
        }
    }

    if (nDegenerate > 0)
        Log(EInfo, "Encountered " SIZE_T_FMT
            " degenerate segments!", nDegenerate);
    if (nSkipped > 0)
        Log(EInfo, "Skipped " SIZE_T_FMT " segments.", nSkipped);
    Log(EInfo, "Done (took %i ms)", timer->getMilliseconds());

    vertexStartsFiber.push_back(true);

    m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
}

HairShape::HairShape(Stream *stream, InstanceManager *manager)
    : Shape(stream, manager) {
    Float radius = stream->readFloat();
    size_t vertexCount = stream->readSize();

    std::vector<Point> vertices(vertexCount);
    std::vector<bool> vertexStartsFiber(vertexCount+1);
    stream->readFloatArray((Float *) &vertices[0], vertexCount * 3);

    for (size_t i=0; i<vertexCount; ++i)
        vertexStartsFiber[i] = stream->readBool();
    vertexStartsFiber[vertexCount] = true;

    m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
}

void HairShape::serialize(Stream *stream, InstanceManager *manager) const {
    Shape::serialize(stream, manager);

    const std::vector<Point> &vertices = m_kdtree->getVertices();
    const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();

    stream->writeFloat(m_kdtree->getRadius());
    stream->writeSize(vertices.size());
    stream->writeFloatArray((Float *) &vertices[0], vertices.size() * 3);
    for (size_t i=0; i<vertices.size(); ++i)
        stream->writeBool(vertexStartsFiber[i]);
}

bool HairShape::rayIntersect(const Ray &ray, Float mint,
        Float maxt, Float &t, void *temp) const {
    return m_kdtree->rayIntersect(ray, mint, maxt, t, temp);
}

bool HairShape::rayIntersect(const Ray &ray, Float mint, Float maxt) const {
    return m_kdtree->rayIntersect(ray, mint, maxt);
}

void HairShape::fillIntersectionRecord(const Ray &ray,
    const void *temp, Intersection &its) const {
    /* No UV coordinates for now */
    its.uv = Point2(0,0);
    its.dpdu = Vector(0,0,0);
    its.dpdv = Vector(0,0,0);

    const HairKDTree::IntersectionStorage *storage =
        static_cast<const HairKDTree::IntersectionStorage *>(temp);
    HairKDTree::IndexType iv = storage->iv;
    its.p = storage->p;

    const Vector axis = m_kdtree->tangent(iv);
    its.shape = this;
    its.geoFrame.s = axis;
    const Vector relHitPoint = its.p - m_kdtree->firstVertex(iv);
    its.geoFrame.n = Normal(normalize(relHitPoint - dot(axis, relHitPoint) * axis));
    its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);

    /* Migitate roundoff error issues by a normal shift of the computed intersection point */
    const Vector local = its.geoFrame.toLocal(relHitPoint);
    its.p += its.geoFrame.n * (m_kdtree->getRadius() - std::sqrt(local.y*local.y+local.z*local.z));

    its.shFrame.n = its.geoFrame.n;
    coordinateSystem(its.shFrame.n, its.dpdu, its.dpdv);
    its.hasUVPartials = false;
    its.instance = this;
    its.time = ray.time;
}

ref<TriMesh> HairShape::createTriMesh() {
    size_t nSegments = m_kdtree->getSegmentCount();
    /// Use very approximate geometry for large hair meshes
    const uint32_t phiSteps = (nSegments > 100000) ? 4 : 10;
    const Float dPhi   = (2*M_PI) / phiSteps;

    ref<TriMesh> mesh = new TriMesh("Hair mesh approximation",
        phiSteps*2*nSegments, phiSteps*2*nSegments, true, false, false);

    Point *vertices = mesh->getVertexPositions();
    Normal *normals = mesh->getVertexNormals();
    Triangle *triangles = mesh->getTriangles();
    size_t triangleIdx = 0, vertexIdx = 0;

    const std::vector<Point> &hairVertices = m_kdtree->getVertices();
    const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();
    const Float radius = m_kdtree->getRadius();
    Float *cosPhi = new Float[phiSteps];
    Float *sinPhi = new Float[phiSteps];
    for (size_t i=0; i<phiSteps; ++i) {
        sinPhi[i] = std::sin(i*dPhi);
        cosPhi[i] = std::cos(i*dPhi);
    }

    uint32_t hairIdx = 0;
    for (HairKDTree::IndexType iv=0; iv<(HairKDTree::IndexType) hairVertices.size()-1; iv++) {
        if (!vertexStartsFiber[iv+1]) {
            for (uint32_t phi=0; phi<phiSteps; ++phi) {
                Vector tangent = m_kdtree->tangent(iv);
                Vector dir = Frame(tangent).toWorld(
                        Vector(cosPhi[phi], sinPhi[phi], 0));
                Normal miterNormal1 = m_kdtree->firstMiterNormal(iv);
                Normal miterNormal2 = m_kdtree->secondMiterNormal(iv);
                Float t1 = dot(miterNormal1, radius*dir) / dot(miterNormal1, tangent);
                Float t2 = dot(miterNormal2, radius*dir) / dot(miterNormal2, tangent);

                Normal normal(normalize(dir));
                normals[vertexIdx] = normal;
                vertices[vertexIdx++] = m_kdtree->firstVertex(iv) + radius*dir - tangent*t1;
                normals[vertexIdx] = normal;
                vertices[vertexIdx++] = m_kdtree->secondVertex(iv) + radius*dir - tangent*t2;

                uint32_t idx0 = 2*(phi + hairIdx*phiSteps), idx1 = idx0+1;
                uint32_t idx2 = (2*phi+2) % (2*phiSteps) + 2*hairIdx*phiSteps, idx3 = idx2+1;
                triangles[triangleIdx].idx[0] = idx0;
                triangles[triangleIdx].idx[1] = idx2;
                triangles[triangleIdx].idx[2] = idx1;
                triangleIdx++;
                triangles[triangleIdx].idx[0] = idx1;
                triangles[triangleIdx].idx[1] = idx2;
                triangles[triangleIdx].idx[2] = idx3;
                triangleIdx++;
            }
            hairIdx++;
        }
    }
    Assert(triangleIdx == phiSteps*2*nSegments);
    Assert(vertexIdx == phiSteps*2*nSegments);

    delete[] cosPhi;
    delete[] sinPhi;

    mesh->copyAttachments(this);
    mesh->configure();

    return mesh.get();
}

const KDTreeBase<AABB> *HairShape::getKDTree() const {
    return m_kdtree.get();
}

const std::vector<Point> &HairShape::getVertices() const {
    return m_kdtree->getVertices();
}

const std::vector<bool> &HairShape::getStartFiber() const {
    return m_kdtree->getStartFiber();
}

AABB HairShape::getAABB() const {
    return m_kdtree->getAABB();
}

size_t HairShape::getPrimitiveCount() const {
    return m_kdtree->getHairCount();
}

size_t HairShape::getEffectivePrimitiveCount() const {
    return m_kdtree->getHairCount();
}

Float HairShape::getSurfaceArea() const {
    Log(EError, "HairShape::getSurfaceArea(): Not implemented.");
    return -1;
}

std::string HairShape::toString() const {
    std::ostringstream oss;
    oss << "Hair[" << endl
        << "   numVertices = " << m_kdtree->getVertexCount() << ","
        << "   numSegments = " << m_kdtree->getSegmentCount() << ","
        << "   numHairs = " << m_kdtree->getHairCount() << ","
        << "   radius = " << m_kdtree->getRadius()
        << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS(HairKDTree, false, KDTreeBase)
MTS_IMPLEMENT_CLASS_S(HairShape, false, Shape)
MTS_EXPORT_PLUGIN(HairShape, "Hair intersection shape");
MTS_NAMESPACE_END
