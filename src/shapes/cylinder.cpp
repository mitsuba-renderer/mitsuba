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
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{cylinder}{Cylinder intersection primitive}
 * \order{2}
 * \parameters{
 *     \parameter{p0}{\Point}{
 *       Object-space starting point of the cylinder's centerline \default{(0, 0, 0)}
 *     }
 *     \parameter{p1}{\Point}{
 *       Object-space endpoint of the cylinder's centerline \default{(0, 0, 1)}
 *     }
 *     \parameter{radius}{\Float}{
 *       Radius of the cylinder in object-space units \default{1}
 *     }
 *     \parameter{flipNormals}{\Boolean}{
 *        Is the cylinder inverted, i.e. should the normal vectors
 *        be flipped? \default{\code{false}, i.e. the normals point outside}
 *     }
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 * }
 * \renderings{
 *     \rendering{Cylinder with the default one-sided shading}
 *         {shape_cylinder_onesided}
 *     \rendering{Cylinder with two-sided shading, see \lstref{cylinder-twosided}}
 *         {shape_cylinder_twosided}
 * }
 * This shape plugin describes a simple cylinder intersection primitive.
 * It should always be preferred over approximations modeled using
 * triangles. Note that the cylinder does not have endcaps -- also,
 * it's interior has inward-facing normals, which most scattering
 * models in Mitsuba will treat as fully absorbing. If this is not
 * desirable, consider using the \pluginref{twosided} plugin.
 *
 * \begin{xml}[caption={A simple example for instantiating a
 * cylinder, whose interior is visible}, label=lst:cylinder-twosided]
 * <shape type="cylinder">
 *     <float name="radius" value="0.3"/>
 *     <bsdf type="twosided">
 *         <bsdf type="diffuse"/>
 *     </bsdf>
 * </shape>
 * \end{xml}
 */
class Cylinder : public Shape {
private:
    Transform m_objectToWorld;
    Transform m_worldToObject;
    Float m_radius, m_length, m_invSurfaceArea;
    bool m_flipNormals;
public:
    Cylinder(const Properties &props) : Shape(props) {
        Float radius = props.getFloat("radius", 1.0f);
        Point p1 = props.getPoint("p0", Point(0.0f, 0.0f, 0.0f));
        Point p2 = props.getPoint("p1", Point(0.0f, 0.0f, 1.0f));
        Vector d = p2 - p1;
        Float length = d.length();
        m_objectToWorld =
            Transform::translate(Vector(p1)) *
            Transform::fromFrame(Frame(d / length)) *
            Transform::scale(Vector(radius, radius, length));

        if (props.hasProperty("toWorld"))
            m_objectToWorld = props.getTransform("toWorld") * m_objectToWorld;

        /// Are the cylinder normals pointing inwards? default: no
        m_flipNormals = props.getBoolean("flipNormals", false);

        // Remove the scale from the object-to-world transform
        m_radius = m_objectToWorld(Vector(1,0,0)).length();
        m_length = m_objectToWorld(Vector(0,0,1)).length();
        m_objectToWorld = m_objectToWorld * Transform::scale(
            Vector(1/m_radius, 1/m_radius, 1/m_length));

        m_worldToObject = m_objectToWorld.inverse();
        m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
        Assert(m_length > 0 && m_radius > 0);
    }

    Cylinder(Stream *stream, InstanceManager *manager)
        : Shape(stream, manager) {
        m_objectToWorld = Transform(stream);
        m_radius = stream->readFloat();
        m_length = stream->readFloat();
        m_flipNormals = stream->readBool();
        m_worldToObject = m_objectToWorld.inverse();
        m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);
        m_objectToWorld.serialize(stream);
        stream->writeFloat(m_radius);
        stream->writeFloat(m_length);
        stream->writeBool(m_flipNormals);
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
        Ray ray;

        /* Transform into the local coordinate system and normalize */
        m_worldToObject(_ray, ray);

        const double
            ox = ray.o.x,
            oy = ray.o.y,
            dx = ray.d.x,
            dy = ray.d.y;

        const double A = dx*dx + dy*dy;
        const double B = 2 * (dx*ox + dy*oy);
        const double C = ox*ox + oy*oy - m_radius*m_radius;

        double nearT, farT;
        if (!solveQuadraticDouble(A, B, C, nearT, farT))
            return false;

        if (!(nearT <= maxt && farT >= mint)) /* NaN-aware conditionals */
            return false;

        const double zPosNear = ray.o.z + ray.d.z * nearT;
        const double zPosFar = ray.o.z + ray.d.z * farT;

        if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
            t = (Float) nearT;
        } else if (zPosFar >= 0 && zPosFar <= m_length) {
            if (farT > maxt)
                return false;
            t = (Float) farT;
        } else {
            return false;
        }

        return true;
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
        Ray ray;

        /* Transform into the local coordinate system and normalize */
        m_worldToObject(_ray, ray);

        const double
            ox = ray.o.x,
            oy = ray.o.y,
            dx = ray.d.x,
            dy = ray.d.y;

        const double A = dx*dx + dy*dy;
        const double B = 2 * (dx*ox + dy*oy);
        const double C = ox*ox + oy*oy - m_radius*m_radius;

        double nearT, farT;
        if (!solveQuadraticDouble(A, B, C, nearT, farT))
            return false;

        if (nearT > maxt || farT < mint)
            return false;

        const double zPosNear = ray.o.z + ray.d.z * nearT;
        const double zPosFar = ray.o.z + ray.d.z * farT;
        if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
            return true;
        } else if (zPosFar >= 0 && zPosFar <= m_length && farT <= maxt) {
            return true;
        } else {
            return false;
        }
    }

    void fillIntersectionRecord(const Ray &ray,
            const void *temp, Intersection &its) const {
        its.p = ray(its.t);
        Point local = m_worldToObject(its.p);

        Float phi = std::atan2(local.y, local.x);
        if (phi < 0)
            phi += 2*M_PI;
        its.uv.x = phi / (2*M_PI);
        its.uv.y = local.z / m_length;

        Vector dpdu = Vector(-local.y, local.x, 0) * (2*M_PI);
        Vector dpdv = Vector(0, 0, m_length);
        its.shape = this;
        its.dpdu = m_objectToWorld(dpdu);
        its.dpdv = m_objectToWorld(dpdv);
        its.geoFrame.s = normalize(its.dpdu);
        its.geoFrame.t = normalize(its.dpdv);
        its.geoFrame.n = Normal(cross(its.geoFrame.s, its.geoFrame.t));

        /* Migitate roundoff error issues by a normal shift of the computed intersection point */
        its.p += its.geoFrame.n * (m_radius - std::sqrt(local.x*local.x+local.y*local.y));

        if (m_flipNormals)
            its.geoFrame.n *= -1;
        its.shFrame.n = its.geoFrame.n;
        its.hasUVPartials = false;
        its.instance = NULL;
        its.time = ray.time;
    }

    void samplePosition(PositionSamplingRecord &pRec, const Point2 &sample) const {
        Float sinTheta, cosTheta;
        math::sincos(sample.y * (2 * M_PI), &sinTheta, &cosTheta);

        Point p(cosTheta*m_radius, sinTheta*m_radius, sample.x * m_length);
        Normal n(cosTheta, sinTheta, 0.0f);

        if (m_flipNormals)
            n *= -1;

        pRec.p = m_objectToWorld(p);
        pRec.n = normalize(m_objectToWorld(n));
        pRec.pdf = m_invSurfaceArea;
        pRec.measure = EArea;
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return m_invSurfaceArea;
    }

    inline AABB getAABB() const {
        Vector x1 = m_objectToWorld(Vector(m_radius, 0, 0));
        Vector x2 = m_objectToWorld(Vector(0, m_radius, 0));
        Point p0 = m_objectToWorld(Point(0, 0, 0));
        Point p1 = m_objectToWorld(Point(0, 0, m_length));
        AABB result;

        /* To bound the cylinder, it is sufficient to find the
           smallest box containing the two circles at the endpoints.
           This can be done component-wise as follows */

        for (int i=0; i<3; ++i) {
            Float range = std::sqrt(x1[i]*x1[i] + x2[i]*x2[i]);

            result.min[i] = std::min(std::min(result.min[i],
                        p0[i]-range), p1[i]-range);
            result.max[i] = std::max(std::max(result.max[i],
                        p0[i]+range), p1[i]+range);
        }

        return result;
    }

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

        Vector B, A = cylD - dot(cylD, planeNrml)*planeNrml;

        Float length = A.length();
        if (length != 0) {
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
        if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius,
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
            Float alpha = ellipseAxes[0][j], beta = ellipseAxes[1][j];
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

    AABB getClippedAABB(const AABB &box) const {
        /* Compute a base bounding box */
        AABB base(getAABB());
        base.clip(box);

        Point cylPt = m_objectToWorld(Point(0, 0, 0));
        Vector cylD(m_objectToWorld(Vector(0, 0, 1)));

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

        clippedAABB.clip(box);
        return clippedAABB;
    }

    ref<TriMesh> createTriMesh() {
        /// Choice of discretization
        const size_t phiSteps = 20;
        const Float dPhi   = (2*M_PI) / phiSteps;

        ref<TriMesh> mesh = new TriMesh("Cylinder approximation",
            phiSteps*2, phiSteps*2, true, false, false);

        Point *vertices = mesh->getVertexPositions();
        Normal *normals = mesh->getVertexNormals();
        Triangle *triangles = mesh->getTriangles();
        size_t triangleIdx = 0, vertexIdx = 0;

        for (size_t phi=0; phi<phiSteps; ++phi) {
            Float sinPhi = std::sin(phi * dPhi);
            Float cosPhi = std::cos(phi * dPhi);
            uint32_t idx0 = (uint32_t) vertexIdx, idx1 = idx0+1;
            uint32_t idx2 = (vertexIdx+2) % (2*phiSteps), idx3 = idx2+1;
            normals[vertexIdx] = m_objectToWorld(Normal(cosPhi, sinPhi, 0) * (m_flipNormals ? (Float) -1 : (Float) 1));
            vertices[vertexIdx++] = m_objectToWorld(Point(cosPhi*m_radius, sinPhi*m_radius, 0));
            normals[vertexIdx] = m_objectToWorld(Normal(cosPhi, sinPhi, 0) * (m_flipNormals ? (Float) -1 : (Float) 1));
            vertices[vertexIdx++] = m_objectToWorld(Point(cosPhi*m_radius, sinPhi*m_radius, m_length));

            triangles[triangleIdx].idx[0] = idx0;
            triangles[triangleIdx].idx[1] = idx2;
            triangles[triangleIdx].idx[2] = idx1;
            triangleIdx++;
            triangles[triangleIdx].idx[0] = idx1;
            triangles[triangleIdx].idx[1] = idx2;
            triangles[triangleIdx].idx[2] = idx3;
            triangleIdx++;
        }

        mesh->copyAttachments(this);
        mesh->configure();

        return mesh.get();
    }

#if 0
    AABB getAABB() const {
        const Point a = m_objectToWorld(Point(0, 0, 0));
        const Point b = m_objectToWorld(Point(0, 0, m_length));

        const Float r = m_radius;
        AABB result;
        result.expandBy(a - Vector(r, r, r));
        result.expandBy(a + Vector(r, r, r));
        result.expandBy(b - Vector(r, r, r));
        result.expandBy(b + Vector(r, r, r));
        return result;
    }
#endif

    Float getSurfaceArea() const {
        return 2*M_PI*m_radius*m_length;
    }

    void getNormalDerivative(const Intersection &its,
            Vector &dndu, Vector &dndv, bool shadingFrame) const {
        dndu = its.dpdu / (m_radius * (m_flipNormals ? -1 : 1));
        dndv = Vector(0.0f);
    }

    size_t getPrimitiveCount() const {
        return 1;
    }

    size_t getEffectivePrimitiveCount() const {
        return 1;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Cylinder[" << endl
            << "  radius = " << m_radius << "," << endl
            << "  length = " << m_length << "," << endl
            << "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
            << "  bsdf = " << indent(m_bsdf.toString()) << "," << endl;
        if (isMediumTransition())
            oss << "  interiorMedium = " << indent(m_interiorMedium.toString()) << "," << endl
                << "  exteriorMedium = " << indent(m_exteriorMedium.toString()) << "," << endl;
        oss << "  emitter = " << indent(m_emitter.toString()) << "," << endl
            << "  sensor = " << indent(m_sensor.toString()) << "," << endl
            << "  subsurface = " << indent(m_subsurface.toString())
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(Cylinder, false, Shape)
MTS_EXPORT_PLUGIN(Cylinder, "Cylinder intersection primitive");
MTS_NAMESPACE_END
