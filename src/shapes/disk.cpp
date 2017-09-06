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
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{disk}{Disk intersection primitive}
 * \order{4}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies a linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 *     \parameter{flipNormals}{\Boolean}{
 *        Is the disk inverted, i.e. should the normal vectors
 *        be flipped? \default{\code{false}}
 *     }
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
 *         <emitter type="area">
 *             <spectrum name="intensity" value="4"/>
 *         </emitter>
 *     </shape>
 * </scene>
 * \end{xml}
 */
class Disk : public Shape {
public:
    Disk(const Properties &props) : Shape(props) {
        m_objectToWorld = new AnimatedTransform(props.getAnimatedTransform("toWorld", Transform()));

        if (props.getBoolean("flipNormals", false))
            m_objectToWorld->prependScale(Vector(1, 1, -1));
    }

    Disk(Stream *stream, InstanceManager *manager)
            : Shape(stream, manager) {
        m_objectToWorld = new AnimatedTransform(stream);
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);
        m_objectToWorld->serialize(stream);
    }

    void configure() {
        Shape::configure();

        const Transform &trafo = m_objectToWorld->eval(0);
        Vector dpdu = trafo(Vector(1, 0, 0));
        Vector dpdv = trafo(Vector(0, 1, 0));

        if (std::abs(dot(normalize(dpdu), normalize(dpdv))) > 1e-3f)
            Log(EError, "Error: 'toWorld' transformation contains shear!");

        if (std::abs(dpdu.length() / dpdv.length() - 1) > 1e-3f)
            Log(EError, "Error: 'toWorld' transformation contains a non-uniform scale!");

        m_invSurfaceArea = 1.0f / (M_PI * dpdu.length() * dpdu.length());
    }

    AABB getAABB() const {
        std::set<Float> times;
        m_objectToWorld->collectKeyframes(times);

        AABB aabb;
        for (std::set<Float>::iterator it = times.begin(); it != times.end(); ++it) {
            const Transform &trafo = m_objectToWorld->eval(*it);
            aabb.expandBy(trafo(Point( 1,  0, 0)));
            aabb.expandBy(trafo(Point(-1,  0, 0)));
            aabb.expandBy(trafo(Point( 0,  1, 0)));
            aabb.expandBy(trafo(Point( 0, -1, 0)));
        }
        return aabb;
    }

    Float getSurfaceArea() const {
        const Transform &trafo = m_objectToWorld->eval(0);
        Vector dpdu = trafo(Vector(1, 0, 0));
        Vector dpdv = trafo(Vector(0, 1, 0));
        return M_PI * dpdu.length() * dpdv.length();
    }

    inline bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
        Ray ray;
        m_objectToWorld->eval(ray.time).inverse().transformAffine(_ray, ray);
        Float hit = -ray.o.z / ray.d.z;

        if (!(hit >= mint && hit <= maxt))
            return false;

        Point local = ray(hit);

        if (local.x * local.x + local.y * local.y <= 1) {
            t = hit;

            if (temp) {
                Float *data = static_cast<Float *>(temp);
                data[0] = local.x;
                data[1] = local.y;
            }

            return true;
        } else {
            return false;
        }
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
        const Transform &trafo = m_objectToWorld->eval(ray.time);

        its.shape = this;
        if (r != 0) {
            its.dpdu = trafo(Vector(cosPhi, sinPhi, 0));
            its.dpdv = trafo(Vector(-sinPhi, cosPhi, 0));
        } else {
            its.dpdu = trafo(Vector(1, 0, 0));
            its.dpdv = trafo(Vector(0, 1, 0));
        }

        its.shFrame.n = normalize(trafo(Normal(0, 0, 1)));
        its.uv = Point2(r, phi * INV_TWOPI);
        its.p = ray(its.t);
        its.hasUVPartials = false;
        its.instance = NULL;
        its.time = ray.time;
    }

    ref<TriMesh> createTriMesh() {
        const uint32_t phiSteps = 40;

        ref<TriMesh> mesh = new TriMesh(getName(),
            phiSteps-1, 2*phiSteps, true, true, false);

        Point *vertices = mesh->getVertexPositions();
        Normal *normals = mesh->getVertexNormals();
        Point2 *texcoords = mesh->getVertexTexcoords();
        Triangle *triangles = mesh->getTriangles();

        Float dphi = (2 * M_PI) / (Float) (phiSteps-1);

        const Transform &trafo = m_objectToWorld->eval(0.0f);
        Point center = trafo(Point(0.0f));
        Normal normal = normalize(trafo(Normal(0, 0, 1)));

        for (uint32_t i=0; i<phiSteps; ++i) {
            Float phi = i*dphi;
            vertices[i] = center;
            vertices[phiSteps+i] = trafo(
                Point(std::cos(phi), std::sin(phi), 0)
            );

            normals[i] = normal;
            normals[phiSteps+i] = normal;
            texcoords[i] = Point2(0.0f, phi * INV_TWOPI);
            texcoords[phiSteps+i] = Point2(1.0f, phi * INV_TWOPI);
        }

        for (uint32_t i=0; i<phiSteps-1; ++i) {
            triangles[i].idx[0] = i;
            triangles[i].idx[1] = i+phiSteps;
            triangles[i].idx[2] = i+phiSteps+1;
        }

        mesh->copyAttachments(this);
        mesh->configure();

        return mesh.get();
    }

    void getNormalDerivative(const Intersection &its,
            Vector &dndu, Vector &dndv, bool shadingFrame) const {
        dndu = dndv = Vector(0.0f);
    }

    void samplePosition(PositionSamplingRecord &pRec, const Point2 &sample) const {
        const Transform &trafo = m_objectToWorld->eval(pRec.time);
        Point2 p = warp::squareToUniformDiskConcentric(sample);

        pRec.p = trafo(Point3(p.x, p.y, 0));
        pRec.n = normalize(trafo(Normal(0, 0, 1)));
        pRec.pdf = m_invSurfaceArea;
        pRec.measure = EArea;
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return m_invSurfaceArea;
    }

    size_t getPrimitiveCount() const {
        return 1;
    }

    size_t getEffectivePrimitiveCount() const {
        return 1;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Disk[" << endl
            << "  objectToWorld = " << indent(m_objectToWorld->toString()) << "," << endl
            << "  bsdf = " << indent(m_bsdf.toString()) << "," << endl;
        if (isMediumTransition()) {
            oss << "  interiorMedium = " << indent(m_interiorMedium.toString()) << "," << endl
                << "  exteriorMedium = " << indent(m_exteriorMedium.toString()) << "," << endl;
        }
        oss << "  emitter = " << indent(m_emitter.toString()) << "," << endl
            << "  sensor = " << indent(m_sensor.toString()) << "," << endl
            << "  subsurface = " << indent(m_subsurface.toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<AnimatedTransform> m_objectToWorld;
    Float m_invSurfaceArea;
};

MTS_IMPLEMENT_CLASS_S(Disk, false, Shape)
MTS_EXPORT_PLUGIN(Disk, "Disk intersection primitive");
MTS_NAMESPACE_END
