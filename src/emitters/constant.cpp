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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/bsphere.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{constant}{Constant environment emitter}
 * \icon{emitter_constant}
 * \order{10}
 * \parameters{
 *     \parameter{radiance}{\Spectrum}{
 *         Specifies the emitted radiance in units of
 *         power per unit area per unit steradian.
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This plugin implements a constant environment emitter, which surrounds
 * the scene and radiates diffuse illumination towards it. This is often
 * a good default light source when the goal is to visualize some loaded
 * geometry that uses basic (e.g. diffuse) materials.
 */
class ConstantBackgroundEmitter : public Emitter {
public:
    ConstantBackgroundEmitter(const Properties &props) : Emitter(props) {
        m_type |= EOnSurface | EEnvironmentEmitter;
        m_radiance = props.getSpectrum("radiance", Spectrum::getD65());
    }

    ConstantBackgroundEmitter(Stream *stream, InstanceManager *manager)
     : Emitter(stream, manager) {
        m_radiance = Spectrum(stream);
        m_sceneBSphere = BSphere(stream);
        m_geoBSphere = BSphere(stream);
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        m_radiance.serialize(stream);
        m_sceneBSphere.serialize(stream);
        m_geoBSphere.serialize(stream);
    }

    ref<Shape> createShape(const Scene *scene) {
        /* Create a bounding sphere that surrounds the scene */
        BSphere sceneBSphere(scene->getAABB().getBSphere());
        sceneBSphere.radius = std::max(Epsilon, sceneBSphere.radius * 1.5f);
        BSphere geoBSphere(scene->getKDTree()->getAABB().getBSphere());

        if (sceneBSphere != m_sceneBSphere || geoBSphere != m_geoBSphere) {
            m_sceneBSphere = sceneBSphere;
            m_geoBSphere = geoBSphere;
            configure();
        }

        Transform trafo =
            Transform::translate(Vector(m_sceneBSphere.center)) *
            Transform::scale(Vector(m_sceneBSphere.radius));

        Properties props("sphere");
        props.setTransform("toWorld", trafo);
        props.setBoolean("flipNormals", true);
        Shape *shape = static_cast<Shape *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Shape), props));
        shape->addChild(this);
        shape->configure();

        return shape;
    }

    void configure() {
        Emitter::configure();
        Float surfaceArea = 4 * M_PI *
            m_sceneBSphere.radius * m_sceneBSphere.radius;
        m_invSurfaceArea = 1 / surfaceArea;
        m_power = m_radiance * surfaceArea * M_PI;
    }

    Spectrum eval(const Intersection &its, const Vector &d) const {
        if (dot(its.shFrame.n, d) <= 0)
            return Spectrum(0.0f);
        else
            return m_radiance;
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
        Vector d = warp::squareToUniformSphere(sample);

        pRec.p = m_sceneBSphere.center + d * m_sceneBSphere.radius;
        pRec.n = -d;
        pRec.measure = EArea;
        pRec.pdf = m_invSurfaceArea;

        return m_power;
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return m_radiance * M_PI;
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return m_invSurfaceArea;
    }

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
            PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
        Vector local = warp::squareToCosineHemisphere(sample);
        dRec.d = Frame(pRec.n).toWorld(local);
        dRec.pdf = warp::squareToCosineHemispherePdf(local);
        dRec.measure = ESolidAngle;
        return Spectrum(1.0f);
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        Float dp = dot(dRec.d, pRec.n);

        if (dRec.measure != ESolidAngle || dp < 0)
            dp = 0.0f;

        return Spectrum(INV_PI * dp);
    }

    Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        Float dp = dot(dRec.d, pRec.n);

        if (dRec.measure != ESolidAngle || dp < 0)
            dp = 0.0f;

        return INV_PI * dp;
    }

    Spectrum sampleRay(Ray &ray,
            const Point2 &spatialSample,
            const Point2 &directionalSample,
            Float time) const {
        Vector v0 = warp::squareToUniformSphere(spatialSample);
        Vector v1 = warp::squareToCosineHemisphere(directionalSample);

        ray.setOrigin(m_geoBSphere.center + v0 * m_geoBSphere.radius);
        ray.setDirection(Frame(-v0).toWorld(v1));
        ray.setTime(time);

        return m_radiance * (4 * M_PI * M_PI * m_geoBSphere.radius * m_geoBSphere.radius);
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec,
            const Point2 &sample) const {
        Vector d;
        Float pdf;

        if (!dRec.refN.isZero()) {
            d = warp::squareToCosineHemisphere(sample);
            pdf = warp::squareToCosineHemispherePdf(d);
            d = Frame(dRec.refN).toWorld(d);
        } else {
            d = warp::squareToUniformSphere(sample);
            pdf = warp::squareToUniformSpherePdf();
        }

        /* Intersect against the bounding sphere. This is not really
           necessary for path tracing and similar integrators. However,
           to make BDPT+MLT work with this emitter, we should return
           positions that are consistent with respect to the other
           sampling techniques in this class. */

        Ray ray(dRec.ref, d, 0);
        Float nearT, farT;
        dRec.pdf = 0.0f;
        if (!m_sceneBSphere.rayIntersect(ray, nearT, farT))
            return Spectrum(0.0f);

        if (!(nearT < 0 && farT > 0))
            return Spectrum(0.0f);

        dRec.p = ray(farT);
        dRec.n = normalize(m_sceneBSphere.center - dRec.p);
        dRec.measure = ESolidAngle;
        dRec.d = ray.d;
        dRec.dist = farT;
        dRec.pdf = pdf;

        if (!dRec.refN.isZero() && dot(dRec.d, dRec.refN) <= 0) {
            /// Ignore the sample if roundoff errors moved it to the backside
            return Spectrum(0.0f);
        }

        return m_radiance / pdf;
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        Float pdfSA;

        if (!dRec.refN.isZero())
            pdfSA = INV_PI * std::max((Float) 0.0f, dot(dRec.d, dRec.refN));
        else
            pdfSA = warp::squareToUniformSpherePdf();

        if (dRec.measure == ESolidAngle)
            return pdfSA;
        else if (dRec.measure == EArea)
            return pdfSA * absDot(dRec.d, dRec.n)
                / (dRec.dist * dRec.dist);
        else
            return 0.0f;
    }

    AABB getAABB() const {
        /* The scene sets its bounding box so that it contains all shapes and
           emitters, but this particular emitter always wants to be *a little*
           bigger than the scene. To avoid a silly recursion, just return a
           point here. */
        return AABB(m_sceneBSphere.center);
    }

    Spectrum evalEnvironment(const RayDifferential &ray) const {
        return m_radiance;
    }

    bool fillDirectSamplingRecord(DirectSamplingRecord &dRec, const Ray &ray) const {
        Float nearT, farT;

        if (!m_sceneBSphere.rayIntersect(ray, nearT, farT) || nearT > 0 || farT < 0) {
            Log(EWarn, "fillDirectSamplingRecord(): internal error!");
            return false;
        }

        dRec.p = ray(farT);
        dRec.n = normalize(m_sceneBSphere.center - dRec.p);
        dRec.measure = ESolidAngle;
        dRec.object = this;
        dRec.d = ray.d;
        dRec.dist = farT;
        return true;
    }

    Spectrum fillDirectionSamplingRecord(const Ray &ray) const {
        return m_radiance;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ConstantBackgroundEmitter[" << endl
            << "  radiance = " << m_radiance.toString() << "," << endl
            << "  samplingWeight = " << m_samplingWeight << "," << endl
            << "  geoBSphere = " << m_geoBSphere.toString() << "," << endl
            << "  sceneBSphere = " << m_sceneBSphere.toString() << "," << endl
            << "  medium = " << indent(m_medium.toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
protected:
    Spectrum m_radiance, m_power;
    BSphere m_geoBSphere, m_sceneBSphere;
    Float m_invSurfaceArea;
};

// ================ Hardware shader implementation ================

class ConstantBackgroundEmitterShader : public Shader {
public:
    ConstantBackgroundEmitterShader(Renderer *renderer, const Spectrum &radiance)
        : Shader(renderer, EEmitterShader), m_radiance(radiance * M_PI) {
    }

    void resolve(const GPUProgram *program, const std::string &evalName,
            std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_radiance", false));
    }

    void generateCode(std::ostringstream &oss, const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_radiance;" << endl
            << endl
            << "vec3 " << evalName << "_dir(vec3 wo) {" << endl
            << "    return vec3(1.0);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_background(vec3 wo) {" << endl
            << "    const float inv_pi = 0.318309886183791;" << endl
            << "    return " << evalName << "_radiance * inv_pi;" << endl
            << "}" << endl;
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
        int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_radiance);
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_radiance;
};

Shader *ConstantBackgroundEmitter::createShader(Renderer *renderer) const {
    return new ConstantBackgroundEmitterShader(renderer, m_radiance);
}

MTS_IMPLEMENT_CLASS(ConstantBackgroundEmitterShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(ConstantBackgroundEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(ConstantBackgroundEmitter, "Constant background emitter");
MTS_NAMESPACE_END
