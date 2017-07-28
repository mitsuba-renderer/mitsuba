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
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{directional}{Directional emitter}
 * \icon{emitter_directional}
 * \order{4}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional emitter-to-world transformation.
 *        \default{none (i.e. emitter space $=$ world space)}
 *     }
 *     \parameter{direction}{\Vector}{
 *        Alternative to \code{toWorld}: explicitly specifies
 *        the illumination direction. Note that only one of the
 *        two parameters can be used.
 *     }
 *     \parameter{irradiance}{\Spectrum}{
 *         Specifies the amount of power per unit area received
 *         by a hypothetical surface normal to the specified direction
 *         \default{1}
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This emitter plugin implements a distant directional source, which
 * radiates a specified power per unit area along a fixed direction.
 * By default, the emitter radiates in the direction of the postive Z axis.
 */

class DirectionalEmitter : public Emitter {
public:
    DirectionalEmitter(const Properties &props) : Emitter(props) {
        m_type |= EDeltaDirection;

        m_normalIrradiance = props.getSpectrum("irradiance", Spectrum::getD65());
        if (props.hasProperty("direction")) {
            if (props.hasProperty("toWorld"))
                Log(EError, "Only one of the parameters 'direction' and 'toWorld'"
                    "can be used at a time!");

            Vector d(normalize(props.getVector("direction"))), u, unused;
            coordinateSystem(d, u, unused);
            m_worldTransform = new AnimatedTransform(
                Transform::lookAt(Point(0.0f), Point(d), u));
        } else {
            if (props.getTransform("toWorld", Transform()).hasScale())
                Log(EError, "Scale factors in the emitter-to-world "
                    "transformation are not allowed!");
        }
    }

    DirectionalEmitter(Stream *stream, InstanceManager *manager)
     : Emitter(stream, manager) {
        m_normalIrradiance = Spectrum(stream);
        m_bsphere = BSphere(stream);
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        m_normalIrradiance.serialize(stream);
        m_bsphere.serialize(stream);
    }

    ref<Shape> createShape(const Scene *scene) {
        /* Create a bounding sphere that surrounds the scene */
        m_bsphere = scene->getKDTree()->getAABB().getBSphere();
        m_bsphere.radius *= 1.1f;
        configure();
        return NULL;
    }

    void configure() {
        Emitter::configure();
        Float surfaceArea = M_PI * m_bsphere.radius * m_bsphere.radius;
        m_invSurfaceArea = 1.0f / surfaceArea;
        m_power = m_normalIrradiance * surfaceArea;
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample, const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);

        Point2 p = warp::squareToUniformDiskConcentric(sample);

        Vector perpOffset = trafo(Vector(p.x, p.y, 0) * m_bsphere.radius);
        Vector d = trafo(Vector(0, 0, 1));

        pRec.p = m_bsphere.center - d*m_bsphere.radius + perpOffset;
        pRec.n = d;
        pRec.pdf = m_invSurfaceArea;
        pRec.measure = EArea;
        return m_power;
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EArea) ? m_normalIrradiance : Spectrum(0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EArea) ? m_invSurfaceArea : 0.0f;
    }

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
            PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
        dRec.d = pRec.n;
        dRec.pdf = 1.0f;
        dRec.measure = EDiscrete;
        return Spectrum(1.0f);
    }

    Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        return (dRec.measure == EDiscrete) ? 1.0f : 0.0f;
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        return Spectrum((dRec.measure == EDiscrete) ? 1.0f : 0.0f);
    }

    Spectrum sampleRay(Ray &ray,
            const Point2 &spatialSample,
            const Point2 &directionalSample,
            Float time) const {
        const Transform &trafo = m_worldTransform->eval(time);
        Point2 p = warp::squareToUniformDiskConcentric(spatialSample);
        Vector perpOffset = trafo(Vector(p.x, p.y, 0) * m_bsphere.radius);
        Vector d = trafo(Vector(0, 0, 1));
        ray.setOrigin(m_bsphere.center - d*m_bsphere.radius + perpOffset);
        ray.setDirection(d);
        ray.setTime(time);
        return m_power;
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);
        Vector d = trafo(Vector(0,0,1));
        Point diskCenter = m_bsphere.center - d*m_bsphere.radius;

        Float distance = dot(dRec.ref - diskCenter, d);
        if (distance < 0) {
            /* This can happen when doing bidirectional renderings
               involving environment maps and directional sources. Just
               return zero */
            return Spectrum(0.0f);
        }

        dRec.p = dRec.ref - distance * d;
        dRec.d = -d;
        dRec.n = Normal(d);
        dRec.dist = distance;

        dRec.pdf = 1.0f;
        dRec.measure = EDiscrete;
        return m_normalIrradiance;
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return dRec.measure == EDiscrete ? 1.0f : 0.0f;
    }

    AABB getAABB() const {
        return AABB();
    }

    Shader *createShader(Renderer *renderer) const;

    std::string toString() const {
        std::ostringstream oss;
        oss << "DirectionalEmitter[" << endl
            << "  normalIrradiance = " << m_normalIrradiance.toString() << "," << endl
            << "  samplingWeight = " << m_samplingWeight << "," << endl
            << "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
            << "  medium = " << indent(m_medium.toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_normalIrradiance, m_power;
    BSphere m_bsphere;
    Float m_invSurfaceArea;
};

// ================ Hardware shader implementation ================

class DirectionalEmitterShader : public Shader {
public:
    DirectionalEmitterShader(Renderer *renderer)
        : Shader(renderer, EEmitterShader) { }

    void generateCode(std::ostringstream &oss, const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "_dir(vec3 wo) {" << endl
            << "    return vec3(1.0);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
};

Shader *DirectionalEmitter::createShader(Renderer *renderer) const {
    return new DirectionalEmitterShader(renderer);
}

MTS_IMPLEMENT_CLASS(DirectionalEmitterShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DirectionalEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(DirectionalEmitter, "Directional emitter");
MTS_NAMESPACE_END
