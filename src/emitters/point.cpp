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

#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/track.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{point}{Point light source}
 * \icon{emitter_point}
 * \order{1}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional sensor-to-world transformation.
 *        \default{none (i.e. sensor space $=$ world space)}
 *     }
 *     \parameter{position}{\Point}{
 *        Alternative parameter for specifying the light source
 *        position. Note that only one of the parameters
 *        \code{toWorld} and \code{position} can be used at a time.
 *     }
 *     \parameter{intensity}{\Spectrum}{
 *         Specifies the radiant intensity in units of
 *         power per unit steradian.
 *         \default{1}
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This emitter plugin implements a simple point light source, which
 * uniformly radiates illumination into all directions.
 */

class PointEmitter : public Emitter {
public:
    PointEmitter(const Properties &props) : Emitter(props) {
        m_type |= EDeltaPosition;

        if (props.hasProperty("position")) {
            if (props.hasProperty("toWorld"))
                Log(EError, "Only one of the parameters 'position'"
                    " and 'toWorld' can be used!'");
            m_worldTransform = new AnimatedTransform(
                Transform::translate(Vector(props.getPoint("position"))));
        }

        m_intensity = props.getSpectrum("intensity", Spectrum::getD65());
    }

    PointEmitter(Stream *stream, InstanceManager *manager)
     : Emitter(stream, manager) {
        configure();
        m_intensity = Spectrum(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        m_intensity.serialize(stream);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        pRec.p = trafo(Point(0.0f));
        pRec.n = Normal(0.0f);
        pRec.pdf = 1.0f;
        pRec.measure = EDiscrete;
        return m_intensity * (4 * M_PI);
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? (m_intensity * 4*M_PI) : Spectrum(0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? 1.0f : 0.0f;
    }

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
            PositionSamplingRecord &pRec,
            const Point2 &sample,
            const Point2 *extra) const {
        dRec.d = warp::squareToUniformSphere(sample);
        dRec.pdf = INV_FOURPI;
        dRec.measure = ESolidAngle;
        return Spectrum(1.0f);
    }

    Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        return (dRec.measure == ESolidAngle) ? INV_FOURPI : 0.0f;
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        return Spectrum((dRec.measure == ESolidAngle) ? INV_FOURPI : 0.0f);
    }

    Spectrum sampleRay(Ray &ray,
            const Point2 &spatialSample,
            const Point2 &directionalSample,
            Float time) const {
        const Transform &trafo = m_worldTransform->eval(time);
        ray.setTime(time);
        ray.setOrigin(trafo(Point(0.0f)));
        ray.setDirection(warp::squareToUniformSphere(directionalSample));
        return m_intensity * (4 * M_PI);
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);

        dRec.p = trafo.transformAffine(Point(0.0f));
        dRec.pdf = 1.0f;
        dRec.measure = EDiscrete;
        dRec.uv = Point2(0.5f);
        dRec.d = dRec.p - dRec.ref;
        dRec.dist = dRec.d.length();
        Float invDist = 1.0f / dRec.dist;
        dRec.d *= invDist;
        dRec.n = Normal(0.0f);
        dRec.pdf = 1;
        dRec.measure = EDiscrete;

        return m_intensity * (invDist * invDist);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return dRec.measure == EDiscrete ? 1.0f : 0.0f;
    }

    AABB getAABB() const {
        return m_worldTransform->getTranslationBounds();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "PointEmitter[" << endl
            << "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
            << "  intensity = " << m_intensity.toString() << "," << endl
            << "  medium = " << indent(m_medium.toString())
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    Spectrum m_intensity;
};


// ================ Hardware shader implementation ================

class PointEmitterShader : public Shader {
public:
    PointEmitterShader(Renderer *renderer, const Spectrum &intensity)
        : Shader(renderer, EEmitterShader), m_intensity(intensity) {
    }

    void resolve(const GPUProgram *program, const std::string &evalName,
            std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_intensity", false));
    }

    void generateCode(std::ostringstream &oss, const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_intensity;" << endl
            << endl
            << "vec3 " << evalName << "_area(vec2 uv) {" << endl
            << "    return " << evalName << "_intensity * (4*pi);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_dir(vec3 wo) {" << endl
            << "    return vec3(inv_fourpi);" << endl
            << "}" << endl;
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
        int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_intensity);
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_intensity;
};

Shader *PointEmitter::createShader(Renderer *renderer) const {
    return new PointEmitterShader(renderer, m_intensity);
}

MTS_IMPLEMENT_CLASS(PointEmitterShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(PointEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(PointEmitter, "Point emitter");
MTS_NAMESPACE_END
