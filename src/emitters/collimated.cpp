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

MTS_NAMESPACE_BEGIN

/*!\plugin{collimated}{Collimated beam emitter}
 * \icon{emitter_collimated}
 * \order{5}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional emitter-to-world transformation.
 *        \default{none (i.e. emitter space $=$ world space)}
 *     }
 *     \parameter{power}{\Spectrum}{
 *         Specifies the amount of power radiated along the beam
 *         \default{1}
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This emitter plugin implements a collimated beam source, which
 * radiates a specified amount of power along a fixed ray.
 * It can be thought of as the limit of a spot light as its field
 * of view tends to zero.
 *
 * Such a emitter is useful for conducting virtual experiments and
 * testing the renderer for correctness.
 *
 * By default, the emitter is located at the origin and radiates
 * into the positive Z direction $(0,0,1)$. This can
 * be changed by providing a custom \code{toWorld} transformation.
 */

class CollimatedBeamEmitter : public Emitter {
public:
    CollimatedBeamEmitter(const Properties &props) : Emitter(props) {
        m_type |= EDeltaPosition | EDeltaDirection;

        m_power = props.getSpectrum("power", Spectrum::getD65());

        if (props.getTransform("toWorld", Transform()).hasScale())
            Log(EError, "Scale factors in the emitter-to-world "
                "transformation are not allowed!");
    }

    CollimatedBeamEmitter(Stream *stream, InstanceManager *manager)
     : Emitter(stream, manager) {
        configure();
        m_power = Spectrum(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        m_power.serialize(stream);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample, const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        pRec.p = trafo(Point(0.0f));
        pRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
        pRec.pdf = 1.0f;
        pRec.measure = EDiscrete;
        return m_power;
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? m_power : Spectrum(0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? 1.0f : 0.0f;
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
        ray.setTime(time);
        ray.setOrigin(trafo.transformAffine(Point(0.0f)));
        ray.setDirection(trafo(Vector(0.0f, 0.0f, 1.0f)));
        return m_power;
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        /* Direct sampling always fails for a response function on a 0D space */
        dRec.pdf = 0.0f;
        return Spectrum(0.0f);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return 0.0f;
    }

    AABB getAABB() const {
        return m_worldTransform->getTranslationBounds();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "CollimatedBeamEmitter[" << endl
            << "  power = " << m_power.toString() << "," << endl
            << "  samplingWeight = " << m_samplingWeight << "," << endl
            << "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
            << "  medium = " << indent(m_medium.toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_power;
};

MTS_IMPLEMENT_CLASS_S(CollimatedBeamEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(CollimatedBeamEmitter, "Collimated beam emitter");
MTS_NAMESPACE_END
