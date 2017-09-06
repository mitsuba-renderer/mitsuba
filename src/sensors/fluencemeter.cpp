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

#include <mitsuba/render/sensor.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/track.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{fluencemeter}{Fluence meter}
 * \order{7}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional sensor-to-world transformation.
 *        \default{none (i.e. sensor space $=$ world space)}
 *     }
 *     \parameter{shutterOpen, shutterClose}{\Float}{
 *         Specifies the time interval of the measurement---this
 *         is only relevant when the scene is in motion.
 *         \default{0}
 *     }
 * }
 *
 * This sensor plugin implements a simple fluence meter, which measures
 * the average radiance passing through a specified position.
 * By default, the sensor is located at the origin.
 *
 * Such a sensor is useful for conducting virtual experiments and
 * testing the renderer for correctness.
 *
 * \vspace{4mm}
 * \begin{xml}
 * <scene version=$\MtsVer$>
 *     <sensor type="fluencemeter">
 *         <!-- Measure the average radiance traveling
 *              through the point (1,2,3) -->
 *         <transform name="toWorld">
 *             <translate x="1" y="2" z="3"/>
 *         </transform>
 *
 *         <!-- Write the output to a MATLAB M-file. The output file will
 *              contain a 1x1 matrix storing the computed estimate -->
 *         <film type="mfilm"/>
 *
 *         <!-- Use 1024 samples for the measurement -->
 *         <sampler type="independent">
 *             <integer name="sampleCount" value="1024"/>
 *         </sampler>
 *     </sensor>
 *
 *     <!-- ... other scene declarations ... -->
 * </scene>
 * \end{xml}
 */

class FluenceMeter : public Sensor {
public:
    FluenceMeter(const Properties &props) : Sensor(props) {
        m_type |= EDeltaPosition;

        if (props.getTransform("toWorld", Transform()).hasScale())
            Log(EError, "Scale factors in the sensor-to-world "
                "transformation are not allowed!");
    }

    FluenceMeter(Stream *stream, InstanceManager *manager)
     : Sensor(stream, manager) {
        configure();
    }

    Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
            const Point2 &otherSample, Float timeSample) const {
        ray.time = sampleTime(timeSample);
        ray.mint = Epsilon;
        ray.maxt = std::numeric_limits<Float>::infinity();

        const Transform &trafo = m_worldTransform->eval(ray.time);
        ray.setOrigin(trafo(Point(0.0f)));
        ray.setDirection(trafo(warp::squareToUniformSphere(pixelSample)));
        return Spectrum(1.0f);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        pRec.p = trafo(Point(0.0f));
        pRec.n = Normal(0.0f);
        pRec.pdf = 1.0f;
        pRec.measure = EDiscrete;
        return Spectrum(1.0f);
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return Spectrum((pRec.measure == EDiscrete) ? 1.0f : 0.0f);
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

        return Spectrum(INV_FOURPI * invDist * invDist);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return dRec.measure == EDiscrete ? 1.0f : 0.0f;
    }

    AABB getAABB() const {
        return m_worldTransform->getTranslationBounds();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "FluenceMeter[" << endl
            << "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
            << "  sampler = " << indent(m_sampler->toString()) << "," << endl
            << "  film = " << indent(m_film->toString()) << "," << endl
            << "  medium = " << indent(m_medium.toString()) << "," << endl
            << "  shutterOpen = " << m_shutterOpen << "," << endl
            << "  shutterOpenTime = " << m_shutterOpenTime << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(FluenceMeter, false, Sensor)
MTS_EXPORT_PLUGIN(FluenceMeter, "Fluence meter");
MTS_NAMESPACE_END
