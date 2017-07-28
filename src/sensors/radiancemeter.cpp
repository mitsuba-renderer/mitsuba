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

MTS_NAMESPACE_BEGIN

/*!\plugin{radiancemeter}{Radiance meter}
 * \order{6}
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
 * This sensor plugin implements a simple radiance meter, which measures
 * the incident power per unit area per unit solid angle along a
 * certain ray. It can be thought of as the limit of a standard
 * perspective camera as its field of view tends to zero.
 * Hence, when this sensor is given a film with multiple pixels, all
 * of them will record the same value.
 *
 * Such a sensor is useful for conducting virtual experiments and
 * testing the renderer for correctness.
 *
 * By default, the sensor is located at the origin and performs
 * a measurement in the positive Z direction $(0,0,1)$. This can
 * be changed by providing a custom \code{toWorld} transformation:
 *
 * \vspace{4mm}
 * \begin{xml}
 * <scene version=$\MtsVer$>
 *     <sensor type="radiancemeter">
 *         <!-- Measure the amount of radiance traveling
 *              from the origin to (1,2,3) -->
 *         <transform name="toWorld">
 *             <lookat origin="1,2,3"
 *                     target="0,0,0"/>
 *         </transform>
 *
 *         <!-- Write the output to a MATLAB M-file. The output file will
 *              contain a 1x1 matrix storing an estimate of the incident
 *              radiance along the specified ray. -->
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

class RadianceMeter : public Sensor {
public:
    RadianceMeter(const Properties &props) : Sensor(props) {
        m_type |= EDeltaPosition | EDeltaDirection;

        if (props.getTransform("toWorld", Transform()).hasScale())
            Log(EError, "Scale factors in the sensor-to-world "
                "transformation are not allowed!");
    }

    RadianceMeter(Stream *stream, InstanceManager *manager)
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
        ray.setDirection(trafo(Vector(0.0f, 0.0f, 1.0f)));
        return Spectrum(1.0f);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        pRec.p = trafo(Point(0.0f));
        pRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
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
        oss << "RadianceMeter[" << endl
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

MTS_IMPLEMENT_CLASS_S(RadianceMeter, false, Sensor)
MTS_EXPORT_PLUGIN(RadianceMeter, "Radiance meter");
MTS_NAMESPACE_END
