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

/*!\plugin{spherical}{Spherical camera}
 * \order{5}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional camera-to-world transformation.
 *        \default{none (i.e. camera space $=$ world space)}
 *     }
 *     \parameter{shutterOpen, shutterClose}{\Float}{
 *         Specifies the time interval of the measurement---this
 *         is only relevant when the scene is in motion.
 *         \default{0}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{A rendering made using a spherical camera}{sensor_spherical_cbox.jpg}
 * }
 *
 * The spherical camera captures the illumination arriving from all
 * directions and turns it into a latitude-longitude environment map.
 * It is best used with a high dynamic range film that has 2:1 aspect ratio,
 * and the resulting output can then be turned into a distant light source
 * using the \pluginref{envmap} plugin.
 * By default, the camera is located at the origin, which can
 * be changed by providing a custom \code{toWorld} transformation.
 */

class SphericalCamera : public Sensor {
public:
    SphericalCamera(const Properties &props) : Sensor(props) {
        m_type |= EDeltaPosition | EDirectionSampleMapsToPixels;

        if (props.getTransform("toWorld", Transform()).hasScale())
            Log(EError, "Scale factors in the sensor-to-world "
                "transformation are not allowed!");
    }

    SphericalCamera(Stream *stream, InstanceManager *manager)
     : Sensor(stream, manager) {
        configure();
    }

    Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
            const Point2 &otherSample, Float timeSample) const {
        ray.time = sampleTime(timeSample);
        ray.mint = Epsilon;
        ray.maxt = std::numeric_limits<Float>::infinity();

        const Transform &trafo = m_worldTransform->eval(ray.time);

        Float sinPhi, cosPhi, sinTheta, cosTheta;
        math::sincos(pixelSample.x * m_invResolution.x * 2 * M_PI, &sinPhi, &cosPhi);
        math::sincos(pixelSample.y * m_invResolution.y * M_PI, &sinTheta, &cosTheta);

        Vector d(sinPhi*sinTheta, cosTheta, -cosPhi*sinTheta);

        ray.setOrigin(trafo(Point(0.0f)));
        ray.setDirection(trafo(d));
        return Spectrum(1.0f);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
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
            const Point2 &sample, const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);

        Point samplePos(sample.x, sample.y, 0.0f);

        if (extra) {
            /* The caller wants to condition on a specific pixel position */
            samplePos.x = (extra->x + sample.x) * m_invResolution.x;
            samplePos.y = (extra->y + sample.y) * m_invResolution.y;
        }

        pRec.uv = Point2(samplePos.x * m_resolution.x,
            samplePos.y * m_resolution.y);

        Float sinPhi, cosPhi, sinTheta, cosTheta;
        math::sincos(samplePos.x * 2 * M_PI, &sinPhi, &cosPhi);
        math::sincos(samplePos.y * M_PI, &sinTheta, &cosTheta);

        dRec.d = trafo(Vector(sinPhi*sinTheta, cosTheta, -cosPhi*sinTheta));
        dRec.measure = ESolidAngle;
        dRec.pdf = 1 / (2 * M_PI * M_PI * std::max(sinTheta, Epsilon));

        return Spectrum(1.0f);
    }

    Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle)
            return 0.0f;

        Vector d = m_worldTransform->eval(pRec.time).inverse()(dRec.d);
        Float sinTheta = math::safe_sqrt(1-d.y*d.y);

        return 1 / (2 * M_PI * M_PI * std::max(sinTheta, Epsilon));
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle)
            return Spectrum(0.0f);

        Vector d = m_worldTransform->eval(pRec.time).inverse()(dRec.d);
        Float sinTheta = math::safe_sqrt(1-d.y*d.y);

        return Spectrum(1 / (2 * M_PI * M_PI * std::max(sinTheta, Epsilon)));
    }

    bool getSamplePosition(const PositionSamplingRecord &pRec,
            const DirectionSamplingRecord &dRec, Point2 &samplePosition) const {
        Vector d = normalize(m_worldTransform->eval(pRec.time).inverse()(dRec.d));

        samplePosition = Point2(
            math::modulo(std::atan2(d.x, -d.z) * INV_TWOPI, (Float) 1) * m_resolution.x,
            math::safe_acos(d.y) * INV_PI * m_resolution.y
        );

        return true;
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);

        /* Transform the reference point into the local coordinate system */
        Point refP = trafo.inverse().transformAffine(dRec.ref);
        Vector d(refP);
        Float dist = d.length(),
              invDist = 1.0f / dist;
        d *= invDist;

        dRec.uv = Point2(
            math::modulo(std::atan2(d.x, -d.z) * INV_TWOPI, (Float) 1) * m_resolution.x,
            math::safe_acos(d.y) * INV_PI * m_resolution.y
        );

        Float sinTheta = math::safe_sqrt(1-d.y*d.y);

        dRec.p = trafo.transformAffine(Point(0.0f));
        dRec.d = (dRec.p - dRec.ref) * invDist;
        dRec.dist = dist;
        dRec.n = Vector(0.0f);
        dRec.pdf = 1;
        dRec.measure = EDiscrete;

        return Spectrum(
            (1/(2 * M_PI * M_PI * std::max(sinTheta, Epsilon))) * invDist * invDist);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return (dRec.measure == EDiscrete) ? 1.0f : 0.0f;
    }

    AABB getAABB() const {
        return m_worldTransform->getTranslationBounds();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SphericalCamera[" << endl
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

MTS_IMPLEMENT_CLASS_S(SphericalCamera, false, Sensor)
MTS_EXPORT_PLUGIN(SphericalCamera, "Spherical camera");
MTS_NAMESPACE_END
