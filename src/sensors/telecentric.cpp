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
#include <mitsuba/core/frame.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{telecentric}{Telecentric lens camera}
 * \order{4}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional sensor-to-world transformation.
 *        \default{none (i.e. camera space $=$ world space)}
 *     }
 *     \parameter{apertureRadius}{\Float}{
 *         Denotes the radius of the camera's aperture in scene units.
 *         \default{\code{0}}
 *     }
 *     \parameter{focusDistance}{\Float}{
 *         Denotes the world-space distance from the camera's aperture to the
 *         focal plane. \default{\code{0}}
 *     }
 *     \parameter{shutterOpen, shutterClose}{\Float}{
 *         Specifies the time interval of the measurement---this
 *         is only relevant when the scene is in motion.
 *         \default{0}
 *     }
 *     \parameter{nearClip, farClip}{\Float}{
 *         Distance to the near/far clip
 *         planes.\default{\code{near\code}-\code{Clip=1e-2} (i.e.
 *         \code{0.01}) and {\code{farClip=1e4} (i.e. \code{10000})}}
 *     }
 * }
 * \renderings{
 * \rendering{The material test ball viewed through an telecentric camera.
 * Note the orthographic view together with a narrow depth of field.}{sensor_telecentric}
 * \medrendering{A rendering of the Cornell box. The red and green walls are partially visible
 * due to the aperture size.}{sensor_telecentric_2}
 * }
 *
 * This plugin implements a simple model of a camera with a \emph{telecentric lens}.
 * This is a type of lens that produces an in-focus orthographic view on a
 * plane at some distance from the sensor. Points away from this plane are out of focus
 * and project onto a circle of confusion. In comparison to idealized orthographic cameras,
 * telecentric lens cameras exist in the real world and find use in some computer
 * vision applications where perspective effects cause problems.
 * This sensor relates to the \pluginref{orthographic} plugin in the same way
 * that \pluginref{thinlens} does to \pluginref{perspective}.
 *
 * The configuration is identical to the \pluginref{orthographic} plugin, except that
 * the additional parameters \code{apertureRadius} and \code{focusDistance} must be provided.
 */
class TelecentricLensCamera : public ProjectiveCamera {
public:
    TelecentricLensCamera(const Properties &props)
            : ProjectiveCamera(props) {
        m_type |= ENeedsApertureSample
            | EOrthographicCamera
            | EPositionSampleMapsToPixels;

        /* World-space aperture radius */
        m_apertureRadius = props.getFloat("apertureRadius", 0.0f);
    }

    TelecentricLensCamera(Stream *stream, InstanceManager *manager)
            : ProjectiveCamera(stream, manager) {
        m_apertureRadius = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        ProjectiveCamera::serialize(stream, manager);
        stream->writeFloat(m_apertureRadius);
    }

    void configure() {
        ProjectiveCamera::configure();

        const Vector2i &filmSize   = m_film->getSize();
        const Vector2i &cropSize   = m_film->getCropSize();
        const Point2i  &cropOffset = m_film->getCropOffset();

        Vector2 relSize((Float) cropSize.x / (Float) filmSize.x,
            (Float) cropSize.y / (Float) filmSize.y);
        Point2 relOffset((Float) cropOffset.x / (Float) filmSize.x,
            (Float) cropOffset.y / (Float) filmSize.y);

        /**
         * These do the following (in reverse order):
         *
         * 1. Create transform from camera space to [-1,1]x[-1,1]x[0,1] clip
         *    coordinates (not taking account of the aspect ratio yet)
         *
         * 2+3. Translate and scale to shift the clip coordinates into the
         *    range from zero to one, and take the aspect ratio into account.
         *
         * 4+5. Translate and scale the coordinates once more to account
         *     for a cropping window (if there is any)
         */
        m_cameraToSample =
              Transform::scale(Vector(1.0f / relSize.x, 1.0f / relSize.y, 1.0f))
            * Transform::translate(Vector(-relOffset.x, -relOffset.y, 0.0f))
            * Transform::scale(Vector(-0.5f, -0.5f*m_aspect, 1.0f))
            * Transform::translate(Vector(-1.0f, -1.0f/m_aspect, 0.0f))
            * Transform::orthographic(m_nearClip, m_farClip);

        m_sampleToCamera = m_cameraToSample.inverse();

        /* Position differentials on the near plane */
        m_dx = m_sampleToCamera(Point(m_invResolution.x, 0.0f, 0.0f))
             - m_sampleToCamera(Point(0.0f));
        m_dy = m_sampleToCamera(Point(0.0f, m_invResolution.y, 0.0f))
             - m_sampleToCamera(Point(0.0f));

        /* Clip-space transformation for OpenGL */
        m_clipTransform = Transform::translate(
            Vector((1-2*relOffset.x)/relSize.x - 1,
                  -(1-2*relOffset.y)/relSize.y + 1, 0.0f)) *
            Transform::scale(Vector(1.0f / relSize.x, 1.0f / relSize.y, 1.0f));

        const Transform &trafo = m_worldTransform->eval(0.0f);

        m_normalization = 1.0f / (
            trafo(m_sampleToCamera(Vector(1, 0, 0))).length() *
            trafo(m_sampleToCamera(Vector(0, 1, 0))).length());

        m_scale = Vector(
            trafo(Vector(1, 0, 0)).length(),
            trafo(Vector(0, 1, 0)).length(),
            trafo(Vector(0, 0, 1)).length());

        m_aperturePdf = 1 / (M_PI * m_apertureRadius * m_apertureRadius);
        m_maxRotation = 360/M_PI * std::atan(m_apertureRadius/(m_focusDistance * m_scale.x));
    }

    /**
     * \brief Compute the percentage of overlap between an arbitrary disk and a
     * rectangle that is centered at the origin
     *
     * \remark To restrict the possible number of cases to a sane amount, the
     * implementation here assumes that <tt>size.{x,y}/2 < radius</tt>.
     *
     * \param center
     *   Center of the disk
     * \param radius
     *    Radius of the disk
     * \param size
     *    Total width and height of the rectangular region
     * \author Wenzel Jakob
     */
    Float convolve(const Point2 &center, Float radius, const Vector2 &size) const {
        /* Coordinates of the lower right rectangle corner in the disc frame */
        Float invRadius = 1.0f / radius,
              x = (0.5f * size.x - std::abs(center.x)) * invRadius,
              y = (0.5f * size.y - std::abs(center.y)) * invRadius;

        if (x < -1 || y < -1 || (x < 0 && y < 0 && x*x+y*y > 1)) {
            return 0.0f; /* No coverage at all */
        } else if (x >= 1 && y >= 1) {
            return 1.0f; /* Full coverage */
        } else if (x >= 0 && x <= 1 && y >= 0 && y <= 1 && x*x+y*y > 1) {
            /* There are four intersections. Area = full circle - 2 circular segments */
            return INV_PI * (x * std::sqrt(1-x*x) + y * std::sqrt(1-y*y) + std::asin(x) + std::asin(y));
        } else if (x*x + y*y <= 1) {
            /* There are two intersections (case 1): area = 1 circular segment + 1 triangle. */
            Float yt = -std::sqrt(1-x*x), xt = -std::sqrt(1-y*y),
                  ct = x*xt + y*yt;
            return INV_TWOPI * ((x-xt) * (y-yt) + math::safe_acos(ct) - math::safe_sqrt(1-ct*ct));
        } else {
            /* There are two intersections (case 2): area = 1 circular segment only */
            Float h = std::min(x, y);
            return INV_PI * (std::acos(-h) + h*std::sqrt(1-h*h));
        }
    }

    Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
            const Point2 &otherSample, Float timeSample) const {
        Point2 diskSample = warp::squareToUniformDiskConcentric(otherSample)
            * (m_apertureRadius / m_scale.x);
        ray.time = sampleTime(timeSample);

        /* Compute the corresponding position on the
           near plane (in local camera space) */
        Point focusP = m_sampleToCamera.transformAffine(Point(
            pixelSample.x * m_invResolution.x,
            pixelSample.y * m_invResolution.y, 0.0f));
        focusP.z = m_focusDistance/m_scale.z;

        /* Compute the ray origin */
        Point orig(diskSample.x+focusP.x,
            diskSample.y+focusP.y, 0.0f);

        /* Turn these into a normalized ray direction, and
           adjust the ray interval accordingly */
        const Transform &trafo = m_worldTransform->eval(ray.time);

        ray.setOrigin(trafo.transformAffine(orig));
        ray.setDirection(normalize(trafo(focusP - orig)));
        ray.mint = m_nearClip;
        ray.maxt = m_farClip;

        return Spectrum(1.0f);
    }

    Spectrum sampleRayDifferential(RayDifferential &ray, const Point2 &pixelSample,
            const Point2 &otherSample, Float timeSample) const {
        Point2 diskSample = warp::squareToUniformDiskConcentric(otherSample)
            * (m_apertureRadius / m_scale.x);
        ray.time = sampleTime(timeSample);

        /* Compute the corresponding position on the
           near plane (in local camera space) */
        Point focusP = m_sampleToCamera.transformAffine(Point(
            pixelSample.x * m_invResolution.x,
            pixelSample.y * m_invResolution.y, 0.0f));
        focusP.z = m_focusDistance/m_scale.z;

        /* Compute the ray origin */
        Point orig(diskSample.x+focusP.x,
            diskSample.y+focusP.y, 0.0f);

        /* Turn these into a normalized ray direction, and
           adjust the ray interval accordingly */
        const Transform &trafo = m_worldTransform->eval(ray.time);

        ray.setOrigin(trafo.transformAffine(orig));
        ray.setDirection(normalize(trafo(focusP - orig)));
        ray.mint = m_nearClip;
        ray.maxt = m_farClip;
        ray.rxOrigin = trafo(orig + m_dx);
        ray.ryOrigin = trafo(orig + m_dy);
        ray.rxDirection = ray.ryDirection = ray.d;
        ray.hasDifferentials = true;

        return Spectrum(1.0f);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");

        /* This function should sample from a rectangle that has been
           convolved with a circle, but using only two random numbers.
           Turns out that this is actually fairly difficult to do.
           Since the precision requirements are modest, we can get away
           with the following *terrible* hack, which makes four random
           numbers out of two and then does convolution the "trivial" way */

        #if defined(SINGLE_PRECISION)
            uint32_t tmp1 = union_cast<uint32_t>(sample.x + 1.0f) & 0x7FFFFF;
            uint32_t tmp2 = union_cast<uint32_t>(sample.y + 1.0f) & 0x7FFFFF;

            float rand1 = (tmp1 >> 11)   * (1.0f / 0xFFF);
            float rand2 = (tmp2 >> 11)   * (1.0f / 0xFFF);
            float rand3 = (tmp1 & 0x7FF) * (1.0f / 0x7FF);
            float rand4 = (tmp2 & 0x7FF) * (1.0f / 0x7FF);
        #else
            uint64_t tmp1 = union_cast<uint64_t>(sample.x + 1.0) & 0xFFFFFFFFFFFFF;
            uint64_t tmp2 = union_cast<uint64_t>(sample.y + 1.0) & 0xFFFFFFFFFFFFF;

            double rand1 = (tmp1 >> 26)       * (1.0 / 0x3FFFFFF);
            double rand2 = (tmp2 >> 26)       * (1.0 / 0x3FFFFFF);
            double rand3 = (tmp1 & 0x3FFFFFF) * (1.0 / 0x3FFFFFF);
            double rand4 = (tmp2 & 0x3FFFFFF) * (1.0 / 0x3FFFFFF);
        #endif

        Point2 aperturePos = warp::squareToUniformDiskConcentric(Point2(rand1, rand2))
            * (m_apertureRadius / m_scale.x);
        Point2 samplePos(rand3, rand4);

        if (extra) {
            /* The caller wants to condition on a specific pixel position */
            pRec.uv = *extra + samplePos;
            samplePos.x = pRec.uv.x * m_invResolution.x;
            samplePos.y = pRec.uv.y * m_invResolution.y;
        }

        Point p = m_sampleToCamera.transformAffine(Point(
            aperturePos.x + samplePos.x, aperturePos.y + samplePos.y, 0.0f));

        pRec.p = trafo.transformAffine(
            Point(p.x, p.y, 0.0f));
        pRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
        pRec.pdf = m_aperturePdf; /// XXX
        pRec.measure = EArea;
        return Spectrum(1.0f);
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");
        return Spectrum((pRec.measure == EArea) ? m_aperturePdf : 0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");
        return (pRec.measure == EArea) ? m_aperturePdf : 0.0f;
    }

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
            PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");
        const Transform &trafo = m_worldTransform->eval(pRec.time);

        /* Compute the corresponding position on the
           near plane (in local camera space) */
        Point nearP = m_sampleToCamera(Point(
            sample.x, sample.y, 0.0f));

        /* Turn that into a normalized ray direction */
        Vector d = normalize(Vector(nearP));
        dRec.d = trafo(d);
        dRec.measure = ESolidAngle;
        dRec.pdf = m_normalization / (d.z * d.z * d.z);

        return Spectrum(1.0f);
    }

    inline Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");
        if (dRec.measure != ESolidAngle)
            return 0.0f;

        return 1.0f;
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");
        if (dRec.measure != ESolidAngle)
            return Spectrum(0.0f);

        return Spectrum(1.0f);
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        Log(EError, "The telecentric lens camera is currently incompatible "
            "with bidirectional rendering algorithms!");
        Transform trafo    = m_worldTransform->eval(dRec.time),
                  invTrafo = trafo.inverse();

        Float f = m_focusDistance / m_scale.z,
              apertureRadius = m_apertureRadius / m_scale.x;

        Point localP = invTrafo.transformAffine(dRec.ref);

        Float dist = localP.z * m_scale.z;
        if (dist < m_nearClip || dist > m_farClip) {
            dRec.pdf = 0.0f;
            return Spectrum(0.0f);
        }

        /* Circle of confusion */
        Float radius = std::abs(localP.z - f) * apertureRadius/f;
        radius += apertureRadius;

        /* Sample the ray origin */
        Point2 disk = warp::squareToUniformDiskConcentric(sample);
        Point diskP(disk.x*radius+localP.x, disk.y*radius+localP.y, 0.0f);

        /* Compute the intersection with the focal plane */
        Vector localD = localP - diskP;
        Point intersection = diskP + localD * (f/localD.z);

        /* Determine the associated sample coordinates */
        Point uv = m_cameraToSample.transformAffine(intersection);
        if (uv.x < 0 || uv.x > 1 || uv.y < 0 || uv.y > 1) {
            dRec.pdf = 0.0f;
            return Spectrum(0.0f);
        }

        dRec.uv = Point2(uv.x, uv.y);
        dRec.p = trafo(diskP);
        dRec.n = normalize(trafo(Vector(0, 0, 1.0f)));
        dRec.d = dRec.p - dRec.ref;
        dRec.dist = dRec.d.length();
        dRec.d /= dRec.dist;
        dRec.measure = ESolidAngle;

        dRec.pdf = dist*dist / (-dot(dRec.n, dRec.d)
            * M_PI * radius*radius);

        return Spectrum(m_normalization);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return 0.0f;
    }

    Transform getProjectionTransform(const Point2 &apertureSample,
            const Point2 &aaSample) const {
        Point2 offset(
            2.0f * m_invResolution.x * (aaSample.x-.5f),
            2.0f * m_invResolution.y * (aaSample.y-.5f));

        Float angle1 = (apertureSample.x-.5f)*m_maxRotation;
        Float angle2 = (apertureSample.y-.5f)*m_maxRotation;

        return m_clipTransform *
            Transform::translate(Vector(offset.x, offset.y, 0.0f)) *
            Transform::scale(Vector(1.0f, m_aspect, 1.0f)) *
            Transform::glOrthographic(m_nearClip, m_farClip)*
            Transform::scale(Vector(1.0f, 1.0f, m_scale.z)) *
            Transform::translate(Vector(0.0f, 0.0f, -m_focusDistance)) *
            Transform::rotate(Vector(1.0f, 0.0f, 0.0f), angle1) *
            Transform::rotate(Vector(0.0f, 1.0f, 0.0f), angle2) *
            Transform::translate(Vector(0.0f, 0.0f, m_focusDistance));
    }

    AABB getAABB() const {
        AABB bounds;
        bounds.expandBy(m_sampleToCamera(Point(0, 0, 0)));
        bounds.expandBy(m_sampleToCamera(Point(1, 1, 0)));

        return m_worldTransform->getSpatialBounds(bounds);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "TelecentricLensCamera[" << endl
            << "  apertureRadius = " << m_apertureRadius << "," << endl
            << "  focusDistance = " << m_focusDistance << "," << endl
            << "  nearClip = " << m_nearClip << "," << endl
            << "  farClip = " << m_farClip << "," << endl
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
private:
    Transform m_cameraToSample;
    Transform m_sampleToCamera;
    Transform m_clipTransform;
    Float m_apertureRadius;
    Float m_aperturePdf;
    Float m_normalization;
    Vector m_scale;
    Float m_maxRotation;
    Vector m_dx, m_dy;
};

MTS_IMPLEMENT_CLASS_S(TelecentricLensCamera, false, ProjectiveCamera)
MTS_EXPORT_PLUGIN(TelecentricLensCamera, "Telecentric lens camera");
MTS_NAMESPACE_END
