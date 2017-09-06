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
#include <mitsuba/render/shape.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{thinlens}{Perspective camera with a thin lens}
 * \order{2}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional camera-to-world transformation.
 *        \default{none (i.e. camera space $=$ world space)}
 *     }
 *     \parameter{apertureRadius}{\Float}{
 *         Denotes the radius of the camera's aperture in scene units.
 *     }
 *     \parameter{focusDistance}{\Float}{
 *         Denotes the world-space distance from the camera's aperture to the
 *         focal plane. \default{\code{0}}
 *     }
 *     \parameter{focalLength}{\String}{
 *         Denotes the camera's focal length specified using
 *         \code{35mm} film equivalent units. See the main description
 *         for further details.\default{\code{50mm}}
 *     }
 *     \parameter{fov}{\Float}{
 *         An alternative to \code{focalLength}:
 *         denotes the camera's field of view in degrees---must be
 *         between 0 and 180, excluding the extremes.
 *     }
 *     \parameter{fovAxis}{\String}{
 *         When the parameter \code{fov} is given (and only then),
 *         this parameter further specifies the image axis, to
 *         which it applies.
 *         \begin{enumerate}[(i)]
 *             \item \code{\textbf{x}}: \code{fov} maps to the
 *                 \code{x}-axis in screen space.\vspace{-1mm}
 *             \item \code{\textbf{y}}: \code{fov} maps to the
 *                 \code{y}-axis in screen space.\vspace{-1mm}
 *             \item \code{\textbf{diagonal}}: \code{fov}
 *                maps to the screen diagonal.\vspace{-1mm}
 *             \item \code{\textbf{smaller}}: \code{fov}
 *                maps to the smaller dimension
 *                (e.g. \code{x} when \code{width}<\code{height})\vspace{-1mm}
 *             \item \code{\textbf{larger}}: \code{fov}
 *                maps to the larger dimension
 *                (e.g. \code{y} when \code{width}<\code{height})
 *                \vspace{-1mm}
 *         \end{enumerate}
 *         The default is \code{\textbf{x}}.
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
 *\vspace{-4mm}
 * \rendering{The material test ball viewed through a perspective thin lens
 * camera. Points away from the focal plane project onto a circle of confusion.}{sensor_thinlens}
 * \medrendering{A rendering of the Cornell box}{sensor_thinlens_2}
 * }
 *
 * This plugin implements a simple perspective camera model with a thin lens
 * at its circular aperture. It is very similar to the \pluginref{perspective} plugin
 * except that the extra lens element permits rendering with a
 * specifiable (i.e. non-infinite) depth of field. To configure this, it has two
 * extra parameters named \code{apertureRadius} and \code{focusDistance}.
 *
 * By default, the camera's field of view is specified using a 35mm film
 * equivalent focal length, which is first converted into a diagonal field
 * of view and subsequently applied to the camera. This assumes that
 * the film's aspect ratio matches that of 35mm film (1.5:1), though the
 * parameter still behaves intuitively when this is not the case.
 * Alternatively, it is also possible to specify a field of view in degrees
 * along a given axis (see the \code{fov} and \code{fovAxis} parameters).
 *
 * The exact camera position and orientation is most easily expressed using the
 * \code{lookat} tag, i.e.:
 * \begin{xml}
 * <sensor type="thinlens">
 *     <transform name="toWorld">
 *         <!-- Move and rotate the camera so that looks from (1, 1, 1) to (1, 2, 1)
 *              and the direction (0, 0, 1) points "up" in the output image -->
 *         <lookat origin="1, 1, 1" target="1, 2, 1" up="0, 0, 1"/>
 *     </transform>
 *
 *     <!-- Focus on the target -->
 *     <float name="focusDistance" value="1"/>
 *     <float name="apertureRadius" value="0.1"/>
 * </sensor>
 * \end{xml}
 */
class ThinLens : public PerspectiveCamera {
public:
    ThinLens(const Properties &props)
            : PerspectiveCamera(props) {
        m_type |= ENeedsApertureSample
               | EPerspectiveCamera
               | EOnSurface
               | EDirectionSampleMapsToPixels;

        /* World-space aperture radius */
        m_apertureRadius = props.getFloat("apertureRadius");

        if (m_apertureRadius == 0) {
            Log(EWarn, "Can't have a zero aperture radius -- "
                "setting to %f", Epsilon);
            m_apertureRadius = Epsilon;
        }

        if (props.getAnimatedTransform("toWorld", Transform())->eval(0).hasScale())
            Log(EError, "Scale factors in the camera-to-world "
                "transformation are not allowed!");
    }

    ThinLens(Stream *stream, InstanceManager *manager)
            : PerspectiveCamera(stream, manager) {
        m_apertureRadius = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        PerspectiveCamera::serialize(stream, manager);
        stream->writeFloat(m_apertureRadius);
    }

    void configure() {
        PerspectiveCamera::configure();

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
            * Transform::perspective(m_xfov, m_nearClip, m_farClip);

        m_sampleToCamera = m_cameraToSample.inverse();

        /* Position differentials on the near plane */
        m_dx = m_sampleToCamera(Point(m_invResolution.x, 0.0f, 0.0f))
             - m_sampleToCamera(Point(0.0f));
        m_dy = m_sampleToCamera(Point(0.0f, m_invResolution.y, 0.0f))
             - m_sampleToCamera(Point(0.0f));

        /* Precompute some data for importance(). Please
           look at that function for further details */
        Point min(m_sampleToCamera(Point(0, 0, 0))),
              max(m_sampleToCamera(Point(1, 1, 0)));

        m_imageRect.reset();
        m_imageRect.expandBy(Point2(min.x, min.y) / min.z);
        m_imageRect.expandBy(Point2(max.x, max.y) / max.z);
        m_normalization = 1.0f / m_imageRect.getVolume();

        /* Clip-space transformation for OpenGL */
        m_clipTransform = Transform::translate(
            Vector((1-2*relOffset.x)/relSize.x - 1,
                  -(1-2*relOffset.y)/relSize.y + 1, 0.0f)) *
            Transform::scale(Vector(1.0f / relSize.x, 1.0f / relSize.y, 1.0f));

        m_aperturePdf = 1 / (M_PI * m_apertureRadius * m_apertureRadius);
    }

    /**
     * \brief Compute the directional sensor response function
     * of the camera multiplied with the cosine foreshortening
     * factor associated with the image plane
     *
     * \param p
     *     A position on the aperture
     *
     * \param d
     *     A normalized direction vector from the aperture position to the
     *     reference point in question (all in local camera space)
     *
     * \param sample
     *     Optional: a pointer to a 2D point data structure, which will be
     *     used to return the fractional pixel position associated with the
     *     reference point
     */
    inline Float importance(const Point &p, const Vector &d, Point2 *sample = NULL) const {
        /* How is this derived? Imagine a hypothetical image plane at a
           distance of d=1 away from the aperture in camera space.

           Then the visible rectangular portion of the plane has the area

              A = (2 * tan(0.5 * xfov in radians))^2 / aspect

           Since we allow crop regions, the actual visible area is
           potentially reduced:

              A' = A * (cropX / filmX) * (cropY / filmY)

           Perspective transformations of such aligned rectangles produce
           an equivalent scaled (but otherwise undistorted) rectangle
           in screen space. This means that a strategy, which uniformly
           generates samples in screen space has an associated area
           density of 1/A' on this rectangle.

           To compute the solid angle density of a sampled point P on
           the rectangle, we can apply the usual measure conversion term:

              d_omega = 1/A' * distance(P, origin)^2 / cos(theta)

           where theta is the angle that the unit direction vector from
           the origin to P makes with the rectangle. Since

              distance(P, origin)^2 = Px^2 + Py^2 + 1

           and

              cos(theta) = 1/sqrt(Px^2 + Py^2 + 1),

           we have

              d_omega = 1 / (A' * cos^3(theta))
        */

        Float cosTheta = Frame::cosTheta(d);

        /* Check if the direction points behind the camera */
        if (cosTheta <= 0)
            return 0.0f;

        Float invCosTheta = 1.0f / cosTheta;

        /* Check if the associated pixel is visible */
        Point scr = m_cameraToSample(p
            + d * (m_focusDistance*invCosTheta));
        if (scr.x < 0 || scr.x > 1 ||
            scr.y < 0 || scr.y > 1)
            return 0.0f;

        if (sample) {
            sample->x = scr.x * m_resolution.x;
            sample->y = scr.y * m_resolution.y;
        }

        return m_normalization * invCosTheta *
            invCosTheta * invCosTheta;
    }

    Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
            const Point2 &otherSample, Float timeSample) const {
        Point2 tmp = warp::squareToUniformDiskConcentric(otherSample)
            * m_apertureRadius;
        ray.time = sampleTime(timeSample);

        /* Compute the corresponding position on the
           near plane (in local camera space) */
        Point nearP = m_sampleToCamera(Point(
            pixelSample.x * m_invResolution.x,
            pixelSample.y * m_invResolution.y, 0.0f));

        /* Aperture position */
        Point apertureP(tmp.x, tmp.y, 0.0f);

        /* Sampled position on the focal plane */
        Point focusP = nearP * (m_focusDistance / nearP.z);

        /* Turn these into a normalized ray direction, and
           adjust the ray interval accordingly */
        Vector d = normalize(focusP - apertureP);
        Float invZ = 1.0f / d.z;
        ray.mint = m_nearClip * invZ;
        ray.maxt = m_farClip * invZ;

        const Transform &trafo = m_worldTransform->eval(ray.time);
        ray.setOrigin(trafo.transformAffine(apertureP));
        ray.setDirection(trafo(d));
        return Spectrum(1.0f);
    }

    Spectrum sampleRayDifferential(RayDifferential &ray, const Point2 &pixelSample,
            const Point2 &otherSample, Float timeSample) const {
        Point2 tmp = warp::squareToUniformDiskConcentric(otherSample)
            * m_apertureRadius;
        ray.time = sampleTime(timeSample);

        /* Compute the corresponding position on the
           near plane (in local camera space) */
        Point nearP = m_sampleToCamera(Point(
            pixelSample.x * m_invResolution.x,
            pixelSample.y * m_invResolution.y, 0.0f));

        /* Aperture position */
        Point apertureP(tmp.x, tmp.y, 0.0f);

        /* Sampled position on the focal plane */
        Float fDist = m_focusDistance / nearP.z;
        Point focusP  =  nearP       * fDist;
        Point focusPx = (nearP+m_dx) * fDist;
        Point focusPy = (nearP+m_dy) * fDist;

        /* Turn that into a normalized ray direction, and
           adjust the ray interval accordingly */
        Vector d = normalize(focusP - apertureP);
        Float invZ = 1.0f / d.z;
        ray.mint = m_nearClip * invZ;
        ray.maxt = m_farClip * invZ;

        const Transform &trafo = m_worldTransform->eval(ray.time);
        ray.setOrigin(trafo.transformAffine(apertureP));
        ray.setDirection(trafo(d));
        ray.rxOrigin = ray.ryOrigin = ray.o;
        ray.rxDirection = trafo(normalize(Vector(focusPx - apertureP)));
        ray.ryDirection = trafo(normalize(Vector(focusPy - apertureP)));
        ray.hasDifferentials = true;

        return Spectrum(1.0f);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec,
            const Point2 &sample, const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);

        Point2 aperturePos = warp::squareToUniformDiskConcentric(sample)
            * m_apertureRadius;

        pRec.p = trafo.transformAffine(
            Point(aperturePos.x, aperturePos.y, 0.0f));
        pRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
        pRec.pdf = m_aperturePdf;
        pRec.measure = EArea;
        return Spectrum(1.0f);
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return Spectrum((pRec.measure == EArea) ? m_aperturePdf : 0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EArea) ? m_aperturePdf : 0.0f;
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

        /* Compute the corresponding position on the
           near plane (in local camera space) */
        Point nearP = m_sampleToCamera(samplePos);
        nearP.x = nearP.x * (m_focusDistance / nearP.z);
        nearP.y = nearP.y * (m_focusDistance / nearP.z);
        nearP.z = m_focusDistance;

        Point apertureP = trafo.inverse().transformAffine(pRec.p);

        /* Turn that into a normalized ray direction */
        Vector d = normalize(nearP - apertureP);
        dRec.d = trafo(d);
        dRec.measure = ESolidAngle;
        dRec.pdf = m_normalization / (d.z * d.z * d.z);

        return Spectrum(1.0f);
    }

    inline Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle)
            return 0.0f;

        Transform invTrafo = m_worldTransform->eval(pRec.time).inverse();
        return importance(invTrafo.transformAffine(pRec.p), invTrafo(dRec.d));
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle)
            return Spectrum(0.0f);

        Transform invTrafo = m_worldTransform->eval(pRec.time).inverse();
        return Spectrum(importance(invTrafo.transformAffine(pRec.p),
                    invTrafo(dRec.d)));
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);

        /* Transform the reference point into the local coordinate system */
        Point refP = trafo.inverse().transformAffine(dRec.ref);

        /* Check if it is outside of the clip range */
        if (refP.z < m_nearClip || refP.z > m_farClip) {
            dRec.pdf = 0.0f;
            return Spectrum(0.0f);
        }

        /* Sample a position on the aperture (in local coordinates) */
        Point2 tmp = warp::squareToUniformDiskConcentric(sample)
            * m_apertureRadius;
        Point apertureP(tmp.x, tmp.y, 0);

        /* Compute the normalized direction vector from the
           aperture position to the reference point */
        Vector localD(refP - apertureP);
        Float dist = localD.length(),
              invDist = 1.0f / dist;
        localD *= invDist;

        Float value = importance(apertureP, localD, &dRec.uv);
        if (value == 0.0f) {
            dRec.pdf = 0.0f;
            return Spectrum(0.0f);
        }

        dRec.p = trafo.transformAffine(apertureP);
        dRec.d = (dRec.p - dRec.ref) * invDist;
        dRec.dist = dist;
        dRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
        dRec.pdf = m_aperturePdf * dist*dist/(Frame::cosTheta(localD));
        dRec.measure = ESolidAngle;

        /* intentionally missing a cosine factor wrt. the aperture
           disk (it is already accounted for in importance()) */
        return Spectrum(value * invDist * invDist);
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        Float dp = -dot(dRec.n, dRec.d);
        if (dp < 0)
            return 0.0f;

        if (dRec.measure == ESolidAngle)
            return m_aperturePdf * dRec.dist*dRec.dist / dp;
        else if (dRec.measure == EArea)
            return m_aperturePdf;
        else
            return 0.0f;
    }

    Transform getProjectionTransform(const Point2 &apertureSample,
            const Point2 &aaSample) const {
        Float right = std::tan(m_xfov * M_PI/360) * m_nearClip, left = -right;
        Float top = right / m_aspect, bottom = -top;
        Point2 apertureP = warp::squareToUniformDiskConcentric(apertureSample)
            * m_apertureRadius;

        Vector2 jitterScale(
            (top-bottom)/m_film->getSize().y);

        Point2 offset = apertureP * (m_nearClip/m_focusDistance);

        offset += Vector2(
            (right-left)/m_film->getSize().x * (aaSample.x-0.5f),
            (top-bottom)/m_film->getSize().y * (aaSample.y-0.5f));

        return m_clipTransform *
            Transform::glFrustum(left+offset.x, right+offset.x,
                bottom+offset.y, top+offset.y, m_nearClip, m_farClip)
          * Transform::translate(Vector(apertureP.x, apertureP.y, 0));
    }

    AABB getAABB() const {
        const Float r = m_apertureRadius;
        AABB bounds(Point(-r, -r, 0), Point(r, r, 0));
        return m_worldTransform->getSpatialBounds(bounds);
    }

    ref<Shape> createShape(const Scene *scene) {
        ref<AnimatedTransform> trafo = new AnimatedTransform(m_worldTransform);
        trafo->prependScale(Vector(m_apertureRadius));

        Properties props("disk");
        props.setAnimatedTransform("toWorld", trafo);
        Shape *shape = static_cast<Shape *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Shape), props));
        shape->addChild(this);
        shape->configure();

        return shape;
    }

    bool getSamplePosition(const PositionSamplingRecord &pRec,
            const DirectionSamplingRecord &dRec, Point2 &samplePosition) const {
        Transform invTrafo = m_worldTransform->eval(pRec.time).inverse();
        Point localP(invTrafo.transformAffine(pRec.p));
        Point localD(invTrafo(dRec.d));

        if (localD.z <= 0)
            return false;

        Point intersection = localP + localD * (m_focusDistance / localD.z);

        Point screenSample = m_cameraToSample(intersection);
        if (screenSample.x < 0 || screenSample.x > 1 ||
            screenSample.y < 0 || screenSample.y > 1)
            return false;

        samplePosition = Point2(
                screenSample.x * m_resolution.x,
                screenSample.y * m_resolution.y);

        return true;
    }

    Spectrum eval(const Intersection &its, const Vector &d, Point2 &samplePos) const {
        Transform invTrafo = m_worldTransform->eval(its.time).inverse();
        Point localP = invTrafo.transformAffine(its.p);
        Vector localD = invTrafo(d);

        Float result = importance(localP, localD, &samplePos);

        if (result == 0)
            return Spectrum(0.0f);

        return Spectrum(result * m_aperturePdf / Frame::cosTheta(localD));
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ThinLens[" << endl
            << "  fov = [" << getXFov() << ", " << getYFov() << "]," << endl
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
    AABB2 m_imageRect;
    Float m_apertureRadius;
    Float m_aperturePdf;
    Float m_normalization;
    Vector m_dx, m_dy;
};

MTS_IMPLEMENT_CLASS_S(ThinLens, false, PerspectiveCamera)
MTS_EXPORT_PLUGIN(ThinLens, "Perspective thin lens camera");
MTS_NAMESPACE_END
