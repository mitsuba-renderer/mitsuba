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
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/qmc.h>
#include "sunsky/sunmodel.h"

MTS_NAMESPACE_BEGIN

/* Apparent radius of the sun as seen from the earth (in degrees).
   This is an approximation--the actual value is somewhere between
   0.526 and 0.545 depending on the time of year */
#define SUN_APP_RADIUS 0.5358

#if SPECTRUM_SAMPLES == 3
# define SUN_PIXELFORMAT Bitmap::ERGB
#else
# define SUN_PIXELFORMAT Bitmap::ESpectrum
#endif

/*!\plugin{sun}{Sun emitter}
 * \icon{emitter_sun}
 * \order{7}
 * \parameters{
 *     \parameter{turbidity}{\Float}{
 *         This parameter determines the amount of aerosol present in the atmosphere.
 *         Valid range: 2-10. \default{3, corresponding to a clear sky in a temperate climate}
 *     }
 *     \parameter{year, month, day}{\Integer}{Denote the date of the
 *      observation \default{2010, 07, 10}}
 *     \parameter{hour,minute,\showbreak second}{\Float}{Local time
 *       at the location of the observer in 24-hour format\default{15, 00, 00,
 *       i.e. 3PM}}
 *     \parameter{latitude, longitude, timezone}{\Float}{
 *       These three parameters specify the oberver's latitude and longitude
 *       in degrees, and the local timezone offset in hours, which are required
 *       to compute the sun's position. \default{35.6894, 139.6917, 9 --- Tokyo, Japan}
 *     }
 *     \parameter{sunDirection}{\Vector}{Allows to manually
 *       override the sun direction in world space. When this value
 *       is provided, parameters pertaining to the computation
 *       of the sun direction (\code{year, hour, latitude,} etc.
 *       are unnecessary. \default{none}
 *     }
 *     \parameter{resolution}{\Integer}{Specifies the horizontal resolution of the precomputed
 *         image that is used to represent the sun environment map \default{512, i.e. 512$\times$256}}
 *     \parameter{scale}{\Float}{
 *         This parameter can be used to scale the amount of illumination
 *         emitted by the sun emitter. \default{1}
 *     }
 *     \parameter{sunRadiusScale}{\Float}{
 *         Scale factor to adjust the radius of the sun, while preserving its power.
 *         Set to \code{0} to turn it into a directional light source.
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 * This plugin implements the physically-based sun model proposed by
 * Preetham et al. \cite{Preetham1999Practical}. Using the provided position
 * and time information (see \pluginref{sky} for details), it can determine the
 * position of the sun as seen from the position of the observer.
 * The radiance arriving at the earth surface is then found based on the spectral
 * emission profile of the sun and the extinction cross-section of the
 * atmosphere (which depends on the \code{turbidity} and the zenith angle of the sun).
 *
 * Like the \code{blackbody} emission profile (Page~\pageref{sec:blackbody}),
 * the sun model introduces physical units into the rendering process.
 * The radiance values computed by this plugin have units of power ($W$) per
 * unit area ($m^{-2}$) per steradian ($sr^{-1}$) per unit wavelength ($nm^{-1}$).
 * If these units are inconsistent with your scene description, you may use the
 * optional \texttt{scale} parameter to adjust them.
 *
 * This plugin supplies proper spectral power distributions when Mitsuba is
 * compiled in spectral rendering mode. Otherwise, they are simply projected onto
 * a linear RGB color space.
 *
 * \remarks{
 *   \item The sun is an intense light source that subtends a tiny solid angle.
 *   This can be a problem for certain rendering techniques (e.g. path
 *   tracing), which produce high variance output (i.e. noise in renderings)
 *   when the scene also contains specular or glossy or materials.
 * }
 */
class SunEmitter : public Emitter {
public:
    SunEmitter(const Properties &props)
            : Emitter(props) {
        m_scale = props.getFloat("scale", 1.0f);
        m_resolution = props.getInteger("resolution", 512);
        m_sun = computeSunCoordinates(props);
        m_sunRadiusScale = props.getFloat("sunRadiusScale", 1.0f);
        m_turbidity = props.getFloat("turbidity", 3.0f);
        m_stretch = props.getFloat("stretch", 1.0f);
    }

    SunEmitter(Stream *stream, InstanceManager *manager)
            : Emitter(stream, manager) {
        m_scale = stream->readFloat();
        m_sunRadiusScale = stream->readFloat();
        m_turbidity = stream->readFloat();
        m_resolution = stream->readInt();
        m_sun = SphericalCoordinates(stream);
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        stream->writeFloat(m_scale);
        stream->writeFloat(m_sunRadiusScale);
        stream->writeFloat(m_turbidity);
        stream->writeInt(m_resolution);
        m_sun.serialize(stream);
    }

    void configure() {
        SphericalCoordinates sun(m_sun);
        sun.elevation *= m_stretch;
        m_sunDir = toSphere(sun);

        /* Solid angle covered by the sun */
        m_theta = degToRad(SUN_APP_RADIUS * 0.5f);
        m_solidAngle = 2 * M_PI * (1 - std::cos(m_theta));
        m_radiance = computeSunRadiance(m_sun.elevation, m_turbidity) * m_scale;
    }

    bool isCompound() const {
        return true;
    }

    Emitter *getElement(size_t i) {
        if (i != 0)
            return NULL;

        if (m_sunRadiusScale == 0) {
            Properties props("directional");
            const Transform &trafo = m_worldTransform->eval(0);
            props.setVector("direction", -trafo(m_sunDir));
            props.setFloat("samplingWeight", m_samplingWeight);

            props.setSpectrum("irradiance", m_radiance * m_solidAngle);

            Emitter *emitter = static_cast<Emitter *>(
                PluginManager::getInstance()->createObject(
                MTS_CLASS(Emitter), props));

            emitter->configure();
            return emitter;
        }

        /* Rasterizing the sphere to an environment map and checking the
           individual pixels for coverage (which is what Mitsuba 0.3.0 did)
           was slow and not very effective; for instance the power varied
           dramatically with resolution changes. Since the sphere generally
           just covers a few pixels, the code below rasterizes it much more
           efficiently by generating a few thousand QMC samples.

           Step 1: compute a *very* rough estimate of how many
           pixel in the output environment map will be covered
           by the sun */
        size_t pixelCount = m_resolution*m_resolution/2;
        Float cosTheta = std::cos(m_theta * m_sunRadiusScale);

        /* Ratio of the sphere that is covered by the sun */
        Float coveredPortion = 0.5f * (1 - cosTheta);

        /* Approx. number of samples that need to be generated,
           be very conservative */
        size_t nSamples = (size_t) std::max((Float) 100,
            (pixelCount * coveredPortion * 1000));

        ref<Bitmap> bitmap = new Bitmap(SUN_PIXELFORMAT, Bitmap::EFloat,
            Vector2i(m_resolution, m_resolution/2));
        bitmap->clear();
        Frame frame(m_sunDir);

        Point2 factor(bitmap->getWidth() / (2*M_PI),
            bitmap->getHeight() / M_PI);

        Spectrum *target = (Spectrum *) bitmap->getFloatData();
        Spectrum value =
            m_radiance * (2 * M_PI * (1-std::cos(m_theta))) *
            static_cast<Float>(bitmap->getWidth() * bitmap->getHeight())
            / (2 * M_PI * M_PI * nSamples);

        for (size_t i=0; i<nSamples; ++i) {
            Vector dir = frame.toWorld(
                warp::squareToUniformCone(cosTheta, sample02(i)));

            Float sinTheta = math::safe_sqrt(1-dir.y*dir.y);
            SphericalCoordinates sphCoords = fromSphere(dir);

            Point2i pos(
                std::min(std::max(0, (int) (sphCoords.azimuth * factor.x)), bitmap->getWidth()-1),
                std::min(std::max(0, (int) (sphCoords.elevation * factor.y)), bitmap->getHeight()-1));

            target[pos.x + pos.y * bitmap->getWidth()] += value / std::max((Float) 1e-3f, sinTheta);
        }

        /* Instantiate a nested envmap plugin */
        Properties props("envmap");
        Properties::Data bitmapData;
        bitmapData.ptr = (uint8_t *) bitmap.get();
        bitmapData.size = sizeof(Bitmap);
        props.setData("bitmap", bitmapData);
        props.setAnimatedTransform("toWorld", m_worldTransform);
        props.setFloat("samplingWeight", m_samplingWeight);
        Emitter *emitter = static_cast<Emitter *>(
            PluginManager::getInstance()->createObject(
            MTS_CLASS(Emitter), props));
        emitter->configure();
        return emitter;
    }

    AABB getAABB() const {
        NotImplementedError("getAABB");
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SunEmitter[" << endl
            << "  sunDir = " << m_sunDir.toString() << "," << endl
            << "  sunRadiusScale = " << m_sunRadiusScale << "," << endl
            << "  turbidity = " << m_turbidity << "," << endl
            << "  scale = " << m_scale << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    /// Environment map resolution
    int m_resolution;
    /// Constant scale factor applied to the model
    Float m_scale;
    /// Scale factor that can be applied to the sun radius
    Float m_sunRadiusScale;
    /// Angle cutoff for the sun disk (w/o scaling)
    Float m_theta;
    /// Solid angle covered by the sun (w/o scaling)
    Float m_solidAngle;
    /// Position of the sun in spherical coordinates
    SphericalCoordinates m_sun;
    /// Direction of the sun (untransformed)
    Vector m_sunDir;
    /// Turbidity of the atmosphere
    Float m_turbidity;
    /// Radiance arriving from the sun disk
    Spectrum m_radiance;
    /// Stretch factor to extend to the bottom hemisphere
    Float m_stretch;
};

MTS_IMPLEMENT_CLASS_S(SunEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(SunEmitter, "Sun emitter");
MTS_NAMESPACE_END

