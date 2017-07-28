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
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/qmc.h>
#include "sunsky/sunmodel.h"

#if SPECTRUM_SAMPLES == 3
# define SUNSKY_PIXELFORMAT Bitmap::ERGB
#else
# define SUNSKY_PIXELFORMAT Bitmap::ESpectrum
#endif

/* Apparent radius of the sun as seen from the earth (in degrees).
   This is an approximation--the actual value is somewhere between
   0.526 and 0.545 depending on the time of year */
#define SUN_APP_RADIUS 0.5358

MTS_NAMESPACE_BEGIN

/*!\plugin{sunsky}{Sun and sky emitter}
 * \icon{emitter_sunsky}
 * \order{8}
 * \parameters{
 *     \parameter{turbidity}{\Float}{
 *         This parameter determines the amount of aerosol present in the atmosphere.
 *         Valid range: 1-10. \default{3, corresponding to a clear sky in a temperate climate}
 *     }
 *     \parameter{albedo}{\Spectrum}{Specifies the ground albedo \default{0.15}}
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
 *     \parameter{stretch}{\Float}{
 *         Stretch factor to extend emitter below the horizon, must be
 *         in [1,2] \default{\code{1}, i.e. not used}
 *     }
 *     \parameter{resolution}{\Integer}{Specifies the horizontal resolution of the precomputed
 *         image that is used to represent the sun environment map \default{512, i.e. 512$\times$256}}
 *     \parameter{sunScale}{\Float}{
 *         This parameter can be used to separately scale the amount of illumination
 *         emitted by the sun. \default{1}
 *     }
 *     \parameter{skyScale}{\Float}{
 *         This parameter can be used to separately scale the amount of illumination
 *         emitted by the sky.\default{1}
 *     }
 *     \parameter{sunRadiusScale}{\Float}{
 *         Scale factor to adjust the radius of the sun, while preserving its power.
 *         Set to \code{0} to turn it into a directional light source.
 *     }
 * }
 * \vspace{-3mm}
 *
 * \renderings{
 *   \medrendering{\pluginref{sky} emitter}{emitter_sunsky_sky}
 *   \medrendering{\pluginref{sun} emitter}{emitter_sunsky_sun}
 *   \medrendering{\pluginref{sunsky} emitter}{emitter_sunsky_sunsky}
 *   \vspace{-2mm}
 *   \caption{A coated rough copper test ball lit with the three
 *   provided daylight illumination models}
 * }
 * \vspace{1mm}
 * This convenience plugin has the sole purpose of instantiating
 * \pluginref{sun} and \pluginref{sky} and merging them into a joint
 * environment map. Please refer to these plugins individually for more
 * details.
 */
class SunSkyEmitter : public Emitter {
public:
    SunSkyEmitter(const Properties &props)
        : Emitter(props) {
        Float scale = props.getFloat("scale", 1.0f),
              sunScale = props.getFloat("sunScale", scale),
              skyScale = props.getFloat("skyScale", scale),
              sunRadiusScale = props.getFloat("sunRadiusScale", 1.0f);

        const Transform &trafo = m_worldTransform->eval(0);

        Properties skyProps(props);
        skyProps.removeProperty("toWorld");
        if (props.hasProperty("sunDirection"))
            skyProps.setVector("sunDirection", trafo.inverse()(props.getVector("sunDirection")));
        skyProps.setPluginName("sky");
        skyProps.setFloat("scale", skyScale, false);

        ref<Emitter> sky = static_cast<Emitter *>(
            PluginManager::getInstance()->createObject(
            MTS_CLASS(Emitter), skyProps));
        sky->configure();
        props.markQueried("albedo");

        int resolution = props.getInteger("resolution", 512);
        ref<Bitmap> bitmap = new Bitmap(SUNSKY_PIXELFORMAT, Bitmap::EFloat,
            Vector2i(resolution, resolution/2));

        Point2 factor((2*M_PI) / bitmap->getWidth(),
            M_PI / bitmap->getHeight());

        ref<Timer> timer = new Timer();
        Log(EDebug, "Rasterizing sun & skylight emitter to an %ix%i environment map ..",
                resolution, resolution/2);

        Spectrum *data = (Spectrum *) bitmap->getFloatData();

        /* First, rasterize the sky */
        #if defined(MTS_OPENMP)
            #pragma omp parallel for
        #endif
        for (int y=0; y<bitmap->getHeight(); ++y) {
            Float theta = (y+.5f) * factor.y;
            Spectrum *target = data + y * bitmap->getWidth();

            for (int x=0; x<bitmap->getWidth(); ++x) {
                Float phi = (x+.5f) * factor.x;

                RayDifferential ray(Point(0.0f),
                    toSphere(SphericalCoordinates(theta, phi)), 0.0f);

                *target++ = sky->evalEnvironment(ray);
            }
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

        SphericalCoordinates sun = computeSunCoordinates(props);
        Spectrum sunRadiance = computeSunRadiance(sun.elevation,
            props.getFloat("turbidity", 3.0f)) * sunScale;
        sun.elevation *= props.getFloat("stretch", 1.0f);
        Frame sunFrame = Frame(toSphere(sun));

        Float theta = degToRad(SUN_APP_RADIUS * 0.5f);

        if (sunRadiusScale == 0) {
            Float solidAngle = 2 * M_PI * (1 - std::cos(theta));
            Properties props("directional");
            props.setVector("direction", -trafo(sunFrame.n));
            props.setFloat("samplingWeight", m_samplingWeight);
            props.setSpectrum("irradiance", sunRadiance * solidAngle);

            m_dirEmitter = static_cast<Emitter *>(
                PluginManager::getInstance()->createObject(
                MTS_CLASS(Emitter), props));
        } else {
            size_t pixelCount = resolution*resolution/2;
            Float cosTheta = std::cos(theta * sunRadiusScale);

            /* Ratio of the sphere that is covered by the sun */
            Float coveredPortion = 0.5f * (1 - cosTheta);

            /* Approx. number of samples that need to be generated,
               be very conservative */
            size_t nSamples = (size_t) std::max((Float) 100,
                (pixelCount * coveredPortion * 1000));

            factor = Point2(bitmap->getWidth() / (2*M_PI),
                bitmap->getHeight() / M_PI);

            Spectrum value =
                sunRadiance * (2 * M_PI * (1-std::cos(theta))) *
                static_cast<Float>(bitmap->getWidth() * bitmap->getHeight())
                / (2 * M_PI * M_PI * nSamples);

            for (size_t i=0; i<nSamples; ++i) {
                Vector dir = sunFrame.toWorld(
                    warp::squareToUniformCone(cosTheta, sample02(i)));

                Float sinTheta = math::safe_sqrt(1-dir.y*dir.y);
                SphericalCoordinates sphCoords = fromSphere(dir);

                Point2i pos(
                    std::min(std::max(0, (int) (sphCoords.azimuth * factor.x)), bitmap->getWidth()-1),
                    std::min(std::max(0, (int) (sphCoords.elevation * factor.y)), bitmap->getHeight()-1));

                data[pos.x + pos.y * bitmap->getWidth()] += value / std::max((Float) 1e-3f, sinTheta);
            }

        }

        Log(EDebug, "Done (took %i ms)", timer->getMilliseconds());

        /* Instantiate a nested envmap plugin */
        Properties envProps("envmap");
        Properties::Data bitmapData;
        bitmapData.ptr = (uint8_t *) bitmap.get();
        bitmapData.size = sizeof(Bitmap);
        envProps.setData("bitmap", bitmapData);
        envProps.setAnimatedTransform("toWorld", m_worldTransform);
        envProps.setFloat("samplingWeight", m_samplingWeight);
        m_envEmitter = static_cast<Emitter *>(
            PluginManager::getInstance()->createObject(
            MTS_CLASS(Emitter), envProps));

        #if 0
            /* For debugging purposes */
            ref<FileStream> fs = new FileStream("debug.exr", FileStream::ETruncReadWrite);
            bitmap->write(Bitmap::EOpenEXR, fs);
        #endif
    }

    SunSkyEmitter(Stream *stream, InstanceManager *manager)
        : Emitter(stream, manager) {
        m_envEmitter = static_cast<Emitter *>(manager->getInstance(stream));
        if (stream->readBool())
            m_dirEmitter = static_cast<Emitter *>(manager->getInstance(stream));
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        manager->serialize(stream, m_envEmitter.get());
        stream->writeBool(m_dirEmitter.get() != NULL);
        if (m_dirEmitter.get())
            manager->serialize(stream, m_dirEmitter.get());
    }

    void configure() {
        Emitter::configure();
        m_envEmitter->configure();
        if (m_dirEmitter)
            m_dirEmitter->configure();
    }

    bool isCompound() const {
        return true;
    }

    AABB getAABB() const {
        NotImplementedError("getAABB");
    }

    Emitter *getElement(size_t i) {
        if (i == 0)
            return m_envEmitter;
        else if (i == 1)
            return m_dirEmitter;
        else
            return NULL;
    }

    MTS_DECLARE_CLASS()
private:
    ref<Emitter> m_dirEmitter;
    ref<Emitter> m_envEmitter;
};

MTS_IMPLEMENT_CLASS_S(SunSkyEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(SunSkyEmitter, "Sun & sky emitter");
MTS_NAMESPACE_END
