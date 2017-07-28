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
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include "sunsky/sunmodel.h"
#include "sunsky/skymodel.h"

MTS_NAMESPACE_BEGIN

#if SPECTRUM_SAMPLES == 3
# define SKY_PIXELFORMAT Bitmap::ERGB
#else
# define SKY_PIXELFORMAT Bitmap::ESpectrum
#endif

/*!\plugin{sky}{Skylight emitter}
 * \icon{emitter_sky}
 * \order{6}
 * \parameters{
 *     \parameter{turbidity}{\Float}{
 *         This parameter determines the amount of aerosol present
 *         in the atmosphere.
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
 *     \parameter{scale}{\Float}{
 *         This parameter can be used to scale the amount of illumination
 *         emitted by the sky emitter. \default{1}
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional sensor-to-world transformation.
 *        \default{none (i.e. sensor space $=$ world space)}
 *     }
 * }
 *
 * \renderings{
 *     \tinyrendering{5AM}{emitter_sky_small_5}
 *     \tinyrendering{7AM}{emitter_sky_small_7}
 *     \tinyrendering{9AM}{emitter_sky_small_9}
 *     \tinyrendering{11AM}{emitter_sky_small_11}
 *     \tinyrendering{1PM}{emitter_sky_small_13}
 *     \tinyrendering{3PM}{emitter_sky_small_15}
 *     \tinyrendering{5PM}{emitter_sky_small_17}
 *     \tinyrendering{6:30 PM}{emitter_sky_small_1830}\hfill
 *     \vspace{-2mm}
 *     \caption{Time series at the default settings (Equidistant fisheye
 *     projection of the sky onto a disk. East is left.)}
 *     \vspace{3mm}
 * }
 *
 * This plugin provides the physically-based skylight model by
 * Ho\v{s}ek and Wilkie \cite{Hosek2012Analytic}. It can be used to
 * create predictive daylight renderings of scenes under clear skies,
 * which is useful for architectural and computer vision applications.
 * The implementation in Mitsuba is based on code that was
 * generously provided by the authors.
 *
 * The model has two main parameters: the turbidity of the atmosphere
 * and the position of the sun.
 * The position of the sun in turn depends on a number of secondary
 * parameters, including the \code{latitude}, \code{longitude},
 * and \code{timezone} at the location of the observer, as well as the
 * current \code{year}, \code{month}, \code{day}, \code{hour},
 * \code{minute}, and \code{second}.
 * Using all of these, the elevation and azimuth of the sun are computed
 * using the PSA algorithm by Blanco et al. \cite{Blanco2001Computing},
 * which is accurate to about 0.5 arcminutes (\nicefrac{1}{120} degrees).
 * Note that this algorithm does not account for daylight
 * savings time where it is used, hence a manual correction of the
 * time may be necessary.
 * For detailed coordinate and timezone information of various cities, see
 * \url{http://www.esrl.noaa.gov/gmd/grad/solcalc}.
 *
 * If desired, the world-space solar vector may also be specified
 * using the \code{sunDirection} parameter, in which case all of the
 * previously mentioned time and location parameters become irrelevant.
 *
 * \renderings{
 *     \tinyrendering{1}{emitter_sky_turb_1}
 *     \tinyrendering{2}{emitter_sky_turb_2}
 *     \tinyrendering{3}{emitter_sky_turb_3}
 *     \tinyrendering{4}{emitter_sky_turb_4}
 *     \tinyrendering{5}{emitter_sky_turb_5}
 *     \tinyrendering{6}{emitter_sky_turb_6}
 *     \tinyrendering{8}{emitter_sky_turb_8}
 *     \tinyrendering{10}{emitter_sky_turb_10}
 *     \caption{Sky light for different turbidity values (default configuration at 5PM)}
 * }
 *
 * \emph{Turbidity}, the other important parameter, specifies the aerosol
 * content of the atmosphere. Aerosol particles cause additional scattering that
 * manifests in a halo around the sun, as well as color fringes near the
 * horizon.
 * Smaller turbidity values ($\sim 1-2$) produce an
 * arctic-like clear blue sky, whereas larger values ($\sim 8-10$)
 * create an atmosphere that is more typical of a warm, humid day.
 * Note that this model does not aim to reproduce overcast, cloudy, or foggy
 * atmospheres with high corresponding turbidity values. An photographic
 * environment map may be more appropriate in such cases.

 * The default coordinate system of the emitter associates the up
 * direction with the $+Y$ axis. The east direction is associated with $+X$
 * and the north direction is equal to $+Z$. To change this coordinate
 * system, rotations can be applied using the \code{toWorld} parameter
 * (see \lstref{sky-up} for an example).
 *
 * By default, the emitter will not emit any light below the
 * horizon, which means that these regions are black when
 * observed directly. By setting the \code{stretch} parameter to values
 * between $1$ and $2$, the sky can be extended to cover these directions
 * as well. This is of course a complete kludge and only meant as a quick
 * workaround for scenes that are not properly set up.
 *
 * Instead of evaluating the full sky model every on every radiance query,
 * the implementation precomputes a low resolution environment map
 * (512$\times$ 256) of the entire sky that is then forwarded to the
 * \pluginref{envmap} plugin---this dramatically improves rendering
 * performance. This resolution is generally plenty since the sky radiance
 * distribution is so smooth, but it can be adjusted manually if
 * necessary using the \code{resolution} parameter.
 *
 * Note that while the model encompasses sunrise and sunset configurations,
 * it does not extend to the night sky, where illumination from stars, galaxies,
 * and the moon dominate. When started with a sun configuration that lies
 * below the horizon, the plugin will fail with an error message.
 *
 * \vspace{5mm}
 * \begin{xml}[caption={Rotating the sky emitter for scenes that use $Z$ as
 * the ``up'' direction}, label=lst:sky-up]
 * <emitter type="sky">
 *     <transform name="toWorld">
 *         <rotate x="1" angle="90"/>
 *     </transform>
 * </emitter>
 * \end{xml}
 *
 * \subsubsection*{Physical units and spectral rendering}
 * Like the \code{blackbody} emission profile (Page~\pageref{sec:blackbody}),
 * the sky model introduces physical units into the rendering process.
 * The radiance values computed by this plugin have units of power ($W$) per
 * unit area ($m^{-2}$) per steradian ($sr^{-1}$) per unit wavelength ($nm^{-1}$).
 * If these units are inconsistent with your scene description, you may use the
 * optional \texttt{scale} parameter to adjust them.
 *
 * When Mitsuba is compiled for spectral rendering, the plugin switches
 * from RGB to a spectral variant of the skylight model, which relies on
 * precomputed data between $320$ and $720 nm$ sampled at $40nm$-increments.
 *
 * \subsubsection*{Ground albedo}
 * The albedo of the ground (e.g. due to rock, snow, or vegetation) can have a
 * noticeable and nonlinear effect on the appearance of the sky.
 * \figref{sky_groundalbedo} shows an example of this effect. By default,
 * the ground albedo is set to a 15% gray.
 *
 * \renderings{
 *    \rendering{3 PM}{emitter_sky_mattest_3pm}
 *    \rendering{6:30 PM}{emitter_sky_mattest_630pm}
 *    \vspace{-3mm}
 *    \caption{Renderings with the \pluginref{plastic} material
 *    under default conditions. Note that these only contain
 *    skylight illumination. For a model that also includes the sun,
 *    refer to \pluginref{sunsky}.}
 *    \vspace{-6mm}
 * }
 * \vspace{5mm}
 * \renderings{
 *   \medrendering{\code{albedo}=0%}{emitter_sky_albedo_0}
 *   \medrendering{\code{albedo}=100%}{emitter_sky_albedo_1}
 *   \medrendering{\code{albedo}=20% green}{emitter_sky_albedo_green}
 *   \caption{\label{fig:sky_groundalbedo}Influence
 *   of the ground albedo on the appearance of the sky}
 * }
 */
class SkyEmitter : public Emitter {
public:
    SkyEmitter(const Properties &props)
            : Emitter(props) {
        m_scale = props.getFloat("scale", 1.0f);
        m_turbidity = props.getFloat("turbidity", 3.0f);
        m_stretch = props.getFloat("stretch", 1.0f);
        m_resolution = props.getInteger("resolution", 512);
        m_albedo = props.getSpectrum("albedo", Spectrum(0.2f));
        m_sun = computeSunCoordinates(props);
        m_extend = props.getBoolean("extend", false);

        if (m_turbidity < 1 || m_turbidity > 10)
            Log(EError, "The turbidity parameter must be in the range [1,10]!");
        if (m_stretch < 1 || m_stretch > 2)
            Log(EError, "The stretch parameter must be in the range [1,2]!");
        for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
            if (m_albedo[i] < 0 || m_albedo[i] > 1)
                Log(EError, "The albedo parameter must be in the range [0,1]!");
        }

        Float sunElevation = 0.5f * M_PI - m_sun.elevation;

        if (sunElevation < 0)
            Log(EError, "The sun is below the horizon -- this is not supported by the sky model.");

        #if SPECTRUM_SAMPLES == 3
            for (int i=0; i<SPECTRUM_SAMPLES; ++i)
                m_state[i] = arhosek_rgb_skymodelstate_alloc_init(
                    m_turbidity, m_albedo[i], sunElevation);
        #else
            for (int i=0; i<SPECTRUM_SAMPLES; ++i)
                m_state[i] = arhosekskymodelstate_alloc_init(
                    m_turbidity, m_albedo[i], sunElevation);
        #endif

        configure();
    }

    SkyEmitter(Stream *stream, InstanceManager *manager)
            : Emitter(stream, manager) {
        m_scale = stream->readFloat();
        m_turbidity = stream->readFloat();
        m_stretch = stream->readFloat();
        m_resolution = stream->readInt();
        m_extend = stream->readBool();
        m_albedo = Spectrum(stream);
        m_sun = SphericalCoordinates(stream);

        Float sunElevation = 0.5f * M_PI - m_sun.elevation;
        #if SPECTRUM_SAMPLES == 3
            for (int i=0; i<SPECTRUM_SAMPLES; ++i)
                m_state[i] = arhosek_rgb_skymodelstate_alloc_init(
                    m_turbidity, m_albedo[i], sunElevation);
        #else
            for (int i=0; i<SPECTRUM_SAMPLES; ++i)
                m_state[i] = arhosekskymodelstate_alloc_init(
                    m_turbidity, m_albedo[i], sunElevation);
        #endif

        configure();
    }

    ~SkyEmitter() {
        #if SPECTRUM_SAMPLES == 3
            for (int i=0; i<SPECTRUM_SAMPLES; ++i)
                arhosek_tristim_skymodelstate_free(m_state[i]);
        #else
            for (int i=0; i<SPECTRUM_SAMPLES; ++i)
                arhosekskymodelstate_free(m_state[i]);
        #endif
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        stream->writeFloat(m_scale);
        stream->writeFloat(m_turbidity);
        stream->writeFloat(m_stretch);
        stream->writeInt(m_resolution);
        stream->writeBool(m_extend);
        m_albedo.serialize(stream);
        m_sun.serialize(stream);
    }

    bool isCompound() const {
        return true;
    }

    Emitter *getElement(size_t i) {
        if (i != 0)
            return NULL;

        ref<Timer> timer = new Timer();
        Log(EDebug, "Rasterizing skylight emitter to an %ix%i environment map ..",
                m_resolution, m_resolution/2);
        ref<Bitmap> bitmap = new Bitmap(SKY_PIXELFORMAT, Bitmap::EFloat,
            Vector2i(m_resolution, m_resolution/2));

        Point2 factor((2*M_PI) / bitmap->getWidth(),
            M_PI / bitmap->getHeight());

        #if defined(MTS_OPENMP)
            #pragma omp parallel for
        #endif
        for (int y=0; y<bitmap->getHeight(); ++y) {
            Float theta = (y+.5f) * factor.y;
            Spectrum *target = (Spectrum *) bitmap->getFloatData()
                + y * bitmap->getWidth();

            for (int x=0; x<bitmap->getWidth(); ++x) {
                Float phi = (x+.5f) * factor.x;

                *target++ = getSkyRadiance(SphericalCoordinates(theta, phi));
            }
        }

        Log(EDebug, "Done (took %i ms)", timer->getMilliseconds());

        #if defined(MTS_DEBUG_SUNSKY)
        /* Write a debug image for inspection */
        {
            int size = 513 /* odd-sized */, border = 2;
            int fsize = size+2*border, hsize = size/2;
            ref<Bitmap> debugBitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, Vector2i(fsize));
            debugBitmap->clear();

            #if defined(MTS_OPENMP)
                #pragma omp parallel for
            #endif
            for (int y=0; y<size; ++y) {
                float *target = debugBitmap->getFloat32Data() + ((y + border) * fsize + border) * 3;

                for (int x=0; x<size; ++x) {
                    Float xp = -(x - hsize) / (Float) hsize;
                    Float yp = -(y - hsize) / (Float) hsize;

                    Float radius = std::sqrt(xp*xp + yp*yp);

                    Spectrum result(0.0f);
                    if (radius < 1) {
                        Float theta = radius * 0.5f * M_PI;
                        Float phi = std::atan2(xp, yp);
                        result = getSkyRadiance(SphericalCoordinates(theta, phi));
                    }

                    Float r, g, b;
                    result.toLinearRGB(r, g, b);

                    *target++ = (float) r;
                    *target++ = (float) g;
                    *target++ = (float) b;
                }
            }

            ref<FileStream> fs = new FileStream("sky.exr", FileStream::ETruncReadWrite);
            debugBitmap->write(Bitmap::EOpenEXR, fs);
        }
        #endif

        /* Instantiate a nested environment map plugin */
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

    Spectrum evalEnvironment(const RayDifferential &ray) const {
        return getSkyRadiance(fromSphere(ray.d));
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SkyEmitter[" << endl
            << "  turbidity = " << m_turbidity << "," << endl
            << "  sunPos = " << m_sun.toString() << endl
            << "  resolution = " << m_resolution << endl
            << "  stretch = " << m_stretch << endl
            << "  scale = " << m_scale << endl
            << "]";
        return oss.str();
    }
protected:
    AABB getAABB() const {
        NotImplementedError("getAABB");
    }

    /// Calculates the spectral radiance of the sky in the specified direction.
    Spectrum getSkyRadiance(const SphericalCoordinates &coords) const {
        Float theta = coords.elevation / m_stretch;

        if (std::cos(theta) <= 0) {
            if (!m_extend)
                return Spectrum(0.0f);
            else
                theta = 0.5f * M_PI - Epsilon; /* super-unrealistic mode */
        }

        /* Compute the angle between the sun and (theta, phi) in radians */
        Float cosGamma = std::cos(theta) * std::cos(m_sun.elevation)
            + std::sin(theta) * std::sin(m_sun.elevation)
            * std::cos(coords.azimuth - m_sun.azimuth);

        Float gamma = math::safe_acos(cosGamma);

        Spectrum result;
        for (int i=0; i<SPECTRUM_SAMPLES; i++) {
            #if SPECTRUM_SAMPLES == 3
                result[i] = (Float) (arhosek_tristim_skymodel_radiance(m_state[i],
                    theta, gamma, i) / 106.856980); // (sum of Spectrum::CIE_Y)
            #else
                std::pair<Float, Float> bin = Spectrum::getBinCoverage(i);
                result[i] = (Float) arhosekskymodel_radiance(m_state[i],
                    theta, gamma, 0.5f * (bin.first + bin.second));
            #endif
        }

        result.clampNegative();

        if (m_extend)
            result *= math::smoothStep((Float) 0, (Float) 1, 2 - 2*coords.elevation*INV_PI);

        return result * m_scale;
    }

    MTS_DECLARE_CLASS()
protected:
    /// Environment map resolution in pixels
    int m_resolution;
    /// Constant scale factor applied to the model
    Float m_scale;
    /// Sky turbidity
    Float m_turbidity;
    /// Position of the sun in spherical coordinates
    SphericalCoordinates m_sun;
    /// Stretch factor to extend to the bottom hemisphere
    Float m_stretch;
    /// Extend to the bottom hemisphere (super-unrealistic mode)
    bool m_extend;
    /// Ground albedo
    Spectrum m_albedo;

    /// State vector for the sky model
    #if SPECTRUM_SAMPLES == 3
        ArHosekTristimSkyModelState *m_state[SPECTRUM_SAMPLES];
    #else
        ArHosekSkyModelState *m_state[SPECTRUM_SAMPLES];
    #endif
};

MTS_IMPLEMENT_CLASS_S(SkyEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(SkyEmitter, "Skylight emitter");
MTS_NAMESPACE_END
