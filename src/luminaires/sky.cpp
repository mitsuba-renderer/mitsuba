/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2010 by Wenzel Jakob and others.

	Mitsuba is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Mitsuba is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include "sun.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{sky}{Skylight luminaire}
 * \parameters{
 *     \parameter{turbidity}{\Float}{
 *         This parameter determines the amount of scattering particles (or
 *         `haze') in the atmosphere. Smaller values ($\sim 2$) produce a 
 *         clear blue sky, larger values ($\sim 8$) lead to an overcast sky, 
 *         and a very high values ($\sim 20$) cause a color shift towards 
 *         orange and red. \default{3}
 *     }
 *     \parameter{day}{\Integer}{Solar day used to compute the sun's position. 
 *       Must be in the range between 1 and 365. \default{180}}
 *     \parameter{time}{\Float}{Fractional time used to compute the sun's
 *       position. A time of 4:15 PM corresponds to 16.25. \default{15.00}}
 *     \parameter{latitude, longitude}{\Float}{
 *       These two parameters specify the oberver's latitude and longitude 
 *       in degrees, which are required to compute the sun's position.
 *       \default{35.6894, 139.6917 --- Tokyo, Japan}
 *     }
 *     \parameter{standardMeridian}{\Integer}{Denotes the
 *       standard meridian of the time zone for finding
 *       the sun's position \default{135 --- Japan standard time}
 *     }
 *     \parameter{sunDirection}{\Vector}{Allows to manually 
 *       override the sun direction in world space. When this value
 *       is provided, parameters pertaining to the computation 
 *       of the sun direction (\code{day, time, latitude, longitude,} 
 *       and \code{standardMeridian}) are unnecessary. \default{none}
 *     }
 *     \parameter{extend}{\Boolean}{
 *         Extend luminaire below the horizon? \default{\code{false}}
 *     }
 *     \parameter{resolution}{\Integer}{Specifies the resolution of the precomputed
 *         image that is used to represent the sky environment map
 *         \default{256}}
 *     \parameter{skyScale}{\Float}{
 *         This parameter can be used to scale the the amount of illumination
 *         emitted by the sky luminaire, for instance to change its units. To
 *         switch from photometric ($\nicefrac{W}{m^2\cdot sr}$) 
 *         to arbitrary but convenient units in the $[0,1]$ range, set 
 *         this parameter to \code{1e-5}. \default{1}
 *     }
 * }
 *
 * \renderings{
 *     \tinyrendering{6AM}{preetham_06}
 *     \tinyrendering{8AM}{preetham_08}
 *     \tinyrendering{10AM}{preetham_10}
 *     \tinyrendering{12PM}{preetham_12}
 *     \tinyrendering{2PM}{preetham_14}
 *     \tinyrendering{4PM}{preetham_16}
 *     \tinyrendering{6PM}{preetham_18}
 *     \tinyrendering{8PM}{preetham_20}\hfill
 *     \vspace{-3mm}
 *     \caption{Time series with the default settings (shown by
 *     projecting the sky onto a disk. East is left.)}
 * }
 *
 * This plugin implements the physically-based skylight model proposed by 
 * Preetham et al. \cite{Preetham1999Practical}. It can be used for realistic 
 * daylight renderings of scenes under clear and overcast skies, assuming
 * that the sky is observed from a position either on or close to the surface 
 * of the earth. 
 *
 * Numerous parameters allow changing the both the position on Earth, as
 * well as the time of observation. These are used to compute the sun 
 * direction which, together with \code{turbidity}, constitutes the main 
 * parameter of the model. If desired, the sun direction can also be 
 * specified manually.
 *
 * \renderings{
 *     \tinyrendering{2}{preetham_turb_2}
 *     \tinyrendering{3}{preetham_turb_3}
 *     \tinyrendering{4}{preetham_turb_4}
 *     \tinyrendering{5}{preetham_turb_5}
 *     \tinyrendering{6}{preetham_turb_6}
 *     \tinyrendering{7}{preetham_turb_7}
 *     \tinyrendering{8}{preetham_turb_8}
 *     \tinyrendering{9}{preetham_turb_9}
 *     \vspace{-3mm}
 *     \caption{Sky light for different turbidity values (fixed time \& location)}
 *     \vspace{-6mm}
 * }
 *
 * \emph{Turbidity}, the other important parameter, specifies the amount of 
 * atmospheric extinction due to larger particles ($t_l$), as opposed to 
 * molecules ($t_m$). Lower values correspond to a clear sky, and higher values
 * produce illumination resembling that of a hazy, overcast sky. Formally, 
 * the turbidity is defined as the ratio between the combined extinction
 * cross-section and the cross-section only due to molecules, i.e.
 * $T=\frac{t_m+t_l}{t_m}$. Values between 1 and 30 are possible, though 
 * the model will be most accurate for values between 2 and 6, to which 
 * it was fit using numerical optimization.
  
 * The default coordinate system of the luminaire associates the up
 * direction with the $+Y$ axis. The east direction is associated with $+X$
 * and the north direction is equal to $+Z$. To change this coordinate
 * system, rotations can be applied using the \code{toWorld} parameter.
 *
 * By default, the luminaire will not emit any light below the
 * horizon, which means that these regions will be black when they
 * are observed directly. By setting the \code{extend} parameter to 
 * \code{true}, the emitted radiance at the horizon will be extended to 
 * the entire bottom hemisphere. Note that this will significantly 
 * increase the amount of illumination present in the scene.
 *
 * For performance reasons, the implementation precomputes an environment 
 * map of the entire sky that is then forwarded to the \pluginref{envmap} 
 * plugin. The resolution of this environment map can affect the quality
 * of the result. Due to the smoothness of the sky illumination,
 * \code{resolution} values of around 256 (the default) are usually
 * more than sufficient.
 *
 * Note that while the model encompasses sunrise and sunset configurations,
 * it does not extend to the night sky, where illumination from stars, galaxies,
 * and the moon dominate. The model also currently does not handle cloudy skies.
 * The implementation in Mitsuba is based on code by Preetham et al. It was
 * ported by Tom Kazimiers.
 *
 * \begin{xml}[caption={Rotating the sky luminaire for scenes that use $Z$ as
 * the ``up'' direction}, label=lst:sky-up]
 * <luminaire type="sky">
 *     <transform name="toWorld">
 *         <rotate x="1" angle="90"/>
 *     </transform>
 * </luminaire>
 * \end{xml}
 *
 * \renderings{
 *    \medrendering{3 PM}{luminaire_sky_3pm}
 *    \medrendering{6 PM}{luminaire_sky_6pm}
 *    \medrendering{8 PM}{luminaire_sky_8pm}
 *    \vspace{-3mm}
 *    \caption{Renderings with the \pluginref{plastic} material
 *    under default conditions}
 *    \vspace{-6mm}
 * }
 */
class SkyLuminaire : public Luminaire {
public:
	SkyLuminaire(const Properties &props)
			: Luminaire(props) {
		m_scale = props.getFloat("skyScale", 1.0f);
		m_turbidity = props.getFloat("turbidity", 3.0f);
		if (m_turbidity < 1 || m_turbidity > 30)
			Log(EError, "The turbidity parameter must be in the range [1,30]!");

		m_extend = props.getBoolean("extend", false);

		m_resolution = props.getInteger("resolution", 256);

		Point2 sunPos = configureSunPosition(props);
		m_thetaS = sunPos.x;
		m_phiS = sunPos.y;

		configure();
	}

	SkyLuminaire(Stream *stream, InstanceManager *manager) 
		    : Luminaire(stream, manager) {
		m_scale = stream->readFloat();
		m_turbidity = stream->readFloat();
		m_thetaS = stream->readFloat();
		m_phiS = stream->readFloat();
		m_extend = stream->readBool();
		m_resolution = stream->readInt();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);
		stream->writeFloat(m_scale);
		stream->writeFloat(m_turbidity);
		stream->writeFloat(m_thetaS);
		stream->writeFloat(m_phiS);
		stream->writeBool(m_extend);
		stream->writeInt(m_resolution);
	}

	/**
	 * Precalculates some inernal varibles. It needs a valid sun position,
	 * hence the configureSunPos... method must have been called previously.
	 */
	void configure() {
		Float theta2 = m_thetaS * m_thetaS;
		Float theta3 = theta2 * m_thetaS;
		Float turb2 = m_turbidity * m_turbidity;

		/* calculate zenith chromaticity */
		m_zenithX = 
			(+0.00165f*theta3 - 0.00374f*theta2 + 0.00208f*m_thetaS + 0.0f)     * turb2 +
			(-0.02902f*theta3 + 0.06377f*theta2 - 0.03202f*m_thetaS + 0.00394f) * m_turbidity  +
			(+0.11693f*theta3 - 0.21196f*theta2 + 0.06052f*m_thetaS + 0.25885f);

		m_zenithY =
			(+0.00275f*theta3 - 0.00610f*theta2 + 0.00316f*m_thetaS  + 0.0f)     * turb2 +
			(-0.04214f*theta3 + 0.08970f*theta2 - 0.04153f*m_thetaS  + 0.00515f) * m_turbidity  +
			(+0.15346f*theta3 - 0.26756f*theta2 + 0.06669f*m_thetaS  + 0.26688f);

		/* calculate zenith luminance */
		Float chi = (4.0f/9.0f - m_turbidity / 120.0f) * (M_PI - 2 * m_thetaS);

		m_zenithL = (4.0453f * m_turbidity - 4.9710f) * std::tan(chi)
			- 0.2155f * m_turbidity + 2.4192f;

		/* Evaluate quadratic polynomials to find the Perez sky
		 * model coefficients for the x, y and luminance components */
		m_perezL[0] =  0.17872f * m_turbidity  - 1.46303f;
		m_perezL[1] = -0.35540f * m_turbidity  + 0.42749f;
		m_perezL[2] = -0.02266f * m_turbidity  + 5.32505f;
		m_perezL[3] =  0.12064f * m_turbidity  - 2.57705f;
		m_perezL[4] = -0.06696f * m_turbidity  + 0.37027f;

		m_perezX[0] = -0.01925f * m_turbidity - 0.25922f;
		m_perezX[1] = -0.06651f * m_turbidity + 0.00081f;
		m_perezX[2] = -0.00041f * m_turbidity + 0.21247f;
		m_perezX[3] = -0.06409f * m_turbidity - 0.89887f;
		m_perezX[4] = -0.00325f * m_turbidity + 0.04517f;

		m_perezY[0] = -0.01669f * m_turbidity - 0.26078f;
		m_perezY[1] = -0.09495f * m_turbidity + 0.00921f;
		m_perezY[2] = -0.00792f * m_turbidity + 0.21023f;
		m_perezY[3] = -0.04405f * m_turbidity - 1.65369f;
		m_perezY[4] = -0.01092f * m_turbidity + 0.05291f;
	}
	
	bool isCompound() const {
		return true;
	}

	Luminaire *getElement(int i) {
		if (i != 0)
			return NULL;
		int thetaBins = m_resolution, phiBins = m_resolution*2;

		ref<Bitmap> bitmap = new Bitmap(phiBins, thetaBins, 128);
		Point2 factor(M_PI / thetaBins, (2*M_PI) / phiBins);
		float *target = bitmap->getFloatData();
		for (int i=0; i<thetaBins; ++i) {
			Float theta = (i+.5f)*factor.x;
			for (int j=0; j<phiBins; ++j) {
				Float phi = (j+.5f)*factor.y;
				Spectrum s = getSkySpectralRadiance(theta, phi) * m_scale;
				Float r, g, b;
				s.toLinearRGB(r, g, b);
				*target++ = r; *target++ = g;
				*target++ = b; *target++ = 1;
			}
		}

		/* Instantiate a nested envmap plugin */
		Properties props("envmap");
		Properties::Data bitmapData;
		bitmapData.ptr = (uint8_t *) bitmap.get();
		bitmapData.size = sizeof(Bitmap);
		props.setData("bitmap", bitmapData);
		props.setTransform("toWorld", m_luminaireToWorld);
		props.setFloat("samplingWeight", m_samplingWeight);
		Luminaire *luminaire = static_cast<Luminaire *>(
			PluginManager::getInstance()->createObject(
			MTS_CLASS(Luminaire), props));
		luminaire->configure();
		return luminaire;
	}

	Spectrum Le(const Ray &ray) const {
		Point2 coords = fromSphere(ray.d);
		return getSkySpectralRadiance(coords.x, coords.y) * m_scale;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SkyLuminaire[" << endl
			<< "  turbidity = " << m_turbidity << "," << endl 
			<< "  sunPos = [theta: " << m_thetaS << ", phi: "<< m_phiS << "]," << endl 
			<< "  zenithL = " << m_zenithL << "," << endl
			<< "  skyScale = " << m_scale << endl
			<< "]";
		return oss.str();
	}
private:
	/**
	 * Calculates the angle between two spherical cooridnates. All
	 * angles in radians, theta angles measured from up/zenith direction.
	 */
	inline Float getAngleBetween(const Float thetav, const Float phiv,
			const Float theta, const Float phi) const {
		const Float cospsi = std::sin(thetav) * std::sin(theta) * std::cos(phi - phiv)
			+ std::cos(thetav) * std::cos(theta);

		if (cospsi > 1.0f)
			return 0.0f;
		if (cospsi < -1.0f)
			return M_PI;
		return std::acos(cospsi);
	}

	/**
	 * Calculates the distribution of the sky radiance with two
	 * Perez functions:
	 *
	 *        Perez(Theta, Gamma)
	 *   d = ---------------------
	 *        Perez(0, ThetaSun)
	 *
	 * From IES Lighting Handbook pg 361
	 */
	inline Float getDistribution(const Float *lam, const Float theta,
			const Float gamma) const {
		const Float cosGamma = std::cos(gamma);
		const Float num = ((1 + lam[0] * std::fastexp(lam[1] / std::cos(theta)))
			* (1 + lam[2] * std::fastexp(lam[3] * gamma)
			+ lam[4] * cosGamma * cosGamma));

		const Float cosTheta = std::cos(m_thetaS);
		const Float den = ( (1 + lam[0] * std::fastexp(lam[1] /* / cos 0 */))
			* (1 + lam[2] * std::fastexp(lam[3] * m_thetaS)
			+ lam[4] * cosTheta * cosTheta));

		return num / den;
	}

	/**
	 * Calculates the spectral radiance of the sky in the specified directiono.
	 */
	Spectrum getSkySpectralRadiance(Float theta, Float phi) const {
		if (!m_extend && std::cos(theta) <= 0)
			return Spectrum(0.0f);
		/* Clip directions that are extremely close to grazing (for numerical
		 * stability) or entirely below the horizon. This effectively extends 
		 * the horizon luminance value to the bottom hemisphere */
		theta = std::min(theta, (M_PI * 0.5f) - 0.001f);
		/* get angle between sun (zenith is 0, 0) and point (theta, phi) */
		const Float gamma = getAngleBetween(theta, phi, m_thetaS, m_phiS);
		/* Compute xyY values by calculating the distribution for the point
		 * point of interest and multiplying it with the the components
		 * zenith value. */
		const Float x = m_zenithX * getDistribution(m_perezX, theta, gamma);
		const Float y = m_zenithY * getDistribution(m_perezY, theta, gamma);
		const Float Y = m_zenithL * getDistribution(m_perezL, theta, gamma);

		/* Convert xyY to XYZ */
		const Float yFrac = Y / y;
		const Float X = yFrac * x;
		/* It seems the following is necassary to stay always above zero */
		const Float z = std::max((Float) 0.0f, 1.0f - x - y);
		const Float Z = yFrac * z;

		/* Create spectrum from XYZ values */
		Spectrum dstSpect;
		dstSpect.fromXYZ(X, Y, Z, Spectrum::EIlluminant);
		/* The produced spectrum might contain out-of-gamut colors.
		 * The common solution is to clamp resulting values to zero. */
		dstSpect.clampNegative();
		return dstSpect;
	}

	MTS_DECLARE_CLASS()
protected:
	/* Environment map resolution */
	int m_resolution;
	/* Constant scale factor applied to the model */
	Float m_scale;
	/* The turbidity of the sky ranges normally from 1 to 30.
	   For clear skies values in range [2,6] are useful. */
	Float m_turbidity;
	/* Position of the sun in spherical coordinates */
	Float m_thetaS, m_phiS;
	/* Radiance at the zenith, in xyY */
	Float m_zenithL, m_zenithX, m_zenithY;
	/* The distribution coefficints are called A, B, C, D and E by
	   Preetham. They exist for x, y and Y (here called L). The
	   following attributes save the precalculated version of each */
	Float m_perezL[5], m_perezX[5], m_perezY[5];
	/* Extend to the bottom hemisphere? */
	bool m_extend;
};

MTS_IMPLEMENT_CLASS_S(SkyLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(SkyLuminaire, "Preetham sky luminaire");
MTS_NAMESPACE_END

