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
#include <mitsuba/core/util.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

#define SAMPLE_UNIFORMLY 1

/* Define this to automatically set the viewing angles phi to zero if
 * its theta is zero. This means that there is only one possible
 * viewing zenith angle. Some implementitions do this, but it is
 * unclear why one should use it. */
// #define ENFORCE_SINGLE_ZENITH_ANGLE

MTS_NAMESPACE_BEGIN

/*
 * A sun and skylight luminaire. In its local coordinate system, the sun
 * is "above" on positive Y direction. So when configuring, keep in mind
 * that south = x, east = y and up = z. All times in decimal form (6.25
 * = 6:15 AM) and all angles in radians.
 *
 * The model behind it is described by Preetham et al. (2002).
 */
class SkyLuminaire : public Luminaire {
public:
	/**
	 * Creates a new Sky luminaire. The defaults values origate
	 * from the sample code of Preetham et al. and create an over
	 * head sun on a clear day.
	 */
	SkyLuminaire(const Properties &props)
			: Luminaire(props) {
		/* Transformation from the luminaire's local coordinates to
		 *  world coordiantes */
		m_luminaireToWorld = Transform::rotate(Vector(1, 0, 0),-90)
			* props.getTransform("toWorld", Transform());
		m_worldToLuminaire = m_luminaireToWorld.inverse();

		m_skyScale = props.getFloat("skyScale", Float(1.0));
		m_turbidity = props.getFloat("turbidity", Float(2.0));
		m_aConst = props.getFloat("aConst", 1.0);
		m_bConst = props.getFloat("bConst", 1.0);
		m_cConst = props.getFloat("cConst", 1.0);
		m_dConst = props.getFloat("dConst", 1.0);
		m_eConst = props.getFloat("eConst", 1.0);
		m_clipBelowHorizon = props.getBoolean("clipBelowHorizon", true);
		m_exposure = props.getFloat("exposure", 1.0/15.0);

		/* Do some input checks for sun position information */
		bool hasSunDir = props.hasProperty("sunDirection");
		bool hasLatitude = props.hasProperty("latitude");
		bool hasLongitude = props.hasProperty("longitude");
		bool hasStdMrd = props.hasProperty("standardMeridian");
		bool hasJulDay = props.hasProperty("julianDay");
		bool hasTimeOfDay = props.hasProperty("timeOfDay");

		bool hasSomeLocInfo = hasLatitude || hasLongitude || hasStdMrd || hasJulDay || hasTimeOfDay;
		bool hasAllLocInfo = hasLatitude && hasLongitude && hasStdMrd && hasJulDay && hasTimeOfDay;

		if (!(hasSunDir || hasSomeLocInfo)) {
			/* If no sun position has been specified, use a default one. */
			hasSomeLocInfo = hasAllLocInfo = true;
		} else if (hasSunDir && hasSomeLocInfo) {
			/* We don't allow the use of both position formats, raise error. */
			throw std::runtime_error("Please decide for either positioning the sun by a direction vector or by some location and time information.");
		} else if (hasSomeLocInfo && !hasAllLocInfo) {
			/* The sun positioning by location and time information is missing
			* some parameters in the input, raise error. */
			throw std::runtime_error("Please give all required parameters for specifing the sun's position by location and time information. At least one is missing.");
		} /* else we have complete positioning information */

		/* The direction should always relate to "up" beeing in positive Z */
		Vector sunDir = props.getVector("sunDirection", Vector(0.0, 0.0, 1.0));
		Float lat = props.getFloat("latitude", 51.050891);
		Float lon  = props.getFloat("longitude", 13.694458);
		Float stdMrd = props.getFloat("standardMeridian", 0);
		Float julianDay = props.getFloat("julianDay", 200);
		Float timeOfDay = props.getFloat("timeOfDay", 15.00);
		
		/* configure position of sun */
		if (hasSunDir)
			configureSunPosition(sunDir);
		else
			configureSunPosition(lat, lon, stdMrd, julianDay, timeOfDay);

		configure();
	}

	SkyLuminaire(Stream *stream, InstanceManager *manager) 
		    : Luminaire(stream, manager) {
		m_skyScale = stream->readFloat();
		m_turbidity = stream->readFloat();
		m_thetaS = stream->readFloat();
		m_phiS = stream->readFloat();
		m_aConst = stream->readFloat();
		m_bConst = stream->readFloat();
		m_cConst = stream->readFloat();
		m_dConst = stream->readFloat();
		m_eConst = stream->readFloat();
		m_clipBelowHorizon = stream->readBool();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);
		stream->writeFloat(m_skyScale);
		stream->writeFloat(m_turbidity);
		stream->writeFloat(m_thetaS);
		stream->writeFloat(m_phiS);
		stream->writeFloat(m_aConst);
		stream->writeFloat(m_bConst);
		stream->writeFloat(m_cConst);
		stream->writeFloat(m_dConst);
		stream->writeFloat(m_eConst);
		stream->writeBool(m_clipBelowHorizon);
	}

	/**
	 * Precalculates some inernal varibles. It needs a valid sun position.
	 * Hence it needs a configureSunPos... method called before itself.
	 */
	void configure() {
		Float theta2 = m_thetaS * m_thetaS;
		Float theta3 = theta2 * m_thetaS;

		Float turb2 = m_turbidity * m_turbidity;

		/* calculate zenith chromaticity */
		m_zenithX =
		(+0.00165*theta3 - 0.00374*theta2 + 0.00208*m_thetaS + 0)       * turb2 +
		(-0.02902*theta3 + 0.06377*theta2 - 0.03202*m_thetaS + 0.00394) * m_turbidity  +
		(+0.11693*theta3 - 0.21196*theta2 + 0.06052*m_thetaS + 0.25885);

		m_zenithY =
		(+0.00275*theta3 - 0.00610*theta2 + 0.00316*m_thetaS  + 0)       * turb2 +
		(-0.04214*theta3 + 0.08970*theta2 - 0.04153*m_thetaS  + 0.00515) * m_turbidity  +
		(+0.15346*theta3 - 0.26756*theta2 + 0.06669*m_thetaS  + 0.26688);

		/* calculate zenith luminance */
		Float chi = (4.0/9.0 - m_turbidity / 120.0) * (M_PI - 2 * m_thetaS);

		m_zenithL = (4.0453 * m_turbidity - 4.9710) * tan(chi)
			- 0.2155 * m_turbidity + 2.4192;

		m_perezL[0] =  ( 0.17872 * m_turbidity  - 1.46303) * m_aConst;
		m_perezL[1] =  (-0.35540 * m_turbidity  + 0.42749) * m_bConst;
		m_perezL[2] =  (-0.02266 * m_turbidity  + 5.32505) * m_cConst;
		m_perezL[3] =  ( 0.12064 * m_turbidity  - 2.57705) * m_dConst;
		m_perezL[4] =  (-0.06696 * m_turbidity  + 0.37027) * m_eConst;

		m_perezX[0] =  (-0.01925 * m_turbidity - 0.25922) * m_aConst;
		m_perezX[1] =  (-0.06651 * m_turbidity + 0.00081) * m_bConst;
		m_perezX[2] =  (-0.00041 * m_turbidity + 0.21247) * m_cConst;
		m_perezX[3] =  (-0.06409 * m_turbidity - 0.89887) * m_dConst;
		m_perezX[4] =  (-0.00325 * m_turbidity + 0.04517) * m_eConst;

		m_perezY[0] =  (-0.01669 * m_turbidity - 0.26078) * m_aConst;
		m_perezY[1] =  (-0.09495 * m_turbidity + 0.00921) * m_bConst;
		m_perezY[2] =  (-0.00792 * m_turbidity + 0.21023) * m_cConst;
		m_perezY[3] =  (-0.04405 * m_turbidity - 1.65369) * m_dConst;
		m_perezY[4] =  (-0.01092 * m_turbidity + 0.05291) * m_eConst;



		int thetaBins = 512, phiBins = thetaBins*2;
		ref<Bitmap> bitmap = new Bitmap(phiBins, thetaBins, 128);
		Point2 factor(M_PI / thetaBins, (2*M_PI) / phiBins);
		for (int i=0; i<thetaBins; ++i) {
			Float theta = (i+.5f)*factor.x;
			for (int j=0; j<phiBins; ++j) {
				Float phi = (j+.5f)*factor.x;
				Spectrum s = Le(sphericalDirection(theta, phi));
				Float r, g, b;
				s.toLinearRGB(r, g, b);
				bitmap->getFloatData()[(j+i*phiBins)*4 + 0] = r;
				bitmap->getFloatData()[(j+i*phiBins)*4 + 1] = g;
				bitmap->getFloatData()[(j+i*phiBins)*4 + 2] = b;
				bitmap->getFloatData()[(j+i*phiBins)*4 + 3] = 1;
			}
		}
		ref<FileStream> fs = new FileStream("out.exr", FileStream::ETruncReadWrite);
		bitmap->save(Bitmap::EEXR, fs);
	}

	/**
	 * Configures the position of the sun. This calculation is based on
	 * your position on the world and time of day.
	 * From IES Lighting Handbook pg 361.
	 */
	void configureSunPosition(const Float lat, const Float lon,
			const int stdMrd, const int julDay, const Float timeOfDay) {
		const Float solarTime = timeOfDay
			+ (0.170 * sin(4.0 * M_PI * (julDay - 80.0) / 373.0)
			- 0.129 * sin(2.0 * M_PI * (julDay - 8.0) / 355.0))
			+ (stdMrd - lon) / 15.0;

		const Float solarDeclination = (0.4093 * sin(2 * M_PI * (julDay - 81) / 368));

		const Float solarAltitude = asin(sin(degToRad(lat))
			* sin(solarDeclination) - cos(degToRad(lat))
			* cos(solarDeclination) * cos(M_PI * solarTime / 12.0));

		const Float opp = -cos(solarDeclination) * sin(M_PI * solarTime / 12.0);
		const Float adj = -(cos(degToRad(lat)) * sin(solarDeclination)
			+ sin(degToRad(lat)) * cos(solarDeclination)
			* cos(M_PI * solarTime / 12.0));

		const Float solarAzimuth = atan2(opp,adj);

		m_phiS = -solarAzimuth;
		m_thetaS = M_PI / 2.0 - solarAltitude;
	}

	/**
	 * Configures the position of the sun by using a vector that points
	 * to the sun. It is expected, that +x = south, +y = east, +z = up.
	 */
	void configureSunPosition(const Vector& sunDir) {
		Vector wh = normalize(sunDir);
		Point2 sunPos = toSphericalCoordinates(wh);
		m_thetaS = sunPos.x;
		m_phiS = sunPos.y;
	}

	void preprocess(const Scene *scene) {
		/* Get the scene's bounding sphere and slightly enlarge it */
		m_bsphere = scene->getBSphere();
		m_bsphere.radius *= 1.01f;
		m_surfaceArea = m_bsphere.radius * m_bsphere.radius * M_PI;
		m_invSurfaceArea = 1/m_surfaceArea;
	}

	Spectrum getPower() const {
		/* TODO */
		return m_average * (M_PI * 4 * M_PI
		* m_bsphere.radius * m_bsphere.radius);
	}

	inline Spectrum Le(const Vector &direction) const {
		/* Compute sky light radiance for direction */
		Vector d = normalize(m_worldToLuminaire(direction));

		if (m_clipBelowHorizon) {
			/* if sun is below horizon, return black */
			if (d.z < 0.0f)
				return Spectrum(0.0f);
		}

		/* make all "zero values" the same */
		if (d.z < 0.001f)
			d = normalize(Vector(d.x, d.y,  0.001f));

		const Point2 dSpherical = toSphericalCoordinates(d);
		const Float theta = dSpherical.x;

#ifndef ENFORCE_SINGLE_ZENITH_ANGLE
		const Float phi = dSpherical.y;
#else
		Float phi;
		if (fabs(theta) < 1e-5)
			phi = 0.0f;
		else
			phi = dSpherical.y;
#endif

		Spectrum L;
		getSkySpectralRadiance(theta, phi, L);
		L *= m_skyScale;
		
		return L;
	}

	inline Spectrum Le(const Ray &ray) const {
		return Le(normalize(ray.d));
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		return Le(-lRec.d);
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
			const Point2 &sample) const {
		lRec.d = sampleDirection(sample, lRec.pdf, lRec.value);
		lRec.sRec.p = p - lRec.d * (2 * m_bsphere.radius);
	}

	void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		SkyLuminaire::sample(its.p, lRec, sample);
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec, bool delta) const {
#if defined(SAMPLE_UNIFORMLY)
		return 1.0f / (4 * M_PI);
#endif
	}

	Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec, bool delta) const {
		return SkyLuminaire::pdf(its.p, lRec, delta);
	}

	/**
	 * This is the tricky bit - we want to sample a ray that
	 * has uniform density over the set of all rays passing
	 * through the scene.
	 * For more detail, see "Using low-discrepancy sequences and 
	 * the Crofton formula to compute surface areas of geometric models"
	 * by Li, X. and Wang, W. and Martin, R.R. and Bowyer, A. 
	 * (Computer-Aided Design vol 35, #9, pp. 771--782)
	 */
	void sampleEmission(EmissionRecord &eRec, 
		const Point2 &sample1, const Point2 &sample2) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		/* Chord model - generate the ray passing through two uniformly
		   distributed points on a sphere containing the scene */
		Vector d = squareToSphere(sample1);
		eRec.sRec.p = m_bsphere.center + d * m_bsphere.radius;
		eRec.sRec.n = Normal(-d);
		Point p2 = m_bsphere.center + squareToSphere(sample2) * m_bsphere.radius;
		eRec.d = p2 - eRec.sRec.p;
		Float length = eRec.d.length();

		if (length == 0) {
			eRec.value = Spectrum(0.0f);
			eRec.pdfArea = eRec.pdfDir = 1.0f;
			return;
		}

		eRec.d /= length;
		eRec.pdfArea = 1.0f / (4 * M_PI * m_bsphere.radius * m_bsphere.radius);
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		eRec.value = Le(-eRec.d);
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		if (eRec.type == EmissionRecord::ENormal) {
			Vector d = squareToSphere(sample);
			eRec.sRec.p = m_bsphere.center + d * m_bsphere.radius;
			eRec.sRec.n = Normal(-d);
			eRec.pdfArea = 1.0f / (4 * M_PI * m_bsphere.radius * m_bsphere.radius);
			eRec.value = Spectrum(M_PI);
		} else {
			/* Preview mode, which is more suitable for VPL-based rendering: approximate 
			   the infinitely far-away source with set of diffuse point sources */
			const Float radius = m_bsphere.radius * 1.5f;
			Vector d = squareToSphere(sample);
			eRec.sRec.p = m_bsphere.center + d * radius;
			eRec.sRec.n = Normal(-d);
			eRec.pdfArea = 1.0f / (4 * M_PI * radius * radius);
			eRec.value = Le(d) * M_PI;
		}
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		Float radius = m_bsphere.radius;
		if (eRec.type == EmissionRecord::EPreview) 
			radius *= 1.5f;
		Point p2 = m_bsphere.center + squareToSphere(sample) * radius;
		eRec.d = p2 - eRec.sRec.p;
		Float length = eRec.d.length();

		if (length == 0.0f) {
			eRec.pdfDir = 1.0f;
			return Spectrum(0.0f);
		}

		eRec.d /= length;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		if (eRec.type == EmissionRecord::ENormal)
			return Le(-eRec.d) * INV_PI;
		else
			return Spectrum(INV_PI);
	}
    
	Spectrum fDirection(const EmissionRecord &eRec) const {
		if (eRec.type == EmissionRecord::ENormal)
			return Le(-eRec.d) * INV_PI;
		else
			return Spectrum(INV_PI);
	}

	Spectrum fArea(const EmissionRecord &eRec) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		return Spectrum(M_PI);
	}

	Spectrum f(const EmissionRecord &eRec) const {
		if (eRec.type == EmissionRecord::ENormal)
			return Le(-eRec.d) * INV_PI;
		else
			return Spectrum(INV_PI);
	}

	void pdfEmission(EmissionRecord &eRec, bool delta) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = delta ? 0.0f : INV_PI * dp;
		else
			eRec.pdfDir = 0;
		eRec.pdfArea = delta ? 0.0f : m_invSurfaceArea;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SkyLuminaire[" << std::endl
			<< "  sky sscale = " << m_skyScale << "," << std::endl
			<< "  power = " << getPower().toString() << "," << std::endl
			<< "  sun pos = theta: " << m_thetaS << ", phi: "<< m_phiS << ","  << std::endl 
			<< "  turbidity = " << m_turbidity << "," << std::endl 
			<< "]";
		return oss.str();
	}

	bool isBackgroundLuminaire() const {
		return true;
	}

	Vector sampleDirection(Point2 sample, Float &pdf, Spectrum &value) const {
#if defined(SAMPLE_UNIFORMLY)
		pdf = 1.0f / (4*M_PI);
		Vector d = squareToSphere(sample);
		value = Le(-d);
		return d;
#endif
	}

private:
	/**
	 * Calculates the anlgle between two spherical cooridnates. All
	 * angles in radians, theta angles measured from up/zenith direction.
	 */
	inline Float getAngleBetween(const Float thetav, const Float phiv,
			const Float theta, const Float phi) const {
		const Float cospsi = sin(thetav) * sin(theta) * cos(phi - phiv)
			+ cos(thetav) * cos(theta);

		if (cospsi > 1.0f)
			return 0.0f;
		if (cospsi < -1.0f)
			return M_PI;
		return acos(cospsi);
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
		const Float cosGamma = cos(gamma);
		const Float num = ( (1 + lam[0] * exp(lam[1] / cos(theta)))
			* (1 + lam[2] * exp(lam[3] * gamma)
			+ lam[4] * cosGamma * cosGamma));

		const Float cosTheta = cos(m_thetaS);
		const Float den = ( (1 + lam[0] * exp(lam[1] /* / cos 0 */))
			* (1 + lam[2] * exp(lam[3] * m_thetaS)
			+ lam[4] * cosTheta * cosTheta));

		return (num / den);
	}

	/**
	 * Calculates the spectral radiance of the sky in the specified directiono.
	 */
	void getSkySpectralRadiance(const Float theta, const Float phi, Spectrum &dstSpect) const {
		/* add bottom half of hemisphere with horizon colour */
		const Float theta_fin = std::min(theta, (M_PI * 0.5f) - 0.001f);
		/* get angle between sun (zenith is 0, 0) and point (theta, phi) */
		const Float gamma = getAngleBetween(theta, phi, m_thetaS, m_phiS);
		/* Compute xyY values by calculating the distribution for the point
		 * point of interest and multiplying it with the the components
		 * zenith value. */
		const Float x = m_zenithX * getDistribution(m_perezX, theta_fin, gamma);
		const Float y = m_zenithY * getDistribution(m_perezY, theta_fin, gamma);
		const Float Y = m_zenithL * getDistribution(m_perezL, theta_fin, gamma);

		/* Apply an exponential exposure function */
		// Y = 1.0 - exp(-m_exposure * Y);
		/* Convert xyY to XYZ */
		const Float yFrac = Y / y;
		const Float X = yFrac * x;
		/* It seems the following is necassary to stay always above zero */
		const Float z = std::max(0.0f, 1.0f - x - y);
		const Float Z = yFrac * z;

		/* Create spectrum from XYZ values */
		dstSpect.fromXYZ(X, Y, Z);
		/* The produced spectrum contains out-of-gamut colors.
		 * It is common to clamp resulting values to zero. */
		dstSpect.clampNegative();
	}

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_average;
	Float m_surfaceArea;
	Float m_invSurfaceArea;
	BSphere m_bsphere;
	Float m_exposure;
	Float m_skyScale;
	/* The turbidity of the sky ranges normally from 0 to 30+.
	 * For clear skies values in range [2,6] are useful. */
	Float m_turbidity;
	Float m_thetaS, m_phiS;
	Float m_zenithL;
	Float m_zenithX;
	Float m_zenithY;
	/* The distribution coefficints are called A, B, C, D and E by
	 * Preedham. Since they exist for x, y and Y (here called L)
	 * we save the precalculated version of it
	 *
	 * distribution coefficients for luminance distribution function */
	Float m_perezL[5];
	/* distribution coefficients for x distribution function */
	Float m_perezX[5];
	/* distribution coefficients for y distribution function */
	Float m_perezY[5];
	/* distribution function tuning coefficients a, b, c, d and e.
	 * They can be used to scale the luminance, x and y distribution
	 * coefficints. */
	Float m_aConst, m_bConst, m_cConst, m_dConst, m_eConst;

	/* To disable clipping of incoming sky light below the horizon, set this
	 * to false. If set to true black is returned for queries. Be awere that
	 * a huge amount of additional radiance is coming in (i.e. it gets a
	 * lot brighter). */
	bool m_clipBelowHorizon;
};

MTS_IMPLEMENT_CLASS_S(SkyLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(SkyLuminaire, "Sky luminaire");
MTS_NAMESPACE_END

