/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__SUN_H)
#define __SUN_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

Vector toSphere(Float theta, Float phi) {
	/* Spherical-to-cartesian coordinate mapping with 
	   theta=0 => Y=1 */
	Float cosTheta = std::cos(theta), 
		  sinTheta = std::sqrt(1-cosTheta*cosTheta),
		  cosPhi = std::cos(phi), sinPhi = std::sin(phi);
	return Vector(
		sinTheta * sinPhi, cosTheta, -sinTheta*cosPhi);
}

Point2 fromSphere(const Vector &d) {
	Float theta = std::acos(std::max((Float) -1.0f, 
			std::min((Float) 1.0f, d.y)));
	Float phi = std::atan2(d.x,-d.z);
	if (phi < 0)
		phi += 2*M_PI;
	return Point2(theta, phi);
}

/**
 * Configures the position of the sun. This calculation is based on
 * your position on the world and time of day.
 * From IES Lighting Handbook pg 361.
 */
Point2 configureSunPosition(Float lat, Float lon, int stdMrd, 
		int day, Float time) {
	const Float solarTime = time
		+ (0.170f * std::sin(4.0f * M_PI * (day - 80.0f) / 373.0f)
		- 0.129f * std::sin(2.0f * M_PI * (day - 8.0f) / 355.0f))
		+ (stdMrd - lon) / 15.0f;

	const Float solarDeclination = (0.4093f * std::sin(2 * M_PI
			* (day - 81.0f) / 368.0f));

	lat = degToRad(lat);
	const Float solarAltitude = std::asin(std::sin(lat)
		* std::sin(solarDeclination) - std::cos(lat)
		* std::cos(solarDeclination) * std::cos(M_PI * solarTime / 12.0f));

	const Float opp = -std::cos(solarDeclination) 
		* std::sin(M_PI * solarTime / 12.0f);
	const Float adj = -(std::cos(lat) * std::sin(solarDeclination)
		+ std::sin(lat) * std::cos(solarDeclination)
		* std::cos(M_PI * solarTime / 12.0f));

	const Float solarAzimuth = std::atan2(opp, adj);

	return Point2(
		M_PI / 2.0f - solarAltitude,
		-solarAzimuth);
}

Point2 configureSunPosition(const Vector& sunDir, const Transform &worldToLuminaire) {
	return fromSphere(normalize(worldToLuminaire(sunDir)));
}

Point2 configureSunPosition(const Properties &props) {
	/* configure position of sun */
	if (props.hasProperty("sunDirection")) {
		if (props.hasProperty("latitude") || props.hasProperty("longitude")
			|| props.hasProperty("standardMeridian") || props.hasProperty("day")
			|| props.hasProperty("time"))
			SLog(EError, "Both the 'sunDirection' parameter and time/location "
					"information were provided -- only one of them can be specified at a time!");

		return configureSunPosition(
			props.getVector("sunDirection"),
			props.getTransform("toWorld", Transform()).inverse());
	} else {
		Float lat = props.getFloat("latitude", 35.6894f);
		Float lon  = props.getFloat("longitude", 139.6917f);
		int stdMrd = props.getInteger("standardMeridian", 135);

		int day = props.getInteger("day", 180);
		if (day < 1 || day > 365)
			SLog(EError, "The day parameter must be in the range [1, 365]!");

		Float time = props.getFloat("time", 15.00f);
		if (time < 0 || time > 24)
			SLog(EError, "The time parameter must be in the range [0, 24]!");

		return configureSunPosition(lat, lon, stdMrd, day, time);
	}
}

/* All data lifted from MI */
/* Units are either [] or cm^-1. refer when in doubt MI */

// k_o Spectrum table from pg 127, MI.
Float k_oWavelengths[64] = {
	300, 305, 310, 315, 320, 325, 330, 335, 340, 345,
	350, 355, 445, 450, 455, 460, 465, 470, 475, 480, 
	485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 
	535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 
	585, 590, 595, 600, 605, 610, 620, 630, 640, 650, 
	660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 
	760, 770, 780, 790,
};

Float k_oAmplitudes[65] = {
	10.0, 4.8, 2.7, 1.35, .8, .380, .160, .075, .04, .019, .007, 
	.0, .003, .003, .004, .006, .008, .009, .012, .014, .017, 
	.021, .025, .03, .035, .04, .045, .048, .057, .063, .07, 
	.075, .08, .085, .095, .103, .110, .12, .122, .12, .118, 
	.115, .12, .125, .130, .12, .105, .09, .079, .067, .057, 
	.048, .036, .028, .023, .018, .014, .011, .010, .009, 
	.007, .004, .0, .0
};

// k_g Spectrum table from pg 130, MI.
Float k_gWavelengths[4] = {
	759, 760, 770, 771
};

Float k_gAmplitudes[4] = {
	0, 3.0, 0.210, 0
};

// k_wa Spectrum table from pg 130, MI.
Float k_waWavelengths[13] = {
	689, 690, 700, 710, 720,
	730, 740, 750, 760, 770,
	780, 790, 800
};

Float k_waAmplitudes[13] = {
	0, 0.160e-1, 0.240e-1, 0.125e-1,
	0.100e+1, 0.870, 0.610e-1, 0.100e-2,
	0.100e-4, 0.100e-4, 0.600e-3,
	0.175e-1, 0.360e-1
};

Float solWavelengths[38] = {
	380, 390, 400, 410, 420, 430, 440, 450, 
	460, 470, 480, 490, 500, 510, 520, 530, 
	540, 550, 560, 570, 580, 590, 600, 610, 
	620, 630, 640, 650, 660, 670, 680, 690, 
	700, 710, 720, 730, 740, 750
};

Float solAmplitudes[38] = {
	165.5, 162.3, 211.2, 258.8, 258.2,
	242.3, 267.6, 296.6, 305.4, 300.6,
	306.6, 288.3, 287.1, 278.2, 271.0,
	272.3, 263.6, 255.0, 250.6, 253.1,
	253.5, 251.3, 246.3, 241.7, 236.8,
	232.1, 228.2, 223.4, 219.7, 215.3,
	211.0, 207.3, 202.4, 198.7, 194.3,
	190.7, 186.3, 182.6
};

Spectrum computeSunRadiance(Float theta, Float turbidity) {
    InterpolatedSpectrum k_oCurve(k_oWavelengths, k_oAmplitudes, 64);
   	InterpolatedSpectrum k_gCurve(k_gWavelengths, k_gAmplitudes, 4);
    InterpolatedSpectrum k_waCurve(k_waWavelengths, k_waAmplitudes, 13);
    InterpolatedSpectrum solCurve(solWavelengths, solAmplitudes, 38);
    Float  data[91], wavelengths[91];  // (800 - 350) / 5  + 1

    Float beta = 0.04608365822050f * turbidity - 0.04586025928522f;

    Float m = 1.0f/(std::cos(theta) + 0.15f*std::pow(93.885f-theta/M_PI*180.0f, (Float) -1.253f));  // Relative Optical Mass

    Float lambda;
	int i;
    for(i = 0, lambda = 350; i < 91; i++, lambda+=5) {
		// Rayleigh Scattering
		// Results agree with the graph (pg 115, MI) */
		Float tauR = std::fastexp(-m * 0.008735f * std::pow(lambda/1000.0f, (Float) -4.08));

		// Aerosal (water + dust) attenuation
		// beta - amount of aerosols present 
		// alpha - ratio of small to large particle sizes. (0:4,usually 1.3)
		// Results agree with the graph (pg 121, MI) 
		const Float alpha = 1.3f;
		Float tauA = exp(-m * beta * std::pow(lambda/1000.0f, -alpha));  // lambda should be in um

		// Attenuation due to ozone absorption  
		// lOzone - amount of ozone in cm(NTP) 
		// Results agree with the graph (pg 128, MI) 
		const Float lOzone = .35f;
		Float tauO = std::fastexp(-m * k_oCurve.eval(lambda) * lOzone);

		// Attenuation due to mixed gases absorption  
		// Results agree with the graph (pg 131, MI)
		Float tauG = std::fastexp(-1.41f * k_gCurve.eval(lambda) * m / std::pow(1 + 118.93f
			* k_gCurve.eval(lambda) * m, (Float) 0.45f));

		// Attenuation due to water vapor absorbtion  
		// w - precipitable water vapor in centimeters (standard = 2) 
		// Results agree with the graph (pg 132, MI)
		const Float w = 2.0;
		Float tauWA = std::fastexp(-0.2385f * k_waCurve.eval(lambda) * w * m /
				std::pow(1 + 20.07f * k_waCurve.eval(lambda) * w * m, (Float) 0.45f));

		data[i] = (Float) 100.0f * solCurve.eval(lambda) * tauR * tauA * tauO * tauG * tauWA;
		wavelengths[i] = lambda;
    }
    InterpolatedSpectrum interpolated(wavelengths, data, 91);
	Spectrum discretized;
	discretized.fromContinuousSpectrum(interpolated);
	discretized *= 300;
	return discretized;
}


MTS_NAMESPACE_END

#endif /* __SUN_H */
