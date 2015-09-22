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

#if !defined(__SUN_H)
#define __SUN_H

#include <mitsuba/mitsuba.h>

#define EARTH_MEAN_RADIUS 6371.01	// In km
#define ASTRONOMICAL_UNIT 149597890	// In km

MTS_NAMESPACE_BEGIN

struct DateTimeRecord {
	int year;
	int month;
	int day;
	Float hour;
	Float minute;
	Float second;

	std::string toString() const {
		std::ostringstream oss;
		oss << "DateTimeRecord[year = " << year
			<< ", month= " << month
			<< ", day = " << day
			<< ", hour = " << hour
			<< ", minute = " << minute
			<< ", second = " << second << "]";
		return oss.str();
	}
};

struct LocationRecord {
	Float longitude;
	Float latitude;
	Float timezone;

	std::string toString() const {
		std::ostringstream oss;
		oss << "LocationRecord[latitude = " << latitude
			<< ", longitude = " << longitude
			<< ", timezone = " << timezone << "]";
		return oss.str();
	}
};

struct SphericalCoordinates {
	Float elevation;
	Float azimuth;

	inline SphericalCoordinates() { }

	inline SphericalCoordinates(Float elevation, Float azimuth)
		: elevation(elevation), azimuth(azimuth) { }

	inline SphericalCoordinates(Stream *stream) {
		elevation = stream->readFloat();
		azimuth = stream->readFloat();
	}

	void serialize(Stream *stream) const {
		stream->writeFloat(elevation);
		stream->writeFloat(azimuth);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SphericalCoordinates[elevation = " << radToDeg(elevation)
			<< ", azimuth = " << radToDeg(azimuth) << "]";
		return oss.str();
	}
};

Vector toSphere(const SphericalCoordinates coords) {
	Float sinTheta, cosTheta, sinPhi, cosPhi;

	math::sincos(coords.elevation, &sinTheta, &cosTheta);
	math::sincos(coords.azimuth, &sinPhi, &cosPhi);

	return Vector(sinPhi*sinTheta, cosTheta, -cosPhi*sinTheta);
}

SphericalCoordinates fromSphere(const Vector &d) {
	Float azimuth = std::atan2(d.x, -d.z);
	Float elevation = math::safe_acos(d.y);
	if (azimuth < 0)
		azimuth += 2*M_PI;
	return SphericalCoordinates(elevation, azimuth);
}

/**
 * \brief Compute the elevation and azimuth of the sun as seen by an observer
 * at \c location at the date and time specified in \c dateTime.
 *
 * Based on "Computing the Solar Vector" by Manuel Blanco-Muriel,
 * Diego C. Alarcon-Padilla, Teodoro Lopez-Moratalla, and Martin Lara-Coira,
 * in "Solar energy", vol 27, number 5, 2001 by Pergamon Press.
 */
SphericalCoordinates computeSunCoordinates(const DateTimeRecord &dateTime, const LocationRecord &location) {
	// Main variables
	double elapsedJulianDays, decHours;
	double eclipticLongitude, eclipticObliquity;
	double rightAscension, declination;
	double elevation, azimuth;

	// Auxiliary variables
	double dY;
	double dX;

	/* Calculate difference in days between the current Julian Day
	   and JD 2451545.0, which is noon 1 January 2000 Universal Time */
	{
		// Calculate time of the day in UT decimal hours
		decHours = dateTime.hour - location.timezone +
			(dateTime.minute + dateTime.second / 60.0 ) / 60.0;

		// Calculate current Julian Day
		int liAux1 = (dateTime.month-14) / 12;
		int liAux2 = (1461*(dateTime.year + 4800 + liAux1)) / 4
			+ (367 * (dateTime.month - 2 - 12 * liAux1)) / 12
			- (3 * ((dateTime.year + 4900 + liAux1) / 100)) / 4
			+ dateTime.day - 32075;
		double dJulianDate = (double) liAux2 - 0.5 + decHours / 24.0;

		// Calculate difference between current Julian Day and JD 2451545.0
		elapsedJulianDays = dJulianDate - 2451545.0;
	}

	/* Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
	   ecliptic in radians but without limiting the angle to be less than 2*Pi
	   (i.e., the result may be greater than 2*Pi) */
	{
		double omega = 2.1429 - 0.0010394594 * elapsedJulianDays;
		double meanLongitude = 4.8950630 + 0.017202791698 * elapsedJulianDays; // Radians
		double anomaly = 6.2400600 + 0.0172019699 * elapsedJulianDays;

		eclipticLongitude = meanLongitude + 0.03341607 * std::sin(anomaly)
			+ 0.00034894 * std::sin(2*anomaly) - 0.0001134
			- 0.0000203 * std::sin(omega);

		eclipticObliquity = 0.4090928 - 6.2140e-9 * elapsedJulianDays
			+ 0.0000396 * std::cos(omega);
	}

	/* Calculate celestial coordinates ( right ascension and declination ) in radians
	   but without limiting the angle to be less than 2*Pi (i.e., the result may be
	   greater than 2*Pi) */
	{
		double sinEclipticLongitude = std::sin(eclipticLongitude);
		dY = std::cos(eclipticObliquity) * sinEclipticLongitude;
		dX = std::cos(eclipticLongitude);
		rightAscension = std::atan2(dY, dX);
		if (rightAscension < 0.0)
			rightAscension += 2*M_PI;
		declination = std::asin(std::sin(eclipticObliquity) * sinEclipticLongitude);
	}

	// Calculate local coordinates (azimuth and zenith angle) in degrees
	{
		double greenwichMeanSiderealTime = 6.6974243242
			+ 0.0657098283 * elapsedJulianDays + decHours;

		double localMeanSiderealTime = degToRad((Float) ((greenwichMeanSiderealTime * 15
			+ location.longitude)));

		double latitudeInRadians = degToRad(location.latitude);
		double cosLatitude = std::cos(latitudeInRadians);
		double sinLatitude = std::sin(latitudeInRadians);

		double hourAngle = localMeanSiderealTime - rightAscension;
		double cosHourAngle = std::cos(hourAngle);

		elevation = std::acos(cosLatitude * cosHourAngle
			* std::cos(declination) + std::sin(declination) * sinLatitude);

		dY = -std::sin(hourAngle);
		dX = std::tan(declination) * cosLatitude - sinLatitude * cosHourAngle;

		azimuth = std::atan2(dY, dX);
		if (azimuth < 0.0)
			azimuth += 2*M_PI;

		// Parallax Correction
		elevation += (EARTH_MEAN_RADIUS / ASTRONOMICAL_UNIT) * std::sin(elevation);
	}

	return SphericalCoordinates((Float) elevation, (Float) azimuth);
}

SphericalCoordinates computeSunCoordinates(const Vector& sunDir, const Transform &worldToLuminaire) {
	return fromSphere(normalize(worldToLuminaire(sunDir)));
}

SphericalCoordinates computeSunCoordinates(const Properties &props) {
	/* configure position of sun */
	if (props.hasProperty("sunDirection")) {
		if (props.hasProperty("latitude") || props.hasProperty("longitude")
			|| props.hasProperty("timezone") || props.hasProperty("day")
			|| props.hasProperty("time"))
			SLog(EError, "Both the 'sunDirection' parameter and time/location "
					"information were provided -- only one of them can be specified at a time!");

		return computeSunCoordinates(
			props.getVector("sunDirection"),
			props.getAnimatedTransform("toWorld", Transform())->eval(0).inverse());
	} else {
		LocationRecord location;
		DateTimeRecord dateTime;

		location.latitude  = props.getFloat("latitude", 35.6894f);
		location.longitude = props.getFloat("longitude", 139.6917f);
		location.timezone  = props.getFloat("timezone", 9);
		dateTime.year      = props.getInteger("year", 2010);
		dateTime.day       = props.getInteger("day", 10);
		dateTime.month     = props.getInteger("month", 7);
		dateTime.hour      = props.getFloat("hour", 15.0f);
		dateTime.minute    = props.getFloat("minute", 0.0f);
		dateTime.second    = props.getFloat("second", 0.0f);

		SphericalCoordinates coords = computeSunCoordinates(dateTime, location);

		SLog(EDebug, "Computed sun position for %s and %s: %s",
			location.toString().c_str(), dateTime.toString().c_str(),
			coords.toString().c_str());

		return coords;
	}
}

/* The following is from the implementation of "A Practical Analytic Model for
   Daylight" by A.J. Preetham, Peter Shirley, and Brian Smits */

/* All data lifted from MI. Units are either [] or cm^-1. refer when in doubt MI */

// k_o Spectrum table from pg 127, MI.
Float k_oWavelengths[64] = {
	300, 305, 310, 315, 320, 325, 330, 335, 340, 345,
	350, 355, 445, 450, 455, 460, 465, 470, 475, 480,
	485, 490, 495, 500, 505, 510, 515, 520, 525, 530,
	535, 540, 545, 550, 555, 560, 565, 570, 575, 580,
	585, 590, 595, 600, 605, 610, 620, 630, 640, 650,
	660, 670, 680, 690, 700, 710, 720, 730, 740, 750,
	760, 770, 780, 790
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

/* Wavelengths corresponding to the table below */
Float solWavelengths[38] = {
	380, 390, 400, 410, 420, 430, 440, 450,
	460, 470, 480, 490, 500, 510, 520, 530,
	540, 550, 560, 570, 580, 590, 600, 610,
	620, 630, 640, 650, 660, 670, 680, 690,
	700, 710, 720, 730, 740, 750
};

/* Solar amplitude in watts / (m^2 * nm * sr) */
Float solAmplitudes[38] = {
	16559.0, 16233.7, 21127.5, 25888.2, 25829.1,
	24232.3, 26760.5, 29658.3, 30545.4, 30057.5,
	30663.7, 28830.4, 28712.1, 27825.0, 27100.6,
	27233.6, 26361.3, 25503.8, 25060.2, 25311.6,
	25355.9, 25134.2, 24631.5, 24173.2, 23685.3,
	23212.1, 22827.7, 22339.8, 21970.2, 21526.7,
	21097.9, 20728.3, 20240.4, 19870.8, 19427.2,
	19072.4, 18628.9, 18259.2
};

Spectrum computeSunRadiance(Float theta, Float turbidity) {
    InterpolatedSpectrum k_oCurve(k_oWavelengths, k_oAmplitudes, 64);
   	InterpolatedSpectrum k_gCurve(k_gWavelengths, k_gAmplitudes, 4);
    InterpolatedSpectrum k_waCurve(k_waWavelengths, k_waAmplitudes, 13);
    InterpolatedSpectrum solCurve(solWavelengths, solAmplitudes, 38);
    Float data[91], wavelengths[91];  // (800 - 350) / 5  + 1

    Float beta = 0.04608365822050f * turbidity - 0.04586025928522f;

	// Relative Optical Mass
    Float m = 1.0f / (std::cos(theta) + 0.15f *
		std::pow(93.885f - theta/M_PI*180.0f, (Float) -1.253f));

    Float lambda;
	int i = 0;
    for(i = 0, lambda = 350; i < 91; i++, lambda += 5) {
		// Rayleigh Scattering
		// Results agree with the graph (pg 115, MI) */
		Float tauR = math::fastexp(-m * 0.008735f * std::pow(lambda/1000.0f, (Float) -4.08));

		// Aerosol (water + dust) attenuation
		// beta - amount of aerosols present
		// alpha - ratio of small to large particle sizes. (0:4,usually 1.3)
		// Results agree with the graph (pg 121, MI)
		const Float alpha = 1.3f;
		Float tauA = math::fastexp(-m * beta * std::pow(lambda/1000.0f, -alpha));  // lambda should be in um

		// Attenuation due to ozone absorption
		// lOzone - amount of ozone in cm(NTP)
		// Results agree with the graph (pg 128, MI)
		const Float lOzone = .35f;
		Float tauO = math::fastexp(-m * k_oCurve.eval(lambda) * lOzone);

		// Attenuation due to mixed gases absorption
		// Results agree with the graph (pg 131, MI)
		Float tauG = math::fastexp(-1.41f * k_gCurve.eval(lambda) * m / std::pow(1 + 118.93f
			* k_gCurve.eval(lambda) * m, (Float) 0.45f));

		// Attenuation due to water vapor absorbtion
		// w - precipitable water vapor in centimeters (standard = 2)
		// Results agree with the graph (pg 132, MI)
		const Float w = 2.0;
		Float tauWA = math::fastexp(-0.2385f * k_waCurve.eval(lambda) * w * m /
				std::pow(1 + 20.07f * k_waCurve.eval(lambda) * w * m, (Float) 0.45f));

		data[i] = solCurve.eval(lambda) * tauR * tauA * tauO * tauG * tauWA;
		wavelengths[i] = lambda;
    }

    InterpolatedSpectrum interpolated(wavelengths, data, 91);
	Spectrum discretized;
	discretized.fromContinuousSpectrum(interpolated);
	discretized.clampNegative();

	return discretized;
}


MTS_NAMESPACE_END

#endif /* __SUN_H */
