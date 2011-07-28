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

Point2 configureSunPosition(const Vector& sunDir, const Transform &luminaireToWorld) {
	return fromSphere(normalize(luminaireToWorld(sunDir)));
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
			props.getTransform("toWorld", Transform()));
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

MTS_NAMESPACE_END

#endif /* __SUN_H */
