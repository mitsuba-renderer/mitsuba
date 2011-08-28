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

#include <mitsuba/render/photon.h>

MTS_NAMESPACE_BEGIN

/* Precompute cosine/sine values for quick conversions
   from quantized spherical coordinates to floating
   point vectors. */
Float Photon::m_cosTheta[256];
Float Photon::m_sinTheta[256];
Float Photon::m_cosPhi[256];
Float Photon::m_sinPhi[256];
Float Photon::m_expTable[256];

bool Photon::m_precompTableReady = Photon::createPrecompTables();

Float NewPhoton::m_cosTheta[256];
Float NewPhoton::m_sinTheta[256];
Float NewPhoton::m_cosPhi[256];
Float NewPhoton::m_sinPhi[256];
Float NewPhoton::m_expTable[256];

bool NewPhoton::m_precompTableReady = NewPhoton::createPrecompTables();


bool Photon::createPrecompTables() {
	for (int i=0; i<256; i++) {
		Float angle = (Float) i * ((Float) M_PI / 256.0f);
		m_cosPhi[i] = std::cos(2.0f * angle);
		m_sinPhi[i] = std::sin(2.0f * angle);
		m_cosTheta[i] = std::cos(angle);
		m_sinTheta[i] = std::sin(angle);
		m_expTable[i] = std::ldexp((Float) 1, i - (128+8));
	}
	m_expTable[0] = 0;

	return true;
}

Photon::Photon(Stream *stream) {
	stream->readSingleArray(pos, 3);
#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
	power = Spectrum(stream);
	phi = stream->readUChar();
	theta = stream->readUChar();
	phiN = stream->readUChar();
	thetaN = stream->readUChar();
#else
	stream->read(power, 8);
#endif
	depth = stream->readUShort();
	axis = stream->readUChar();
	unused = 0;
}

Photon::Photon(const Point &p, const Normal &normal,
			   const Vector &dir, const Spectrum &P,
			   uint16_t _depth) {
	if (P.isNaN()) 
		SLog(EWarn, "Creating an invalid photon with power: %s", P.toString().c_str());

	/* Possibly convert to single precision floating point
	   (if Mitsuba is configured to use double precision) */
	pos[0] = (float) p.x;
	pos[1] = (float) p.y;
	pos[2] = (float) p.z;
	depth = _depth;
	unused = 0;
	axis = -1;

	/* Convert the direction into an approximate spherical 
	   coordinate format to reduce storage requirements */
	theta = (uint8_t) std::min(255,
		(int) (std::acos(dir.z) * (256.0 / M_PI)));

	int tmp = std::min(255,
		(int) (std::atan2(dir.y, dir.x) * (256.0 / (2.0 * M_PI))));
	if (tmp < 0)
		phi = (uint8_t) (tmp + 256);
	else
		phi = (uint8_t) tmp;
	
	if (normal.isZero()) {
		thetaN = phiN = 0;
	} else {
		thetaN = (uint8_t) std::min(255,
			(int) (std::acos(normal.z) * (256.0 / M_PI)));
		tmp = std::min(255,
			(int) (std::atan2(normal.y, normal.x) * (256.0 / (2.0 * M_PI))));
		if (tmp < 0)
			phiN = (uint8_t) (tmp + 256);
		else
			phiN = (uint8_t) tmp;
	}

#if defined(DOUBLE_PRECISION) || SPECTRUM_SAMPLES > 3
	power = P;
#else
	/* Pack the photon power into Greg Ward's RGBE format */
	P.toRGBE(power);
#endif
}

bool NewPhoton::createPrecompTables() {
	for (int i=0; i<256; i++) {
		Float angle = (Float) i * ((Float) M_PI / 256.0f);
		m_cosPhi[i] = std::cos(2.0f * angle);
		m_sinPhi[i] = std::sin(2.0f * angle);
		m_cosTheta[i] = std::cos(angle);
		m_sinTheta[i] = std::sin(angle);
		m_expTable[i] = std::ldexp((Float) 1, i - (128+8));
	}
	m_expTable[0] = 0;

	return true;
}

NewPhoton::NewPhoton(Stream *stream) {
	position = Point(stream);
#if defined(SINGLE_PRECISION) && SPECTRUM_SAMPLES == 3
	stream->read(data.power, 8);
#else
	data.power = Spectrum(stream);
	data.phi = stream->readUChar();
	data.theta = stream->readUChar();
	data.phiN = stream->readUChar();
	data.thetaN = stream->readUChar();
#endif
	data.depth = stream->readUShort();
	flags = stream->readUChar();
}

NewPhoton::NewPhoton(const Point &p, const Normal &normal,
			   const Vector &dir, const Spectrum &P,
			   uint16_t _depth) {
	if (!P.isValid()) 
		SLog(EWarn, "Creating an invalid photon with power: %s", P.toString().c_str());

	/* Possibly convert to single precision floating point
	   (if Mitsuba is configured to use double precision) */
	position = p;
	data.depth = _depth;
	flags = 0;

	/* Convert the direction into an approximate spherical 
	   coordinate format to reduce storage requirements */
	data.theta = (uint8_t) std::min(255,
		(int) (std::acos(dir.z) * (256.0 / M_PI)));

	int tmp = std::min(255,
		(int) (std::atan2(dir.y, dir.x) * (256.0 / (2.0 * M_PI))));
	if (tmp < 0)
		data.phi = (uint8_t) (tmp + 256);
	else
		data.phi = (uint8_t) tmp;
	
	if (normal.isZero()) {
		data.thetaN = data.phiN = 0;
	} else {
		data.thetaN = (uint8_t) std::min(255,
			(int) (std::acos(normal.z) * (256.0 / M_PI)));
		tmp = std::min(255,
			(int) (std::atan2(normal.y, normal.x) * (256.0 / (2.0 * M_PI))));
		if (tmp < 0)
			data.phiN = (uint8_t) (tmp + 256);
		else
			data.phiN = (uint8_t) tmp;
	}

#if defined(SINGLE_PRECISION) && SPECTRUM_SAMPLES == 3
	/* Pack the photon power into Greg Ward's RGBE format */
	P.toRGBE(data.power);
#else
	data.power = P;
#endif
}


MTS_NAMESPACE_END
