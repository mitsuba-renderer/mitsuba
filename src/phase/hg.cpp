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

#include <mitsuba/render/medium.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

/**
 * Phase function by Henyey and Greenstein (1941). Parameterizable
 * from backward- through isotropic- to forward scattering.
 */
class HGPhaseFunction : public PhaseFunction {
public:
	HGPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
		/* Asymmetry parameter: must
		   lie in [-1, 1] where >0 is forward scattering and <0 is backward
		   scattering. */
		m_g = props.getFloat("g", 0.8f);
	}

	HGPhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
		m_g = stream->readFloat();
	}

	virtual ~HGPhaseFunction() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);

		stream->writeFloat(m_g);
	}

	inline Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
			ESampledType &sampledType, const Point2 &sample) const {

		Float cosTheta;
		if (std::abs(m_g) < Epsilon) {
			cosTheta = 1 - 2*sample.x;
		} else {
			Float sqrTerm = (1 - m_g * m_g) / (1 - m_g + 2 * m_g * sample.x);
			cosTheta = (1 + m_g * m_g - sqrTerm * sqrTerm) / (2 * m_g);
		}

		Float sinTheta = std::sqrt(std::max((Float) 0, 1.0f-cosTheta*cosTheta));
		Float phi = 2*M_PI*sample.y, cosPhi = std::cos(phi), sinPhi = std::sin(phi);

		Vector dir(
			sinTheta * cosPhi,
			sinTheta * sinPhi,
			cosTheta);
		wo = Frame(-wi).toWorld(dir);
		
		sampledType = ENormal;
		return Spectrum(1.0f);
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
			ESampledType &sampledType, Float &pdf, const Point2 &sample) const {
		HGPhaseFunction::sample(mRec, wi, wo, sampledType, sample);

		pdf = f(mRec, wi, wo)[0];
		return Spectrum(pdf);
	}

	Spectrum f(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const {
		return Spectrum(1/(4*M_PI) * (1 - m_g*m_g) /
			std::pow(1.f + m_g*m_g - 2.f * m_g * dot(-wi, wo), (Float) 1.5f));
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HGPhaseFunction[g=" << m_g << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_g;
};

MTS_IMPLEMENT_CLASS_S(HGPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(HGPhaseFunction, "Henyey-Greenstein phase function");
MTS_NAMESPACE_END
