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

MTS_NAMESPACE_BEGIN

/**
 * Basic isotropic phase function
 */
class IsotropicPhaseFunction : public PhaseFunction {
public:
	IsotropicPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
	}

	IsotropicPhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
	}

	virtual ~IsotropicPhaseFunction() { }


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
			ESampledType &sampledType, const Point2 &sample) const {
		wo = squareToSphere(sample);
		sampledType = ENormal;
		return Spectrum(1.0f);
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
			ESampledType &sampledType, Float &pdf, const Point2 &sample) const {
		wo = squareToSphere(sample);
		sampledType = ENormal;
		pdf = 1/(4 * (Float) M_PI);
		return Spectrum(pdf);
	}

	Spectrum f(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const {
		return Spectrum(1/(4 * (Float) M_PI));
	}

	std::string toString() const {
		return "IsotropicPhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
};


MTS_IMPLEMENT_CLASS_S(IsotropicPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(IsotropicPhaseFunction, "Isotropic phase function");
MTS_NAMESPACE_END
