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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include "microflake.h"

MTS_NAMESPACE_BEGIN

class MicroflakePhaseFunction : public PhaseFunction {
public:
	MicroflakePhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
	}

	MicroflakePhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
	}

	virtual ~MicroflakePhaseFunction() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
	}

	void configure() {
	}

	Spectrum sample(PhaseFunctionQueryRecord &pRec, 
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());
		pRec.wo = squareToSphere(sample);
		return Spectrum(1.0f);
	}

	Spectrum sample(PhaseFunctionQueryRecord &pRec, 
			Float &pdf, Sampler *sampler) const {
		pRec.wo = squareToSphere(sampler->next2D());
		pdf = 1/(4 * (Float) M_PI);
		return Spectrum(pdf);
	}

	Spectrum f(const PhaseFunctionQueryRecord &pRec) const {
		return Spectrum(1/(4 * (Float) M_PI));
	}

	std::string toString() const {
		return "MicroflakePhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
};


MTS_IMPLEMENT_CLASS_S(MicroflakePhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(MicroflakePhaseFunction, "Microflake phase function");
MTS_NAMESPACE_END
