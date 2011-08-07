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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{isotropic}{Isotropic phase function}
 * \order{1}
 
 * \renderings{
 *     \rendering{Isotropic}{phase_isotropic}
 *     \rendering{Anisotropic micro-flakes}{phase_microflakes_005}
 *     \caption{Heterogeneous volume renderings of a scarf model
 *     with isotropic and anisotropic phase functions.}
 * }
 *
 *
 * This phase function simulates completely uniform scattering,
 * where all directionality is lost after a single scattering 
 * interaction. It does not have any parameters.
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

	Float sample(PhaseFunctionQueryRecord &pRec, 
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());
		pRec.wo = squareToSphere(sample);
		return 1.0f;
	}

	Float sample(PhaseFunctionQueryRecord &pRec, 
			Float &pdf, Sampler *sampler) const {
		pRec.wo = squareToSphere(sampler->next2D());
		pdf = 1/(4 * (Float) M_PI);
		return 1.0f;
	}

	Float eval(const PhaseFunctionQueryRecord &pRec) const {
		return 1/(4 * (Float) M_PI);
	}

	std::string toString() const {
		return "IsotropicPhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
};


MTS_IMPLEMENT_CLASS_S(IsotropicPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(IsotropicPhaseFunction, "Isotropic phase function");
MTS_NAMESPACE_END
