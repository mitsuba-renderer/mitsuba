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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/warp.h>

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
        configure();
    }

    virtual ~IsotropicPhaseFunction() { }

    void configure() {
        PhaseFunction::configure();
        m_type = EIsotropic | EAngleDependence;
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        PhaseFunction::serialize(stream, manager);
    }

    Float sample(PhaseFunctionSamplingRecord &pRec,
            Sampler *sampler) const {
        Point2 sample(sampler->next2D());
        pRec.wo = warp::squareToUniformSphere(sample);
        return 1.0f;
    }

    Float sample(PhaseFunctionSamplingRecord &pRec,
            Float &pdf, Sampler *sampler) const {
        pRec.wo = warp::squareToUniformSphere(sampler->next2D());
        pdf = warp::squareToUniformSpherePdf();
        return 1.0f;
    }

    Float eval(const PhaseFunctionSamplingRecord &pRec) const {
        return warp::squareToUniformSpherePdf();
    }

    Float getMeanCosine() const {
        return 0.0f;
    }

    std::string toString() const {
        return "IsotropicPhaseFunction[]";
    }

    MTS_DECLARE_CLASS()
};


MTS_IMPLEMENT_CLASS_S(IsotropicPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(IsotropicPhaseFunction, "Isotropic phase function");
MTS_NAMESPACE_END
