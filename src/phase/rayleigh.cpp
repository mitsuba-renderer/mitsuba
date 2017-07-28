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
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

inline Float cuberoot(Float x) {
    return (x < 0.0f) ?
            -pow(-x, (Float) (1.0f/3.0f)) :
            pow(x, (Float) (1.0f/3.0f));
}

/*!\plugin{rayleigh}{Rayleigh phase function}
 * \order{3}
 *
 * Scattering by particles that are much smaller than the wavelength
 * of light (e.g. individual molecules in the atmosphere) is well-approximated
 * by the Rayleigh phase function. This plugin implements an unpolarized
 * version of this scattering model (i.e the effects of polarization are ignored).
 * This plugin is useful for simulating scattering in planetary atmospheres.
 *
 * This model has no parameters.
 */
class RayleighPhaseFunction : public PhaseFunction {
public:
    RayleighPhaseFunction(const Properties &props)
        : PhaseFunction(props) { }

    RayleighPhaseFunction(Stream *stream, InstanceManager *manager)
        : PhaseFunction(stream, manager) { }

    virtual ~RayleighPhaseFunction() { }

    void serialize(Stream *stream, InstanceManager *manager) const {
        PhaseFunction::serialize(stream, manager);
    }

    inline Float sample(PhaseFunctionSamplingRecord &pRec,
            Sampler *sampler) const {
        Point2 sample(sampler->next2D());

        Float z = 2 * (2*sample.x - 1),
              tmp = std::sqrt(z*z+1),
              A = cuberoot(z+tmp),
              B = cuberoot(z-tmp),
              cosTheta = A + B,
              sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
              phi = 2*M_PI*sample.y,
              cosPhi = std::cos(phi),
              sinPhi = std::sin(phi);

        Vector dir(
            sinTheta * cosPhi,
            sinTheta * sinPhi,
            cosTheta);

        pRec.wo = Frame(-pRec.wi).toWorld(dir);
        return 1.0f;
    }

    Float sample(PhaseFunctionSamplingRecord &pRec,
            Float &pdf, Sampler *sampler) const {
        RayleighPhaseFunction::sample(pRec, sampler);
        pdf = RayleighPhaseFunction::eval(pRec);
        return 1.0f;
    }


    Float eval(const PhaseFunctionSamplingRecord &pRec) const {
        Float mu = dot(pRec.wi, pRec.wo);
        return (3.0f/(16.0f*M_PI)) * (1+mu*mu);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "RayleighPhaseFunction[]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(RayleighPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(RayleighPhaseFunction, "Rayleigh phase function");
MTS_NAMESPACE_END
