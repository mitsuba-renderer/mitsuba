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
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{hg}{Henyey-Greenstein phase function}
 * \order{2}
 * \parameters{
 *     \parameter{g}{\Float}{
 *       This parameter must be somewhere in the range $-1$ to $1$
 *       (but not equal to $-1$ or $1$). It denotes the \emph{mean cosine} 
 *       of scattering interactions. A value greater than zero indicates that
 *       medium interactions predominantly scatter incident light into a similar
 *       direction (i.e. the medium is \emph{forward-scattering}), whereas
 *       values smaller than zero cause the medium to be
 *       scatter more light in the opposite direction.
 *     }
 * }
 * This plugin implements the phase function model proposed by 
 * Henyey and Greenstein \cite{Henyey1941Diffuse}. It is 
 * parameterizable from backward- ($g<0$) through 
 * isotropic- ($g=0$) to forward ($g>0$) scattering.
 */
class HGPhaseFunction : public PhaseFunction {
public:
	HGPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
		/* Asymmetry parameter: must
		   lie in [-1, 1] where >0 is forward scattering and <0 is backward
		   scattering. */
		m_g = props.getFloat("g", 0.8f);
		if (m_g >= 1 || m_g <= -1)
			Log(EError, "Asymmetry parameter must be in the interval (-1, 1)!");
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

	inline Float sample(PhaseFunctionQueryRecord &pRec,
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());

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
		pRec.wo = Frame(-pRec.wi).toWorld(dir);

		return 1.0f;
	}

	Float sample(PhaseFunctionQueryRecord &pRec,
			Float &pdf, Sampler *sampler) const {
		HGPhaseFunction::sample(pRec, sampler);
		pdf = HGPhaseFunction::eval(pRec);
		return 1.0f;
	}


	Float eval(const PhaseFunctionQueryRecord &pRec) const {
		return 1/(4*M_PI) * (1 - m_g*m_g) /
			std::pow(1.f + m_g*m_g - 2.f * m_g * dot(-pRec.wi, pRec.wo), (Float) 1.5f);
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
