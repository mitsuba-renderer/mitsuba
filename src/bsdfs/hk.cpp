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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{hk}{Hanrahan-Krueger BRDF}
 *
 * This plugin is an implementation of the Hanrahan-Krueger scattering model 
 * \cite{Hanrahan1993Reflection}. This model provides an analytic BSDF that 
 * accounts for single scattering within thin slabs of scattering media having
 * index matched boundaries. 
 *
 * It also accounts for $0^{\mathrm{th}}$-order scattered light (i.e. attenuated
 * light) that did not scatter within the medium. When used in conjuction with 
 * the \pluginref{coating} plugin, it can also acount for refraction and reflection
 * at the boundaries of the medium when the indices of refraction are mismatched.
 *
 * This BSDF needs a nested phase function to model the scattering within the
 * medium. When no phase function is given, it will use an isotropic one ($g=0$)
 * by default. A sample usage is given below: 
 *
 * \begin{xml}
 * <bsdf type="hk">
 *     <spectrum name="sigmaS" value="4"/>
 *     <spectrum name="sigmaA" value="0.1"/>
 *     <float name="thickness" value="0.5"/>
 *     <phase type="hg">
 *         <float name="g" value="0.98"/>
 *     </phase>
 * </bsdf>
 * \end{xml}
 *
 * When \texttt{sigmaS} = \texttt{sigmaA}$\ = 0$, any associated geometry 
 * will be invisible.
 *
 * The sampling of this BSDF is either with respect to the phase function PDF or
 * with the delta transmission function. The weighting between the two is decided
 * based on the probability of an event within the medium.
 *
 * The implementation is based on code by Tom Kamzimier and Marios Papas.
*/

class HanrahanKrueger : public BSDF {
public:
	HanrahanKrueger(const Properties &props) : BSDF(props) {
		/* Scattering events per scene unit */
		m_sigmaS = props.getSpectrum("sigmaS", Spectrum(4.0f)); 
		
		/* Absorption events per scene unit */
		m_sigmaA = props.getSpectrum("sigmaA", Spectrum(0.1f));  

		/* Slab Thickness in scene units*/
		m_d = props.getFloat("thickness",0.5); 

	}

	HanrahanKrueger(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_sigmaS = Spectrum(stream);
		m_sigmaA = Spectrum(stream);
		m_d      = stream->readFloat();

		configure();
	}

	void configure() {
		if (m_phase == NULL)
			m_phase = static_cast<PhaseFunction *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(PhaseFunction), Properties("isotropic")));

		m_components.clear();
		m_components.push_back(EGlossyReflection   | EFrontSide | EBackSide | ECanUseSampler);
		m_components.push_back(EGlossyTransmission | EFrontSide | EBackSide | ECanUseSampler);
		m_components.push_back(EDeltaTransmission  | EFrontSide | EBackSide | ECanUseSampler);

		/* Precompute some helper quantities to save rendering cycles */
		m_sigmaT = m_sigmaS + m_sigmaA;
		m_tauD   = m_sigmaT * m_d;
		/* Avoid divisions by 0 */
		for (int i = 0; i < SPECTRUM_SAMPLES; i++)
			m_albedo[i] = m_sigmaT[i] > 0.0f ? (m_sigmaS[i]/m_sigmaT[i]) : 0.0f;

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_albedo; /* Very approximate .. */
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		Spectrum result(0.0f);

		if (measure == EDiscrete) {
			/* Figure out if the specular transmission is specifically requested */
			bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 2);

			/* Return the attenuated light if requested */
			if (hasSpecularTransmission &&
				std::abs(1+dot(bRec.wi, bRec.wo)) < Epsilon)
				result = (-m_tauD/std::abs(Frame::cosTheta(bRec.wi))).exp();
		} else if (measure == ESolidAngle) {
			/* Sample single scattering events */
			bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
			bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

			const Float cosThetaI = Frame::cosTheta(bRec.wi);
			const Float cosThetaO = Frame::cosTheta(bRec.wo);
			bool reflection = cosThetaI * cosThetaO > 0;

			/* ==================================================================== */
			/*                        Reflection component                          */
			/* ==================================================================== */

			if (hasGlossyReflection && reflection) {
				MediumSamplingRecord dummy;
				PhaseFunctionQueryRecord pRec(dummy,bRec.wi,bRec.wo); 
				const Float phaseVal = m_phase->eval(pRec);

				result = m_albedo * (phaseVal*cosThetaI/(cosThetaI+cosThetaO)) *
					(Spectrum(1.0f)-((-m_tauD/std::abs(cosThetaI))+(-m_tauD/std::abs(cosThetaO))).exp());
			}

			/* ==================================================================== */
			/*                       Transmission component                         */
			/* ==================================================================== */

			if (hasGlossyTransmission && !reflection) {
				MediumSamplingRecord dummy;
				PhaseFunctionQueryRecord pRec(dummy,bRec.wi,bRec.wo);
				const Float phaseVal = m_phase->eval(pRec);

				/* Hanrahan etal 93 Single Scattering transmission term */
				if (cosThetaI + cosThetaO == 0.0f) {
					/* avoid division by zero */
					result += m_albedo * phaseVal*m_tauD/std::abs(cosThetaO) * 
									((-m_tauD/std::abs(cosThetaO)).exp());
				} else {
					/* Guaranteed to be positive even if |cosThetaO| > |cosThetaI| */
					result += m_albedo * phaseVal*std::abs(cosThetaI)/(std::abs(cosThetaI)-std::abs(cosThetaO)) * 
						((-m_tauD/std::abs(cosThetaI)).exp() - (-m_tauD/std::abs(cosThetaO)).exp());
				}
			}
		}


		return result;
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool hasSingleScattering = (bRec.typeMask & EGlossy)
			&& (bRec.component == -1 || bRec.component == 0 || bRec.component == 1);
		bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
			&& (bRec.component == -1 || bRec.component == 2);

		Float probSpecularTransmission = (-m_tauD/std::abs(Frame::cosTheta(bRec.wi))).exp().average();

		if (measure == EDiscrete) {
			bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 2);
			/* Return the attenuated light if requested */
			if (hasSpecularTransmission &&
				std::abs(1+dot(bRec.wi, bRec.wo)) < Epsilon)
				return hasSingleScattering ? probSpecularTransmission : 1.0f;
		} else if (hasSingleScattering && measure == ESolidAngle) {
			bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
			bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
			bool reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0;

			if ((!hasGlossyReflection && reflection) ||
				(!hasGlossyTransmission && !reflection))
				return 0.0f;

			/* Sampled according to the phase function lobe(s) */
			MediumSamplingRecord dummy;
			PhaseFunctionQueryRecord pRec(dummy, bRec.wi, bRec.wo);
			Float pdf = m_phase->pdf(pRec);
			if (hasSpecularTransmission)
				pdf *= 1-probSpecularTransmission;
			return pdf;
		}
		return 0.0f;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, Float &_pdf, const Point2 &_sample) const {
		AssertEx(bRec.sampler != NULL, "The BSDFQueryRecord needs to have a sampler!");

		bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
			&& (bRec.component == -1 || bRec.component == 2);
		bool hasSingleScattering = (bRec.typeMask & EGlossy)
			&& (bRec.component == -1 || bRec.component == 0 || bRec.component == 1);

		/* Probability for a specular transmission is approximated by the average (per wavelength) 
		 * probability of a photon exiting without a scattering event or an absorption event */
		Float probSpecularTransmission = (-m_tauD/std::abs(Frame::cosTheta(bRec.wi))).exp().average();

		bool choseSpecularTransmission = hasSpecularTransmission;

		Point2 sample(_sample);
		if (hasSpecularTransmission && hasSingleScattering) {
			if (sample.x > probSpecularTransmission) {
				sample.x = (sample.x - probSpecularTransmission) / (1 - probSpecularTransmission);
				choseSpecularTransmission = false;
			}
		}

		if (choseSpecularTransmission) {
			/* The specular transmission component was sampled */
			bRec.sampledComponent = 2;
			bRec.sampledType = EDeltaTransmission;

			bRec.wo = -bRec.wi;

			_pdf = hasSingleScattering ? probSpecularTransmission : 1.0f;
			return eval(bRec, EDiscrete);
		} else {
			/* The glossy transmission/scattering component should be sampled */
			bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
			bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

			/* Sample According to the phase function lobes */
			PhaseFunctionQueryRecord pRec(MediumSamplingRecord(), bRec.wi, bRec.wo);
			m_phase->sample(pRec, _pdf, bRec.sampler);

			/* Store the sampled direction */
			bRec.wo = pRec.wo;
			
			bool reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0;
			if ((!hasGlossyReflection && reflection) ||
				(!hasGlossyTransmission && !reflection))
				return Spectrum(0.0f);

			/* Notify that the scattering component was sampled */
			bRec.sampledComponent = reflection ? 0 : 1;
			bRec.sampledType = EGlossy;

			_pdf *= (hasSpecularTransmission ? (1 - probSpecularTransmission) : 1.0f);

			/* Guard against numerical imprecisions */
			if (_pdf == 0) 
				return Spectrum(0.0f);
			else
				return eval(bRec, ESolidAngle);

		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf = 0;
		Spectrum result = HanrahanKrueger::sample(bRec, pdf, sample);

		if (result.isZero())
			return Spectrum(0.0f);
		else
			return result / pdf;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		m_sigmaS.serialize(stream);
		m_sigmaA.serialize(stream);
		stream->writeFloat(m_d);

	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
			Assert(m_phase == NULL);
			m_phase = static_cast<PhaseFunction *>(child);
		} else {
			Log(EError, "Invalid child node! (\"%s\")",
				cClass->getName().c_str());
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HanrahanKrueger[" << endl
   			<< "  sigmaS = " << m_sigmaS.toString() << "," << std::endl
   			<< "  sigmaA = " << m_sigmaA.toString() << "," << std::endl
   			<< "  sigmaT = " << m_sigmaT.toString() << "," << std::endl
   			<< "  albedo = " << m_albedo.toString() << "," << std::endl
   			<< "  tauD = "   << m_tauD.toString()   << "," << std::endl
   			<< "  d = "      << m_d                 << "," << std::endl
   			<< "  phase = "  << m_phase->toString() << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<PhaseFunction> m_phase;
	Spectrum m_sigmaS, m_sigmaA;
	Float m_d;

	/* These values should not be serialized, they are evaluated 
	 * in configure() 
	 */
	Spectrum m_sigmaT, m_tauD, m_albedo; 

};

MTS_IMPLEMENT_CLASS_S(HanrahanKrueger, false, BSDF)
MTS_EXPORT_PLUGIN(HanrahanKrueger, "Hanrahan-Krueger BSDF");
MTS_NAMESPACE_END
