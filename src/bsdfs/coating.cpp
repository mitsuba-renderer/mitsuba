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
#include <mitsuba/hw/basicshader.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*! \plugin{coating}{Smooth dielectric coating}
 * \order{9}
 * \icon{bsdf_coating}
 *
 * \parameters{
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{bk7} / 1.5046}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{thickness}{\Float}{Denotes the thickness of the layer (to 
 *      model absorption --- should be specified in inverse units of \code{sigmaA})\default{1}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{The absorption coefficient of the 
 *      coating layer. \default{0, i.e. there is no absorption}}
 *     \parameter{specular\showbreak Transmittance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular transmission component. Note 
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 *     \parameter{\Unnamed}{\BSDF}{A nested BSDF model that should be coated.}
 * }
 * 
 * \renderings{
 *     \rendering{Rough copper}
 *         {bsdf_coating_uncoated}
 *     \rendering{The same material coated with a single layer of 
 *         clear varnish (see \lstref{coating-roughcopper})}
 *         {bsdf_coating_roughconductor}
 * }
 *
 * This plugin implements a smooth dielectric coating (e.g. a layer of varnish) 
 * in the style of the paper ``Arbitrarily Layered Micro-Facet Surfaces'' by 
 * Weidlich and Wilkie \cite{Weidlich2007Arbitrarily}. Any BSDF in Mitsuba 
 * can be coated using this plugin, and multiple coating layers can even 
 * be applied in sequence. This allows designing interesting custom materials 
 * like car paint or glazed metal foil. The coating layer can optionally be 
 * tinted (i.e. filled with an absorbing medium), in which case this model also 
 * accounts for the directionally dependent absorption within the layer.
 *
 * Note that the plugin discards illumination that undergoes internal
 * reflection within the coating. This can lead to a noticeable energy
 * loss for materials that reflect much of their energy near or below the critical
 * angle (i.e. diffuse or very rough materials). 
 * Therefore, users are discouraged to use this plugin to coat smooth
 * diffuse materials, since there is a separately available plugin
 * named \pluginref{plastic}, which covers the same case and does not
 * suffer from energy loss.\newpage
 *
 * \renderings{
 *     \smallrendering{$\code{thickness}=0$}{bsdf_coating_0}
 *     \smallrendering{$\code{thickness}=1$}{bsdf_coating_1}
 *     \smallrendering{$\code{thickness}=5$}{bsdf_coating_5}
 *     \smallrendering{$\code{thickness}=15$}{bsdf_coating_15}
 *     \caption{The effect of the layer thickness parameter on
 *        a tinted coating ($\code{sigmaT}=(0.1, 0.2, 0.5)$)}
 * }
 *
 * \vspace{4mm}
 *
 * \begin{xml}[caption=Rough copper coated with a transparent layer of 
 *     varnish, label=lst:coating-roughcopper]
 * <bsdf type="coating">
 *     <float name="intIOR" value="1.7"/>
 *     
 *     <bsdf type="roughconductor">
 *         <string name="material" value="Cu"/>
 *         <float name="alpha" value="0.1"/>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 *
 * \renderings{
 *     \rendering{Coated rough copper with a bump map applied on top}{bsdf_coating_coatedbump}
 *     \rendering{Bump mapped rough copper with a coating on top}{bsdf_coating_bumpcoating}
 *     \caption{Some interesting materials can be created simply by applying
 *     Mitsuba's material modifiers in different orders.}
 * }
 *
 * \subsubsection*{Technical details}
 * Evaluating the internal component of this model entails refracting the 
 * incident and exitant rays through the dielectric interface, followed by
 * querying the nested material with this modified direction pair. The result 
 * is attenuated by the two Fresnel transmittances and the absorption, if
 * any.
 */
class SmoothCoating : public BSDF {
public:
	/// \sa refractTo()
	enum EDestination {
		EInterior = 0,
		EExterior = 1
	};

	SmoothCoating(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = lookupIOR(props, "intIOR", "bk7");

		/* Specifies the external index of refraction at the interface */
		m_extIOR = lookupIOR(props, "extIOR", "air");

		/* Specifies the layer's thickness using the inverse units of sigmaA */
		m_thickness = props.getFloat("thickness", 1);

		/* Specifies the absorption within the layer */
		m_sigmaA = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaA", Spectrum(0.0f)));

		/* Specifies a multiplier for the specular reflectance component */
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		if (m_intIOR < 0 || m_extIOR < 0 || m_intIOR == m_extIOR)
			Log(EError, "The interior and exterior indices of "
				"refraction must be positive and differ!");
	}

	SmoothCoating(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_thickness = stream->readFloat();
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_sigmaA = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");

		unsigned int extraFlags = 0;
		if (!m_sigmaA->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		for (int i=0; i<m_nested->getComponentCount(); ++i) 
			m_components.push_back(m_nested->getType(i) | extraFlags);

		m_components.push_back(EDeltaReflection | EFrontSide | EBackSide
			| (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));

		m_usesRayDifferentials = m_nested->usesRayDifferentials()
			|| m_sigmaA->usesRayDifferentials()
			|| m_specularReflectance->usesRayDifferentials();

		/* Compute weights that further steer samples towards
		   the specular or nested components */
		Float avgAbsorption = (m_sigmaA->getAverage()
			 *(-2*m_thickness)).exp().average();

		m_specularSamplingWeight = 1.0f / (avgAbsorption + 1.0f);

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);

		BSDF::configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		stream->writeFloat(m_thickness);
		manager->serialize(stream, m_nested.get());
		manager->serialize(stream, m_sigmaA.get());
		manager->serialize(stream, m_specularReflectance.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (m_nested != NULL)
				Log(EError, "Only a single nested BRDF can be added!");
			m_nested = static_cast<BSDF *>(child);
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "sigmaA") {
			m_sigmaA = static_cast<Texture *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	/// Refraction in local coordinates 
	Vector refractTo(EDestination dest,
			const Vector &wi, Float &F) const {
		Float etaI, etaT;
		if (dest == EInterior) {
			etaI = m_extIOR;
			etaT = m_intIOR;
		} else {
			etaI = m_intIOR;
			etaT = m_extIOR;
		}

		Float cosThetaI = Frame::cosTheta(wi);
		bool entering = cosThetaI > 0.0f;

		/* Using Snell's law, calculate the squared sine of the
		   angle between the normal and the transmitted ray */
		Float eta = etaI / etaT,
			  sinThetaTSqr = eta*eta * Frame::sinTheta2(wi);

		if (sinThetaTSqr >= 1.0f) {
			/* Total internal reflection */
			F = 1.0f;

			return Vector(0.0f);
		} else {
			Float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

			/* Compute the Fresnel transmittance */
			F = fresnelDielectric(std::abs(Frame::cosTheta(wi)),
				cosThetaT, etaI, etaT);

			/* Retain the directionality of the vector */
			return Vector(eta*wi.x, eta*wi.y, 
				entering ? cosThetaT : -cosThetaT);
		}
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		if (measure == EDiscrete && sampleSpecular &&
			std::abs(1-dot(reflect(bRec.wi), bRec.wo)) < Epsilon) {
			return m_specularReflectance->getValue(bRec.its) *
				fresnel(std::abs(Frame::cosTheta(bRec.wi)), m_extIOR, m_intIOR);
		} else if (sampleNested) {
			Float R12, R21;
			BSDFQueryRecord bRecInt(bRec);
			bRecInt.wi = refractTo(EInterior, bRec.wi, R12);
			bRecInt.wo = refractTo(EInterior, bRec.wo, R21);

			if (R12 == 1 || R21 == 1) /* Total internal reflection */
				return Spectrum(0.0f);

			Spectrum result = m_nested->eval(bRecInt, measure)
				* (1-R12) * (1-R21);

			Spectrum sigmaA = m_sigmaA->getValue(bRec.its) * m_thickness;
			if (!sigmaA.isZero()) 
				result *= (-sigmaA *
					(1/std::abs(Frame::cosTheta(bRecInt.wi)) +
					 1/std::abs(Frame::cosTheta(bRecInt.wo)))).exp();

			if (measure == ESolidAngle) {
				Float eta = m_extIOR / m_intIOR;
				/* Solid angle compression & irradiance conversion factors */
				result *= eta * eta * 
					  Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)
				   / (Frame::cosTheta(bRecInt.wi) * Frame::cosTheta(bRecInt.wo));
			}

			return result;
		}
	
		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);
		
		Float R12;
		Vector wiPrime = refractTo(EInterior, bRec.wi, R12);

		/* Reallocate samples */
		Float probSpecular = (R12*m_specularSamplingWeight) /
			(R12*m_specularSamplingWeight + 
			(1-R12) * (1-m_specularSamplingWeight));

		if (measure == EDiscrete && sampleSpecular &&
			std::abs(1-dot(reflect(bRec.wi), bRec.wo)) < Epsilon) {
			return sampleNested ? probSpecular : 1.0f;
		} else if (sampleNested) {
			Float R21;
			BSDFQueryRecord bRecInt(bRec);
			bRecInt.wi = wiPrime;
			bRecInt.wo = refractTo(EInterior, bRec.wo, R21);

			if (R12 == 1 || R21 == 1) /* Total internal reflection */
				return 0.0f;

			Float pdf = m_nested->pdf(bRecInt, measure);

			if (measure == ESolidAngle) {
				Float eta = m_extIOR / m_intIOR;
				pdf *= eta * eta * Frame::cosTheta(bRec.wo)
			          / Frame::cosTheta(bRecInt.wo);
			}

			return sampleSpecular ? (pdf * (1 - probSpecular)) : pdf;
		} else {
			return 0.0f;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &_sample) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		if ((!sampleSpecular && !sampleNested))
			return Spectrum(0.0f);

		Float R12;
		Vector wiPrime = refractTo(EInterior, bRec.wi, R12);

		/* Reallocate samples */
		Float probSpecular = (R12*m_specularSamplingWeight) /
			(R12*m_specularSamplingWeight + 
			(1-R12) * (1-m_specularSamplingWeight));

		bool choseSpecular = sampleSpecular;

		Point2 sample(_sample);
		if (sampleSpecular && sampleNested) {
			if (sample.x > probSpecular) {
				sample.x = (sample.x - probSpecular) / (1 - probSpecular);
				choseSpecular = false;
			} else {
				sample.x /= probSpecular;
			}
		}

		if (choseSpecular) {
			bRec.sampledComponent = m_components.size()-1;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = sampleNested ? probSpecular : 1.0f;
			return m_specularReflectance->getValue(bRec.its) * (R12/pdf);
		} else {
			if (R12 == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			Vector wiBackup = bRec.wi;
			bRec.wi = wiPrime;
			Spectrum result = m_nested->sample(bRec, pdf, sample);
			bRec.wi = wiBackup;
			if (result.isZero()) 
				return Spectrum(0.0f);

			Vector woPrime = bRec.wo;

			Spectrum sigmaA = m_sigmaA->getValue(bRec.its) * m_thickness;
			if (!sigmaA.isZero()) 
				result *= (-sigmaA *
					(1/std::abs(Frame::cosTheta(wiPrime)) +
					 1/std::abs(Frame::cosTheta(woPrime)))).exp();

			Float R21;
			bRec.wo = refractTo(EExterior, woPrime, R21);
			if (R21 == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			if (sampleSpecular) {
				pdf *= 1.0f - probSpecular;
				result /= 1.0f - probSpecular;
			}

			result *= (1 - R12) * (1 - R21);

			if (BSDF::getMeasure(bRec.sampledType) == ESolidAngle) {
				/* Solid angle compression & irradiance conversion factors */
				Float eta = m_extIOR / m_intIOR, etaSqr = eta*eta;

				result *= Frame::cosTheta(bRec.wi) / Frame::cosTheta(wiPrime);
				pdf *= etaSqr * Frame::cosTheta(bRec.wo) / Frame::cosTheta(woPrime);
			}

			return result;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf;
		return SmoothCoating::sample(bRec, pdf, sample);
	}

	Shader *createShader(Renderer *renderer) const;

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothCoating[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
			<< "  sigmaA = " << indent(m_sigmaA->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  thickness = " << m_thickness << "," << endl
			<< "  nested = " << indent(m_nested.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	Float m_specularSamplingWeight;
	Float m_intIOR, m_extIOR;
	ref<Texture> m_sigmaA;
	ref<Texture> m_specularReflectance;
	ref<BSDF> m_nested;
	Float m_thickness;
};

// ================ Hardware shader implementation ================ 

/**
 * Simple GLSL version -- uses Schlick's approximation and approximates the
 * ideally specular reflection with a somewhat smoothed out reflection lobe
 */
class SmoothCoatingShader : public Shader {
public:
	SmoothCoatingShader(Renderer *renderer, Float extIOR, Float intIOR, const BSDF *nested,
			const Texture *sigmaA) : Shader(renderer, EBSDFShader), m_nested(nested), m_sigmaA(sigmaA) {
		m_nestedShader = renderer->registerShaderForResource(m_nested.get());
		m_sigmaAShader = renderer->registerShaderForResource(m_sigmaA.get());
		m_R0 = fresnel(1.0f, extIOR, intIOR);
		m_eta = extIOR / intIOR;
	}

	bool isComplete() const {
		return m_nestedShader.get() != NULL
			&& m_sigmaAShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nested.get());
		renderer->unregisterShaderForResource(m_sigmaA.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedShader.get());
		deps.push_back(m_sigmaAShader.get());
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_eta", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_alpha", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_R0);
		program->setParameter(parameterIDs[1], m_eta);
		program->setParameter(parameterIDs[2], 0.4f);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_R0;" << endl
			<< "uniform float " << evalName << "_eta;" << endl
			<< "uniform float " << evalName << "_alpha;" << endl
			<< endl
			<< "float " << evalName << "_schlick(float ct) {" << endl
			<< "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
			<< "    return " << evalName << "_R0 + (1.0 - " << evalName << "_R0) * ct5;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_refract(vec3 wi, out float T) {" << endl
			<< "    float cosThetaI = cosTheta(wi);" << endl
			<< "    bool entering = cosThetaI > 0.0;" << endl
			<< "    float eta = " << evalName << "_eta;" << endl
			<< "    float sinThetaTSqr =  eta * eta * sinTheta2(wi);" << endl
			<< "    if (sinThetaTSqr >= 1.0) {" << endl
			<< "        T = 0.0; /* Total internal reflection */" << endl
			<< "        return vec3(0.0);" << endl
			<< "    } else {" << endl
			<< "        float cosThetaT = sqrt(1.0 - sinThetaTSqr);" << endl
			<< "        T = 1.0 - " << evalName << "_schlick(1.0 - abs(cosThetaI));" << endl
			<< "        return vec3(eta*wi.x, eta*wi.y, entering ? cosThetaT : -cosThetaT);" << endl
			<< "    }" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_D(vec3 m) {" << endl
			<< "    float ct = cosTheta(m);" << endl
			<< "    if (cosTheta(m) <= 0.0)" << endl
			<< "        return 0.0;" << endl
			<< "    float ex = tanTheta(m) / " << evalName << "_alpha;" << endl
			<< "    return exp(-(ex*ex)) / (pi * " << evalName << "_alpha" << endl
			<< "        * " << evalName << "_alpha * pow(cosTheta(m), 4.0));" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
			<< "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
			<< "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
			<< "        return 0.0;" << endl
			<< "    float nDotM = cosTheta(m);" << endl
			<< "    return min(1.0, min(" << endl
			<< "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
			<< "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    float T12, T21;" << endl
			<< "    vec3 wiPrime = " << evalName << "_refract(wi, T12);" << endl
			<< "    vec3 woPrime = " << evalName << "_refract(wo, T21);" << endl
			<< "    vec3 nested = " << depNames[0] << "(uv, wiPrime, woPrime);" << endl
			<< "    vec3 sigmaA = " << depNames[1] << "(uv);" << endl
			<< "    vec3 result = nested * " << evalName << "_eta * " << evalName << "_eta" << endl
			<< "                  * T12 * T21 * (cosTheta(wi)*cosTheta(wo)) /" << endl
			<< "                  (cosTheta(wiPrime)*cosTheta(woPrime));" << endl
			<< "    if (sigmaA != vec3(0.0))" << endl
			<< "        result *= exp(-sigmaA * (1/abs(cosTheta(wiPrime)) + " << endl
			<< "                                 1/abs(cosTheta(woPrime))));" << endl
			<< "    if (cosTheta(wi)*cosTheta(wo) > 0) {" << endl
			<< "        vec3 H = normalize(wi + wo);" << endl
			<< "        float D = " << evalName << "_D(H)" << ";" << endl
			<< "        float G = " << evalName << "_G(H, wi, wo);" << endl
			<< "        float F = " << evalName << "_schlick(1-dot(wi, H));" << endl
			<< "        result += vec3(F * D * G / (4*cosTheta(wi)));" << endl
			<< "    }" << endl
			<< "    return result;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const BSDF> m_nested;
	ref<Shader> m_nestedShader;
	ref<const Texture> m_sigmaA;
	ref<Shader> m_sigmaAShader;
	Float m_R0, m_eta;
};

Shader *SmoothCoating::createShader(Renderer *renderer) const { 
	return new SmoothCoatingShader(renderer, m_extIOR, m_intIOR, 
		m_nested.get(), m_sigmaA.get());
}

MTS_IMPLEMENT_CLASS(SmoothCoatingShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothCoating, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothCoating, "Smooth dielectric coating");
MTS_NAMESPACE_END
