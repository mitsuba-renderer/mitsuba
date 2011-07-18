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
 * Weidlich and Wilkie \cite{Weidlich2007Arbitrarily}. Any non-transmissive
 * BSDF in Mitsuba can be coated using this plugin, and multiple coating layers
 * can be applied in sequence. This allows designing interesting custom materials 
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
 * suffer from energy loss.
 *
 * Evaluating the internal component of this model entails refracting the 
 * incident and exitant rays through the dielectric interface, followed by
 * querying the nested material with this modified direction pair. The result 
 * is attenuated by the two Fresnel transmittances and the absorption, if
 * any.\newpage
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
 */
class SmoothCoating : public BSDF {
public:
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
	}

	SmoothCoating(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_thickness = stream->readFloat();
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_sigmaA = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");
		if (m_nested->getType() & BSDF::ETransmission)
			Log(EError, "Tried to put a smooth coating layer on top of a BSDF "
				"with a transmission component -- this is currently not allowed!");

		unsigned int extraFlags = 0;
		if (!m_sigmaA->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		for (int i=0; i<m_nested->getComponentCount(); ++i) 
			m_components.push_back(m_nested->getType(i) | extraFlags);

		m_components.push_back(EDeltaReflection | EFrontSide);

		m_usesRayDifferentials = m_nested->usesRayDifferentials()
			|| m_sigmaA->usesRayDifferentials();

		/* Compute weights that further steer samples towards
		   the specular or nested components */
		Float avgAbsorption = (m_sigmaA->getAverage()
			 *(-2*m_thickness)).exp().average();

		m_specularSamplingWeight = 1.0f / (avgAbsorption + 1.0f);

		BSDF::configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		stream->writeFloat(m_thickness);
		manager->serialize(stream, m_nested.get());
		manager->serialize(stream, m_sigmaA.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (m_nested != NULL)
				Log(EError, "Only a single nested BRDF can be added!");
			m_nested = static_cast<BSDF *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	/**
	 * \brief Refraction in local coordinates 
	 *
	 * To be used when some of the data is already available
	 */
	inline Vector refract(const Vector &wi, Float eta, Float cosThetaT) const {
		return Vector(-eta*wi.x, -eta*wi.y, cosThetaT);
	}

	/// Refraction in local coordinates (full version)
	inline Vector refract(const Vector &wi, Float &F) const {
		Float cosThetaI = Frame::cosTheta(wi),
			  etaI = m_extIOR, etaT = m_intIOR;

		bool entering = cosThetaI > 0.0f;

		/* Determine the respective indices of refraction */
		if (!entering)
			std::swap(etaI, etaT);

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
				cosThetaT, m_extIOR, m_intIOR);

			return Vector(-eta*wi.x, -eta*wi.y, 
				entering ? -cosThetaT : cosThetaT);
		}
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		if (measure == EDiscrete && sampleSpecular &&
			std::abs(1-dot(reflect(bRec.wi), bRec.wo)) < Epsilon) {
			return Spectrum(fresnel(
				Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR));
		} else if (sampleNested) {
			Float R12, R21;
			BSDFQueryRecord bRec2(bRec);
			bRec2.wi = -refract(bRec.wi, R12);
			bRec2.wo = -refract(bRec.wo, R21);

			if (R12 == 1 || R21 == 1) /* Total internal reflection */
				return Spectrum(0.0f);

			Float eta = m_extIOR / m_intIOR;
			Spectrum result = m_nested->eval(bRec2, measure)
				* ((1-R12) * (1-R21) * eta * eta);

			Spectrum sigmaA = m_sigmaA->getValue(bRec.its) * m_thickness;
			if (!sigmaA.isZero()) 
				result *= (-sigmaA *
					(1/std::abs(Frame::cosTheta(bRec2.wi)) +
					 1/std::abs(Frame::cosTheta(bRec2.wo)))).exp();

			if (measure == ESolidAngle)
				result *= Frame::cosTheta(bRec.wo)
					    / Frame::cosTheta(bRec2.wo);

			return result;
		}
	
		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);
		
		Float R12;
		Vector wiPrime = -refract(bRec.wi, R12);

		/* Reallocate samples */
		Float probSpecular = (R12*m_specularSamplingWeight) /
			(R12*m_specularSamplingWeight + 
			(1-R12) * (1-m_specularSamplingWeight));

		if (measure == EDiscrete && sampleSpecular &&
			std::abs(1-dot(reflect(bRec.wi), bRec.wo)) < Epsilon) {
			return sampleNested ? probSpecular : 1.0f;
		} else if (sampleNested) {
			Float R21;
			BSDFQueryRecord bRec2(bRec);
			bRec2.wi = wiPrime;
			bRec2.wo = -refract(bRec.wo, R21);

			if (R12 == 1 || R21 == 1) /* Total internal reflection */
				return 0.0f;
		
			Float pdf = m_nested->pdf(bRec2, measure);

			if (measure == ESolidAngle)
				pdf *= Frame::cosTheta(bRec.wo)
					 / Frame::cosTheta(bRec2.wo);

			Float eta = m_extIOR / m_intIOR;
			pdf *= eta * eta;

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

		if ((!sampleSpecular && !sampleNested) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		/* Refract the incident direction and compute the Fresnel reflectance */
		Float eta = m_extIOR / m_intIOR,
			  sinThetaTSqr = eta*eta * Frame::sinTheta2(bRec.wi),
			  R12, cosThetaT = 0;

		if (sinThetaTSqr >= 1.0f) {
			R12 = 1.0f; /* Total internal reflection */
		} else {
			cosThetaT = -std::sqrt(1.0f - sinThetaTSqr);
			R12 = fresnelDielectric(Frame::cosTheta(bRec.wi),
					-cosThetaT, m_extIOR, m_intIOR);
		}

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
			}
		}
		
		if (choseSpecular) {
			bRec.sampledComponent = m_components.size()-1;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = sampleNested ? probSpecular : 1.0f;
			return Spectrum(R12);
		} else {
			if (R12 == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			Vector wiBackup = bRec.wi;
			bRec.wi = -refract(bRec.wi, eta, cosThetaT);

			Spectrum result = m_nested->sample(bRec, pdf, sample);
			if (result.isZero()) 
				return Spectrum(0.0f);

			Spectrum sigmaA = m_sigmaA->getValue(bRec.its) * m_thickness;
			if (!sigmaA.isZero()) 
				result *= (-sigmaA *
					(1/std::abs(Frame::cosTheta(bRec.wi)) +
					 1/std::abs(Frame::cosTheta(bRec.wo)))).exp();

			Float R21, cosThetaWoPrime = Frame::cosTheta(bRec.wo);
			bRec.wo = refract(-bRec.wo, R21);
			bRec.wi = wiBackup;

			if (R21 == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);
			bool sampledSA = (BSDF::getMeasure(bRec.sampledType) == ESolidAngle);
			Float cosRatio = Frame::cosTheta(bRec.wo) / cosThetaWoPrime,
				  commonTerms = (sampledSA ? cosRatio : 1.0f)* eta * eta;

			pdf *= (sampleSpecular ? (1 - probSpecular) : 1.0f) * commonTerms;
			result *= (1 - R12) * (1 - R21) * commonTerms;

			return result;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf;
		Spectrum result = SmoothCoating::sample(bRec, pdf, sample);

		if (result.isZero())
			return Spectrum(0.0f);
		else
			return result / pdf;
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
			<< "  thickness = " << m_thickness << "," << endl
			<< "  nested = " << indent(m_nested->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	Float m_specularSamplingWeight;
	Float m_intIOR, m_extIOR;
	ref<Texture> m_sigmaA;
	ref<BSDF> m_nested;
	Float m_thickness;
};

// ================ Hardware shader implementation ================ 

/* Crude GLSL approximation -- just forwards to the nested model */
class SmoothCoatingShader : public Shader {
public:
	SmoothCoatingShader(Renderer *renderer, const BSDF *nested) 
		: Shader(renderer, EBSDFShader), m_nested(nested) {
		m_nestedShader = renderer->registerShaderForResource(m_nested.get());
	}

	bool isComplete() const {
		return m_nestedShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nested.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << depNames[0] << "(uv, wi, wo);" << endl
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
};

Shader *SmoothCoating::createShader(Renderer *renderer) const { 
	return new SmoothCoatingShader(renderer, m_nested.get());
}

MTS_IMPLEMENT_CLASS(SmoothCoatingShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothCoating, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothCoating, "Smooth dielectric coating");
MTS_NAMESPACE_END
