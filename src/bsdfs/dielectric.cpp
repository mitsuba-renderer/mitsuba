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
#include <mitsuba/render/consttexture.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{dielectric}{Smooth dielectric material}
 *
 * \parameters{
 *     \parameter{intIOR}{\Float}{Interior index of refraction \default{1.5046}}
 *     \parameter{extIOR}{\Float}{Exterior index of refraction \default{1.0}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the reflectance component\default{1.0}}
 *     \lastparameter{specular\showbreak Transmittance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the transmittance component\default{1.0}}
 * }
 *
 * \renderings{
 *     \medrendering{Air$\leftrightarrow$Water (IOR: 1.33) interface. 
 *         See \lstref{dielectric-water}.}{bsdf_dielectric_water}
 *     \medrendering{Air$\leftrightarrow$Diamond (IOR: 2.419)}{bsdf_dielectric_diamond}
 *     \medrendering{Air$\leftrightarrow$Glass (IOR: 1.504) interface and absorption within. 
 *         See \lstref{dielectric-glass}.}{bsdf_dielectric_glass}
 * }
 *
 * This plugin models an interface between two dielectric materials having mismatched 
 * indices of refraction (for instance, water and air). Exterior and interior IOR values
 * can each be independently specified, where ``exterior'' refers to the side that contains 
 * the surface normal. When no parameters are given, the plugin activates the defaults, which
 * describe a borosilicate glass BK7/air interface.
 *
 * In this model, the microscopic surface structure of the surface is assumed to be perfectly 
 * smooth, resulting in a degenerate\footnote{Meaning that for any given incoming ray of light,
 * the model always scatters into a discrete set of directions, as opposed to a continuum.} 
 * BSDF described by a Dirac delta function. For a similar model that describes a rough 
 * surface microstructure, take a look at the \pluginref{roughdielectric} plugin. 
 *
 * \begin{xml}[caption=A simple air-to-water interface, label=lst:dielectric-water]
 * <shape type="...">
 *     <bsdf type="dielectric">
 *         <float name="intIOR" value="1.33"/>
 *         <float name="extIOR" value="1.0"/>
 *     </bsdf>
 * <shape>
 * \end{xml}
 *
 * When using this model, it is crucial that the scene contains
 * meaningful and mutally compatible indices of refraction changes---see
 * \figref{glass-explanation} for a description of what this entails.
 * 
 * In many cases, we will want to additionally describe the \emph{medium} within a
 * dielectric material. This requires the use of a rendering technique that is
 * aware of media (e.g. the volumetric path tracer). An example of how one might
 * describe a slightly absorbing piece of glass is given below:
 *
 * \begin{xml}[caption=A glass material with absorption, label=lst:dielectric-glass]
 * <shape type="...">
 *     <bsdf type="dielectric">
 *         <float name="intIOR" value="1.504"/>
 *         <float name="extIOR" value="1.0"/>
 *     </bsdf>
 *     <medium type="homogeneous" name="interior">
 *			<rgb name="sigmaS" value="0, 0, 0"/>
 *			<rgb name="sigmaA" value="4, 4, 2"/>
 *     </medium>
 * <shape>
 * \end{xml}
 */
class SmoothDielectric : public BSDF {
public:
	SmoothDielectric(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = props.getFloat("intIOR", 1.5046f);
		/* Specifies the external index of refraction at the interface */
		m_extIOR = props.getFloat("extIOR", 1);

		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_specularTransmittance = new ConstantSpectrumTexture(
			props.getSpectrum("specularTransmittance", Spectrum(1.0f)));

		m_components.push_back(EDeltaReflection | EFrontSide | EBackSide);
		m_components.push_back(EDeltaTransmission | EFrontSide | EBackSide);
		m_usesRayDifferentials = false;
	}

	SmoothDielectric(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));
		m_components.push_back(EDeltaReflection | EFrontSide | EBackSide);
		m_components.push_back(EDeltaTransmission | EFrontSide | EBackSide);
		m_usesRayDifferentials = 
			m_specularReflectance->usesRayDifferentials() ||
			m_specularTransmittance->usesRayDifferentials();
	}

	virtual ~SmoothDielectric() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularReflectance") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularTransmittance") {
			m_specularTransmittance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularTransmittance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	/// Refraction in local coordinates
	inline Vector refract(const Vector &wi, Float eta, Float cosThetaT) const {
		return Vector(-eta*wi.x, -eta*wi.y, cosThetaT);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		
		if (!sampleTransmission && !sampleReflection)
			return Spectrum(0.0f);

		Float cosThetaI = Frame::cosTheta(bRec.wi),
			  etaI = m_extIOR,
			  etaT = m_intIOR;

		bool entering = cosThetaI > 0.0f;

		/* Determine the respective indices of refraction */
		if (!entering)
			std::swap(etaI, etaT);

		/* Using Snell's law, calculate the squared sine of the
		   angle between the normal and the transmitted ray */
		Float eta = etaI / etaT,
			  sinThetaTSqr = eta*eta * Frame::sinTheta2(bRec.wi);

		Float Fr, cosThetaT = 0;
		if (sinThetaTSqr >= 1.0f) {
			/* Total internal reflection */
			Fr = 1.0f;
		} else {
			cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

			/* Compute the Fresnel refletance */
			Fr = fresnelDielectric(std::abs(cosThetaI),
				cosThetaT, etaI, etaT);

			if (entering)
				cosThetaT = -cosThetaT;
		}

		/* Calculate the refracted/reflected vectors+coefficients */
		if (sampleTransmission && sampleReflection) {
			/* Importance sample according to the reflectance/transmittance */
			if (sample.x <= Fr) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);

				return m_specularReflectance->getValue(bRec.its);
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;
;
				/* Given cos(N, transmittedRay), compute the 
				   transmitted direction */
				bRec.wo = refract(bRec.wi, eta, cosThetaT);

				/* When transporting radiance, account for the solid angle
				   change at boundaries with different indices of refraction. */
				return m_specularTransmittance->getValue(bRec.its) 
					* (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1);
			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			return m_specularReflectance->getValue(bRec.its) * Fr;
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;

			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			bRec.wo = refract(bRec.wi, eta, cosThetaT);

			/* When transporting radiance, account for the solid angle
			   change at boundaries with different indices of refraction. */
			return m_specularTransmittance->getValue(bRec.its) 
				* ((1-Fr) * (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1));
		}
	}

	Spectrum eval(const BSDFQueryRecord &bRec) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		bool reflection = BSDF::cosTheta(bRec.wo) * BSDF::cosTheta(bRec.wi) > 0;
		
		if ((reflection && !sampleReflection) || 
		   (!reflection && !sampleTransmission))
			return Spectrum(0.0f);

		Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		if (reflection) {
			return m_specularReflectance->getValue(bRec.its) * fr;
		} else {
			Float etaI = m_extIOR, etaT = m_intIOR;
			bool entering = Frame::cosTheta(bRec.wi) > 0.0f;
			if (!entering)
				std::swap(etaI, etaT);

			Float factor = (bRec.quantity == ERadiance) 
				? (etaI*etaI) / (etaT*etaT) : 1.0f;

			return m_specularTransmittance->getValue(bRec.its)  * factor * (1 - fr);
		}
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		if (measure != EDiscrete || (!sampleReflection && !sampleTransmission))
			return 0.0f;

		bool reflection = Frame::costheta(bRec.wo)
			* Frame::cosTheta(bRec.wi) > 0;

		if (sampleTransmission && sampleReflection) {
			Float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
			return reflection ? Fr : (1 - Fr);
		} else {
			return sampleReflection == reflection ? 1.0f : 0.0f;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		
		if (!sampleTransmission && !sampleReflection)
			return Spectrum(0.0f);

		Float cosThetaI = Frame::cosTheta(bRec.wi),
			  etaI = m_extIOR,
			  etaT = m_intIOR;

		bool entering = cosThetaI > 0.0f;

		/* Determine the respective indices of refraction */
		if (!entering)
			std::swap(etaI, etaT);

		/* Using Snell's law, calculate the squared sine of the
		   angle between the normal and the transmitted ray */
		Float eta = etaI / etaT,
			  sinThetaTSqr = eta*eta * Frame::sinTheta2(bRec.wi);

		Float Fr, cosThetaT = 0;
		if (sinThetaTSqr >= 1.0f) {
			/* Total internal reflection */
			Fr = 1.0f;
		} else {
			cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

			/* Compute the Fresnel refletance */
			Fr = fresnelDielectric(std::abs(cosThetaI),
				cosThetaT, etaI, etaT);

			if (entering)
				cosThetaT = -cosThetaT;
		}

		/* Calculate the refracted/reflected vectors+coefficients */
		if (sampleTransmission && sampleReflection) {
			/* Importance sample according to the reflectance/transmittance */
			if (sample.x <= Fr) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);

				pdf = Fr;
				return m_specularReflectance->getValue(bRec.its) * Fr;
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;

				/* Given cos(N, transmittedRay), compute the 
				   transmitted direction */
				bRec.wo = refract(bRec.wi, eta, cosThetaT);
					
				pdf = 1-Fr;

				/* When transporting radiance, account for the solid angle
				   change at boundaries with different indices of refraction. */
				return m_specularTransmittance->getValue(bRec.its) 
					* (1-Fr) * (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1);
			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = 1;
			return m_specularReflectance->getValue(bRec.its) * Fr;
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;

			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			bRec.wo = refract(bRec.wi, eta, cosThetaT);
			pdf = 1;

			/* When transporting radiance, account for the solid angle
			   change at boundaries with different indices of refraction. */
			return m_specularTransmittance->getValue(bRec.its) 
				* ((1-Fr) * (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1));
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothDielectric[" << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR;
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
};

/* Fake glass shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least 
   something that suggests the presence of a transparent boundary */
class SmoothDielectricShader : public Shader {
public:
	SmoothDielectricShader(Renderer *renderer) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl;
		oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
};

Shader *SmoothDielectric::createShader(Renderer *renderer) const { 
	return new SmoothDielectricShader(renderer);
}

MTS_IMPLEMENT_CLASS(SmoothDielectricShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothDielectric, "Smooth dielectric BSDF");
MTS_NAMESPACE_END
