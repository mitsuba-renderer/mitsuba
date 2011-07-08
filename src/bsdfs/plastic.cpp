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
#include <mitsuba/render/texture.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{plastic}{Smooth plastic material}
 * \order{7}
 * \parameters{
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{polypropylene} / 1.49}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the specular component\default{1.0}}
 *     \lastparameter{specular\showbreak Transmittance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the diffuse component\default{0.5}}
 * }
 *
 * \renderings{
 *     \medrendering{Air$\leftrightarrow$Water (IOR: 1.33) interface. 
 *         See \lstref{dielectric-water}.}{bsdf_dielectric_water}
 *     \medrendering{Air$\leftrightarrow$Diamond (IOR: 2.419)}{bsdf_dielectric_diamond}
 *     \medrendering{Air$\leftrightarrow$Glass (IOR: 1.504) interface and absorption within. 
 *         See \lstref{dielectric-glass}.}{bsdf_dielectric_glass}
 * }
 */
class SmoothPlastic : public BSDF {
public:
	SmoothPlastic(const Properties &props) : BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = lookupIOR(props, "intIOR", "polypropylene");

		/* Specifies the external index of refraction at the interface */
		m_extIOR = lookupIOR(props, "extIOR", "air");

		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_diffuseReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));

		m_components.push_back(EDeltaReflection | EFrontSide);
		m_components.push_back(EDiffuseReflection | EFrontSide);
		m_usesRayDifferentials = false;
	}

	SmoothPlastic(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_components.push_back(EDeltaReflection | EFrontSide);
		m_components.push_back(EDiffuseReflection | EFrontSide);
		m_usesRayDifferentials = 
			m_specularReflectance->usesRayDifferentials() ||
			m_diffuseReflectance->usesRayDifferentials();
	}

	virtual ~SmoothPlastic() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_diffuseReflectance.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularReflectance") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "diffuseReflectance") {
			m_diffuseReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_diffuseReflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void configure() {
		BSDF::configure();

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
		m_diffuseReflectance = ensureEnergyConservation(
			m_diffuseReflectance, "diffuseReflectance", 1.0f);
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleSpecular   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (Frame::cosTheta(bRec.wo) <= 0 || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
				  
		Float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		if (measure == EDiscrete && sampleSpecular) {
			/* Check if the provided direction pair matches an ideal
			   specular reflection; tolerate some roundoff errors */
			bool reflection = std::abs(1 - dot(reflect(bRec.wi), bRec.wo)) < Epsilon;
			if (reflection)
				return m_specularReflectance->getValue(bRec.its) * Fr;
		} else if (measure == ESolidAngle && sampleDiffuse) {
			if (sampleDiffuse)
				return m_diffuseReflectance->getValue(bRec.its) 
					* (INV_PI * Frame::cosTheta(bRec.wo) * (1-Fr));
		}

		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleSpecular   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (Frame::cosTheta(bRec.wo) <= 0 || Frame::cosTheta(bRec.wi) <= 0)
			return 0.0f;

		if (measure == EDiscrete && sampleSpecular) {
			/* Check if the provided direction pair matches an ideal
			   specular reflection; tolerate some roundoff errors */
			bool reflection = std::abs(1 - dot(reflect(bRec.wi), bRec.wo)) < Epsilon;
			if (reflection) 
				return sampleDiffuse ? 
					fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR) : 1.0f;
		} else if (measure == ESolidAngle && sampleDiffuse) {
			return Frame::cosTheta(bRec.wo) * INV_PI *
				sampleSpecular ? (1 - fresnel(Frame::cosTheta(bRec.wi),
				m_extIOR, m_intIOR)) : 1.0f;
		}

		return 0.0f;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		bool sampleSpecular   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);
		
		if ((!sampleDiffuse && !sampleSpecular) || Frame::cosTheta(bRec.wi) < 0)
			return Spectrum(0.0f);

		Float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		if (sampleDiffuse && sampleSpecular) {
			/* Importance sample wrt. the Fresnel reflectance */
			if (sample.x <= Fr) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);

				return m_specularReflectance->getValue(bRec.its);
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDiffuseReflection;
				bRec.wo = squareToHemispherePSA(Point2(
					(sample.x - Fr) / (1 - Fr),
					sample.y
				));

				return m_diffuseReflectance->getValue(bRec.its);
			}
		} else if (sampleSpecular) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			return m_specularReflectance->getValue(bRec.its) * Fr;
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDiffuseReflection;

			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);
				
			bRec.wo = squareToSphere(sample);
				
			return m_diffuseReflectance->getValue(bRec.its) * (1-Fr);
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		bool sampleSpecular   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);
		
		if ((!sampleDiffuse && !sampleSpecular) || Frame::cosTheta(bRec.wi) < 0)
			return Spectrum(0.0f);

		Float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		if (sampleDiffuse && sampleSpecular) {
			/* Importance sample wrt. the Fresnel reflectance */
			if (sample.x <= Fr) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);

				pdf = Fr;
				return m_specularReflectance->getValue(bRec.its) * Fr;
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDiffuseReflection;
				bRec.wo = squareToHemispherePSA(Point2(
					(sample.x - Fr) / (1 - Fr),
					sample.y
				));

				pdf = (1-Fr) * Frame::cosTheta(bRec.wo) * INV_PI;

				return m_diffuseReflectance->getValue(bRec.its) 
					* (INV_PI * Frame::cosTheta(bRec.wo) * (1-Fr));
			}
		} else if (sampleSpecular) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = 1;
			return m_specularReflectance->getValue(bRec.its) * Fr;
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDiffuseReflection;

			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);
				
			bRec.wo = squareToSphere(sample);
				
			pdf = Frame::cosTheta(bRec.wo) * INV_PI;

			return m_diffuseReflectance->getValue(bRec.its) 
				* (INV_PI * Frame::cosTheta(bRec.wo) * (1-Fr));
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothPlastic[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR;
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
};

MTS_IMPLEMENT_CLASS_S(SmoothPlastic, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothPlastic, "Smooth plastic BSDF");
MTS_NAMESPACE_END
