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

/*! \plugin{dielectric}{Ideal dielectric/glass material}
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
 * This plugin models an interface between two materials having mismatched 
 * indices of refraction (e.g. a boundary between water and air).
 * The microscopic surface structure of the surface is assumed to be perfectly 
 * smooth, resulting in a degenerate BSDF described by a Dirac delta function.
 * For non-smooth interfaces, please take a look at the \pluginref{roughglass}
 * model. When using this plugin, it is crucial that the scene contains
 * meaningful and mutally compatible index of refraction change -- see
 * \figref{glass-explanation} for an example.
 * 
 * The default settings of this plugin are set to a borosilicate glass BK7/air 
 * interface.
 */
class Dielectric : public BSDF {
public:
	Dielectric(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = props.getFloat("intIOR", 1.5046f);
		/* Specifies the external index of refraction at the interface */
		m_extIOR = props.getFloat("extIOR", 1);

		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_specularTransmittance = new ConstantSpectrumTexture(
			props.getSpectrum("specularTransmittance", Spectrum(1.0f)));

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection | EFrontSide | EBackSide;
		m_type[1] = EDeltaTransmission | EFrontSide | EBackSide;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	Dielectric(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection | EFrontSide | EBackSide;
		m_type[1] = EDeltaTransmission | EFrontSide | EBackSide;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	virtual ~Dielectric() {
		delete[] m_type;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		return 0.0f;
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

				return m_specularReflectance->getValue(bRec.its) 
					/ std::abs(Frame::cosTheta(bRec.wo));
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;

				/* Given cos(N, transmittedRay), compute the 
				   transmitted direction */
				bRec.wo = refract(bRec.wi, eta, cosThetaT);

				/* When transporting radiance, account for the solid angle
				   change at boundaries with different indices of refraction. */
				return m_specularTransmittance->getValue(bRec.its) 
					* (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1)
					/ std::abs(Frame::cosTheta(bRec.wo));
			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			return m_specularReflectance->getValue(bRec.its) * (Fr
				/ std::abs(Frame::cosTheta(bRec.wo)));
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;

			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			bRec.wo = refract(bRec.wi, eta, cosThetaT);

			/* When transporting radiance, account for the solid angle
			   change at boundaries with different indices of refraction. */
			return m_specularTransmittance->getValue(bRec.its) 
				* ((1-Fr) * (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1))
				/ std::abs(Frame::cosTheta(bRec.wo));
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

				pdf = Fr * std::abs(Frame::cosTheta(bRec.wo));
				return m_specularReflectance->getValue(bRec.its) * Fr;
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;

				/* Given cos(N, transmittedRay), compute the 
				   transmitted direction */
				bRec.wo = refract(bRec.wi, eta, cosThetaT);
					
				pdf = (1-Fr) * std::abs(Frame::cosTheta(bRec.wo));

				/* When transporting radiance, account for the solid angle
				   change at boundaries with different indices of refraction. */
				return m_specularTransmittance->getValue(bRec.its) 
					* (1-Fr) * (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1);
			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = std::abs(Frame::cosTheta(bRec.wo));
			return m_specularReflectance->getValue(bRec.its) * Fr;
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;

			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			bRec.wo = refract(bRec.wi, eta, cosThetaT);
			pdf = std::abs(Frame::cosTheta(bRec.wo));

			/* When transporting radiance, account for the solid angle
			   change at boundaries with different indices of refraction. */
			return m_specularTransmittance->getValue(bRec.its) 
				* ((1-Fr) * (bRec.quantity == ERadiance ? (eta*eta) : (Float) 1));
		}
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		bool reflection = bRec.wo.z * bRec.wi.z > 0;

		Float result = 0.0f;
		if (sampleTransmission && sampleReflection) {
			Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
			result = reflection ? fr : (1-fr);
		} else if (sampleReflection) {
			result = reflection ? 1.0f : 0.0f;
		} else if (sampleTransmission) {
			result = reflection ? 0.0f : 1.0f;
		}
		return result * std::abs(Frame::cosTheta(bRec.wo));
	}

	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		bool reflection = bRec.wo.z * bRec.wi.z > 0;
		Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
		
		if (sampleReflection && !sampleTransmission && !reflection) 
			return Spectrum(0.0f);
		else if (!sampleReflection && sampleTransmission && reflection)
			return Spectrum(0.0f);

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

	std::string toString() const {
		std::ostringstream oss;
		oss << "Dielectric[" << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR;
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
};


MTS_IMPLEMENT_CLASS_S(Dielectric, false, BSDF)
MTS_EXPORT_PLUGIN(Dielectric, "Dielectric BSDF");
MTS_NAMESPACE_END
