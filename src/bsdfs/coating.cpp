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
 * }
 *
 * This class implements a smooth dielectric coating in the style of 
 * Weidlich and Wilkie \cite{Weidlich2007Arbitrarily}.
 */
class SmoothCoating : public BSDF {
public:
	SmoothCoating(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = lookupIOR(props, "intIOR", "bk7");

		/* Specifies the external index of refraction at the interface */
		m_extIOR = lookupIOR(props, "extIOR", "air");
		
		/* Specifies the layer's thickness using the inverse units of sigmaT */
		m_thickness = props.getFloat("thickness", 1);

		/* Specifies the attenuation within the varnish layer */
		m_sigmaT = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaT", Spectrum(0.0f)));
	}

	SmoothCoating(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_thickness = stream->readFloat();
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_sigmaT = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");
		if (m_nested->getType() & BSDF::ETransmission)
			Log(EError, "Tried to put a smooth coating layer on top of a BSDF "
				"with a transmission component -- this is currently not allowed!");

		m_components.clear();
		for (int i=0; i<m_nested->getComponentCount(); ++i)
			m_components.push_back(m_nested->getType(i));
		m_components.push_back(EDeltaReflection | EFrontSide);
		
		m_usesRayDifferentials = m_nested->usesRayDifferentials()
			|| m_sigmaT->usesRayDifferentials();

		BSDF::configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		stream->writeFloat(m_thickness);
		manager->serialize(stream, m_nested.get());
		manager->serialize(stream, m_sigmaT.get());
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
			  etaI = m_extIOR,
			  etaT = m_intIOR;

		bool entering = cosThetaI > 0.0f;

		/* Determine the respective indices of refraction */
		if (!entering)
			std::swap(etaI, etaT);

		/* Using Snell's law, calculate the squared sine of the
		   angle between the normal and the transmitted ray */
		Float eta = etaI / etaT,
			  sinThetaTSqr = eta*eta * Frame::sinTheta2(wi);

		Float cosThetaT = 0;
		if (sinThetaTSqr >= 1.0f) {
			/* Total internal reflection */
			F = 1.0f;

			return Vector(0.0f);
		} else {
			cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

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
			Assert(bRec2.wi.z >= 0);
			Assert(bRec2.wo.z >= 0);

			if (R12 == 1 || R21 == 1) /* Total internal reflection */
				return Spectrum(0.0f);

			Spectrum result = m_nested->eval(bRec2, measure)
				* ((1-R12) * (1-R21));

			Spectrum sigmaT = m_sigmaT->getValue(bRec.its) * m_thickness;
			if (!sigmaT.isZero()) 
				result *= (-sigmaT *
					(1/std::abs(Frame::cosTheta(bRec2.wi)) +
					 1/std::abs(Frame::cosTheta(bRec2.wo)))).exp();

			if (measure == ESolidAngle)
				result *= Frame::cosTheta(bRec2.wo);
				
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

		if (measure == EDiscrete && sampleSpecular &&
			std::abs(1-dot(reflect(bRec.wi), bRec.wo)) < Epsilon) {
			return sampleNested ? fresnel(
				Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR) : 1.0f;
		} else if (sampleNested) {
			Float R12, R21;
			BSDFQueryRecord bRec2(bRec);
			bRec2.wi = -refract(bRec.wi, R12);
			bRec2.wo = -refract(bRec.wo, R21);

			if (R12 == 1 || R21 == 1) /* Total internal reflection */
				return 0.0f;

			Float pdf = m_nested->pdf(bRec2, measure);
			if (measure == ESolidAngle) {
				Float eta = m_extIOR / m_intIOR;
				pdf /= eta * eta;
			}

			return sampleSpecular ? (pdf * (1-R12)) : pdf;
		} else {
			return 0.0f;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &_sample) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		if ((!sampleNested && !sampleNested) || Frame::cosTheta(bRec.wi) < 0)
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

		Point2 sample(_sample);
		if (sampleNested && sampleNested) {
			if (sample.x <= R12) {
				bRec.sampledComponent = m_components.size()-1;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);
				pdf = R12;
				return Spectrum(R12);
			} else {
				Vector wiBackup = bRec.wi;
			//	bRec.wi = -refract(bRec.wi, eta, cosThetaT);
				bRec.wi = -refract(bRec.wi, R12);
				sample.x = (sample.x - R12) / (1 - R12);

				Spectrum result = m_nested->sample(bRec, pdf, sample);
				if (result.isZero())
					return Spectrum(0.0f);

				Spectrum sigmaT = m_sigmaT->getValue(bRec.its) * m_thickness;
				if (!sigmaT.isZero()) 
					result *= (-sigmaT *
						(1/std::abs(Frame::cosTheta(bRec.wi)) +
						 1/std::abs(Frame::cosTheta(bRec.wo)))).exp();

				Float R21, cosThetaWoPrime = Frame::cosTheta(bRec.wo);
				bRec.wi = wiBackup;
				bRec.wo = refract(-bRec.wo, R21);

				if (R21 == 1.0f) /* Total internal reflection */
					return Spectrum(0.0f);
	
				pdf *= 1 - R12;
				if (BSDF::getMeasure(bRec.sampledType) == ESolidAngle)
					pdf /= eta * eta;

				result *= (1 - R12) * (1 - R21) * cosThetaWoPrime;

				return result;
			}
		} else if (sampleSpecular) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = 1.0f;
			return Spectrum(R12);
		} else {
			// XXX not implemented
			return Spectrum(0.0f);
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

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothCoating[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  sigmaT = " << indent(m_sigmaT->toString()) << "," << endl
			<< "  thickness = " << m_thickness << "," << endl
			<< "  nested = " << indent(m_nested->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR;
	ref<Texture> m_sigmaT;
	ref<BSDF> m_nested;
	Float m_thickness;
};

MTS_IMPLEMENT_CLASS_S(SmoothCoating, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothCoating, "Smooth varnish layer");
MTS_NAMESPACE_END
