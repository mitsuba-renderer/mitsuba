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

MTS_NAMESPACE_BEGIN

/*! \plugin{coating}{Smooth dieletric coating}
 *
 * \parameters{
 *     \parameter{intIOR}{\Float}{Interior index of refraction \default{1.5046}}
 *     \parameter{extIOR}{\Float}{Exterior index of refraction \default{1.0}}
 * }
 *
 *  XXX cancel out cosine factors?
 *  XXX did I get the measure conversion terms right?
 *  XXX allow testing interface to verify delta components
 */
class SmoothVarnish : public BSDF {
public:
	SmoothVarnish(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = props.getFloat("intIOR", 1.5046f);
		/* Specifies the external index of refraction at the interface */
		m_extIOR = props.getFloat("extIOR", 1);
		/* Specifies the layer's thickness using the inverse units of sigmaT */
		m_thickness = props.getFloat("thickness", 1);
		/* Specifies the attenuation within the varnish layer */
		m_sigmaT = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaT", Spectrum(0.0f)));
	}

	SmoothVarnish(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_thickness = stream->readFloat();
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_sigmaT = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	virtual ~SmoothVarnish() {
		delete[] m_type;
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");
		if (m_nested->getType() & BSDF::ETransmission)
			Log(EError, "Tried to put a smooth varnish layer on top of a BSDF "
				"with a transmission component -- this is currently not allowed!");
		if (m_nested->getType() & BSDF::EDelta)
			Log(EError, "Tried to put a smooth varnish layer on top of a material with a "
				"Dirac delta distribution -- this is currently not allowed!");
		if (m_type)
			delete[] m_type;

		m_componentCount = 1 + m_nested->getComponentCount();
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection | EFrontSide;
		m_combinedType = m_type[0];
		for (int i=0; i<m_nested->getComponentCount(); ++i) {
			m_type[i+1] = m_nested->getType(i);
			m_combinedType |= m_type[i+1];
		}
		m_usesRayDifferentials = m_nested->usesRayDifferentials()
			|| m_sigmaT->usesRayDifferentials();
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

	/// Fully complete local coordinate refraction routine
	inline Vector refract(const Vector &wi, Float &Fr) const {
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
			Fr = 1.0f;
			return Vector(0.0f);
		} else {
			cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

			/* Compute the Fresnel refletance */
			Fr = fresnelDielectric(std::abs(cosThetaI),
				cosThetaT, etaI, etaT);

			if (entering)
				cosThetaT = -cosThetaT;
			return Vector(-eta*wi.x, -eta*wi.y, cosThetaT);
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &_sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
				&& (bRec.component == -1 || (bRec.component > 0 
				&& bRec.component < m_nested->getComponentCount() + 1));

		if ((!sampleNested && !sampleReflection) || Frame::cosTheta(bRec.wi) < 0)
			return Spectrum(0.0f);

		Float cosThetaI = Frame::cosTheta(bRec.wi),
			  etaI = m_extIOR,
			  etaT = m_intIOR;

		/* Using Snell's law, calculate the squared sine of the
		   angle between the normal and the transmitted ray */
		Float eta = etaI / etaT,
			  sinThetaTSqr = eta*eta * Frame::sinTheta2(bRec.wi);

		Float Fr, FrOut, cosThetaT = 0;
		if (sinThetaTSqr >= 1.0f) {
			/* Total internal reflection */
			Fr = 1.0f;
		} else {
			cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

			/* Compute the Fresnel refletance */
			Fr = fresnelDielectric(cosThetaI,
				cosThetaT, etaI, etaT);

			cosThetaT = -cosThetaT;
		}

		Point2 sample(_sample);
		if (sampleNested && sampleReflection) {
			if (sample.x <= Fr) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);

				pdf = Fr * std::abs(Frame::cosTheta(bRec.wo));
				return Spectrum(Fr);
			} else {
				Vector wiBackup = bRec.wi;
				bRec.wi = -refract(bRec.wi, eta, cosThetaT);
				sample.x = (sample.x - Fr) / (1 - Fr);

				Spectrum result = m_nested->sample(bRec, pdf, sample);
				if (result.isZero())
					return Spectrum(0.0f);
				bRec.sampledComponent++;

				Spectrum sigmaT = m_sigmaT->getValue(bRec.its) * m_thickness;
				if (!sigmaT.isZero()) 
					result *= (-sigmaT *
						(1/std::abs(Frame::cosTheta(bRec.wi)) +
						 1/std::abs(Frame::cosTheta(bRec.wo)))).exp();

				Float cosThetaWoPrime = Frame::cosTheta(bRec.wo);
				bRec.wi = wiBackup;
				bRec.wo = refract(-bRec.wo, FrOut);
			
				if (FrOut == 1)
					return Spectrum(0.0f);
	
				pdf *= (1 - Fr) * eta * eta;

				result *= 
					(1 - Fr) * (1 - FrOut)
					* std::abs(cosThetaWoPrime *
					/ Frame::cosTheta(bRec.wo));
				return result;

			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = std::abs(Frame::cosTheta(bRec.wo));
			return Spectrum(Fr);
		} else {
			if (Fr == 1.0f) /* Total internal reflection */
				return Spectrum(0.0f);

			Vector wiBackup = bRec.wi;
			bRec.wi = -refract(bRec.wi, eta, cosThetaT);
			sample.x = (sample.x - Fr) / (1 - Fr);

			Spectrum result = m_nested->sample(bRec, pdf, sample);
			if (result.isZero())
				return Spectrum(0.0f);
			bRec.sampledComponent++;

			Spectrum sigmaT = m_sigmaT->getValue(bRec.its) * m_thickness;
			if (!sigmaT.isZero()) 
				result *= (-sigmaT *
					(1/std::abs(Frame::cosTheta(bRec.wi)) +
					 1/std::abs(Frame::cosTheta(bRec.wo)))).exp();

			Float cosThetaWoPrime = Frame::cosTheta(bRec.wo);
			bRec.wi = wiBackup;
			bRec.wo = refract(-bRec.wo, FrOut);
			
			if (FrOut == 1)
				return Spectrum(0.0f);
	
			pdf *= (1 - Fr) * eta * eta;

			result *= 
				(1 - Fr) * (1 - FrOut)
				* std::abs(cosThetaWoPrime *
				/ Frame::cosTheta(bRec.wo));

			return result;
		}
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 0);

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleSpecular)
			return 0.0f;

		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
				&& (bRec.component == -1 || (bRec.component > 0 
				&& bRec.component < m_nested->getComponentCount() + 1));

		Float pdf = std::abs(Frame::cosTheta(bRec.wo));
		if (sampleNested)
			pdf *= fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		return pdf;
	}

	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 0);

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleSpecular)
			return Spectrum(0.0f);

		return Spectrum(fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR));
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
				&& (bRec.component == -1 || (bRec.component > 0 
				&& bRec.component < m_nested->getComponentCount() + 1));

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleNested)
			return 0.0f;

		Float T12, T21;
		Vector wiPrime = -refract(bRec.wi, T12);
		Vector woPrime = -refract(bRec.wo, T21);

		if (T12 == 1 || T21 == 1) /* Total internal reflection */
			return 0.0f;

		BSDFQueryRecord bRec2(bRec);
		if (bRec2.component != -1)
			bRec2.component++;
		bRec2.wi = wiPrime;
		bRec2.wo = woPrime;

		Float eta = m_extIOR / m_intIOR;
		return m_nested->pdf(bRec2) * T12 * eta * eta;
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
				&& (bRec.component == -1 || (bRec.component > 0 
				&& bRec.component < m_nested->getComponentCount() + 1));

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleNested)
			return Spectrum(0.0f);

		Float T12, T21;
		Vector wiPrime = -refract(bRec.wi, T12);
		Vector woPrime = -refract(bRec.wo, T21);

		if (T12 == 1 || T21 == 1) /* Total internal reflection */
			return Spectrum(0.0f);

		BSDFQueryRecord bRec2(bRec);
		if (bRec2.component != -1)
			bRec2.component++;
		bRec2.wi = wiPrime;
		bRec2.wo = woPrime;
		return m_nested->f(bRec2) * T12 * T21;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothVarnish[" << endl
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
	ref<BSDF> m_nested;
	ref<Texture> m_sigmaT;
	Float m_thickness;
};

MTS_IMPLEMENT_CLASS_S(SmoothVarnish, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothVarnish, "Smooth varnish layer");
MTS_NAMESPACE_END
