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

MTS_NAMESPACE_BEGIN

const bool importanceSampleComponents = true;

/**
 * Models an interface between two materials with non-matched indices of refraction.
 * The microscopic surface structure is assumed to be perfectly flat, resulting 
 * in a BSDF equal to a Dirac delta function. 
 *
 * The default settings are set to a borosilicate glass BK7/air interface.
 */
class Dielectric : public BSDF {
public:
	Dielectric(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = props.getFloat("intIOR", 1.5046f);
		/* Specifies the external index of refraction at the interface */
		m_extIOR = props.getFloat("extIOR", 1);
		/* Reflectance modulation term */
		m_reflectance = props.getSpectrum("specularReflectance", Spectrum(1.0f));
		/* Transmittance modulation term */
		m_transmittance = props.getSpectrum("specularTransmittance", Spectrum(1.0f));

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
		m_transmittance = Spectrum(stream);
		m_reflectance = Spectrum(stream);

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
		m_transmittance.serialize(stream);
		m_reflectance.serialize(stream);
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

	inline void reflect(const Vector &wo, Vector &wi) const {
		wi = Vector(-wo.x, -wo.y, wo.z);
	}

	Float refract(Float intIOR, Float extIOR, 
			const Vector &wi, Vector &wo, ETransportQuantity quantity) const {
		Float cosTheta1 = Frame::cosTheta(wi);
		bool entering = cosTheta1 > 0.0f;

		/* Swap the indices of refraction if the interaction starts
		   at the inside of the object */
		if (!entering)
			std::swap(intIOR, extIOR);

		Float eta = extIOR/intIOR;

		/* Using Snell's law, calculate the squared sine of the
		   angle between the normal and the transmitted ray */
		Float sinTheta2Sqr = eta*eta * Frame::sinTheta2(wi);

		if (sinTheta2Sqr > 1.0f) /* Total internal reflection! */
			return 0.0f;

		/* Compute the cosine, but guard against numerical imprecision */
		Float cosTheta2 = std::sqrt(std::max((Float) 0.0f, 1.0f - sinTheta2Sqr));
		if (entering)
			cosTheta2 = -cosTheta2;

		/* With cos(N, transmittedRay) on tap, calculating the 
		   transmission direction is straightforward */
		wo = Vector(-eta*wi.x, -eta*wi.y, cosTheta2);

		/* Finally compute transmission coefficient. When transporting
		   radiance, account for the change at boundaries with different 
		   indices of refraction. */
		if (quantity == ERadiance)
			return (extIOR*extIOR)/(intIOR*intIOR);
		else
			return 1.0f;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf=1;
		Spectrum spec = Dielectric::sample(bRec, pdf, sample);
		if (pdf == 0 || spec.isZero())
			return Spectrum(0.0f);
		return spec/pdf;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		/* Calculate the refracted/reflected vectors+coefficients */
		if (sampleTransmission && sampleReflection) {
			/* Importance sample according to the reflectance/transmittance */
			if (sample.x < (importanceSampleComponents ? fr : 0.5f)) {
				reflect(bRec.wi, bRec.wo);
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				pdf = (importanceSampleComponents ? fr : 0.5f) * std::abs(Frame::cosTheta(bRec.wo));
				/* Cancel out the cosine term */
				return m_reflectance * fr;
			} else {
				pdf = importanceSampleComponents ? (1-fr) : 0.5f;
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;

				Float result = refract(m_intIOR, m_extIOR, bRec.wi, bRec.wo, bRec.quantity);
				if (result == 0)
					return Spectrum(0.0f);
				pdf *= std::abs(Frame::cosTheta(bRec.wo));

				return m_transmittance * result * (1-fr);
			}
		} else if (sampleReflection) {
			reflect(bRec.wi, bRec.wo);
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			pdf = std::abs(Frame::cosTheta(bRec.wo));
			return m_reflectance * fr;
		} else if (sampleTransmission) {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;

			Float result = refract(m_intIOR, m_extIOR, bRec.wi, bRec.wo, bRec.quantity);
			pdf = std::abs(Frame::cosTheta(bRec.wo));
			if (result == 0)
				return Spectrum(0.0f);

			return m_transmittance * result * (1-fr);
		}
		return Spectrum(0.0f);
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		bool reflection = bRec.wo.z * bRec.wi.z > 0;

		Float result = 0.0f;
		if (sampleTransmission && sampleReflection) {
			if (!importanceSampleComponents) {
				result = 0.5f;
			} else {
				Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
				if (reflection)
					result = fr;
				else
					result = 1-fr;
			}
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
		Float intIOR = m_intIOR, extIOR = m_extIOR;
		Float fr = fresnel(Frame::cosTheta(bRec.wi), extIOR, intIOR);
		if (sampleReflection && !sampleTransmission && !reflection) 
			return Spectrum(0.0f);
		else if (!sampleReflection && sampleTransmission && reflection)
			return Spectrum(0.0f);
		if (reflection)
			return m_reflectance * fr;
		else {
			bool entering = Frame::cosTheta(bRec.wi) > 0.0f;
			if (!entering)
				std::swap(intIOR, extIOR);

			Float factor = (bRec.quantity == ERadiance) 
				? (extIOR*extIOR)/(intIOR*intIOR) : 1.0f;

			return m_transmittance  * factor * (1-fr);
		}

	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Dielectric[" << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  reflectance = " << m_reflectance.toString() << "," << endl
			<< "  transmittance = " << m_transmittance.toString() << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR;
	Spectrum m_reflectance;
	Spectrum m_transmittance;
};


MTS_IMPLEMENT_CLASS_S(Dielectric, false, BSDF)
MTS_EXPORT_PLUGIN(Dielectric, "Dielectric BSDF");
MTS_NAMESPACE_END
