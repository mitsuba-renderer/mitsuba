/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>

MTS_NAMESPACE_BEGIN

/**
 * Rough glass BSDF model based on
 * "Microfacet Models for Refraction through Rough Surfaces"
 * by Bruce Walter, Stephen R. Marschner, Hongsong Li
 * and Kenneth E. Torrance
 */
class RoughGlass : public BSDF {
public:
	RoughGlass(const Properties &props) 
		: BSDF(props) {
		m_specularReflectance = props.getSpectrum("specularReflectance", 
			Spectrum(1.0f));
		m_specularTransmittance = props.getSpectrum("specularTransmittance", 
			Spectrum(1.0f));
		m_alphaB = props.getFloat("alphaB", .1f);
		m_intIOR = props.getFloat("intIOR", 1.5);
		m_extIOR = props.getFloat("extIOR", 1.0f);

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EGlossyReflection;
		m_type[1] = EGlossyTransmission;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	RoughGlass(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_specularReflectance = Spectrum(stream);
		m_specularTransmittance = Spectrum(stream);
		m_alphaB = stream->readFloat();
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EGlossyReflection;
		m_type[1] = EGlossyTransmission;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	virtual ~RoughGlass() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}

	/**
	 * Beckmann distribution function for gaussian random surfaces
	 * @param thetaM Tangent of the angle between M and N.
	 */
	inline Float beckmannD(const Vector &m, Float alphaB) const {
		Float ex = Frame::tanTheta(m) / alphaB;
		Float value = std::exp(-(ex*ex)) / (M_PI * alphaB*alphaB * 
			std::pow(Frame::cosTheta(m), (Float) 4.0f));
		if (value < Epsilon)
			return 0;
		return value;
	}

	/**
	 * Sample microsurface normals according to 
	 * the Beckmann distribution
	 */
	Normal sampleBeckmannD(Point2 sample, Float alphaB) const {
		Float thetaM = std::atan(std::sqrt(-alphaB*alphaB 
			* std::log(1.0f - sample.x)));
		Float phiM = (2.0f * M_PI) * sample.y;
		return Normal(sphericalDirection(thetaM, phiM));
	}

	/**
	 * Smith's shadow-masking function G1 for the Beckmann distribution
	 * @param m The microsurface normal
	 * @param v An arbitrary direction
	 */
	Float smithBeckmannG1(const Vector &v, const Vector &m) const {
		if (dot(v, m)*Frame::cosTheta(v) <= 0)
			return 0.0;

		const Float tanTheta = std::abs(Frame::tanTheta(v));

		if (tanTheta == 0.0f)
			return 1.0f;

		const Float a = 1.0f / (m_alphaB * tanTheta);
		const Float aSqr = a * a;

		if (a >= 1.6f)
			return 1.0f;

		return (3.535f * a + 2.181f * aSqr)/(1.0f + 2.276f * a + 2.577f * aSqr);
	}

	inline Vector reflect(const Vector &wi, const Normal &n) const {
		return Vector(n*(2.0f*dot(n, wi))) - wi;
	}

	inline Float signum(Float value) const {
		return (value < 0) ? -1.0f : 1.0f;
	}

	Float refract(const Vector &wi, Vector &wo, ETransportQuantity quantity) const {
		Float cosTheta1 = Frame::cosTheta(wi);
		Float intIOR = m_intIOR, extIOR = m_extIOR;
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

		/* Use the sin^2+cos^2=1 identity - max() guards against
		   numerical imprecision*/
		Float cosTheta2 = std::sqrt(std::max((Float) 0.0f, 1.0f - sinTheta2Sqr));
		if (entering)
			cosTheta2 = -cosTheta2;

		/* Having cos(N, transmittedRay), calculating the actual
		   direction becomes easy. */
		wo = Vector(-eta*wi.x, -eta*wi.y, cosTheta2);

		/* Finally compute transmission coefficient. When transporting
		   radiance, account for the change at boundaries with different 
		   indices of refraction. */

		if (quantity == ERadiance)
			return (extIOR*extIOR)/(intIOR*intIOR) * 
				(1.0f - fresnel(Frame::cosTheta(wi), extIOR, intIOR));
		else
			return 1.0f - fresnel(Frame::cosTheta(wi), extIOR, intIOR);
	}

	inline Spectrum fReflection(const BSDFQueryRecord &bRec) const {
		Float intIOR = m_intIOR, extIOR = m_extIOR;

		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(intIOR, extIOR);

		/* Calculate the reflection half-vector (and possibly flip it
		   so that it lies inside the hemisphere around the normal) */
		Vector Hr = normalize(bRec.wo+bRec.wi) 
			* signum(Frame::cosTheta(bRec.wo));

		/* Fresnel factor */
		Float F = fresnel(dot(bRec.wi, Hr), m_extIOR, m_intIOR);

		/* Microsurface normal distribution */
		Float D = beckmannD(Hr, m_alphaB);

		/* Smith's shadow-masking function for the Beckmann distribution */
		Float G = smithBeckmannG1(bRec.wi, Hr) * smithBeckmannG1(bRec.wo, Hr);

		/* Calculate the total amount of reflection */
		Float value = F * D * G / 
			(4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));
		
		return m_specularReflectance * value; 
	}

	Spectrum fTransmission(const BSDFQueryRecord &bRec) const {
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return Spectrum(0.0f);

		Float etaI = m_extIOR, etaO = m_intIOR;
		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(etaI, etaO);

		/* Calculate the transmission half-vector */
		Vector Ht = -normalize(bRec.wi*etaI+bRec.wo*etaO);

		/* Fresnel factor */
		Float F = 1.0f - fresnel(dot(bRec.wi, Ht), m_extIOR, m_intIOR);

		/* Microsurface normal distribution */
		Float D = beckmannD(Ht, m_alphaB);

		/* Smith's shadow-masking function for the Beckmann distribution */
		Float G;
		if (Ht.z > 0) 
			G = smithBeckmannG1(bRec.wi, Ht) * smithBeckmannG1(bRec.wo, Ht);
		else
			G = smithBeckmannG1(bRec.wi, -Ht) * smithBeckmannG1(bRec.wo, -Ht);

		/* Calculate the total amount of transmission */
		Float value = F * D * G * std::abs((dot(bRec.wi, Ht)*dot(bRec.wo, Ht)) / 
			(Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)));

		if (bRec.quantity == ERadiance)
			value *= (etaI*etaI)/(etaO*etaO);

		/* Half-angle Jacobian for refraction */
		Float sqrtDenom = etaI * dot(bRec.wi, Ht) + etaO * dot(bRec.wo, Ht);
		value *= (etaO*etaO)/(sqrtDenom*sqrtDenom);

		return m_specularTransmittance * value;
	}


	inline Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (!bRec.typeMask & m_combinedType)
			return result;

		bool hasReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		if (hasReflection)
			result += fReflection(bRec);
		if (hasTransmission)
			result += fTransmission(bRec);

		return result;
	}
	
	inline Float pdfReflection(const BSDFQueryRecord &bRec, Float alphaB) const {
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) < 0)
			return 0.0f;

		Vector Hr = normalize(bRec.wo+bRec.wi) 
			* signum(Frame::cosTheta(bRec.wi));

		/* Jacobian of the half-direction transform */
		Float dwhr_dwo = 1.0f / (4.0f * absDot(bRec.wo, Hr));

		return beckmannD(Hr, alphaB) * std::abs(Frame::cosTheta(Hr)) * dwhr_dwo;
	}

	inline Float pdfTransmission(const BSDFQueryRecord &bRec, Float alphaB) const {
		Float etaI = m_extIOR, etaO = m_intIOR;

		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return 0.0f;

		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(etaI, etaO);

		Vector Ht = -normalize(bRec.wi*etaI+bRec.wo*etaO);

		/* Jacobian of the half-direction transform. */
		Float sqrtDenom = etaI * dot(bRec.wi, Ht) + etaO * dot(bRec.wo, Ht);
		Float dwht_dwo = (etaO*etaO * absDot(bRec.wo, Ht)) / (sqrtDenom*sqrtDenom);

		return beckmannD(Ht, alphaB) * std::abs(Frame::cosTheta(Ht)) * dwht_dwo;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool hasReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alphaB. This in practice limits the weights to 
		   values <= 4. See also \ref sample() */
		Float alphaB = m_alphaB * (1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		if (hasReflection && hasTransmission) {
			Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
			fr = std::min(std::max(fr, (Float) 0.05f), (Float) 0.95f);
			return fr * pdfReflection(bRec, alphaB) +
				   (1-fr) * pdfTransmission(bRec, alphaB);
		} else if (hasReflection) {
			return pdfReflection(bRec, alphaB);
		} else if (hasTransmission) {
			return pdfTransmission(bRec, alphaB);
		}

		return 0.0f;
	}

	inline Spectrum sampleReflection(BSDFQueryRecord &bRec, Float alphaB) const {
		/* Sample M, the microsurface normal */
		Normal m = sampleBeckmannD(bRec.sample, alphaB);

		/* Perfect specular reflection along the microsurface normal */
		bRec.wo = reflect(bRec.wi, m);

		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Sidedness test */
		if (Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		return f(bRec) / pdf(bRec);
	}

	inline Spectrum sampleTransmission(BSDFQueryRecord &bRec, Float alphaB) const {
		/* Sample M, the microsurface normal */
		Frame mFrame(sampleBeckmannD(bRec.sample, alphaB));

		/* Perfect specular reflection along the microsurface normal */
		if (refract(mFrame.toLocal(bRec.wi), bRec.wo, bRec.quantity) == 0)
			return Spectrum(0.0f);

		bRec.wo = mFrame.toWorld(bRec.wo);

		bRec.sampledComponent = 1;
		bRec.sampledType = EGlossyTransmission;

		if (Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo) >= 0)
			return Spectrum(0.0f);
			
		return f(bRec) / pdf(bRec);
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		bool hasReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alphaB. This in practice limits the weights to 
		   values <= 4. The change is of course also accounted for 
		   in \ref pdf(), hence no error is introduced. */
		Float alphaB = m_alphaB * (1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		if (hasReflection && hasTransmission) {
			Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
			fr = std::min(std::max(fr, (Float) 0.05f), (Float) 0.95f);
			if (bRec.sample.x < fr) {
				bRec.sample.x /= fr;
				return sampleReflection(bRec, alphaB);
			} else {
				bRec.sample.x = (bRec.sample.x - fr) / (1-fr);
				return sampleTransmission(bRec, alphaB);
			}
		} else if (hasReflection) {
			return sampleReflection(bRec, alphaB);
		} else if (hasTransmission) {
			return sampleTransmission(bRec, alphaB);
		}

		return Spectrum(0.0f);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		m_specularReflectance.serialize(stream);
		m_specularTransmittance.serialize(stream);
		stream->writeFloat(m_alphaB);
		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughGlass[" << std::endl
			<< "  specularReflectance=" << m_specularReflectance.toString() << "," << std::endl
			<< "  specularTransmittance=" << m_specularTransmittance.toString() << "," << std::endl
			<< "  intIOR=" << m_intIOR << "," << std::endl
			<< "  extIOR=" << m_extIOR << "," << std::endl
			<< "  alphaB=" << m_alphaB << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_specularReflectance;
	Spectrum m_specularTransmittance;
	Float m_alphaB, m_intIOR, m_extIOR;
};

MTS_IMPLEMENT_CLASS_S(RoughGlass, false, BSDF)
MTS_EXPORT_PLUGIN(RoughGlass, "Rough glass BSDF");
MTS_NAMESPACE_END
