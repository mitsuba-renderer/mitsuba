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

/**
 * Rough metal BRDF model based on
 * "Microfacet Models for Refraction through Rough Surfaces"
 * by Bruce Walter, Stephen R. Marschner, Hongsong Li
 * and Kenneth E. Torrance.
 *
 * This is similar to the 'microfacet' implementation, but
 * the Fresnel term is now that of a conductor.
 */
class RoughMetal : public BSDF {
public:
	RoughMetal(const Properties &props) 
		: BSDF(props) {
		m_specularReflectance = new ConstantTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_alphaB = props.getFloat("alphaB", .1f);
		m_ior = props.getSpectrum("ior", Spectrum(0.370f));  /* Gold */
		m_k = props.getSpectrum("k", Spectrum(2.820f));

		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EGlossyReflection | EFrontSide;
		m_usesRayDifferentials = false;
	}

	RoughMetal(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaB = stream->readFloat();
		m_ior = Spectrum(stream);
		m_k = Spectrum(stream);

		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EGlossyReflection | EFrontSide;
		m_usesRayDifferentials = 
			m_specularReflectance->usesRayDifferentials();
	}

	virtual ~RoughMetal() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}

	/**
	 * Beckmann distribution function for gaussian random surfaces
	 * \param thetaM Tangent of the angle between M and N.
	 */
	Float beckmannD(const Vector &m) const {
		Float ex = Frame::tanTheta(m) / m_alphaB;
		return std::exp(-(ex*ex)) / (M_PI * m_alphaB*m_alphaB * 
			std::pow(Frame::cosTheta(m), (Float) 4.0f));
	}

	/**
	 * Sample microsurface normals according to 
	 * the Beckmann distribution
	 */
	Normal sampleBeckmannD(Point2 sample) const {
		Float thetaM = std::atan(std::sqrt(-m_alphaB*m_alphaB 
			* std::log(1.0f - sample.x)));
		Float phiM = (2.0f * M_PI) * sample.y;
		return Normal(sphericalDirection(thetaM, phiM));
	}

	/**
	 * Smith's shadow-masking function G1 for the Beckmann distribution
	 * \param m The microsurface normal
	 * \param v An arbitrary direction
	 */
	Float smithBeckmannG1(const Vector &v, const Vector &m) const {
		if (dot(v, m) * Frame::cosTheta(v) <= 0)
			return 0.0;

		const Float tanTheta = Frame::tanTheta(v);

		if (tanTheta == 0.0f)
			return 1.0f;

		const Float a = 1.0f / (m_alphaB * tanTheta);
		const Float aSqr = a * a;

		if (a >= 1.6f)
			return 1.0f;

		return (3.535f * a + 2.181f * aSqr) / 
			(1.0f + 2.276f * a + 2.577f * aSqr);
	}

	inline Vector reflect(const Vector &wi, const Normal &n) const {
		return Vector(n*(2.0f*dot(n, wi))) - wi;
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType)
			|| bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return Spectrum(0.0f);\

		Vector Hr = normalize(bRec.wi+bRec.wo);

		/* Fresnel factor */
		Spectrum F = fresnelConductor(dot(bRec.wi, Hr), m_ior, m_k);

		/* Microsurface normal distribution */
		Float D = beckmannD(Hr);
		/* Smith's shadow-masking function for the Beckmann distribution */
		Float G = smithBeckmannG1(bRec.wi, Hr) * smithBeckmannG1(bRec.wo, Hr);
		/* Calculate the total amount of specular reflection */
		Spectrum specRef = F * (D * G / 
			(4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)));

		return m_specularReflectance->getValue(bRec.its) * specRef; 
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;

		Vector Hr = normalize(bRec.wi + bRec.wo);
		/* Jacobian of the half-direction transform. */
		Float dwhr_dwo = 1.0f / (4.0f * absDot(bRec.wo, Hr));
		return beckmannD(Hr) * Frame::cosTheta(Hr) * dwhr_dwo;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (bRec.wi.z <= 0)
			return Spectrum(0.0f);

		/* Sample M, the microsurface normal */
		Normal m = sampleBeckmannD(sample);
		/* Perfect specular reflection along the microsurface normal */
		bRec.wo = reflect(bRec.wi, m);

		bRec.sampledComponent = 1;
		bRec.sampledType = EGlossyReflection;

		if (bRec.wo.z <= 0)
			return Spectrum(0.0f);

		return f(bRec) / pdf(bRec);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_specularReflectance.get());
		stream->writeFloat(m_alphaB);
		m_ior.serialize(stream);
		m_k.serialize(stream);
	}
	
	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularReflectance") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}


	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughMetal[" << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << std::endl
			<< "  ior = " << m_ior.toString() << "," << std::endl
			<< "  k = " << m_k.toString() << "," << std::endl
			<< "  alphaB = " << m_alphaB << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_specularReflectance;
	Float m_alphaB;
	Spectrum m_ior, m_k;
};

MTS_IMPLEMENT_CLASS_S(RoughMetal, false, BSDF)
MTS_EXPORT_PLUGIN(RoughMetal, "Rough metal BRDF");
MTS_NAMESPACE_END
