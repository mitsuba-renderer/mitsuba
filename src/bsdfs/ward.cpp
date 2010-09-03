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
#include <mitsuba/render/texture.h>

MTS_NAMESPACE_BEGIN

/**
 * Anisotropic Ward BRDF model based on
 *   "Measuring and Modeling Anisotropic Reflection" by 
 *   Gregory J. Ward, SIGGRAPH 1992
 * and
 *   "Notes on the Ward BRDF" by Bruce Walter, Technical Report 
 *   PCG-05-06, Cornell University
 */
class Ward : public BSDF {
public:
	Ward(const Properties &props) 
		: BSDF(props) {
		m_diffuseReflectance = new ConstantTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
		m_specularReflectance = new ConstantTexture(
			props.getSpectrum("specularReflectance", Spectrum(0.2f)));
		
		m_kd = props.getFloat("diffuseAmount", 1.0f);
		m_ks = props.getFloat("specularAmount", 1.0f);

		bool verifyEnergyConservation = props.getBoolean("verifyEnergyConservation", true);

		if (verifyEnergyConservation && (m_kd * m_diffuseReflectance->getMaximum().max() 
				+ m_ks * m_specularReflectance->getMaximum().max() > 1.0f)) {
			Log(EWarn, "%s: Energy conservation is potentially violated!", props.getID().c_str());
			Log(EWarn, "Max. diffuse reflectance = %f * %f = %f", m_kd, m_diffuseReflectance->getMaximum().max(), m_kd*m_diffuseReflectance->getMaximum().max());
			Log(EWarn, "Max. specular reflectance = %f * %f = %f", m_ks, m_specularReflectance->getMaximum().max(), m_ks*m_specularReflectance->getMaximum().max());
			Float normalization = 1/(m_kd * m_diffuseReflectance->getMaximum().max() + m_ks * m_specularReflectance->getMaximum().max());
			Log(EWarn, "Reducing the albedo to %.1f%% of the original value to be on the safe side. "
				"Specify verifyEnergyConservation=false to prevent this.", normalization * 100);
			m_kd *= normalization; m_ks *= normalization;
		}

		Float avgDiffReflectance = m_diffuseReflectance->getAverage().average() * m_kd;
		Float avgSpecularReflectance = m_specularReflectance->getAverage().average() * m_ks;

		m_specularSamplingWeight = props.getFloat("specularSamplingWeight", 
			avgSpecularReflectance / (avgDiffReflectance + avgSpecularReflectance));
		m_diffuseSamplingWeight = 1.0f - m_specularSamplingWeight;

		m_alphaX = props.getFloat("alphaX", .1f);
		m_alphaY = props.getFloat("alphaY", .1f);
		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	Ward(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaX = stream->readFloat();
		m_alphaY = stream->readFloat();
		m_kd = stream->readFloat();
		m_ks = stream->readFloat();
		m_specularSamplingWeight = stream->readFloat();
		m_diffuseSamplingWeight = stream->readFloat();

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = 
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials();
	}

	virtual ~Ward() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_diffuseReflectance->getValue(its) * m_kd;
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return result;

		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (hasGlossy) {
			Vector H = bRec.wi+bRec.wo;
			Float factor1 = 1.0f / (4.0f * M_PI * m_alphaX * m_alphaY * 
				std::sqrt(Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo)));
			Float factor2 = H.x / m_alphaX, factor3 = H.y / m_alphaY;
			Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
			Float specRef = factor1 * std::exp(exponent) * m_ks;
			result += m_specularReflectance->getValue(bRec.its) * specRef;
		}

		if (hasDiffuse) 
			result += m_diffuseReflectance->getValue(bRec.its) * (INV_PI * m_kd);

		return result;
	}
		
	inline Float pdfSpec(const BSDFQueryRecord &bRec) const {
		Vector H = normalize(bRec.wi+bRec.wo);
		Float factor1 = 1.0f / (4.0f * M_PI * m_alphaX * m_alphaY * 
			dot(H, bRec.wi) * std::pow(Frame::cosTheta(H), 3));
		Float factor2 = H.x / m_alphaX, factor3 = H.y / m_alphaY;
		Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
		Float specPdf = factor1 * std::exp(exponent);
		return specPdf;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;

		if (hasDiffuse && hasGlossy) {
			return m_specularSamplingWeight * pdfSpec(bRec) +
				   m_diffuseSamplingWeight * pdfLambertian(bRec);
		} else if (hasDiffuse) {
			return pdfLambertian(bRec);
		} else if (hasGlossy) {
			return pdfSpec(bRec);
		}

		return 0.0f;
	}

	inline Spectrum sampleSpecular(BSDFQueryRecord &bRec) const {
		Float phiH = std::atan(m_alphaY/m_alphaX 
			* std::tan(2.0f * M_PI * bRec.sample.y));
		if (bRec.sample.y > 0.5f)
			phiH += M_PI;
		Float cosPhiH = std::cos(phiH);
		Float sinPhiH = std::sqrt(std::max((Float) 0.0f, 
			1.0f-cosPhiH*cosPhiH));

		Float thetaH = std::atan(std::sqrt(std::max((Float) 0.0f, 
			-std::log(bRec.sample.x) / (
				(cosPhiH*cosPhiH)/(m_alphaX*m_alphaX) +
				(sinPhiH*sinPhiH)/(m_alphaY*m_alphaY)
		))));
		Vector H = sphericalDirection(thetaH, phiH);
		bRec.wo = H * (2.0f * dot(bRec.wi, H)) - bRec.wi;

		bRec.sampledComponent = 1;
		bRec.sampledType = EGlossyReflection;

		if (Frame::cosTheta(bRec.wo) <= 0.0f)
			return Spectrum(0.0f);

		return f(bRec) / pdf(bRec);
	}

	inline Float pdfLambertian(const BSDFQueryRecord &bRec) const {
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	inline Spectrum sampleLambertian(BSDFQueryRecord &bRec) const {
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		return f(bRec) / pdf(bRec);
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0)
			return Spectrum(0.0f);

		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (sampleDiffuse && sampleGlossy) {
			if (bRec.sample.x <= m_specularSamplingWeight) {
				bRec.sample.x = bRec.sample.x / m_specularSamplingWeight;
				return sampleSpecular(bRec);
			} else {
				bRec.sample.x = (bRec.sample.x - m_specularSamplingWeight)
					/ m_diffuseSamplingWeight;
				return sampleLambertian(bRec);
			}
		} else if (sampleDiffuse) {
			return sampleLambertian(bRec);
		} else if (sampleGlossy) {
			return sampleSpecular(bRec);
		}

		return Spectrum(0.0f);
	}
	
	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "diffuseColor") {
			m_diffuseReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_diffuseReflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "specularColor") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}


	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_specularReflectance.get());
		stream->writeFloat(m_alphaX);
		stream->writeFloat(m_alphaY);
		stream->writeFloat(m_kd);
		stream->writeFloat(m_ks);
		stream->writeFloat(m_specularSamplingWeight);
		stream->writeFloat(m_diffuseSamplingWeight);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Ward[" << endl
			<< "  diffuseColor = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularColor = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  alphaX = " << m_alphaX << "," << endl
			<< "  alphaY = " << m_alphaY << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_diffuseReflectance;
	ref<const Texture> m_specularReflectance;
	Float m_alphaX, m_alphaY;
	Float m_kd, m_ks;
	Float m_specularSamplingWeight;
	Float m_diffuseSamplingWeight;
};

MTS_IMPLEMENT_CLASS_S(Ward, false, BSDF);
MTS_EXPORT_PLUGIN(Ward, "Anisotropic Ward BRDF");
MTS_NAMESPACE_END
