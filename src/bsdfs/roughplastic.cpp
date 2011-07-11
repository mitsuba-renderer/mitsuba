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
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

#define TRANSMITTANCE_PRECOMP_NODES 200

/*!\plugin{roughplastic}{Rough plastic material}
 * \order{8}
 * \parameters{
 *     \parameter{distribution}{\String}{
 *          Specifies the type of microfacet normal distribution 
 *          used to model the surface roughness.
 *       \begin{enumerate}[(i)]
 *           \item \code{beckmann}: Physically-based distribution derived from
 *               Gaussian random surfaces. This is the default.
 *           \item \code{ggx}: New distribution proposed by
 *              Walter et al. \cite{Walter07Microfacet}, which is meant to better handle 
 *              the long tails observed in measurements of ground surfaces. 
 *              Renderings with this distribution may converge slowly.
 *           \item \code{phong}: Classical $\cos^p\theta$ distribution.
 *              Due to the underlying microfacet theory, 
 *              the use of this distribution here leads to more realistic 
 *              behavior than the separately available \pluginref{phong} plugin.
 *              \vspace{-4mm}
 *       \end{enumerate}
 *     }
 *     \parameter{alpha}{\Float}{
 *         Specifies the roughness of the unresolved surface microgeometry. 
 *         When the Beckmann distribution is used, this parameter is equal to the 
 *         \emph{root mean square} (RMS) slope of the microfacets. 
 *         \default{0.1}. 
 *     }
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{bk7} / 1.5046}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the specular reflectance component\default{1.0}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the diffuse reflectance component\default{0.5}}
 * }
 * \renderings{
 *     \rendering{Beckmann, $\alpha=0.1$}{bsdf_roughplastic_ggx}
 *     \rendering{GGX, $\alpha=0.3$}{bsdf_roughplastic_ggx}
 * }
 *
 */
class RoughPlastic : public BSDF {
public:
	RoughPlastic(const Properties &props) : BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_diffuseReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));

		/* Specifies the internal index of refraction at the interface */
		m_intIOR = lookupIOR(props, "intIOR", "bk7");

		/* Specifies the external index of refraction at the interface */
		m_extIOR = lookupIOR(props, "extIOR", "air");

		if (m_intIOR < 0 || m_extIOR < 0 || m_intIOR == m_extIOR)
			Log(EError, "The interior and exterior indices of "
				"refraction must be positive and differ!");

		m_distribution = MicrofacetDistribution(
			props.getString("distribution", "beckmann")
		);

		if (m_distribution.isAnisotropic())
			Log(EError, "The 'roughplastic' plugin currently does not support "
				"anisotropic microfacet distributions!");

		m_alpha = m_distribution.transformRoughness(
			props.getFloat("alpha", 0.1f));

		m_usesRayDifferentials = false;
	}

	RoughPlastic(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_distribution = MicrofacetDistribution(
			(MicrofacetDistribution::EType) stream->readUInt()
		);
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_roughTransmittance = static_cast<CubicSpline *>(manager->getInstance(stream));
		m_alpha = stream->readFloat();
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();

		m_usesRayDifferentials = 
			m_specularReflectance->usesRayDifferentials() ||
			m_diffuseReflectance->usesRayDifferentials();

		configure();
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide);
		m_components.push_back(EDiffuseReflection | EFrontSide);

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
		m_diffuseReflectance = ensureEnergyConservation(
			m_diffuseReflectance, "diffuseReflectance", 1.0f);

		/* Compute weights that further steer samples towards
		   the specular or diffuse components */
		Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
			  sAvg = m_specularReflectance->getAverage().getLuminance();
		m_specularSamplingWeight = sAvg / (dAvg + sAvg);

		/* Precompute the rough transmittance through the interface */
		m_roughTransmittance = m_distribution.computeRoughTransmittance(
				m_extIOR, m_intIOR, m_alpha, TRANSMITTANCE_PRECOMP_NODES);

		BSDF::configure();
	}

	virtual ~RoughPlastic() { }

	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleSpecular = (bRec.typeMask & EGlossyReflection) &&
			(bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection) &&
			(bRec.component == -1 || bRec.component == 1);
			
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			(!sampleSpecular && !sampleDiffuse))
			return Spectrum(0.0f);
	
		Spectrum result(0.0f);
		if (sampleSpecular) {
			/* Calculate the reflection half-vector */
			const Vector H = normalize(bRec.wo+bRec.wi);

			/* Evaluate the microsurface normal distribution */
			const Float D = m_distribution.eval(H, m_alpha);

			/* Fresnel term */
			const Float F = fresnel(dot(bRec.wi, H), m_extIOR, m_intIOR);

			/* Smith's shadow-masking function */
			const Float G = m_distribution.G(bRec.wi, bRec.wo, H, m_alpha);

			/* Calculate the specular reflection component */
			Float value = F * D * G / 
				(4.0f * Frame::cosTheta(bRec.wi));

			result += m_specularReflectance->getValue(bRec.its) * value; 
		}

		if (sampleDiffuse) 
			result += m_diffuseReflectance->getValue(bRec.its) * (INV_PI
				* m_roughTransmittance->eval(Frame::cosTheta(bRec.wi))
				* Frame::cosTheta(bRec.wo));

		return result;
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleSpecular = (bRec.typeMask & EGlossyReflection) &&
			(bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection) &&
			(bRec.component == -1 || bRec.component == 1);

		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			(!sampleSpecular && !sampleDiffuse))
			return 0.0f;

		/* Calculate the reflection half-vector */
		const Vector H = normalize(bRec.wo+bRec.wi);

		Float probDiffuse, probSpecular;
		if (sampleSpecular && sampleDiffuse) {
			/* Find the probability of sampling the specular component */
			probSpecular = 1-m_roughTransmittance->eval(Frame::cosTheta(bRec.wi));

			/* Reallocate samples */
			probSpecular = (probSpecular*m_specularSamplingWeight) /
				(probSpecular*m_specularSamplingWeight + 
				(1-probSpecular) * (1-m_specularSamplingWeight));

			probDiffuse = 1 - probSpecular;
		} else {
			probDiffuse = probSpecular = 1.0f;
		}

		Float result = 0.0f;
		if (sampleSpecular) {
			/* Jacobian of the half-direction transform */
			const Float dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));

			/* Evaluate the microsurface normal distribution */
			const Float prob = m_distribution.pdf(H, m_alpha);

			result = prob * dwh_dwo * probSpecular;
		}

		if (sampleDiffuse) 
			result += Frame::cosTheta(bRec.wo) * INV_PI * probDiffuse;

		return result;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, Float &_pdf, const Point2 &_sample) const {
		bool sampleSpecular = (bRec.typeMask & EGlossyReflection) &&
			(bRec.component == -1 || bRec.component == 0);
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection) &&
			(bRec.component == -1 || bRec.component == 1);
		
		if (Frame::cosTheta(bRec.wi) <= 0 || (!sampleSpecular && !sampleDiffuse))
			return Spectrum(0.0f);


		bool choseSpecular = sampleSpecular;
		Point2 sample(_sample);

		Float probSpecular, probDiffuse;
		if (sampleSpecular && sampleDiffuse) {
			/* Find the probability of sampling the diffuse component */
			probSpecular = 1 - m_roughTransmittance->eval(Frame::cosTheta(bRec.wi));

			/* Reallocate samples */
			probSpecular = (probSpecular*m_specularSamplingWeight) /
				(probSpecular*m_specularSamplingWeight + 
				(1-probSpecular) * (1-m_specularSamplingWeight));

			probDiffuse = 1 - probSpecular;

			if (sample.x <= probSpecular) {
				sample.x /= probSpecular;
			} else {
				sample.x = (sample.x - probSpecular) / (1 - probSpecular);
				choseSpecular = false;
			}
		}

		if (choseSpecular) {
			/* Perfect specular reflection based on the microsurface normal */
			Normal m = m_distribution.sample(sample, m_alpha);
			bRec.wo = reflect(bRec.wi, m);
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDiffuseReflection;
			bRec.wo = squareToHemispherePSA(sample);
			m = normalize(bRec.wo+bRec.wi);
		}

		/* Guard against numerical imprecisions */
		_pdf = pdf(bRec, ESolidAngle);

		if (_pdf == 0) 
			return Spectrum(0.0f);
		else
			return eval(bRec, ESolidAngle);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf;
		Spectrum result = RoughPlastic::sample(bRec, pdf, sample);

		if (result.isZero())
			return Spectrum(0.0f);
		else
			return result / pdf;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t) m_distribution.getType());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_roughTransmittance.get());
		stream->writeFloat(m_alpha);
		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
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

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughPlastic[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  distribution = " << m_distribution.toString() << "," << endl
			<< "  alpha = " << m_alpha << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
			<< "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
			<< "  intIOR = " << m_intIOR << "," << endl
			<< "  extIOR = " << m_extIOR << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution m_distribution;
	ref<CubicSpline> m_roughTransmittance;
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
	Float m_alpha, m_intIOR, m_extIOR;
	Float m_specularSamplingWeight;
};

/* Fake plastic shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least 
   something that suggests the presence of a translucent boundary */
class RoughPlasticShader : public Shader {
public:
	RoughPlasticShader(Renderer *renderer) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
};

Shader *RoughPlastic::createShader(Renderer *renderer) const { 
	return new RoughPlasticShader(renderer);
}

MTS_IMPLEMENT_CLASS(RoughPlasticShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughPlastic, false, BSDF)
MTS_EXPORT_PLUGIN(RoughPlastic, "Rough plastic BSDF");
MTS_NAMESPACE_END
