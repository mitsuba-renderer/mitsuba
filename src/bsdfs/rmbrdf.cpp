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
#include <mitsuba/render/phase.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/basicshader.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{rmbrdf}{Random medium BRDF}
 *
 * \parameters{
 *     \parameter{material}{\String}{Name of a material preset, see 
 *           \tblref{medium-coefficients}. \default{\texttt{skin1}}}
 *     \parameter{sigmaS}{\Spectrum\Or\Texture}{Specifies the scattering coefficient 
 *      of the layer. \default{based on \code{material}}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{Specifies the absorption coefficient 
 *      of the layer. \default{based on \code{material}}}
 *     \parameter{sigmaT \& albedo}{\Spectrum\Or\Texture}{
 *      Optional: Alternatively, the scattering and absorption coefficients may also be
 *      specified using the extinction coefficient \code{sigmaT} and the 
 *      single-scattering albedo. Note that only one of the parameter passing 
 *      conventions can be used at a time (i.e. use either \code{sigmaS\&sigmaA} 
 *      \emph{or} \code{sigmaT\&albedo})}
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{based on \code{material}}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{g}{\Float\Or\String}{Specifies the phase function anisotropy
 *     --- see the \pluginref{hg} plugin for details\default{0, i.e. isotropic}}
 *     \parameter{alpha}{\Float\Or\Texture}{
 *         Specifies the roughness of the unresolved surface micro-geometry.
 *         \default{0.1, i.e. the surface has a slightly rough finish}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Rendering using the whole milk material preset}{bsdf_sssbrdf}
 * }
 *
 * This plugin implements a BRDF scattering model that emulates interactions
 * with a random medium embedded inside a dielectric layer. By
 * approximating these events using a BRDF, any scattered illumination
 * is assumed to exit the material \emph{directly} at the original point of incidence.
 * To simulate actual subsurface scattering, refer to Sections~\ref{sec:media} 
 * and \ref{sec:subsurface}.
 *
 * Note that renderings with this BRDF will usually look very similar to what might
 * also be obtained using \pluginref{plastic}. The plugin's reason for existance
 * is that can be configured using parameters that are traditionally reserved
 * for participating media.
 *
 * \subsection*{Implementation details}
 * Internally, the model is implemented by instantiating a Hanrahan-Krueger
 * BSDF for single scattering in an infinitely thick layer together with 
 * an approximate multiple scattering component based on Jensen's 
 * \cite{Jensen2001Practical} integrated dipole BRDF. These are then 
 * embedded into a dielectric layer using either the \pluginref{coating} 
 * or \pluginref{roughcoating} plugins depending on whether or not
 * \code{alpha}=0.
 * This yields a very convenient parameterization of a scattering model
 * that behaves similarly to a coated diffuse material, but expressed
 * in terms of the scattering and absorption coefficients \code{sigmaS} 
 * and \code{sigmaA}.
 */
class RandomMediumBRDF : public BSDF {
public:
	RandomMediumBRDF(const Properties &props)
			: BSDF(props), m_configured(false) {

		Spectrum sigmaS, sigmaA; // ignored here
		Float eta;
		lookupMaterial(props, sigmaS, sigmaA, &eta, false);

		Float g = props.getFloat("g", 0.0f);
		Properties hgProps("hg");
		hgProps.setFloat("g", g);

		Float alpha = props.getFloat("alpha", 0.1f);

		ref<PhaseFunction> hg = static_cast<PhaseFunction *> (
			PluginManager::getInstance()->createObject(
			MTS_CLASS(PhaseFunction), hgProps));

		Properties hkProps(props);
		hkProps.setPluginName("hk");
		hkProps.setFloat("thickness", std::numeric_limits<Float>::infinity());
		m_hk = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), hkProps));
		m_hk->addChild(hg);

		Properties coatingProps(props);
		coatingProps.setPluginName(alpha == 0 ? "coating" : "roughcoating");
		if (!props.hasProperty("intIOR"))
			coatingProps.setFloat("intIOR", eta);
		if (alpha != 0)
			coatingProps.setFloat("alpha", alpha);

		m_coating = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), coatingProps));

		Properties dipoleProps(props);
		dipoleProps.setPluginName("dipolebrdf");
		m_dipole = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), dipoleProps));

		Properties mixtureProps("mixturebsdf");
		mixtureProps.setString("weights", "1.0, 1.0");
		mixtureProps.setBoolean("ensureEnergyConservation", false);
		m_mixture = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), mixtureProps));

		props.markQueried("material");
		props.markQueried("sigmaS");
		props.markQueried("sigmaA");
		props.markQueried("coating");
		props.markQueried("sigmaT");
		props.markQueried("intIOR");
		props.markQueried("extIOR");
		props.markQueried("alpha");
	}

	RandomMediumBRDF(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager), m_configured(true) {
		m_coating = static_cast<BSDF *>(manager->getInstance(stream));
		m_hk = static_cast<BSDF *>(manager->getInstance(stream));
		m_dipole = static_cast<BSDF *>(manager->getInstance(stream));
	}

	void configure() {
		if (!m_configured) {
			m_configured = true;
			m_hk->configure();
			m_dipole->configure();
			m_mixture->addChild(m_hk);
			m_mixture->addChild(m_dipole);
			m_mixture->configure();
			m_coating->addChild(m_mixture);
			m_coating->configure();

			m_components.clear();
			for (int i=0; i<m_coating->getComponentCount(); ++i) {
				unsigned int type = m_coating->getType(i);
				type &= ~BSDF::EBackSide;
				m_components.push_back(type);
			}
			BSDF::configure();
		}
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_coating->getDiffuseReflectance(its);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);
		return m_coating->eval(bRec, measure);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;
		return m_coating->pdf(bRec, measure);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		return m_coating->sample(bRec, pdf, sample);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		return m_coating->sample(bRec, sample);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
		manager->serialize(stream, m_coating.get());
		manager->serialize(stream, m_hk.get());
		manager->serialize(stream, m_dipole.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "sigmaS" || name == "sigmaA" || name == "sigmaT" || name == "albedo") {
				m_hk->addChild(name, child);
				m_dipole->addChild(name, child);
			} else if (name == "alpha") {
				m_coating->addChild(name, child);
			} else {
				BSDF::addChild(name, child);
			}
		} else {
			BSDF::addChild(name, child);
		}
	}

	Shader *createShader(Renderer *renderer) const {
		return m_coating->createShader(renderer);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RandomMediumBRDF[" << endl
   			<< "  name = \"" << m_name << "\"" << endl
   			<< "  nested = " << indent(m_coating->toString()) << endl
			<< "]";
		return oss.str();
	}
	
	MTS_DECLARE_CLASS()
private:
	ref<BSDF> m_coating;
	ref<BSDF> m_dipole;
	ref<BSDF> m_hk;
	ref<BSDF> m_mixture;
	bool m_configured;
};

MTS_IMPLEMENT_CLASS_S(RandomMediumBRDF, false, BSDF)
MTS_EXPORT_PLUGIN(RandomMediumBRDF, "Random medium BRDF");
MTS_NAMESPACE_END
