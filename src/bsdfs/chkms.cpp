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
#include <mitsuba/render/medium.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/basicshader.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{hk}{Coated Hanrahan-Krueger BSDF with multiple scattering}
 *
 * \parameters{
 *     \parameter{material}{\String}{Name of a material preset, see 
 *           \tblref{medium-coefficients}. \default{\texttt{skin1}}}
 *     \parameter{sigmaS}{\Spectrum\Or\Texture}{Specifies the scattering coefficient 
 *      of the scattering layer. \default{based on \code{material}}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{Specifies the absorption coefficient 
 *      of the scattering layer. \default{based on \code{material}}}
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{bk7} / 1.5046}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{\Unnamed}{\Phase}{A nested phase function instance that represents 
 *      the type of scattering interactions occurring within the layer}
 * }
 *
 */
class CoatedHanrahanKruegerMS : public BSDF {
public:
	CoatedHanrahanKruegerMS(const Properties &props)
			: BSDF(props), m_configured(false) {
		Properties hkProps(props);
		hkProps.setPluginName("hk");
		hkProps.setFloat("thickness", std::numeric_limits<Float>::infinity());
		m_hk = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), hkProps));

		Properties coatingProps(props);
		coatingProps.setPluginName("coating");
		m_coating = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), coatingProps));

		props.markQueried("material");
		props.markQueried("sigmaS");
		props.markQueried("sigmaA");
		props.markQueried("intIOR");
		props.markQueried("extIOR");
	}

	CoatedHanrahanKruegerMS(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager), m_configured(true) {
		m_hk = static_cast<BSDF *>(manager->getInstance(stream));
		m_coating = static_cast<BSDF *>(manager->getInstance(stream));
	}

	void configure() {
		if (!m_configured) {
			m_configured = true;
			m_hk->configure();
			m_coating->addChild("", m_hk);
			m_coating->configure();

			m_components.clear();
			for (size_t i=0; i<m_coating->getComponentCount(); ++i)
				m_components.push_back(m_coating->getType(i));
			BSDF::configure();
		}
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_coating->getDiffuseReflectance(its);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		return m_coating->eval(bRec, measure);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		return m_coating->pdf(bRec, measure);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		return m_coating->sample(bRec, pdf, sample);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		return m_coating->sample(bRec, sample);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_hk.get());
		manager->serialize(stream, m_coating.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
			m_hk->addChild(name, child);
		} else {
			Log(EError, "Invalid child node! (\"%s\")",
				cClass->getName().c_str());
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "CoatedHanrahanKruegerMS[" << endl
   			<< "  coating = " << indent(m_coating->toString()) << endl
			<< "]";
		return oss.str();
	}
	
	MTS_DECLARE_CLASS()
private:
	ref<BSDF> m_coating, m_hk;
	bool m_configured;
};

MTS_IMPLEMENT_CLASS_S(CoatedHanrahanKruegerMS, false, BSDF)
MTS_EXPORT_PLUGIN(CoatedHanrahanKruegerMS, "Hanrahan-Krueger BSDF");
MTS_NAMESPACE_END
