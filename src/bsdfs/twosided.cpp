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
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{twosided}{Two-sided BRDF adapter}
 * 
 * Turns a nested one-sided BRDF onto a two-sided version that
 * can be used to render meshes where the back-side is visible.
 *
 * \begin{xml}
 * <bsdf type="twosided">
 *   <bsdf type="lambertian"/>
 * </bsdf>
 * \end{xml}
 */
class TwoSidedBRDF : public BSDF {
public:
	TwoSidedBRDF(const Properties &props) 
		: BSDF(props) { }

	TwoSidedBRDF(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_nestedBRDF = static_cast<BSDF *>(manager->getInstance(stream));
		configure();
	}

	virtual ~TwoSidedBRDF() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_nestedBRDF.get());
	}

	void configure() {
		if (!m_nestedBRDF)
			Log(EError, "A nested one-sided material is required!");
		m_usesRayDifferentials = m_nestedBRDF->usesRayDifferentials();
		m_components.clear();
		for (int i=0; i<m_nestedBRDF->getComponentCount(); ++i) 
			m_components.push_back(m_nestedBRDF->getType(i) | EFrontSide | EBackSide);
		BSDF::configure();
		if (m_combinedType & BSDF::ETransmission)
			Log(EError, "Only materials without "
				"a transmission component can be nested!");
	}

	Spectrum eval(const BSDFQueryRecord &bRec) const {
		BSDFQueryRecord b(bRec);
		if (Frame::cosTheta(b.wi) < 0) {
			b.wi.z *= -1;
			b.wo.z *= -1;
		}
		return m_nestedBRDF->eval(b);
	}


	Float pdf(const BSDFQueryRecord &bRec) const {
		BSDFQueryRecord b(bRec);
		if (b.wi.z < 0) {
			b.wi.z *= -1;
			b.wo.z *= -1;
		}
		return m_nestedBRDF->pdf(b);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		bool flipped = false;
		if (Frame::cosTheta(bRec.wi) < 0) {
			bRec.wi.z *= -1;
			flipped = true;
		}
	
		Spectrum result = m_nestedBRDF->sample(bRec, sample);

		if (flipped) {
			bRec.wi.z *= -1;
			if (!result.isZero()) 
				bRec.wo.z *= -1;
		}
		return result;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		bool flipped = false;
		if (Frame::cosTheta(bRec.wi) < 0) {
			bRec.wi.z *= -1;
			flipped = true;
		}

		Spectrum result = m_nestedBRDF->sample(bRec, pdf, sample);

		if (flipped) {
			bRec.wi.z *= -1;
	
			if (!result.isZero() && pdf != 0)
				bRec.wo.z *= -1;
		}
		return result;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(BSDF::m_theClass)) {
			if (m_nestedBRDF != NULL)
				Log(EError, "Only a single nested BRDF can be added!");
			m_nestedBRDF = static_cast<BSDF *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "TwoSided[" << endl
			<< "  nestedBRDF = " << indent(m_nestedBRDF->toString()) << endl
			<< "]";
		return oss.str();
	}


	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<BSDF> m_nestedBRDF;
};


// ================ Hardware shader implementation ================ 

class TwoSidedShader : public Shader {
public:
	TwoSidedShader(Renderer *renderer, 
			const BSDF *nestedBRDF) : Shader(renderer, EBSDFShader), 
			m_nestedBRDF(nestedBRDF) {
		m_nestedBRDFShader = renderer->registerShaderForResource(nestedBRDF);
	}

	bool isComplete() const {
		return m_nestedBRDFShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedBRDFShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nestedBRDF);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) <= 0.0) {" << endl
			<< "    	wi.z *= -1; wo.z *= -1;" << endl
			<< "    }" << endl
			<< "    return " << depNames[0] << "(uv, wi, wo);" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) <= 0.0) {" << endl
			<< "    	wi.z *= -1; wo.z *= -1;" << endl
			<< "    }" << endl
			<< "    return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	const BSDF *m_nestedBRDF;
	ref<Shader> m_nestedBRDFShader;
	bool m_complete;
};

Shader *TwoSidedBRDF::createShader(Renderer *renderer) const { 
	return new TwoSidedShader(renderer, m_nestedBRDF.get());
}

MTS_IMPLEMENT_CLASS(TwoSidedShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(TwoSidedBRDF, false, BSDF)
MTS_EXPORT_PLUGIN(TwoSidedBRDF, "Two-sided BRDF adapter");
MTS_NAMESPACE_END
