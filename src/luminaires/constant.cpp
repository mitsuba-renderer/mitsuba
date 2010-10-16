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

#include <mitsuba/render/scene.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Constant background light source
 */
class ConstantLuminaire : public Luminaire {
public:
	ConstantLuminaire(const Properties &props) : Luminaire(props) {
		m_intensity = props.getSpectrum("intensity", Spectrum(1.0f));
		m_type = EOnSurface | EDiffuseDirection;
	}

	ConstantLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_intensity = Spectrum(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);

		m_intensity.serialize(stream);
	}

	void preprocess(const Scene *scene) {
		/* Get the scene's bounding sphere and slightly enlarge it */
		m_bsphere = scene->getBSphere();
		m_bsphere.radius *= 1.01f;
		m_surfaceArea = m_bsphere.radius * m_bsphere.radius * M_PI;
	}

	Spectrum getPower() const {
		return m_intensity * m_surfaceArea * M_PI;
	}

	Spectrum Le(const Ray &ray) const {
		return m_intensity;
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		return m_intensity;
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		lRec.d = squareToSphere(sample);
		lRec.sRec.p = p - lRec.d * (2 * m_bsphere.radius);
		lRec.pdf = 1.0f / (4*M_PI);
		lRec.Le = m_intensity;
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec) const {
		return 1.0f / (4*M_PI);
	}

	/* Sampling routine for surfaces - just do BSDF sampling */
	void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		const BSDF *bsdf = its.shape->getBSDF();
		BSDFQueryRecord bRec(its, sample);
		Spectrum val = bsdf->sample(bRec);
		if (!val.isBlack()) {
			lRec.pdf = bsdf->pdf(bRec);
			lRec.Le = m_intensity;
			lRec.d = -its.toWorld(bRec.wo);
			lRec.sRec.p = its.p - lRec.d * (2 * m_bsphere.radius);
		} else {
			lRec.pdf = 0;
		}
	}

	inline Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec) const {
		const BSDF *bsdf = its.shape->getBSDF();
		BSDFQueryRecord bRec(its, its.toLocal(-lRec.d));
		return bsdf->pdf(bRec);
	}

	/**
	 * This is the tricky bit - we want to sample a ray that
	 * has uniform density over the set of all rays passing
	 * through the scene.
	 * For more detail, see "Using low-discrepancy sequences and 
	 * the Crofton formula to compute surface areas of geometric models"
	 * by Li, X. and Wang, W. and Martin, R.R. and Bowyer, A. 
	 * (Computer-Aided Design vol 35, #9, pp. 771--782)
	 */
	void sampleEmission(EmissionRecord &eRec, 
		const Point2 &sample1, const Point2 &sample2) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		/* Chord model - generate the ray passing through two uniformly
		   distributed points on a sphere containing the scene */
		Vector d = squareToSphere(sample1);
		eRec.sRec.p = m_bsphere.center + d * m_bsphere.radius;
		eRec.sRec.n = Normal(-d);
		Point p2 = m_bsphere.center + squareToSphere(sample2) * m_bsphere.radius;
		eRec.d = p2 - eRec.sRec.p;
		Float length = eRec.d.length();

		if (length == 0) {
			eRec.P = Spectrum(0.0f);
			eRec.pdfArea = eRec.pdfDir = 1.0f;
			return;
		}

		eRec.d /= length;
		eRec.pdfArea = 1.0f / (4 * M_PI * m_bsphere.radius * m_bsphere.radius);
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		eRec.P = m_intensity;
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		Float radius = m_bsphere.radius;
		if (eRec.type == EmissionRecord::EPreview) {
			/* This is more suitable for VPL-based rendering */
			radius *= 1.5;
		}
		Vector d = squareToSphere(sample);
		eRec.sRec.p = m_bsphere.center + d * radius;
		eRec.sRec.n = Normal(-d);
		eRec.pdfArea = 1.0f / (4 * M_PI * radius * radius);
		eRec.P = m_intensity * M_PI;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		Float radius = m_bsphere.radius;
		if (eRec.type == EmissionRecord::EPreview) 
			radius *= 1.5f;
		Point p2 = m_bsphere.center + squareToSphere(sample) * radius;
		eRec.d = p2 - eRec.sRec.p;
		Float length = eRec.d.length();

		if (length == 0.0f) {
			eRec.pdfDir = 1.0f;
			return Spectrum(0.0f);
		}
		
		eRec.d /= length;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		return Spectrum(INV_PI);
	}

	void pdfEmission(EmissionRecord &eRec) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = INV_PI * dp;
		else
			eRec.pdfDir = 0;
		eRec.pdfArea = 1.0f / (4 * M_PI * m_bsphere.radius * m_bsphere.radius);
	}

	Spectrum f(const EmissionRecord &eRec) const {
		return Spectrum(INV_PI);
	}

	Spectrum fArea(const EmissionRecord &eRec) const {
		return m_intensity * M_PI;
	}

	bool isBackgroundLuminaire() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "ConstantLuminaire[" << std::endl
			<< "  intensity = " << m_intensity.toString() << "," << std::endl
			<< "  power = " << getPower().toString() << std::endl
			<< "]";
		return oss.str();
	}
	
	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	BSphere m_bsphere;
};

// ================ Hardware shader implementation ================ 

class ConstantLuminaireShader : public Shader {
public:
	ConstantLuminaireShader(Renderer *renderer, const Spectrum &intensity) 
		: Shader(renderer, ELuminaireShader), m_emittance(intensity * M_PI) {
	}
	
	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_emittance", false));
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_emittance;" << endl
			<< endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "    return vec3(0.31831);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_background(vec3 wo) {" << endl
			<< "    return " << evalName << "_emittance * 0.31831;" << endl
			<< "}" << endl;
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_emittance);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_emittance;
};

Shader *ConstantLuminaire::createShader(Renderer *renderer) const { 
	return new ConstantLuminaireShader(renderer, m_intensity);
}

MTS_IMPLEMENT_CLASS_S(ConstantLuminaire, false, Luminaire)
MTS_IMPLEMENT_CLASS(ConstantLuminaireShader, false, Shader)
MTS_EXPORT_PLUGIN(ConstantLuminaire, "Constant background luminaire");
MTS_NAMESPACE_END
