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
		m_bsphere = BSphere(stream);
		m_surfaceArea = 4 * m_bsphere.radius * m_bsphere.radius * M_PI;
		m_invSurfaceArea = 1/m_surfaceArea;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);

		m_intensity.serialize(stream);
		m_bsphere.serialize(stream);
	}

	void preprocess(const Scene *scene) {
		/* Get the scene's bounding sphere and slightly enlarge it */
		if (m_bsphere.isEmpty()) {
			/* Get the scene's bounding sphere and slightly enlarge it */
			m_bsphere = scene->getBSphere();
			m_bsphere.radius *= 1.01f;
		}
		if (scene->getCamera()) {
			BSphere old = m_bsphere;
			m_bsphere.expandBy(scene->getCamera()->getPosition());
			if (old != m_bsphere)
				m_bsphere.radius *= 1.01f;
		}
		m_surfaceArea = 4 * m_bsphere.radius * m_bsphere.radius * M_PI;
		m_invSurfaceArea = 1/m_surfaceArea;
	}

	Spectrum getPower() const {
		return m_intensity * m_surfaceArea * M_PI;
	}

	Spectrum Le(const Ray &ray) const {
		return m_intensity;
	}

	void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		Vector d = squareToSphere(sample);
	
		Float nearHit, farHit;
		if (m_bsphere.contains(p) && m_bsphere.rayIntersect(Ray(p, d, 0.0f), nearHit, farHit)) {
			lRec.sRec.p = p + d * nearHit;
			lRec.pdf = 1.0f / (4*M_PI);
			lRec.sRec.n = normalize(m_bsphere.center - lRec.sRec.p);
			lRec.d = -d;
			lRec.value = m_intensity;
		} else {
			lRec.pdf = 0.0f;
		}
	}
	
	Float pdf(const Point &p, const LuminaireSamplingRecord &lRec, bool delta) const {
		return delta ? 0.0f : 1.0f / (4*M_PI);
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
			eRec.value = Spectrum(0.0f);
			eRec.pdfArea = eRec.pdfDir = 1.0f;
			return;
		}

		eRec.d /= length;
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		eRec.value = m_intensity;
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
		eRec.value = m_intensity * M_PI;
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

	void pdfEmission(EmissionRecord &eRec, bool delta) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = delta ? 0.0f : INV_PI * dp;
		else
			eRec.pdfDir = 0.0f;
		eRec.pdfArea = delta ? 0.0f : m_invSurfaceArea;
	}

	bool createEmissionRecord(EmissionRecord &eRec, const Ray &ray) const {
		Float nearHit, farHit;
		if (!m_bsphere.contains(ray.o) || !m_bsphere.rayIntersect(ray, nearHit, farHit)) {
			Log(EWarn, "Could not create an emission record -- the ray "
				"in question appears to be outside of the scene bounds!");
			return false;
		}

		eRec.type = EmissionRecord::ENormal;
		eRec.sRec.p = ray(nearHit);
		eRec.sRec.n = normalize(m_bsphere.center - eRec.sRec.p);
		eRec.d = -ray.d;
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		eRec.value = m_intensity;
		eRec.luminaire = this;
		return true;
	}

	Spectrum evalDirection(const EmissionRecord &eRec) const {
		return Spectrum(INV_PI);
	}

	Spectrum evalArea(const EmissionRecord &eRec) const {
		return m_intensity * M_PI;
	}

	bool isBackgroundLuminaire() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "ConstantLuminaire[" << std::endl
			<< "  name = \"" << m_name << "\"," << std::endl
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
	Float m_surfaceArea;
	Float m_invSurfaceArea;
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
