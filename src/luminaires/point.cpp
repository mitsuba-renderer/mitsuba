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

MTS_NAMESPACE_BEGIN

/**
 * Simple point light source
 */
class PointLuminaire : public Luminaire {
public:
	PointLuminaire(const Properties &props) : Luminaire(props) {
		m_intensity = props.getSpectrum("intensity", Spectrum(1));
		m_position = m_luminaireToWorld(Point(0,0,0));
		m_surfaceArea = 0.0f;
		m_type = EDeltaPosition | EDiffuseDirection;
	}

	PointLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_intensity = Spectrum(stream);
		m_position = m_luminaireToWorld(Point(0,0,0));
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);

		m_intensity.serialize(stream);
	}

	Spectrum getPower() const {
		return m_intensity * 4 * M_PI;
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		return Spectrum(0.0f);
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec) const {
		/* PDF is a delta function - zero probability when a sample point was not
		   generated using sample() */
		return 0.0f;
	}

	Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec) const {
		return PointLuminaire::pdf(its.p, lRec);
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		Vector lumToP = p - m_position;
		Float invDist = 1.0f / lumToP.length();
		lRec.sRec.p = m_position;
		lRec.d = lumToP * invDist;
		lRec.pdf = 1.0f;
		lRec.Le = m_intensity * (invDist*invDist);
	}
	
	void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		PointLuminaire::sample(its.p, lRec, sample);
	}
	
	void sampleEmission(EmissionRecord &eRec, 
		const Point2 &sample1, const Point2 &sample2) const {
		eRec.sRec.p = m_position;
		eRec.d = squareToSphere(sample2);
		eRec.pdfDir = 1.0f / (4 * M_PI);
		eRec.pdfArea = 1;
		eRec.P = m_intensity;
	}
	
	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.sRec.p = m_position;
		eRec.pdfArea = 1;
		eRec.P = m_intensity * (4*M_PI);
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.d = squareToSphere(sample);
		eRec.pdfDir = 1.0f / (4 * M_PI);
		return Spectrum(1.0f / (4*M_PI));
	}

	void pdfEmission(EmissionRecord &eRec) const {
		eRec.pdfDir = 1.0f / (4 * M_PI);
		eRec.pdfArea = 0;
	}
	
	Spectrum fArea(const EmissionRecord &eRec) const {
		return Spectrum(0.0f);
	}

	Spectrum f(const EmissionRecord &eRec) const {
		return Spectrum(1/(4*M_PI));
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PointLuminaire[" << std::endl
			<< "  intensity = " << m_intensity.toString() << "," << std::endl
			<< "  position = " << m_position.toString() << std::endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	Point m_position;
};

// ================ Hardware shader implementation ================ 

class PointLuminaireShader : public Shader {
public:
	PointLuminaireShader(Renderer *renderer) 
		: Shader(renderer, ELuminaireShader) {
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "    return vec3(0.079577);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
};

Shader *PointLuminaire::createShader(Renderer *renderer) const { 
	return new PointLuminaireShader(renderer);
}

MTS_IMPLEMENT_CLASS(PointLuminaireShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(PointLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(PointLuminaire, "Point luminaire");
MTS_NAMESPACE_END
