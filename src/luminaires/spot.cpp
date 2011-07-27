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
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/**
 * VRML SpotLight-equivalent light source. In its local coordinate system,
 * the spot light is positioned at the origin and points into the positive Z
 * direction. Its intensity linearly ramps up between <tt>cutoffAngle</tt>
 * and <tt>beamWidth</tt>, after which it remains at the maximum value.
 * A projection texture may optionally be supplied.
 */
class SpotLuminaire : public Luminaire {
public:
	SpotLuminaire(const Properties &props) : Luminaire(props) {
		m_intensity = props.getSpectrum("intensity", Spectrum(1.0f));
		m_cutoffAngle = props.getFloat("cutoffAngle", 20);
		m_beamWidth = props.getFloat("beamWidth", m_cutoffAngle * 3.0f/4.0f);
		m_beamWidth = degToRad(m_beamWidth);
		m_cutoffAngle = degToRad(m_cutoffAngle);
		Assert(m_cutoffAngle >= m_beamWidth);
		m_type = EDeltaPosition;
		m_texture = new ConstantSpectrumTexture(
			props.getSpectrum("texture", Spectrum(1.0f)));
	}

	SpotLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_texture = static_cast<Texture *>(manager->getInstance(stream));
		m_intensity = Spectrum(stream);
		m_beamWidth = stream->readFloat();
		m_cutoffAngle = stream->readFloat();
		configure();
	}

	void configure() {
		m_cosBeamWidth = std::cos(m_beamWidth);
		m_cosCutoffAngle = std::cos(m_cutoffAngle);
		m_position = m_luminaireToWorld(Point(0, 0, 0));
		m_uvFactor = std::tan(m_cutoffAngle);
		m_invTransitionWidth = 1.0f / (m_cutoffAngle - m_beamWidth);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);

		manager->serialize(stream, m_texture.get());
		m_intensity.serialize(stream);
		stream->writeFloat(m_beamWidth);
		stream->writeFloat(m_cutoffAngle);
	}

	Spectrum getPower() const {
		/* The getAverage() part here is technically not really correct, 
		   but this only has to be an approximation after all.. */
		if (m_beamWidth == m_cutoffAngle)
			return m_intensity * m_texture->getAverage() * 2 * (Float) M_PI * (1-m_cosCutoffAngle);
		else
			return m_intensity * m_texture->getAverage() * (2 * M_PI * (1 - 
				(std::sin(m_cutoffAngle) - std::sin(m_beamWidth))/
				(m_cutoffAngle - m_beamWidth)));
	}

	inline Spectrum falloffCurve(const Vector &d, bool throughputOnly = false) const {
		Spectrum result(throughputOnly ? Spectrum(1.0f) : m_intensity);
		Vector localDir = m_worldToLuminaire(d);
		const Float cosTheta = localDir.z;
		
		if (cosTheta <= m_cosCutoffAngle)
			return Spectrum(0.0f);

		if (m_texture->getClass() != MTS_CLASS(ConstantSpectrumTexture)) {
			Intersection its;
			its.hasUVPartials = false;
			its.uv.x = 0.5f + 0.5f * localDir.x / (localDir.z * m_uvFactor);
			its.uv.y = 0.5f + 0.5f * localDir.y / (localDir.z * m_uvFactor);
			result *= m_texture->getValue(its);
		}

		if (cosTheta >= m_cosBeamWidth)
			return result;

		return result * ((m_cutoffAngle - std::acos(cosTheta))
				* m_invTransitionWidth);
	}

	Float pdf(const Point &p, const LuminaireSamplingRecord &lRec, bool delta) const {
		/* PDF is a delta function - zero probability when a sample point was not
		   generated using sample() */
		return delta ? 1.0f : 0.0f;
	}
	
	void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		Vector lumToP = p - m_position;
		Float invDist = 1.0f / lumToP.length();
		lRec.sRec.p = m_position;
		lRec.d = lumToP * invDist;
		lRec.pdf = 1.0f;
		lRec.value = falloffCurve(lRec.d) * (invDist*invDist);
	}

	void sampleEmission(EmissionRecord &eRec, 
		const Point2 &sample1, const Point2 &sample2) const {
		eRec.sRec.p = m_position;
		m_luminaireToWorld(squareToCone(m_cosCutoffAngle, sample2), eRec.d);
		eRec.pdfDir = squareToConePdf(m_cosCutoffAngle);
		eRec.pdfArea = 1;
		eRec.value = falloffCurve(eRec.d);
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.sRec.p = m_position;
		eRec.pdfArea = 1;
		eRec.value = m_intensity;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		m_luminaireToWorld(squareToCone(m_cosCutoffAngle, sample), eRec.d);
		eRec.pdfDir = squareToConePdf(m_cosCutoffAngle);
		return falloffCurve(eRec.d, true);
	}

	void pdfEmission(EmissionRecord &eRec, bool delta) const {
		if (Frame::cosTheta(m_worldToLuminaire(eRec.d))>m_cosCutoffAngle)
			eRec.pdfDir = 0;
		else
			eRec.pdfDir = delta ? 0.0f : squareToConePdf(m_cosCutoffAngle);
		eRec.pdfArea = delta ? 1.0f : 0.0f;
	}

	Spectrum evalDirection(const EmissionRecord &eRec) const {
		return falloffCurve(eRec.d, true);
	}

	Spectrum evalArea(const EmissionRecord &eRec) const {
		return Spectrum(0.0f);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "texture") {
			m_texture = static_cast<Texture *>(child);
		} else {
			Luminaire::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SpotLuminaire[" << std::endl
			<< "  name = \"" << m_name << "\"," << std::endl
			<< "  intensity = " << m_intensity.toString() << "," << std::endl
			<< "  texture = " << m_texture.toString() << "," << std::endl
			<< "  position = " << m_position.toString() << "," << std::endl
			<< "  direction = " << normalize(m_luminaireToWorld(Point(0,0,1))-m_luminaireToWorld(Point(0,0,0))).toString() << "," << std::endl
			<< "  beamWidth = " << (m_beamWidth * 180/M_PI) << "," << std::endl
			<< "  cutoffAngle = " << (m_cutoffAngle * 180/M_PI) << std::endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	ref<Texture> m_texture;
	Float m_beamWidth, m_cutoffAngle, m_uvFactor;
	Float m_cosBeamWidth, m_cosCutoffAngle, m_invTransitionWidth;
	Point m_position;
};

// ================ Hardware shader implementation ================ 

class SpotLuminaireShader : public Shader {
public:
	SpotLuminaireShader(Renderer *renderer, Transform worldToLuminaire, 
		Float invTransitionWidth, Float cutoffAngle, Float cosCutoffAngle, 
		Float cosBeamWidth, Float uvFactor, const Texture *texture) 
		: Shader(renderer, ELuminaireShader), m_worldToLuminaire(worldToLuminaire),
		  m_invTransitionWidth(invTransitionWidth), m_cutoffAngle(cutoffAngle), 
		  m_cosCutoffAngle(cosCutoffAngle), m_cosBeamWidth(cosBeamWidth), 
		  m_uvFactor(uvFactor), m_texture(texture) {
		m_textureShader = renderer->registerShaderForResource(m_texture.get());
	}

	bool isComplete() const {
		return m_textureShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_texture.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_textureShader.get());
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_invTransitionWidth;" << endl
			<< "uniform float " << evalName << "_cutoffAngle;" << endl
			<< "uniform float " << evalName << "_cosCutoffAngle;" << endl
			<< "uniform float " << evalName << "_cosBeamWidth;" << endl
			<< "uniform float " << evalName << "_uvFactor;" << endl
			<< "uniform mat4 " << evalName << "_worldToLuminaire;" << endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "    vec3 localDir = (" << evalName << "_worldToLuminaire * vec4(wo, 0)).xyz;" << endl
			<< "    float cosTheta = localDir.z;" << endl
			<< "    if (cosTheta < " << evalName << "_cosCutoffAngle)" << endl
			<< "        return vec3(0.0);" << endl
			<< "    vec2 uv = 0.5 + 0.5 * (localDir.xy / localDir.z * " << evalName << "_uvFactor);" << endl
			<< "    vec3 color = " << depNames[0] << "(uv);" << endl
			<< "    if (cosTheta > " << evalName << "_cosBeamWidth)" << endl
			<< "        return color;" << endl
			<< "    return color * ((" << evalName << "_cutoffAngle - acos(cosTheta))" << endl
			<< "           * " << evalName << "_invTransitionWidth);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_worldToLuminaire", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_invTransitionWidth", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_cutoffAngle", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_cosCutoffAngle", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_cosBeamWidth", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvFactor", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_worldToLuminaire);
		program->setParameter(parameterIDs[1], m_invTransitionWidth);
		program->setParameter(parameterIDs[2], m_cutoffAngle);
		program->setParameter(parameterIDs[3], m_cosCutoffAngle);
		program->setParameter(parameterIDs[4], m_cosBeamWidth);
		program->setParameter(parameterIDs[5], m_uvFactor);
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_worldToLuminaire;
	Float m_invTransitionWidth;
	Float m_cutoffAngle, m_cosCutoffAngle;
	Float m_cosBeamWidth, m_uvFactor;
	ref<const Texture> m_texture;
	ref<Shader> m_textureShader;
};
	
Shader *SpotLuminaire::createShader(Renderer *renderer) const { 
	return new SpotLuminaireShader(renderer, m_worldToLuminaire,
		m_invTransitionWidth, m_cutoffAngle, m_cosCutoffAngle,
		m_cosBeamWidth, m_uvFactor, m_texture.get());
}

MTS_IMPLEMENT_CLASS(SpotLuminaireShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SpotLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(SpotLuminaire, "Spot light");
MTS_NAMESPACE_END
