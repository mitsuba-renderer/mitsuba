#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>

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
		m_beamWidth = props.getFloat("beamWidth", m_cutoffAngle * 3.0/4.0);
		m_beamWidth = degToRad(m_beamWidth);
		m_cutoffAngle = degToRad(m_cutoffAngle);
		Assert(m_cutoffAngle >= m_beamWidth);
		m_cosBeamWidth = std::cos(m_beamWidth);
		m_cosCutoffAngle = std::cos(m_cutoffAngle);
		m_surfaceArea = 0.0f;
		m_position = m_luminaireToWorld(Point(0, 0, 0));
		m_type = EDeltaPosition;
		m_texture = new ConstantTexture(
			props.getSpectrum("texture", Spectrum(1.0f)));
		m_uvFactor = std::tan(m_beamWidth/2);
	}

	SpotLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_texture = static_cast<Texture *>(manager->getInstance(stream));
		m_intensity = Spectrum(stream);
		m_beamWidth = stream->readFloat();
		m_cutoffAngle = stream->readFloat();
		m_cosBeamWidth = std::cos(m_beamWidth);
		m_cosCutoffAngle = std::cos(m_cutoffAngle);
		m_position = m_luminaireToWorld(Point(0, 0, 0));
		m_uvFactor = std::tan(m_beamWidth/2);
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

		if (m_texture->getClass() != ConstantTexture::m_theClass) {
			Intersection its;
			its.hasUVPartials = false;
			its.uv.x = .5+localDir.x / (localDir.z / m_uvFactor);
			its.uv.y = .5+localDir.y / (localDir.z / m_uvFactor);
			result *= m_texture->getValue(its);
		}

		if (cosTheta < m_cosCutoffAngle)
			return Spectrum(0.0f);
		if (cosTheta > m_cosBeamWidth)
			return result;
		return result * ((m_cutoffAngle - std::acos(cosTheta))
			/ (m_cutoffAngle - m_beamWidth));
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		if (lRec.sRec.p == m_position)
			return m_intensity * falloffCurve(lRec.d);
		return Spectrum(0.0f);
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec) const {
		/* PDF is a delta function - zero probability when a sample point was not
		   generated using sample() */
		return 0.0f;
	}
	
	Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec) const {
		return SpotLuminaire::pdf(its.p, lRec);
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		Vector lumToP = p - m_position;
		Float invDist = 1.0f / lumToP.length();
		lRec.sRec.p = m_position;
		lRec.d = lumToP * invDist;
		lRec.pdf = 1.0f;
		lRec.Le = falloffCurve(lRec.d) * (invDist*invDist);
	}

	void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		SpotLuminaire::sample(its.p, lRec, sample);
	}

	void sampleEmission(EmissionRecord &eRec, 
		const Point2 &sample1, const Point2 &sample2) const {
		eRec.sRec.p = m_position;
		m_luminaireToWorld(squareToCone(m_cosCutoffAngle, sample2), eRec.d);
		eRec.pdfDir = squareToConePdf(m_cosCutoffAngle);
		eRec.pdfArea = 1;
		eRec.P = falloffCurve(eRec.d);
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.sRec.p = m_position;
		eRec.pdfArea = 1;
		eRec.P = m_intensity;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		m_luminaireToWorld(squareToCone(m_cosCutoffAngle, sample), eRec.d);
		eRec.pdfDir = squareToConePdf(m_cosCutoffAngle);
		return Spectrum(falloffCurve(eRec.d, true));
	}

	void pdfEmission(EmissionRecord &eRec) const {
		if (Frame::cosTheta(m_worldToLuminaire(eRec.d))>m_cosCutoffAngle)
			eRec.pdfDir = 0;
		else
			eRec.pdfDir = squareToConePdf(m_cosCutoffAngle);
		eRec.pdfArea = 0;
	}

	Spectrum f(const EmissionRecord &eRec) const {
		return falloffCurve(eRec.d, true);
	}
	
	Spectrum fArea(const EmissionRecord &eRec) const {
		return Spectrum(0.0f);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "texture") {
			m_texture = static_cast<Texture *>(child);
		} else {
			Luminaire::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SpotLuminaire[" << std::endl
			<< "  intensity = " << m_intensity.toString() << "," << std::endl
			<< "  texture = " << m_texture.toString() << "," << std::endl
			<< "  position = " << m_position.toString() << "," << std::endl
			<< "  direction = " << normalize(m_luminaireToWorld(Point(0,0,1))-m_luminaireToWorld(Point(0,0,0))).toString() << "," << std::endl
			<< "  beamWidth = " << (m_beamWidth * 180/M_PI) << "," << std::endl
			<< "  cutoffAngle = " << (m_cutoffAngle * 180/M_PI) << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	ref<Texture> m_texture;
	Float m_beamWidth, m_cutoffAngle, m_uvFactor;
	Float m_cosBeamWidth, m_cosCutoffAngle;
	Point m_position;
};

MTS_IMPLEMENT_CLASS_S(SpotLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(SpotLuminaire, "Spot light");
MTS_NAMESPACE_END
