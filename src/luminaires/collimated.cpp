#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Collimated beam with a configurable thickness
 */
class CollimatedBeamLuminaire : public Luminaire {
public:
	CollimatedBeamLuminaire(const Properties &props) : Luminaire(props) {
		m_radius = props.getFloat("radius", 0.01f);
		m_surfaceArea = m_radius * m_radius * M_PI;
		Spectrum power = props.getSpectrum("power", 1);
		m_intensity = props.getSpectrum("intensity", power / m_surfaceArea);
		m_invSurfaceArea = 1 / m_surfaceArea;
		m_direction = m_luminaireToWorld(Vector(0, 0, 1));
		m_type = EDeltaDirection;
	}

	CollimatedBeamLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_intensity = Spectrum(stream);
		m_radius = stream->readFloat();
		m_invSurfaceArea = 1 / m_surfaceArea;
		m_direction = m_luminaireToWorld(Vector(0, 0, 1));
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);
		m_intensity.serialize(stream);
		stream->writeFloat(m_radius);
	}

	Spectrum getPower() const {
		return m_intensity * m_surfaceArea;
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		/* Collimated beam is not part of the scene */
		Log(EWarn, "This function should never be called.");
		return Spectrum(0.0f);
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec) const {
		/* PDF is a delta function - zero probability when a sample point was not
		   generated using sample() */
		return 0.0f;
	}

	Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec) const {
		return CollimatedBeamLuminaire::pdf(its.p, lRec);
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		Point local = m_worldToLuminaire(p);
		Vector2 planeProjection = Vector2(local.x, local.y);

		if (planeProjection.length() > m_radius || local.z < 0) {
			lRec.pdf = 0.0f;
		} else {
			lRec.sRec.p = m_luminaireToWorld(Point(local.x, local.y, 0));
			lRec.d = m_direction;
			lRec.luminaire = this;
			lRec.pdf = 1.0f;
			lRec.Le = m_intensity;
		}
	}

	void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		CollimatedBeamLuminaire::sample(its.p, lRec, sample);
	}

	void sampleEmission(EmissionRecord &eRec, const Point2 &sample1, const Point2 &sample2) const {
		Point2 posOnDisk = squareToDiskConcentric(sample1) * m_radius;
		eRec.sRec.p = m_luminaireToWorld(Point(posOnDisk.x, posOnDisk.y, 0));
		eRec.d = m_direction;
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = 1;
		eRec.P = m_intensity;
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		Point2 posOnDisk = squareToDiskConcentric(sample) * m_radius;
		eRec.sRec.p = m_luminaireToWorld(Point(posOnDisk.x, posOnDisk.y, 0));
		eRec.pdfArea = m_invSurfaceArea;
		eRec.P = m_intensity;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.d = m_direction;
		eRec.pdfDir = 1;
		return Spectrum(1.0f);
	}

	Spectrum fArea(const EmissionRecord &eRec) const {
		return m_intensity;
	}

	Spectrum f(const EmissionRecord &eRec) const {
		/* Collimated beam is not part of the scene */
		Log(EWarn, "This function should never be called.");
		return Spectrum(0.0f);
	}

	void pdfEmission(EmissionRecord &eRec) const {
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = 0;
		eRec.P = m_intensity;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "CollimatedBeamLuminaire[" << std::endl
			<< "  intensity = " << m_intensity.toString() << "," << std::endl
			<< "  power = " << getPower().toString() << "," << std::endl
			<< "  position = " << m_luminaireToWorld(Point(0, 0, 0)).toString() << "," << std::endl
			<< "  direction = " << m_direction.toString() << "," << std::endl
			<< "  radius = " << m_radius << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	Vector m_direction;
	Float m_invSurfaceArea;
	Float m_radius;
};

MTS_IMPLEMENT_CLASS_S(CollimatedBeamLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(CollimatedBeamLuminaire, "Collimated beam luminaire");
MTS_NAMESPACE_END
