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

MTS_NAMESPACE_BEGIN

/**
 * Simple directional luminaire. In an untransformed state, the luminaire 
 * sends light into the positive Z direction. To change the direction,
 * either specify a "toWorld" transformation or set the "direction" parameter.
 */
class DirectionalLuminaire : public Luminaire {
public:
	DirectionalLuminaire(const Properties &props) : Luminaire(props) {
		m_intensity = props.getSpectrum("intensity", Spectrum(1.0f));
		m_diskRadius = 0;
		m_type = EDeltaDirection;
		if (props.hasProperty("toWorld") && props.hasProperty("direction"))
			Log(EError, "Please specify either a direction or a luminaire "
				"to world space transformation");
		else if (props.hasProperty("direction"))
			m_direction = normalize(props.getVector("direction"));
		else 
			m_direction = normalize(m_luminaireToWorld(Vector(0, 0, 1)));
	}

	DirectionalLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_direction = Vector(stream);
		m_intensity = Spectrum(stream);
		m_diskOrigin = Point(stream);
		m_diskRadius = stream->readFloat();
		configure();
	}

	void configure() {
		m_surfaceArea = m_diskRadius * m_diskRadius * M_PI;
		m_invSurfaceArea = 1.0f / m_surfaceArea;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);
		m_direction.serialize(stream);
		m_intensity.serialize(stream);
		m_diskOrigin.serialize(stream);
		stream->writeFloat(m_diskRadius);
	}

	void preprocess(const Scene *scene) {
		/* Get the scene's bounding sphere and slightly enlarge it */
		BSphere bsphere = scene->getBSphere();
		m_diskRadius = bsphere.radius;
		m_diskOrigin = bsphere.center - m_direction * bsphere.radius;

		configure();
	}

	Spectrum getPower() const {
		return m_intensity * m_surfaceArea;
	}

	Float pdf(const Point &p, const LuminaireSamplingRecord &lRec, bool delta) const {
		/* PDF is a delta function - zero probability when a sample point was not
		   generated using sample() */
		return delta ? 1.0f : 0.0f;
	}

	void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		lRec.sRec.p = p - m_direction * (2 * m_diskRadius);
		lRec.d = m_direction;
		lRec.luminaire = this;
		lRec.pdf = 1.0f;
		lRec.value = m_intensity;
	}

	void sampleEmission(EmissionRecord &eRec, const Point2 &sample1, const Point2 &sample2) const {
		Point2 posOnDisk = squareToDiskConcentric(sample1) * m_diskRadius;
		eRec.sRec.p = m_diskOrigin + Frame(m_direction).toWorld(Vector(posOnDisk.x, posOnDisk.y, 0));
		eRec.d = m_direction;
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = 1;
		eRec.value = m_intensity;
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		Point2 posOnDisk = squareToDiskConcentric(sample) * m_diskRadius;
		eRec.sRec.p = m_diskOrigin + Frame(m_direction).toWorld(Vector(posOnDisk.x, posOnDisk.y, 0));
		eRec.pdfArea = m_invSurfaceArea;
		eRec.value = m_intensity;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.d = m_direction;
		eRec.pdfDir = 1;
		return Spectrum(1.0f);
	}

	Spectrum evalArea(const EmissionRecord &eRec) const {
		return m_intensity;
	}

	Spectrum evalDirection(const EmissionRecord &eRec) const {
		/* Directional luminaire beam is not part of the scene */
		return Spectrum(0.0f);
	}

	void pdfEmission(EmissionRecord &eRec, bool delta) const {
		eRec.pdfArea = delta ? 0.0f : m_invSurfaceArea;
		eRec.pdfDir = delta ? 1.0f : 0.0f;
		eRec.value = m_intensity;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "DirectionalLuminaire[" << std::endl
			<< "  name = \"" << m_name << "\"," << std::endl
			<< "  intensity = " << m_intensity.toString() << "," << std::endl
			<< "  power = " << getPower().toString() << "," << std::endl
			<< "  direction = " << m_direction.toString() 
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	Vector m_direction;
	Float m_surfaceArea;
	Float m_invSurfaceArea;
	Float m_diskRadius;
	Point m_diskOrigin;
};

MTS_IMPLEMENT_CLASS_S(DirectionalLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(DirectionalLuminaire, "Directional luminaire");
MTS_NAMESPACE_END
