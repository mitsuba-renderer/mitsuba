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
 * Lambertian area light source - can be attached to an arbitrary shape
 * contained inside the scene. Shadow rays are generally sampled
 * uniformly with respect to surface area, which may lead to high
 * variance (e.g. many of the generated samples are facing away
 * from the point to be shaded). 
 * When the shape in question is a sphere, rays are sampled uniformly 
 * wrt. solid angle, which significantly reduces the variance.
 * Thus, spheres are recommended whenever there is some flexibility 
 * in choosing the luminaire shape.
 */
class AreaLuminaire : public Luminaire {
public:
	AreaLuminaire(const Properties &props) : Luminaire(props), m_shape(NULL) {
		AssertEx(m_luminaireToWorld.isIdentity(), "Error: non-identity transformation found. "
			"Area luminaires inherit their transformation from their associated shape!");
		m_intensity = props.getSpectrum("intensity", Spectrum(1));
		m_type = EDiffuseDirection | EOnSurface;
		m_intersectable = true;
	}

	AreaLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_intensity = Spectrum(stream);
		m_shape = static_cast<Shape *>(manager->getInstance(stream));
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);

		m_intensity.serialize(stream);
		manager->serialize(stream, m_shape);
	}

	Spectrum getPower() const {
		return m_intensity * m_shape->getSurfaceArea() * M_PI;
	}

	Spectrum Le(const ShapeSamplingRecord &sRec, const Vector &d) const {
		if (dot(d, sRec.n) <= 0)
			return Spectrum(0.0f);
		return m_intensity;
	}

	void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		lRec.pdf = m_shape->sampleSolidAngle(lRec.sRec, p, sample);
		lRec.d = p - lRec.sRec.p;

		if (EXPECT_TAKEN(lRec.pdf > 0 && dot(lRec.d, lRec.sRec.n) > 0)) {
			lRec.value = m_intensity;
			lRec.d = normalize(lRec.d);
		} else {
			lRec.pdf = 0;
		}
	}

	Float pdf(const Point &p, const LuminaireSamplingRecord &lRec, bool delta) const {
		return m_shape->pdfSolidAngle(lRec.sRec, p);
	}

	void sampleEmission(EmissionRecord &eRec,
		const Point2 &sample1, const Point2 &sample2) const {
		eRec.pdfArea = m_shape->sampleArea(eRec.sRec, sample1);
		Vector wo = squareToHemispherePSA(sample2);
		eRec.pdfDir = Frame::cosTheta(wo) * INV_PI;
		eRec.d = Frame(eRec.sRec.n).toWorld(wo);
		eRec.value = m_intensity;
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.pdfArea = m_shape->sampleArea(eRec.sRec, sample);
		eRec.value = m_intensity * M_PI;
	}

	Spectrum evalArea(const EmissionRecord &eRec) const {
		return m_intensity * M_PI;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		Vector wo = squareToHemispherePSA(sample);
		eRec.d = Frame(eRec.sRec.n).toWorld(wo);
		eRec.pdfDir = Frame::cosTheta(wo) * INV_PI;
		return Spectrum(INV_PI);
	}

	Spectrum evalDirection(const EmissionRecord &eRec) const {
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			return Spectrum(INV_PI);
		else
			return Spectrum(0.0f);
	}

	void pdfEmission(EmissionRecord &eRec, bool delta) const {
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = delta ? 0.0f : dp * INV_PI;
		else {
			eRec.pdfDir = 0;
		}
		eRec.pdfArea = delta ? 0.0f : m_shape->pdfArea(eRec.sRec);
	}

	void setParent(ConfigurableObject *parent) {
		if (parent->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Shape *shape = static_cast<Shape *>(parent);
			if (parent == m_parent || shape->isCompound())
				return;

			if (m_parent)
				Log(EError, "An area light source cannot be parent of multiple shapes");

			ConfigurableObject::setParent(shape);

			m_shape = shape;
			parent->configure();
		} else {
			Log(EError, "An area light source must be child of a shape instance");
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AreaLuminaire[" << std::endl
			<< "  name = \"" << m_name << "\"," << std::endl
			<< "  intensity = " << m_intensity.toString() << std::endl
			<< "]";
		return oss.str();
	}
	
	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
	const Shape *m_shape;
};

// ================ Hardware shader implementation ================ 

class AreaLuminaireShader : public Shader {
public:
	AreaLuminaireShader(Renderer *renderer, const Spectrum &intensity) 
		: Shader(renderer, ELuminaireShader), m_intensity(intensity) {
	}
	
	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_intensity", false));
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_intensity;" << endl
			<< endl
			<< "vec3 " << evalName << "_area(vec2 uv) {" << endl
			<< "    return " << evalName << "_intensity * 3.1415;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "    if (wo.z < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return vec3(0.31831);" << endl
			<< "}" << endl;
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_intensity);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_intensity;
};

Shader *AreaLuminaire::createShader(Renderer *renderer) const { 
	return new AreaLuminaireShader(renderer, m_intensity);
}

MTS_IMPLEMENT_CLASS(AreaLuminaireShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(AreaLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(AreaLuminaire, "Area luminaire");
MTS_NAMESPACE_END
