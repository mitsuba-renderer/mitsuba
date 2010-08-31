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
			"Area luminaires inherit their transformation from the associated shape!");
		m_intensity = props.getSpectrum("intensity", 1);
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
		return m_intensity * m_surfaceArea * M_PI;
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		if (dot(lRec.d, lRec.sRec.n) <= 0)
			return Spectrum(0.0f);
		return m_intensity;
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		lRec.pdf = m_shape->sampleSolidAngle(lRec.sRec, p, sample);
		lRec.d = p - lRec.sRec.p;

		if (EXPECT_TAKEN(lRec.pdf > 0 && dot(lRec.d, lRec.sRec.n) > 0)) {
			lRec.Le = m_intensity;
			lRec.d = normalize(lRec.d);
		} else {
			lRec.pdf = 0;
		}
	}

	inline void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		AreaLuminaire::sample(its.p, lRec, sample);
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec) const {
		return m_shape->pdfSolidAngle(lRec.sRec, p);
	}
	
	Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec) const {
		return pdf(its.p, lRec);
	}

	void sampleEmission(EmissionRecord &eRec,
		const Point2 &sample1, const Point2 &sample2) const {
		eRec.pdfArea = m_shape->sampleArea(eRec.sRec, sample1);
		Vector wo = squareToHemispherePSA(sample2);
		eRec.pdfDir = Frame::cosTheta(wo) * INV_PI;
		eRec.d = Frame(eRec.sRec.n).toWorld(wo);
		eRec.P = m_intensity;
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.pdfArea = m_shape->sampleArea(eRec.sRec, sample);
		eRec.P = m_intensity * M_PI;
	}

	Spectrum fArea(const EmissionRecord &eRec) const {
		return m_intensity * M_PI;
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		Vector wo = squareToHemispherePSA(sample);
		eRec.d = Frame(eRec.sRec.n).toWorld(wo);
		eRec.pdfDir = Frame::cosTheta(wo) * INV_PI;
		return Spectrum(INV_PI);
	}

	Spectrum f(const EmissionRecord &eRec) const {
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			return Spectrum(INV_PI);
		else
			return Spectrum(0.0f);
	}

	void pdfEmission(EmissionRecord &eRec) const {
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = dp * INV_PI;
		else {
			eRec.pdfDir = 0;
		}
		eRec.pdfArea = m_shape->pdfArea(eRec.sRec);
	}

	void setParent(ConfigurableObject *parent) {
		if (parent->getClass()->derivesFrom(Shape::m_theClass)) {
			Shape *shape = static_cast<Shape *>(parent);
			if (parent == m_parent || shape->isCompound())
				return;

			if (m_parent)
				Log(EError, "An area light source cannot be parent of multiple shapes");

			ConfigurableObject::setParent(shape);

			m_shape = shape;
			parent->configure();
			m_surfaceArea = m_shape->getSurfaceArea();
		} else {
			Log(EError, "An area light source must be child of a shape instance");
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AreaLuminaire[" << std::endl
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
			<< "    return " << evalName << "_intensity;" << endl
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
