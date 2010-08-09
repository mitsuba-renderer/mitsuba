#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

/**
 * Basic isotropic phase function
 */
class IsotropicPhaseFunction : public PhaseFunction {
public:
	IsotropicPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
	}

	IsotropicPhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
	}

	virtual ~IsotropicPhaseFunction() { }


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
			ESampledType &sampledType, const Point2 &sample) const {
		wo = squareToSphere(sample);
		sampledType = ENormal;
		return Spectrum(1.0f);
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
			ESampledType &sampledType, Float &pdf, const Point2 &sample) const {
		wo = squareToSphere(sample);
		sampledType = ENormal;
		pdf = 1/(4 * (Float) M_PI);
		return Spectrum(pdf);
	}

	Spectrum f(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const {
		return Spectrum(1/(4 * (Float) M_PI));
	}

	std::string toString() const {
		return "IsotropicPhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
};


MTS_IMPLEMENT_CLASS_S(IsotropicPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(IsotropicPhaseFunction, "Isotropic phase function");
MTS_NAMESPACE_END
