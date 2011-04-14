#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

std::string PhaseFunctionQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "PhaseFunctionQueryRecord[" << std::endl
		<< "  mRec = " << indent(mRec.toString()) << "," << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  wo = " << wo.toString() << "," << std::endl
		<< "  quantity = " << quantity << std::endl
		<< "]";
	return oss.str();
}

Float PhaseFunction::pdf(const PhaseFunctionQueryRecord &pRec) const {
	return f(pRec);
}
	
bool PhaseFunction::needsDirectionallyVaryingCoefficients() const {
	return false;
}
	
Float PhaseFunction::coeffMultiplier(Float cosTheta) const {
	Log(EError, "coeffMultiplier(): Not implemented! (this is not"
		" an anisotropic medium)");
	return 0.0f;
}

MTS_IMPLEMENT_CLASS(PhaseFunction, true, ConfigurableObject)
MTS_NAMESPACE_END
