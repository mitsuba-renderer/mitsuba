#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

std::string PhaseFunctionSamplingRecord::toString() const {
    std::ostringstream oss;
    oss << "PhaseFunctionSamplingRecord[" << endl
        << "  mRec = " << indent(mRec.toString()) << "," << endl
        << "  wi = " << wi.toString() << "," << endl
        << "  wo = " << wo.toString() << "," << endl
        << "  mode = " << mode << endl
        << "]";
    return oss.str();
}

void PhaseFunction::configure() {
    m_type = 0;
}

Float PhaseFunction::pdf(const PhaseFunctionSamplingRecord &pRec) const {
    return eval(pRec);
}

bool PhaseFunction::needsDirectionallyVaryingCoefficients() const {
    return false;
}

Float PhaseFunction::sigmaDir(Float cosTheta) const {
    Log(EError, "%s::sigmaDir(Float) is not implemented (this is not "
        "an anisotropic medium!)", getClass()->getName().c_str());
    return 0.0f;
}

Float PhaseFunction::sigmaDirMax() const {
    Log(EError, "%s::sigmaDirMax() is not implemented (this is not "
        "an anisotropic medium!)", getClass()->getName().c_str());
    return 0.0f;
}

Float PhaseFunction::getMeanCosine() const {
    Log(EError, "%s::getMeanCosine() is not implemented!",
        getClass()->getName().c_str());
    return 0.0f;
}

MTS_IMPLEMENT_CLASS(PhaseFunction, true, ConfigurableObject)
MTS_NAMESPACE_END
