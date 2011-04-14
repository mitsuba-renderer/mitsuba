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

#include <mitsuba/core/chisquare.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>

#define MICROFLAKE_STATISTICS 1
#include "microflake_fiber.h"

#include <mitsuba/core/plugin.h>///XXX

MTS_NAMESPACE_BEGIN

#if defined(MICROFLAKE_STATISTICS)
static StatsCounter avgSampleIterations("Micro-flake model",
		"Average rejection sampling iterations", EAverage);
#endif

class MicroflakePhaseFunction : public PhaseFunction {
public:
	MicroflakePhaseFunction(const Properties &props) : PhaseFunction(props) {
		m_fiberDistr = GaussianFiberDistribution(props.getFloat("stddev"));
		ChiSquareTest test(7);
		Sampler *sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));
		test.fill(
			boost::bind(&MicroflakePhaseFunction::testSample, this, sampler),
			boost::bind(&MicroflakePhaseFunction::testF, this, _1)
		);
		test.dumpTables("test.m");
		test.runTest(1);
	}

	MicroflakePhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
	}

	virtual ~MicroflakePhaseFunction() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
	}

	void configure() {
	}

	Float f(const PhaseFunctionQueryRecord &pRec) const {
		if (pRec.mRec.orientation.isZero())
			return 0.0f;

		Frame frame(pRec.mRec.orientation);
		Vector wi = frame.toLocal(pRec.wi);
		Vector wo = frame.toLocal(pRec.wo);
		Vector H = wi + wo;
		Float length = H.length();

		if (length == 0)
			return 0.0f;

		Float cosThetaH = H.z/length;
		return 0.5 * m_fiberDistr.pdfCosTheta(cosThetaH)
				/ m_fiberDistr.sigmaT(Frame::cosTheta(wi));
	}

	inline Float sample(PhaseFunctionQueryRecord &pRec, Sampler *sampler) const {
		if (pRec.mRec.orientation.isZero())
			return 0.0f;
		Frame frame(pRec.mRec.orientation);
		Vector wi = frame.toLocal(pRec.wi);

		#if defined(MICROFLAKE_STATISTICS)
			avgSampleIterations.incrementBase();
		#endif

		int iterations = 0, maxIterations = 1000;
		while (true) {
			Vector H = m_fiberDistr.sample(sampler->next2D());
			#if defined(MICROFLAKE_STATISTICS)
				++avgSampleIterations;
			#endif
			++iterations;

			if (sampler->next1D() < absDot(wi, H)) {
				Vector wo = H*(2*dot(wi, H)) - wi;
				pRec.wo = frame.toWorld(wo);
				break;
			}

			if (iterations >= maxIterations) {
				Log(EWarn, "Sample generation unsuccessful after %i iterations"
					" (dp=%f, fiberOrientation=%s, wi=%s)", iterations,
					absDot(pRec.wi, pRec.mRec.orientation),
					pRec.mRec.orientation.toString().c_str(),
					pRec.wi.toString().c_str());
				return 0.0f;
			}
		}

		return 1.0f;
	}

	Float sample(PhaseFunctionQueryRecord &pRec, 
			Float &pdf, Sampler *sampler) const {
		if (sample(pRec, sampler) == 0) {
			pdf = 0;
			return 0.0f;
		}
		pdf = f(pRec);
		return pdf;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MicroflakePhaseFunction[" << endl
			<< "   fiberDistr = " << indent(m_fiberDistr.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	GaussianFiberDistribution m_fiberDistr;
};


MTS_IMPLEMENT_CLASS_S(MicroflakePhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(MicroflakePhaseFunction, "Microflake phase function");
MTS_NAMESPACE_END
