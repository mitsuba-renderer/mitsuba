/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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

#include <mitsuba/core/frame.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>

/// Generate a few statistics related to the implementation?
// #define MICROFLAKE_STATISTICS 1

/// The following file implements the micro-flake distribution
/// for rough fibers
#include "microflake_fiber.h"

MTS_NAMESPACE_BEGIN

#if defined(MICROFLAKE_STATISTICS)
static StatsCounter avgSampleIterations("Micro-flake model",
        "Average rejection sampling iterations", EAverage);
#endif

/*!\plugin{microflake}{Micro-flake phase function}
 * \parameters{
 *     \parameter{stddev}{\Float}{
 *       Standard deviation of the micro-flake normals. This
 *       specifies the roughness of the fibers in the medium.
 *     }
 * }
 *
 * \renderings{
 *     \rendering{\code{stddev}=0.2}{phase_microflakes_02}
 *     \rendering{\code{stddev}=0.05}{phase_microflakes_005}
 * }
 *
 * This plugin implements the anisotropic micro-flake phase function
 * described in
 * ``A radiative transfer framework for rendering materials with
 * anisotropic structure'' by Wenzel Jakob, Adam Arbree,
 * Jonathan T. Moon, Kavita Bala, and Steve Marschner
 * \cite{Jakob2010Radiative}.
 *
 * The implementation in this plugin is specific to rough fibers
 * and uses a Gaussian-type flake distribution. It is much faster
 * than the spherical harmonics approach proposed in the original
 * paper. This distribution, as well as the implemented sampling
 * method, are described in the paper
 * ``Building Volumetric Appearance Models of Fabric using
 * Micro CT Imaging'' by Shuang Zhao, Wenzel Jakob, Steve Marschner,
 * and Kavita Bala \cite{Zhao2011Building}.
 *
 * Note: this phase function must be used with a medium that specifies
 * the local fiber orientation at different points in space. Please
 * refer to \pluginref{heterogeneous} for details.
 *
 * \vspace{1cm}
 *
 * \renderings{
 *   \fbox{\includegraphics[width=\textwidth]{images/medium_knitpatterns}}
 *   \caption{A range of different knit patterns, rendered using the
 *   \pluginref{heterogeneous} and \pluginref{microflake} plugins. Courtesy
 *   of Yuksel et al. \cite{Yuksel2012Stitch}.}
 * }
 */
class MicroflakePhaseFunction : public PhaseFunction {
public:
    MicroflakePhaseFunction(const Properties &props) : PhaseFunction(props) {
        /// Standard deviation of the flake distribution
        m_fiberDistr = GaussianFiberDistribution(props.getFloat("stddev"));
    }

    MicroflakePhaseFunction(Stream *stream, InstanceManager *manager)
        : PhaseFunction(stream, manager) {
        m_fiberDistr = GaussianFiberDistribution(stream->readFloat());
        configure();
    }

    virtual ~MicroflakePhaseFunction() { }

    void configure() {
        PhaseFunction::configure();
        m_type = EAnisotropic | ENonSymmetric;
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        PhaseFunction::serialize(stream, manager);
        stream->writeFloat(m_fiberDistr.getStdDev());
    }

    Float eval(const PhaseFunctionSamplingRecord &pRec) const {
        if (pRec.mRec.orientation.isZero()) {
            /* What to do when the local orientation is undefined */
            #if 0
                return 1.0f / (4 * M_PI);
            #else
                return 0.0f;
            #endif
        }

        Frame frame(pRec.mRec.orientation);
        Vector wi = frame.toLocal(pRec.wi);
        Vector wo = frame.toLocal(pRec.wo);
        Vector H = wi + wo;
        Float length = H.length();

        if (length == 0)
            return 0.0f;

        return 0.5f * m_fiberDistr.pdfCosTheta(Frame::cosTheta(H)/length)
                / m_fiberDistr.sigmaT(Frame::cosTheta(wi));
    }

    inline Float sample(PhaseFunctionSamplingRecord &pRec, Sampler *sampler) const {
        if (pRec.mRec.orientation.isZero()) {
            /* What to do when the local orientation is undefined */
            #if 0
                pRec.wo = warp::squareToUniformSphere(sampler->next2D());
                return 1.0f;
            #else
                return 0.0f;
            #endif
        }

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

    Float sample(PhaseFunctionSamplingRecord &pRec,
            Float &pdf, Sampler *sampler) const {
        if (sample(pRec, sampler) == 0) {
            pdf = 0; return 0.0f;
        }
        pdf = eval(pRec);
        return 1.0f;
    }

    bool needsDirectionallyVaryingCoefficients() const { return true; }

    Float sigmaDir(Float cosTheta) const {
        // Scaled such that replacing an isotropic phase function with an
        // isotropic microflake distribution does not cause changes
        return 2 * m_fiberDistr.sigmaT(cosTheta);
    }

    Float sigmaDirMax() const {
        return sigmaDir(0);
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
