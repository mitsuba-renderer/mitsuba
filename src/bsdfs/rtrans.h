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

#if !defined(__ROUGH_TRANSMITTANCE_H)
#define __ROUGH_TRANSMITTANCE_H

#include <mitsuba/core/fstream.h>
#include <mitsuba/core/spline.h>
#include <mitsuba/core/fresolver.h>
#include "microfacet.h"

#if defined(_MSC_VER)
/// Don't warn about potential divide by zero errors
#pragma warning(disable : 4723)
#endif

MTS_NAMESPACE_BEGIN

/**
 * \brief Utility class for evaluating the transmittance through rough
 * dielectric surfaces modeled using microfacet distributions.
 *
 * The transmittance through a rough dielectric boundary based on
 * a microfacet model depends on several quantities:
 *
 * 1. the relative index of refraction
 * 2. the angle of incidence (of incoming illumination)
 * 3. the used microfacet distribution (beckmann, phong, ggx, ..)
 * 4. the roughness parameter of the microfacet distribution
 *
 * Since a numerical integration over the microfacet distribution is involved
 * in every transmittance evaluation, this function is usually prohibitively
 * expensive. That is the motivation for this class, which instead performs
 * lookups into an extensive precomputed three-dimensional table containing
 * appropriately spaced evaluations of this function. The lookups are combined
 * using tricubic interpolation (more specifically, Catmull-Rom splines).
 *
 * The 3D array is parameterized in a way so that the transmittance
 * function is well-behaved throughout the domain, hence lookups will be
 * quite accurate over a large parameter range (avg. abs. error < 1e-4,
 * for eta in [1, 4] and RMS roughness in [0, 4])
 *
 * In many cases, the flexibility to evaluate the rough transmittance
 * function for any (ior, roughness, angle)-triple is not actually
 * needed, since  the index of refraction parameter might always be
 * constant (and potentially also the roughness).
 *
 * This class therefore provides the operations \ref setEta()
 * and \ref setAlpha(), which reduce the heavy 3D table to
 * an (again spline-interpolated) 1D or 2D slice, which accelerates
 * subsequent lookups.
 *
 * As a final bonus, this class also has support for evaluating the \a diffuse
 * rough transmittance, which is defined as a cosine-weighted integral
 * of the rough transmittance over the incident hemisphere.
 */
class RoughTransmittance : public Object {
public:
    /**
     * \brief Load a rough transmittance data file from disk
     *
     * \param type
     *     Denotes the type of a microfacet distribution,
     *     i.e. Beckmann or GGX
     */
    RoughTransmittance(MicrofacetDistribution::EType type) : m_trans(NULL), m_diffTrans(NULL) {
        std::string name;

        switch (type) {
            case MicrofacetDistribution::EBeckmann: name = "beckmann"; break;
            case MicrofacetDistribution::EPhong: name = "phong"; break;
            case MicrofacetDistribution::EGGX: name = "ggx"; break;
            default:
                SLog(EError, "RoughTransmittance: unsupported distribution type!");
        }

        /* Resolve the precomputed data file */
        fs::path sourceFile = Thread::getThread()->getFileResolver()->resolve(
            formatString("data/microfacet/%s.dat", name.c_str()));

        ref<FileStream> fstream = new FileStream(sourceFile,
                FileStream::EReadOnly);
        fstream->setByteOrder(Stream::ELittleEndian);

        const char header[] = "MTS_TRANSMITTANCE";
        char *fileHeader = (char *) alloca(strlen(header));

        fstream->read(fileHeader, strlen(header));
        if (memcmp(fileHeader, header, strlen(header)) != 0)
            SLog(EError, "Encountered an invalid transmittance data file!");

        m_etaSamples = fstream->readSize();
        m_alphaSamples = fstream->readSize();
        m_thetaSamples = fstream->readSize();

        m_transSize = 2 * m_etaSamples * m_alphaSamples * m_thetaSamples;
        m_diffTransSize = 2 * m_etaSamples * m_alphaSamples;

        SLog(EDebug, "Loading " SIZE_T_FMT "x" SIZE_T_FMT "x" SIZE_T_FMT
            " (%s) rough transmittance samples from \"%s\"", 2*m_etaSamples,
            m_alphaSamples, m_thetaSamples,
            memString((m_transSize + m_diffTransSize) * sizeof(float)).c_str(),
            sourceFile.string().c_str());

        m_trans = new Float[m_transSize];
        m_diffTrans = new Float[m_diffTransSize];
        m_etaFixed = false;
        m_alphaFixed = false;

        m_etaMin = (Float) fstream->readSingle();
        m_etaMax = (Float) fstream->readSingle();
        m_alphaMin = (Float) fstream->readSingle();
        m_alphaMax = (Float) fstream->readSingle();

        SLog(EDebug, "Precomputed data is available for the IOR range "
            "[%.4f, %.1f] and roughness range [%.4f, %.1f]",  m_etaMin,
            m_etaMax, m_alphaMin, m_alphaMax);

        float *temp = new float[m_transSize + m_diffTransSize];
        fstream->readSingleArray(temp, m_transSize + m_diffTransSize);

        float *ptr = temp;
        size_t fdrEntry = 0, dataEntry = 0;
        for (size_t i=0; i<2*m_etaSamples; ++i) {
            for (size_t j=0; j<m_alphaSamples; ++j) {
                for (size_t k=0; k<m_thetaSamples; ++k)
                    m_trans[dataEntry++] = (Float) *ptr++;
                m_diffTrans[fdrEntry++] = (Float) *ptr++;
            }
        }
        delete[] temp;

        SAssert(fstream->getPos() == fstream->getSize());
    }

    /// Release all memory
    virtual ~RoughTransmittance() {
        if (m_trans)
            delete[] m_trans;
        if (m_diffTrans)
            delete[] m_diffTrans;

    }

    /// Return the minimum roughness value that is available in the precomputed data
    inline Float getAlphaMin() { return m_alphaMin; }

    /// Return the maximum roughness value that is available in the precomputed data
    inline Float getAlphaMax() { return m_alphaMax; }

    /// Return the minimum index of refraction that is available in the precomputed data
    inline Float getEtaMin() { return m_etaMin; }

    /// Return the maximum index of refraction that is available in the precomputed data
    inline Float getEtaMax() { return m_etaMax; }

    /**
     * \brief Evaluate the rough transmittance for a given index of refraction,
     * roughness, and angle of incidence.
     *
     * \param cosTheta
     *     Cosine of the angle of incidence
     * \param alpha
     *     Roughness parameter
     * \param eta
     *     Relative index of refraction
     */
    Float eval(Float cosTheta, Float alpha = 0, Float eta = 0) const {
        Float warpedCosTheta = std::pow(std::abs(cosTheta), (Float) 0.25f),
              result;

        if (m_alphaFixed && m_etaFixed) {
            if (!(cosTheta >= 0))
                return 0.f;
            //SAssert(cosTheta >= 0);

            result = evalCubicInterp1D(warpedCosTheta,
                m_trans, m_thetaSamples, 0.0f, 1.0f);
        } else if (m_etaFixed) {
            if (!(cosTheta >= 0))
                return 0.f;
            //SAssert(cosTheta >= 0);

            Float warpedAlpha = std::pow((alpha - m_alphaMin)
                    / (m_alphaMax-m_alphaMin), (Float) 0.25f);

            result = evalCubicInterp2D(Point2(warpedCosTheta, warpedAlpha),
                m_trans, Size2(m_thetaSamples, m_alphaSamples),
                Point2(0.0f), Point2(1.0f));
        } else {
            if (cosTheta < 0) {
                cosTheta = -cosTheta;
                eta = 1.0f / eta;
            }

            Float *data = m_trans;
            if (eta < 1) {
                /* Entering a less dense medium -- skip ahead to the
                   second data block */
                data += m_etaSamples * m_alphaSamples * m_thetaSamples;
                eta = 1.0f / eta;
            }

            if (eta < m_etaMin)
                eta = m_etaMin;

            /* Transform the roughness and IOR values into the warped parameter space */
            Float warpedAlpha = std::pow((alpha - m_alphaMin)
                    / (m_alphaMax-m_alphaMin), (Float) 0.25f);
            Float warpedEta = std::pow((eta - m_etaMin)
                    / (m_etaMax-m_etaMin), (Float) 0.25f);

            result = evalCubicInterp3D(Point3(warpedCosTheta, warpedAlpha, warpedEta),
                data, Size3(m_thetaSamples, m_alphaSamples, m_etaSamples),
                Point3(0.0f), Point3(1.0f));
        }

        return std::min((Float) 1.0f, std::max((Float) 0.0f, result));
    }


    /**
     * \brief Evaluate the \a diffuse rough transmittance for a given
     * index of refraction, roughness, and angle of incidence.
     *
     * The diffuse rough transmittance is cosine-weighted integral
     * of the rough transmittance over the incident hemisphere.
     *
     * \param eta
     *     Relative index of refraction
     * \param alpha
     *     Roughness parameter
     */
    Float evalDiffuse(Float alpha = 0, Float eta = 0) const {
        Float result;

        if (m_alphaFixed && m_etaFixed) {
            result = m_diffTrans[0];
        } else if (m_etaFixed) {
            Float warpedAlpha = std::pow((alpha - m_alphaMin)
                    / (m_alphaMax-m_alphaMin), (Float) 0.25f);

            result = evalCubicInterp1D(warpedAlpha,
                m_diffTrans, m_alphaSamples, 0.0f, 1.0f);
        } else {
            Float *data = m_diffTrans;
            if (eta < 1) {
                /* Entering a less dense medium -- skip ahead to the
                   second data block */
                data += m_etaSamples * m_alphaSamples;
                eta = 1.0f / eta;
            }

            if (eta < m_etaMin)
                eta = m_etaMin;

            /* Transform the roughness and IOR values into the warped parameter space */
            Float warpedAlpha = std::pow((alpha - m_alphaMin)
                    / (m_alphaMax-m_alphaMin), (Float) 0.25f);
            Float warpedEta = std::pow((eta - m_etaMin)
                    / (m_etaMax-m_etaMin), (Float) 0.25f);

            result = evalCubicInterp2D(Point2(warpedAlpha, warpedEta), data,
                Size2(m_alphaSamples, m_etaSamples), Point2(0.0f), Point2(1.0f));

        }

        return std::min((Float) 1.0f, std::max((Float) 0.0f,  result));
    }

    /**
     * \brief Reduce the internal 3D table to 2D by specializing
     * to a constant relative index of refraction
     *
     * Should only be called once!
     */
    void setEta(Float eta) {
        if (m_etaFixed)
            return;

        m_transSize = m_alphaSamples * m_thetaSamples;
        m_diffTransSize = m_alphaSamples;

        SLog(EDebug, "Reducing dimension from 3D to 2D (%s), eta = %f",
            memString((m_transSize + m_diffTransSize) * sizeof(Float)).c_str(), eta);

        Float *trans = m_trans,
              *diffTrans = m_diffTrans;

        if (eta < 1) {
            /* Entering a less dense medium -- skip ahead to the
               second data block */
            trans += m_etaSamples * m_alphaSamples * m_thetaSamples;
            diffTrans += m_etaSamples * m_alphaSamples;
            eta = 1.0f / eta;
        }

        if (eta < m_etaMin)
            eta = m_etaMin;

        Float warpedEta = std::pow((eta - m_etaMin)
                / (m_etaMax-m_etaMin), (Float) 0.25f);

        Float *newTrans = new Float[m_transSize];
        Float *newDiffTrans = new Float[m_diffTransSize];

        Float dAlpha = 1.0f / (m_alphaSamples - 1),
              dTheta = 1.0f / (m_thetaSamples - 1);

        for (size_t i=0; i<m_alphaSamples; ++i) {
            for (size_t j=0; j<m_thetaSamples; ++j) {
                newTrans[i*m_thetaSamples + j] = evalCubicInterp3D(
                    Point3(j*dTheta, i*dAlpha, warpedEta),
                    trans, Size3(m_thetaSamples, m_alphaSamples, m_etaSamples),
                    Point3(0.0f), Point3(1.0f));
            }

            newDiffTrans[i] = evalCubicInterp2D(
                    Point2(i*dAlpha, warpedEta),
                    diffTrans, Size2(m_alphaSamples, m_etaSamples),
                    Point2(0.0f), Point2(1.0f));
        }

        delete[] m_trans;
        delete[] m_diffTrans;

        m_trans = newTrans;
        m_diffTrans = newDiffTrans;
        m_etaFixed = true;
    }

    /**
     * \brief Reduce the internal 2D table (after a preceding call to \ref
     * setEta) to 1D by specializing to a constant roughness
     *
     * Should only be called once!
     */
    void setAlpha(Float alpha) {
        if (!m_etaFixed)
            SLog(EError, "setAlpha(): needs a preceding call to setEta()!");
        if (m_alphaFixed)
            return;

        m_transSize = m_thetaSamples;
        m_diffTransSize = 1;

        SLog(EDebug, "Reducing dimension from 2D to 1D (%s), alpha = %f",
            memString((m_transSize + m_diffTransSize) * sizeof(Float)).c_str(), alpha);

        Float warpedAlpha = std::pow((alpha - m_alphaMin)
                / (m_alphaMax-m_alphaMin), (Float) 0.25f);

        Float *newTrans = new Float[m_transSize];
        Float *newDiffTrans = new Float[m_diffTransSize];

        Float dTheta = 1.0f / (m_thetaSamples - 1);

        for (size_t i=0; i<m_thetaSamples; ++i)
            newTrans[i] = evalCubicInterp2D(
                Point2(i*dTheta, warpedAlpha),
                m_trans, Size2(m_thetaSamples, m_alphaSamples),
                Point2(0.0f), Point2(1.0f));

        newDiffTrans[0] = evalCubicInterp1D(warpedAlpha,
                m_diffTrans, m_alphaSamples, 0.0f, 1.0f);

        delete[] m_trans;
        delete[] m_diffTrans;

        m_trans = newTrans;
        m_diffTrans = newDiffTrans;
        m_alphaFixed = true;
    }

    void checkAlpha(Float alpha) {
        if (alpha < m_alphaMin || alpha > m_alphaMax) {
            SLog(EError, "Error: the requested roughness value alpha=%f is"
                " outside of the supported range [%f, %f]! Please scale "
                " your roughness value/texture to lie within this range.",
                alpha, m_alphaMin, m_alphaMax);
        }
    }

    void checkEta(Float eta) {
        if (eta < 1)
            eta = 1/eta;
        if (eta < m_etaMin || eta > m_etaMax)
            SLog(EError, "Error: the requested relative index of refraction "
                "eta=%f is outside of the supported range [%f, %f]! Please "
                "update your  scene so that it uses realistic IOR values.",
                eta, m_etaMin, m_etaMax);
    }

    /// Create a deep copy of the current instance
    ref<RoughTransmittance> clone() const {
        RoughTransmittance *result = new RoughTransmittance();
        result->m_name = m_name;
        result->m_etaSamples = m_etaSamples;
        result->m_alphaSamples = m_alphaSamples;
        result->m_thetaSamples = m_thetaSamples;
        result->m_etaFixed = m_etaFixed;
        result->m_alphaFixed = m_alphaFixed;
        result->m_etaMin = m_etaMin;
        result->m_etaMax = m_etaMax;
        result->m_alphaMin = m_alphaMin;
        result->m_alphaMax = m_alphaMax;
        result->m_transSize = m_transSize;
        result->m_diffTransSize = m_diffTransSize;
        result->m_trans = new Float[m_transSize];
        result->m_diffTrans = new Float[m_diffTransSize];
        memcpy(result->m_trans, m_trans, m_transSize * sizeof(Float));
        memcpy(result->m_diffTrans, m_diffTrans, m_diffTransSize * sizeof(Float));
        return result;
    }
protected:
    inline RoughTransmittance() { }
protected:
    std::string m_name;
    size_t m_etaSamples;
    size_t m_alphaSamples;
    size_t m_thetaSamples;
    bool m_etaFixed;
    bool m_alphaFixed;
    Float m_etaMin, m_etaMax;
    Float m_alphaMin, m_alphaMax;
    size_t m_transSize;
    size_t m_diffTransSize;
    Float *m_trans, *m_diffTrans;
};

MTS_NAMESPACE_END

#endif /* __ROUGH_TRANSMITTANCE_H */
