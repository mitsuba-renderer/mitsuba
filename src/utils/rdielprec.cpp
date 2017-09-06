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

#include <mitsuba/render/util.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <boost/bind.hpp>
#include <fstream>
#include <iomanip>
#include <omp.h>

#define RESOLUTION_IOR          50
#define RESOLUTION_ROUGHNESS    30
#define RESOLUTION_THETA        100

#define ROUGHNESS_START         0
#define ROUGHNESS_END           0.5

#define IOR_START               (1 + 1e-4)
#define IOR_END                 4

MTS_NAMESPACE_BEGIN

void transmittanceIntegrand(const BSDF *bsdf, const Vector &wi, size_t nPts, const Float *in, Float *out) {
    Intersection its;

    #pragma omp parallel for
    for (int i=0; i<(int) nPts; ++i) {
        BSDFSamplingRecord bRec(its, wi, Vector(), EImportance);
        bRec.typeMask = BSDF::ETransmission;
        Point2 sample(in[2*i], in[2*i+1]);
        if (sample.x == 1)
            sample.x = 1-Epsilon;
        if (sample.y == 1)
            sample.y = 1-Epsilon;
        out[i] = bsdf->sample(bRec, sample)[0];
        if (std::isnan(out[i]))
            SLog(EError, "%s\n\nNaN!", bRec.toString().c_str());
    }
}

void diffTransmittanceIntegrand(Float *data, size_t resolution, size_t nPts, const Float *in, Float *out) {
    #pragma omp parallel for
    for (int i=0; i<(int) nPts; ++i)
        out[i] = 2 * in[i] * interpCubic1D(std::pow(in[i], (Float) 0.25f), data, 0, 1, resolution);
}

class PrecomputeTransmittance : public Utility {
public:
    Float *computeTransmittance(const char *name, Float ior, Float alpha,
            size_t resolution, Float &diffTrans, int inverted) {
        Properties bsdfProps(alpha == 0 ? "dielectric" : "roughdielectric");
        if (inverted) {
            bsdfProps.setFloat("intIOR", 1.00);
            bsdfProps.setFloat("extIOR", ior);
        } else {
            bsdfProps.setFloat("extIOR", 1.00);
            bsdfProps.setFloat("intIOR", ior);
        }
        bsdfProps.setFloat("alpha", alpha);
        bsdfProps.setString("distribution", name);
        ref<BSDF> bsdf = static_cast<BSDF *>(
                PluginManager::getInstance()->createObject(bsdfProps));

        Float stepSize = 1.0f / (resolution-1);
        Float error;

        NDIntegrator intTransmittance(1, 2, 50000, 0, 1e-6f);
        NDIntegrator intDiffTransmittance(1, 1, 50000, 0, 1e-6f);
        Float *transmittances = new Float[resolution];

        for (size_t i=0; i<resolution; ++i) {
            Float t = i * stepSize;
            if (i == 0) /* Don't go all the way to zero */
                t = stepSize/10;

            Float cosTheta = std::pow(t, (Float) 4.0f);

            Vector wi(math::safe_sqrt(1-cosTheta*cosTheta), 0, cosTheta);

            Float min[2] = {0, 0}, max[2] = {1, 1};
            intTransmittance.integrateVectorized(
                boost::bind(&transmittanceIntegrand, bsdf, wi, _1, _2, _3),
                min, max, &transmittances[i], &error, NULL);
        }

        Float min[1] = { 0 }, max[1] = { 1 };
        intDiffTransmittance.integrateVectorized(
            boost::bind(&diffTransmittanceIntegrand, transmittances, resolution, _1, _2, _3),
            min, max, &diffTrans, &error, NULL);

        if (alpha == 0.0f)
            cout << diffTrans << " vs " << 1-fresnelDiffuseReflectance(inverted ? (1.0f / ior) : ior) << endl;

        return transmittances;
    }

    void fit(const char *name, int inverted) {
        size_t resolutionIOR = RESOLUTION_IOR,
               resolutionAlpha = RESOLUTION_ROUGHNESS,
               resolutionTheta = RESOLUTION_THETA;

        Float alphaStart = ROUGHNESS_START,
              alphaEnd   = ROUGHNESS_END,
              iorStart   = IOR_START,
              iorEnd     = IOR_END;

        std::ofstream os(formatString("%s%i.m", name, inverted).c_str());
        os << std::fixed << std::setprecision(8);

        os << "alphaStart=" << alphaStart << ";" << endl;
        os << "alphaEnd=" << alphaEnd << ";" << endl;
        os << "iorStart=" << iorStart << ";" << endl;
        os << "iorEnd=" << iorEnd << ";" << endl;
        os << "alphaSteps=" << resolutionAlpha << ";" << endl;
        os << "thetaSteps=" << resolutionTheta << ";" << endl;
        os << "iorSteps=" << resolutionIOR << ";" << endl;

        os << "transmittance={" << endl;
        ref<FileStream> fstream = new FileStream(formatString("data/microfacet/%s.dat", name).c_str(),
                inverted ? FileStream::EAppendReadWrite : FileStream::ETruncReadWrite);

        fstream->setByteOrder(Stream::ELittleEndian);

        if (!inverted) {
            fstream->write("MTS_TRANSMITTANCE", 17);
            fstream->writeSize(resolutionIOR);
            fstream->writeSize(resolutionAlpha);
            fstream->writeSize(resolutionTheta);
            fstream->writeSingle(iorStart);
            fstream->writeSingle(iorEnd);
            fstream->writeSingle(alphaStart);
            fstream->writeSingle(alphaEnd);
        }

        Float iorStepSize   = 1.0f / (resolutionIOR-1),
              alphaStepSize = 1.0f / (resolutionAlpha-1);

        for (size_t i=0; i<resolutionIOR; ++i) {
            Float t = i * iorStepSize;
            Float ior = iorStart + (iorEnd-iorStart) * std::pow(t, (Float) 4.0f);
            cout << "ior = " << ior << endl;
            os << "\t{" << endl;
            for (size_t j=0; j<resolutionAlpha; ++j) {
                Float t = j * alphaStepSize;
                Float alpha = alphaStart + (alphaEnd-alphaStart) * std::pow(t, (Float) 4.0f);
                cout << "alpha = " << alpha << endl;
                Float diffTrans;
                Float *transmittance = computeTransmittance(name,
                        ior, alpha, resolutionTheta, diffTrans, inverted);
                os << "\t\t{";
                for (size_t k=0; k<resolutionTheta; ++k) {
                    fstream->writeSingle((float) transmittance[k]);
                    os << transmittance[k];
                    if (k+1 < resolutionTheta)
                        os << ", ";
                }
                os << "}";
                if (j + 1 < resolutionAlpha)
                    os << ",";
                os << endl;
                delete[] transmittance;
                fstream->writeSingle((float) diffTrans);
            }
            os << "\t}";
            if (i + 1 < resolutionIOR)
                os << ",";
            os << endl;
            cout << endl;
        }
        os << "};" << endl;
        os.close();
        fstream->close();
    }

    int run(int argc, char **argv) {
        cout << endl;
        cout << "rdielprec: precompute transmittances through a rough dielectric material" << endl;
        cout << "========================================================================" << endl;
        cout << "The transmittance through a rough dielectric boundary modeled using " << endl;
        cout << "a microfacet model depends on several quantities:" << endl;
        cout << endl;
        cout << "  1. the relative index of refraction" << endl;
        cout << "  2. the angle of incidence (of incoming illumination)" << endl;
        cout << "  3. the used microfacet distribution (beckmann, phong, ggx)" << endl;
        cout << "  4. the roughness parameter of the microfacet distribution" << endl;
        cout << endl;
        cout << "Since a 2D integration is involved, the transmittance is quite expensive" << endl;
        cout << "to evaluate at runtime. This utility therefore precomputes the values on" << endl;
        cout << "a densely spaced grid and stores them on disk so that they can be " << endl;
        cout << "interpolated at runtime (e.g. using cubic splines). The output will be " << endl;
        cout << "written to 'data/microfacet'. This process will take 1-3 hours." << endl;
        cout << "It is recommended that you compile Mitsuba in double precision when running" << endl;
        cout << "this utility." << endl;
        cout << endl;
        cout << "(Press any key to proceed)" << endl;
        std::getchar();
        omp_set_num_threads(4);
        ref<Timer> timer = new Timer();

//      fit("beckmann", 0);
//      fit("beckmann", 1);
        fit("phong", 0);
        fit("phong", 1);
//      fit("ggx", 0);
//      fit("ggx", 1);
        cout << endl;
        cout << "Done, took " << timer->getMilliseconds() / 1000.0f << " seconds" << endl;
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(PrecomputeTransmittance, "Precompute transmittance data for rough dielectrics")
MTS_NAMESPACE_END
