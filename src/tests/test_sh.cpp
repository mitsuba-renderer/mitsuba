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

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/shvector.h>

MTS_NAMESPACE_BEGIN

class TestSphericalHarmonics : public TestCase {
public:
    MTS_BEGIN_TESTCASE()
    MTS_DECLARE_TEST(test01_shRotation)
    MTS_DECLARE_TEST(test02_shSampler)
    MTS_END_TESTCASE()

    void test01_shRotation() {
        /* Generate a random SH expansion, rotate it and
           spot-check 100 times against the original evaluated
           at appropriately rotated positions */

        ref<Random> random = new Random();
        int bands = 8;

        SHVector vec1(bands);
        for (int l=0; l<bands; ++l)
            for (int m=-l; m<=l; ++m)
                vec1(l, m) = random->nextFloat();

        Vector axis(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat())));
        Transform trafo = Transform::rotate(axis, random->nextFloat()*360);
        SHRotation rot(vec1.getBands());

        SHVector::rotation(trafo, rot);
        SHVector vec2(bands);

        rot(vec1, vec2);

        for (int i=0; i<100; ++i) {
            Vector dir1(warp::squareToUniformSphere(Point2(random->nextFloat(), random->nextFloat()))), dir2;
            trafo(dir1, dir2);

            Float value1 = vec1.eval(dir2);
            Float value2 = vec2.eval(dir1);
            assertEqualsEpsilon(value1, value2, Epsilon);
        }
    }

    struct ClampedCos {
        Vector axis;
        ClampedCos(Vector axis) : axis(axis) { }
        Float operator()(const Vector &w) const { return std::max((Float) 0, dot(w, axis)); }
    };

    void test02_shSampler() {
        /* Draw 100 samples from a SH expansion of a clamped cosine-shaped
           distribution and verify the returned probabilities */
        int bands = 13, numSamples = 100, depth = 12;

        Vector v = normalize(Vector(1, 2, 3));
        ref<Random> random = new Random();
        SHVector clampedCos = SHVector(bands);
        clampedCos.project(ClampedCos(v), numSamples);
        //Float clampedCosError = clampedCos.l2Error(ClampedCos(v), numSamples);
        clampedCos.normalize();

        //cout << "Projection error = " << clampedCosError << endl;
        //cout << "Precomputing mip-maps" << endl;
        ref<SHSampler> sampler = new SHSampler(bands, depth);
        //cout << "Done: "<< sampler->toString() << endl;
        Float accum = 0;
        int nsamples = 100, nInAvg = 0;
        for (int i=0; i<=nsamples; ++i) {
            Point2 sample(random->nextFloat(), random->nextFloat());
            Float pdf1 = sampler->warp(clampedCos, sample);
            Float pdf2 = dot(v, sphericalDirection(sample.x, sample.y))/M_PI;
            Float relerr = std::abs(pdf1-pdf2)/pdf2;
            if (pdf2 > 0.01) {
                accum += relerr; ++nInAvg;
                assertTrue(relerr < 0.08);
            }
        }
        assertTrue(accum / nInAvg < 0.01);
    }
};

MTS_EXPORT_TESTCASE(TestSphericalHarmonics, "Testcase for Spherical Harmonics code")
MTS_NAMESPACE_END
