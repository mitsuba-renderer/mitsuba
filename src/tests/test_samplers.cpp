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
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/qmc.h>

MTS_NAMESPACE_BEGIN

class TestSamplers : public TestCase {
public:
    MTS_BEGIN_TESTCASE()
    MTS_DECLARE_TEST(test01_Halton)
    MTS_DECLARE_TEST(test02_Hammersley)
    MTS_DECLARE_TEST(test03_radicalInverseIncr)
    MTS_END_TESTCASE()

    void test01_Halton() {
        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(Sampler), Properties("halton")));

        /* MATLAB: p = haltonset(5); net(p,5) */
        Float comparison[] = {
                            0,                 0,                 0,                 0,                 0,
            0.500000000000000, 0.333333333333333, 0.200000000000000, 0.142857142857143, 0.090909090909091,
            0.250000000000000, 0.666666666666667, 0.400000000000000, 0.285714285714286, 0.181818181818182,
            0.750000000000000, 0.111111111111111, 0.600000000000000, 0.428571428571429, 0.272727272727273,
            0.125000000000000, 0.444444444444444, 0.800000000000000, 0.571428571428571, 0.363636363636364
        };

        int pos = 0;
        sampler->generate(Point2i(0));
        for (int i=0; i<5; ++i) {
            for (int j=0; j<5; ++j)
                assertEqualsEpsilon(sampler->next1D(), comparison[pos++], 1e-7);
            sampler->advance();
        }
    }

    void test02_Hammersley() {
        Properties props("hammersley");
        props.setInteger("sampleCount", 5);

        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(Sampler), props));

        Float comparison[] = {
            0.0,                        0,                 0,                 0,                 0,                 0,
            1.0/5.0,    0.500000000000000, 0.333333333333333, 0.200000000000000, 0.142857142857143, 0.090909090909091,
            2.0/5.0,    0.250000000000000, 0.666666666666667, 0.400000000000000, 0.285714285714286, 0.181818181818182,
            3.0/5.0,    0.750000000000000, 0.111111111111111, 0.600000000000000, 0.428571428571429, 0.272727272727273,
            4.0/5.0,    0.125000000000000, 0.444444444444444, 0.800000000000000, 0.571428571428571, 0.363636363636364
        };

        int pos = 0;
        sampler->generate(Point2i(0));
        for (int i=0; i<5; ++i) {
            for (int j=0; j<6; ++j)
                assertEqualsEpsilon(sampler->next1D(), comparison[pos++], 1e-7);
            sampler->advance();
        }
    }

    void test03_radicalInverseIncr() {
        Float x = 0.0f;

        for (int i=0; i<20; ++i) {
            assertEqualsEpsilon(x, radicalInverse(2, i), 0);
            x = radicalInverseIncremental(2, x);
        }
    }
};

MTS_EXPORT_TESTCASE(TestSamplers, "Testcase for sampling-related code")

MTS_NAMESPACE_END
