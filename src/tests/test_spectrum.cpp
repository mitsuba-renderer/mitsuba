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
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

MTS_NAMESPACE_BEGIN

class TestSpectrum : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_spectrum)
	MTS_DECLARE_TEST(test02_interpolatedSpectrum)
	MTS_DECLARE_TEST(test03_blackBody)
	MTS_END_TESTCASE()

	void test01_spectrum() {
		#if SPECTRUM_SAMPLES > 3
		Spectrum spec(0.0f);
		Float binStart = Spectrum::getBinCoverage(2).first;
		Float binEnd = Spectrum::getBinCoverage(2).second;

		spec[2] = 1.0f;
		assertEquals(spec.eval(binStart+Epsilon), 1.0f);
		assertEquals(spec.eval(binEnd-Epsilon), 1.0f);
		assertEquals(spec.eval(binStart-Epsilon), 0.0f);
		assertEquals(spec.eval(binEnd+Epsilon), 0.0f);
		#endif

		assertEqualsEpsilon(Spectrum(1.0f).getLuminance(), 1.0f, 1e-5f);

		#if SPECTRUM_SAMPLES == 30
		/* Verify against the implementation in PBRT */
		Spectrum test;
		Float r, g, b;
		test.fromLinearRGB(0.1f, 0.2f, 0.3f, Spectrum::EReflectance);
		test.toLinearRGB(r, g, b);
		assertEqualsEpsilon(r, 0.124274f, 1e-5f);
		assertEqualsEpsilon(g, 0.190133f, 1e-5f);
		assertEqualsEpsilon(b, 0.270029f, 1e-5f);
		test.fromLinearRGB(0.1f, 0.2f, 0.3f, Spectrum::EIlluminant);
		test.toLinearRGB(r, g, b);
		assertEqualsEpsilon(r, 0.0961602f, 1e-5f);
		assertEqualsEpsilon(g, 0.184446f, 1e-5f);
		assertEqualsEpsilon(b, 0.273519f, 1e-5f);
		#elif SPECTRUM_SAMPLES == 3
		Spectrum test;
		Float x, y, z;
		test.fromXYZ(0.1f, 0.2f, 0.3f);
		test.toXYZ(x, y, z);
		assertEqualsEpsilon(x, 0.1, 1e-5f);
		assertEqualsEpsilon(y, 0.2, 1e-5f);
		assertEqualsEpsilon(z, 0.3, 1e-5f);
		#endif
	}

	void test02_interpolatedSpectrum() {
		/* Test the interpolated spectrum
		   evaluation and integration routines */
		InterpolatedSpectrum spec;

		spec.append(500, 1);
		spec.append(700, 2);
		spec.append(800, 3);

		InterpolatedSpectrum spec2 = spec;

		assertEquals(spec.eval(499), 0.0f);
		assertEquals(spec.eval(500), 1.0f);
		assertEquals(spec.eval(600), 1.5f);
		assertEquals(spec.eval(700), 2.0f);
		assertEquals(spec.eval(750), 2.5f);
		assertEquals(spec.eval(800), 3.0f);
		assertEquals(spec.eval(801), 0.0f);
		assertEqualsEpsilon(spec.average(0, 1000), 550.0f/1000.0f, 1e-6f);
		assertEqualsEpsilon(spec.average(550, 551), 1.2525f, 1e-6f);

		spec2.zeroExtend();
		assertEquals(spec2.eval(500.0f-150.0f), 0.0f);
		assertEquals(spec2.eval(500.0f-150.0f/2.0f), 0.5f);

		spec2.clear();
		spec2.append(0, 1);
		spec2.append(1000, 1);
		Spectrum spec3;
		spec3.fromContinuousSpectrum(spec2);
		#if SPECTRUM_SAMPLES > 3
		assertEqualsEpsilon(spec3, Spectrum(1.0f), 1e-5f);
		#endif
		Float x, y, z;
		spec3.toXYZ(x, y, z);
		assertEqualsEpsilon(x, 1.0f, 1e-3f);
		assertEqualsEpsilon(y, 1.0f, 1e-3f);
		assertEqualsEpsilon(z, 1.0f, 1e-3f);

		#if SPECTRUM_SAMPLES > 3
		spec2.clear();
		spec2.append(Spectrum::getBinCoverage(0).first, 0);
		spec2.append(Spectrum::getBinCoverage(1).first, 1);
		spec2.append(Spectrum::getBinCoverage(2).first, 1);
		spec2.append(Spectrum::getBinCoverage(2).second, 0);
		spec3.fromContinuousSpectrum(spec2);
		Spectrum spec4(0.0f);

		spec4[0] = 0.5f;
		spec4[1] = 1.0f;
		spec4[2] = 0.5f;
		assertEqualsEpsilon(spec3, spec4, 1e-5f);
		#endif
	}

	void test03_blackBody() {
		/* Spot-check a 5000 Kelvin spectrum against values obtained at
		   https://www.sensiac.org/external/resources/calculators/infrared_radiance_calculator.jsf

		   Mitsuba uses units of W * m^-2 * sr^-1 * nm^-1, whereas that
		   page uses W * cm^-2 * sr^-1 * um^-1  --- this is the reason for
		   the factor of 10-difference.
		*/
		BlackBodySpectrum spec(5000);
		assertEqualsEpsilon(spec.eval(400)/10, 874.f, .5f);
		assertEqualsEpsilon(spec.eval(500)/10, 1211.f, .5f);
		assertEqualsEpsilon(spec.eval(600)/10, 1276.f, .5f);
		assertEqualsEpsilon(spec.eval(2000)/10, 115.8f, .5f);
		assertEqualsEpsilon(spec.average(100, 1000) * .09f, 715.f, 1);
	}
};

MTS_EXPORT_TESTCASE(TestSpectrum, "Testcase for manipulating spectral data")
MTS_NAMESPACE_END
