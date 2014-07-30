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

#ifndef MTS_TESTCASE
# define MTS_TESTCASE 1
#endif

#include <mitsuba/core/platform.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/testcase.h>

#if MTS_SSE
#include "simdtonemap.h"

#include <mitsuba/core/matrix.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/sse.h>

using mitsuba::Float;
using mitsuba::Thread;
using mitsuba::EError;


namespace
{

struct RGBA8 {
	uint8_t r;
	uint8_t g;
	uint8_t b;
	uint8_t a;

	operator mitsuba::Vector4i() const {
		return mitsuba::Vector4i(r, g, b, a);
	}
};

struct RGBA32F {
	float r;
	float g;
	float b;
	float a;

	operator RGBA8() const {
		RGBA8 p;
		float cr = std::max(0.0f, std::min(1.0f, r));
		float cg = std::max(0.0f, std::min(1.0f, g));
		float cb = std::max(0.0f, std::min(1.0f, b));
		float ca = std::max(0.0f, std::min(1.0f, a));
		p.r = static_cast<uint8_t>(floor(255.f * cr + 0.5f));
		p.g = static_cast<uint8_t>(floor(255.f * cg + 0.5f));
		p.b = static_cast<uint8_t>(floor(255.f * cb + 0.5f));
		p.a = static_cast<uint8_t>(floor(255.f * ca + 0.5f));
		return p;
	}
};


inline float toSRGB(float value) {
	if (value < 0.0031308f)
		return 12.92f * value;
	return 1.055f * pow(value, 1.0f/2.4f) - 0.055f;
}


struct TonemapGamma
{
	float invWhitePoint, invGamma;
	bool sRGB;

	RGBA32F operator() (const RGBA32F& src) const {
		RGBA32F color;
		color.r = std::max(0.0f, src.r * invWhitePoint);
		color.g = std::max(0.0f, src.g * invWhitePoint);
		color.b = std::max(0.0f, src.b * invWhitePoint);
		color.a = src.a;
		if (sRGB) {
			color.r = toSRGB(color.r);
			color.g = toSRGB(color.g);
			color.b = toSRGB(color.b);
		} else {
			color.r = pow(color.r, invGamma);
			color.g = pow(color.g, invGamma);
			color.b = pow(color.b, invGamma);
		}
		return color;
	}
};


struct TonemapReinhard
{
	typedef mitsuba::Matrix3x3 mat3;
	typedef mitsuba::Vector3f vec3;

	float scale, invWp2, invGamma, multiplier;
	bool sRGB;

	RGBA32F operator() (const RGBA32F& src) const {
		using mitsuba::Matrix3x3;
		using mitsuba::Vector3f;

		const mat3 RGB2XYZ(0.412453f, 0.212671f, 0.019334f,
						   0.357580f, 0.715160f, 0.119193f,
						   0.180423f, 0.072169f, 0.950227f);

		const mat3 XYZ2RGB(3.240479f, -0.969256f,  0.055648f,
						  -1.537150f,  1.875991f, -0.204043f,
						  -0.498535f,  0.041556f,  1.057311f);

		vec3 XYZ = RGB2XYZ * vec3(src.r, src.g, src.b) * multiplier;
		float normalization = 1.0f / (XYZ.x + XYZ.y + XYZ.z);
		float x = XYZ.x*normalization;
		float y = XYZ.y*normalization;
		float Lp = XYZ.y*scale;
		XYZ.y = Lp * (1.0f + Lp*invWp2) / (1.0f + Lp);
		float ratio = XYZ.y / y;
		XYZ.x = ratio * x;
		XYZ.z = ratio * (1.0f - x - y);
		Vector3f color = XYZ2RGB * XYZ;

		color.x = std::max(0.0f, color.x);
		color.y = std::max(0.0f, color.y);
		color.z = std::max(0.0f, color.z);

		RGBA32F frag;
		if (sRGB) {
			frag.r = toSRGB(color.x);
			frag.g = toSRGB(color.y);
			frag.b = toSRGB(color.z);
		} else {
			frag.r = pow(color.x, invGamma);
			frag.g = pow(color.y, invGamma);
			frag.b = pow(color.z, invGamma);
		}
		frag.a = src.a;
		return frag;
	}
};


template <class TMO>
void tonemap(const mitsuba::Bitmap* src, mitsuba::Bitmap* dest, const TMO& tmo){
	SAssert(src->getSize() == dest->getSize());
	SAssert(src->getComponentFormat()  == mitsuba::Bitmap::EFloat32);
	SAssert(src->getPixelFormat()      == mitsuba::Bitmap::ERGBA);
	SAssert(dest->getComponentFormat() == mitsuba::Bitmap::EUInt8);
	SAssert(dest->getPixelFormat()     == mitsuba::Bitmap::ERGBA);

	const RGBA32F* pixels = static_cast<const RGBA32F*>(src->getData());
	RGBA8* out = static_cast<RGBA8*>(dest->getData());
	for (size_t i = 0; i != src->getPixelCount(); ++i) {
		RGBA32F result = tmo(pixels[i]);
		result.r = std::max(0.0f, std::min(1.0f, result.r));
		result.g = std::max(0.0f, std::min(1.0f, result.g));
		result.b = std::max(0.0f, std::min(1.0f, result.b));
		result.a = std::max(0.0f, std::min(1.0f, result.a));
		out[i] = result;
	}
}

void luminanceInfo(const mitsuba::Bitmap* src, const float multiplier,
	float& outMaxLuminance, float& avgLogLuminance) {
	SAssert(src->getComponentFormat()  == mitsuba::Bitmap::EFloat32);
	SAssert(src->getPixelFormat()      == mitsuba::Bitmap::ERGBA);
	const RGBA32F* pixels = static_cast<const RGBA32F*>(src->getData());

	outMaxLuminance = -1.0f;
	double sumLogLuminance = 0.0;
	for (size_t i = 0; i != src->getPixelCount(); ++i) {
		const float& r = pixels[i].r;
		const float& g = pixels[i].g;
		const float& b = pixels[i].b;
		float Y = multiplier * (r * 0.212671f + g * 0.715160f + b * 0.072169f);
		// catch NaNs, negatives, and the Mitsuba banner
		if (Y < 0.0 || Y != Y || Y == 1024.0f) {
			Y = 0.0f;
		}
		float logLuminance = log(1e-3f + Y);
		outMaxLuminance = std::max(Y, outMaxLuminance);
		sumLogLuminance += logLuminance;
	}
	avgLogLuminance = (float) exp(sumLogLuminance / src->getPixelCount());
}


// Online variance calulation by Knuth, referenced by Wikipedia [August 2012]
//   http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
//   Donald E. Knuth (1998). The Art of Computer Programming, volume 2:
//   Seminumerical Algorithms, 3rd edn., p. 232. Boston: Addison-Wesley.
class VarianceFunctor
{
public:
	VarianceFunctor() : m_n(0), m_mean(0.0), m_M2(0.0),
		m_max(-std::numeric_limits<double>::infinity()) {}

	inline void update(double x) {
		++m_n;
		double delta = x - m_mean;
		m_mean += delta / m_n;
		m_M2   += delta * (x - m_mean);
		m_max = fmaxd(m_max, x);
	}

	inline double mean() const {
		return m_mean;
	}

	inline double variance() const {
		assert (m_n > 1);
		return m_M2 / (m_n - 1);
	}

	inline double stddev() const {
		return sqrt(variance());
	}

	inline double max() const {
		return m_max;
	}

	void reset() {
		m_n    = 0;
		m_mean = 0.0;
		m_M2   = 0.0;
		m_max  = -std::numeric_limits<double>::infinity();
	}

private:

	inline static double fmaxd(double x, double y) {
		__m128d tmp = _mm_max_sd(_mm_load_sd(&x), _mm_load_sd(&y));
		return _mm_cvtsd_f64(tmp);
	}

	size_t m_n;
	double m_mean;
	double m_M2;
	double m_max;
};

} // namespace


MTS_NAMESPACE_BEGIN


class TestTonemapperSSE : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(testBasic);
	MTS_DECLARE_TEST(testLuminance);
	MTS_DECLARE_TEST(testGamma);
	MTS_DECLARE_TEST(testReinhard);
	MTS_END_TESTCASE()

	virtual void init();

	void testBasic();

	void testLuminance(int numImageSizes = 5, int runsPerImage = 2);

	void testGamma(int numImageSizes = 5, int runsPerImage = 2,
		int paramsPerImage = 10);

	void testReinhard(int numImageSizes = 5, int runsPerImage = 2,
		int paramsPerImage = 10);

private:

	static inline float pow2(float x) {
		return math::fastexp(x * 0.69314718f);
	}

	/// Generate a new size for a bitmap, including odd-sizes once every 4 times
	inline Vector2i nextSize() {
		uint32_t w = m_rnd->nextUInt(2048) + 1;
		uint32_t h = m_rnd->nextUInt(2048) + 1;
		if (m_rnd->nextFloat() > 0.75f) {
			w |= 1;
			h |= 1;
		}
		Vector2i s(static_cast<int>(w), static_cast<int>(h));
		return s;
	}

	inline float getHDRValue() {
		// Interpret as an f-stop on an image with average luminance 0.18
		float fstop = static_cast<float>(m_rnd->nextStandardNormal());
		// Sometimes add really bright or dim stuff
		if (m_rnd->nextFloat() > 0.98f) {
			fstop *= 8.0f;
		}
		float factor = pow2(fstop);
		return 0.18f * factor;
	}

	/// Fill a HDR bitmap with plausible values. Returns the maximum value
	inline float fill(Bitmap *hdr) {
		SAssert(hdr->getComponentFormat() == Bitmap::EFloat32);
		SAssert(hdr->getPixelFormat()     == Bitmap::ERGBA);
		float *d = hdr->getFloat32Data();
		float maxValue = 0;
		for(size_t i = 0; i != hdr->getPixelCount(); ++i, d += 4) {
			d[0] = getHDRValue();
			d[1] = getHDRValue();
			d[2] = getHDRValue();
			d[3] = m_rnd->nextFloat();
			maxValue = std::max(maxValue, std::max(d[0], std::max(d[1], d[2])));
		}
		return maxValue;
	}

	/// Update the variance functor by the per-component difference
	void updateVariance(const Bitmap *actual, const Bitmap *expected);

	/// Basic single pixel test
	void testBasic(const RGBA32F &pixel, const TonemapCPU::Params &params,
		const Vector4f &delta);

	/**
	 * \brief Fill the HDR image with random values \c numImages times. Each
	 * of those times choose \c numRuns random tonemapping parameters.
	 */
	void testGamma(Bitmap *hdr, int numImages, int numRuns);

	/**
	 * \brief Fill the HDR image with random values \c numImages times. Each
	 * of those times choose \c numRuns random tonemapping parameters.
	 */
	void testReinhard(Bitmap *hdr, int numImages, int numRuns);

	/// Single instance of the luminance test
	void testLuminance(Bitmap *hdr, TonemapCPU *tmo);

	ref<Random> m_rnd;
	ref<Timer> m_timer;
	ref<Timer> m_timerRef;
	ref<Timer> m_timerSIMD;
	VarianceFunctor m_variance;
	VarianceFunctor m_varMaxLum;
	VarianceFunctor m_varAvgLogLum;
};


void TestTonemapperSSE::init() {
	m_rnd = new mitsuba::Random(0xdafa06719c2c08e5ULL);
	m_timer = new mitsuba::Timer(/*autostart*/ false);
	m_timerRef  = new mitsuba::Timer(/*autostart*/ false);
	m_timerSIMD = new mitsuba::Timer(/*autostart*/ false);
}

void TestTonemapperSSE::testGamma(Bitmap *hdr, int numImages, int numRuns) {
	const mitsuba::Vector2i& size = hdr->getSize();
	ref<Bitmap> ldr    = new Bitmap(Bitmap::ERGBA, Bitmap::EUInt8, size);
	ref<Bitmap> ldrRef = new Bitmap(Bitmap::ERGBA, Bitmap::EUInt8, size);

	ref<TonemapCPU> tmoSIMD = new TonemapCPU;
	TonemapGamma tmo;

	for (int imgIdx = 0; imgIdx < numImages; ++imgIdx) {
		const float whitePoint = fill(hdr);
		for (int runIdx = 0; runIdx < numRuns; ++runIdx) {
			// Random set of parameters for the current run
			const bool sRGB = m_rnd->nextFloat() < 0.5f/8.0f;
			const Float invGamma = 1.0f / std::max(0.1f,
				m_rnd->nextStandardNormal() + 2.2f);
			float wFStop = m_rnd->nextStandardNormal();
			float wFactor = pow2(wFStop);
			const Float invWhitePoint = 1 / (whitePoint * wFactor);

			tmoSIMD->setInvGamma(invGamma);
			tmoSIMD->setSRGB(sRGB);
			tmoSIMD->setInvWhitePoint(invWhitePoint);

			tmo.invGamma = invGamma;
			tmo.sRGB = sRGB;
			tmo.invWhitePoint = invWhitePoint;

			m_timerRef->start();
			tonemap(hdr, ldrRef, tmo);
			m_timerRef->stop();

			m_timerSIMD->start();
			tmoSIMD->gammaTonemap(hdr, ldr);
			m_timerSIMD->stop();

			updateVariance(ldr, ldrRef);
		}
	}
}

void TestTonemapperSSE::testReinhard(Bitmap *hdr, int numImages, int numRuns) {
	const mitsuba::Vector2i& size = hdr->getSize();
	ref<Bitmap> ldr    = new Bitmap(Bitmap::ERGBA, Bitmap::EUInt8, size);
	ref<Bitmap> ldrRef = new Bitmap(Bitmap::ERGBA, Bitmap::EUInt8, size);

	ref<TonemapCPU> tmoSIMD = new TonemapCPU;
	TonemapReinhard tmo;

	for (int imgIdx = 0; imgIdx < numImages; ++imgIdx) {
		const float whitePoint = fill(hdr);
		for (int runIdx = 0; runIdx < numRuns; ++runIdx) {
			// Random set of parameters for the current run
			const bool sRGB = m_rnd->nextFloat() < 0.5f/8.0f;
			const Float invGamma = 1.0f / std::max(0.1f,
				m_rnd->nextStandardNormal() + 2.2f);
			float wFStop = m_rnd->nextStandardNormal();
			float wFactor = pow2(wFStop);
			const Float invWhitePoint = 1 / (whitePoint * wFactor);
			const Float multiplier = m_rnd->nextFloat() + 0.5f;
			const Float scale = m_rnd->nextFloat() * 4.0f + 0.01f;

			tmoSIMD->setInvGamma(invGamma);
			tmoSIMD->setSRGB(sRGB);
			tmoSIMD->setInvWhitePoint(invWhitePoint);
			tmoSIMD->setMultiplier(multiplier);
			tmoSIMD->setScale(scale);

			tmo.invGamma = invGamma;
			tmo.sRGB = sRGB;
			tmo.invWp2 = invWhitePoint * invWhitePoint;
			tmo.multiplier = multiplier;
			tmo.scale = scale;

			m_timerRef->start();
			tonemap(hdr, ldrRef, tmo);
			m_timerRef->stop();

			m_timerSIMD->start();
			tmoSIMD->reinhardTonemap(hdr, ldr);
			m_timerSIMD->stop();

			updateVariance(ldr, ldrRef);
		}
	}
}

void TestTonemapperSSE::updateVariance(const Bitmap *actual,
	const Bitmap *expected) {
	SAssert(actual->getSize() == expected->getSize());

	const RGBA8* dActual   = static_cast<const RGBA8*>(actual->getData());
	const RGBA8* dExpected = static_cast<const RGBA8*>(expected->getData());
	const size_t numPixels = actual->getPixelCount();

	for (size_t i = 0; i != numPixels; ++i) {
		Vector4i a = static_cast<Vector4i>(dActual[i]);
		Vector4i b = static_cast<Vector4i>(dExpected[i]);
		Vector4i diff = a - b;
		for (int k = 0; k < 4; ++k) {
			int delta = std::abs(diff[k]);
			m_variance.update(delta);
		}
	}
}

void TestTonemapperSSE::testLuminance(Bitmap *hdr, TonemapCPU *tmo) {
	float maxAvgLum, avgLogLum;

	m_timerRef->start();
	luminanceInfo(hdr, tmo->multiplier(), maxAvgLum, avgLogLum);
	m_timerRef->stop();

	m_timerSIMD->start();
	tmo->setLuminanceInfo(hdr, tmo->multiplier());
	m_timerSIMD->stop();

	float errMax = std::abs(maxAvgLum - tmo->maxLuminance());
	float errAvg = std::abs(avgLogLum - tmo->logAvgLuminance());
	m_varMaxLum.update(errMax);
	m_varAvgLogLum.update(errAvg);
}

void TestTonemapperSSE::testGamma(int numImageSizes, int runsPerImage,
	int paramsPerImage) {
	m_timerRef->reset(false);
	m_timerSIMD->reset(false);
	m_variance.reset();

	m_timer->reset(true);
	for (int i = 0; i < numImageSizes; ++i) {
		mitsuba::Vector2i size = nextSize();
		ref<Bitmap> hdr = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, size);
		testGamma(hdr, runsPerImage, paramsPerImage);
	}
	Float seconds = m_timer->stop();
	Log(EInfo, "Gamma test done: %g s, reference: %g s, SIMD: %g s",
		seconds, m_timerRef->getSeconds(), m_timerSIMD->getSeconds());
	Log(EInfo, "Error - mean: %g, stddev: %g, max: %g",
		m_variance.mean(), m_variance.stddev(), m_variance.max());

	assertTrue(m_variance.mean()   < 1e-5);
	assertTrue(m_variance.stddev() < 1e-2);
	assertTrue(m_variance.max()    < 2);
}

void TestTonemapperSSE::testReinhard(int numImageSizes, int runsPerImage,
	int paramsPerImage) {
	m_timerRef->reset(false);
	m_timerSIMD->reset(false);
	m_variance.reset();

	m_timer->reset(true);
	for (int i = 0; i < numImageSizes; ++i) {
		mitsuba::Vector2i size = nextSize();
		ref<Bitmap> hdr = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, size);
		testReinhard(hdr, runsPerImage, paramsPerImage);
	}
	Float seconds = m_timer->stop();
	Log(EInfo, "Reinhard test done: %g s, reference: %g s, SIMD: %g s",
		seconds, m_timerRef->getSeconds(), m_timerSIMD->getSeconds());
	Log(EInfo, " Error - mean: %g, stddev: %g, max: %g",
		m_variance.mean(), m_variance.stddev(), m_variance.max());

	assertTrue(m_variance.mean()   < 1e-2);
	assertTrue(m_variance.stddev() < 1e-1);
}

void TestTonemapperSSE::testLuminance(int numImageSizes,int runsPerImage) {
	m_timerRef->reset(false);
	m_timerSIMD->reset(false);
	m_varMaxLum.reset();
	m_varAvgLogLum.reset();

	ref<TonemapCPU> tmo = new TonemapCPU;
	m_timer->reset(true);
	for (int i = 0; i < numImageSizes; ++i) {
		mitsuba::Vector2i size = nextSize();
		ref<Bitmap> hdr = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, size);
		for (int j = 0; j < runsPerImage; ++j) {
			fill(hdr);
			tmo->setMultiplier(m_rnd->nextFloat() + 0.5f);
			testLuminance(hdr, tmo);
		}
	}
	Float seconds = m_timer->stop();
	Log(EInfo, "Luminance info test done: %g s, reference: %g s, SIMD: %g s",
		seconds, m_timerRef->getSeconds(), m_timerSIMD->getSeconds());
	Log(EInfo, "Max lum error - mean: %g, stddev: %g, max: %g",
		m_varMaxLum.mean(), m_varMaxLum.stddev(), m_varMaxLum.max());
	Log(EInfo, "Avg log-lum error - mean: %g, stddev: %g, max: %g",
		m_varAvgLogLum.mean(), m_varAvgLogLum.stddev(), m_varAvgLogLum.max());

	assertEquals(static_cast<Float>(m_varMaxLum.max()), static_cast<Float>(0));
	assertTrue(m_varAvgLogLum.mean()   < 1e-5);
	assertTrue(m_varAvgLogLum.stddev() < 1e-6);
}

void TestTonemapperSSE::testBasic(const RGBA32F &pixel,
	const TonemapCPU::Params &params, const Vector4f &delta) {
	// Setup the reference tonemapper
	TonemapReinhard tmo;
	tmo.invWp2     = params.invWhitePoint * params.invWhitePoint;
	tmo.invGamma   = params.invGamma;
	tmo.multiplier = params.multiplier;
	tmo.scale = params.scale;
	tmo.sRGB  = params.isSRGB;

	// Setup the SIMD tonemapper
	ref<TonemapCPU> tmoSIMD = new TonemapCPU;
	tmoSIMD->setInvWhitePoint(params.invWhitePoint);
	tmoSIMD->setInvGamma(params.invGamma);
	tmoSIMD->setMultiplier(params.multiplier);
	tmoSIMD->setScale(params.scale);
	tmoSIMD->setSRGB(params.isSRGB);

	// Basic sanity check of the reference implementation
	const RGBA32F expectedFloat = tmo(pixel);
	const RGBA8 expected = expectedFloat;
	const Vector4i v4i = expected;
	assertEquals(v4i[0], expected.r);
	assertEquals(v4i[1], expected.g);
	assertEquals(v4i[2], expected.b);
	assertEquals(v4i[3], expected.a);

	// Populate and execute the SIMD tonemapper
	ref<Bitmap> hdr = new Bitmap(Bitmap::ERGBA,Bitmap::EFloat32, Vector2i(1,1));
	ref<Bitmap> ldr = new Bitmap(Bitmap::ERGBA,Bitmap::EUInt8,   Vector2i(1,1));
	RGBA32F* data = static_cast<RGBA32F*>(hdr->getData());
	assertEquals(reinterpret_cast<uintptr_t>(data) % 16, 0);
	*data = pixel;
	tmoSIMD->reinhardTonemap(hdr, ldr);

	// Compare the results
	const RGBA8* dataLDR = static_cast<const RGBA8*>(ldr->getData());
	assertEquals(reinterpret_cast<uintptr_t>(dataLDR) % 16, 0);
	assertEqualsEpsilon(dataLDR->r, expected.r, delta[0]);
	assertEqualsEpsilon(dataLDR->g, expected.g, delta[1]);
	assertEqualsEpsilon(dataLDR->b, expected.b, delta[2]);
	assertEqualsEpsilon(dataLDR->a, expected.a, delta[3]);

	// Luminance info
	float Y = params.multiplier * (pixel.r * 0.212671f +
								   pixel.g * 0.715160f +
								   pixel.b * 0.072169f);
	float avgLogLuminance = math::fastexp(math::fastlog(1e-3f + Y));
	tmoSIMD->setLuminanceInfo(hdr, tmo.multiplier);
	assertEquals(tmoSIMD->maxLuminance(), Y);
	// Relative epsilon based on qFuzzyCompare
	const Float epsilon = 0.0001f * std::min(std::abs(avgLogLuminance),
		std::abs(tmoSIMD->logAvgLuminance()));
	assertEqualsEpsilon(tmoSIMD->logAvgLuminance(), avgLogLuminance, epsilon);
}

void TestTonemapperSSE::testBasic() {

	TonemapCPU::Params params;
	RGBA32F pixel;
	Vector4f delta(static_cast<Float>(0));

	params.invWhitePoint = 8.7875755e-011f;
	params.invGamma      = 0.69960159f;
	params.multiplier    = 1.2126911f;
	params.scale         = 0.69272733f;
	params.isSRGB        = true;
	pixel.r = 1.8057191f;
	pixel.g = 0.14349695f;
	pixel.b = 0.12698999f;
	pixel.a = 0.37581384f;
	testBasic(pixel, params, delta);

	params.invWhitePoint = 8.7875755e-011f;
	params.invGamma      = 0.90543908f;
	params.multiplier    = 1.3749540f;
	params.scale         = 0.45611405f;
	params.isSRGB        = false;
	pixel.r = 0.030551272f;
	pixel.g = 0.20581213f;
	pixel.b = 595144.75f;
	pixel.a = 0.48818684f;
	testBasic(pixel, params, delta);

	// In this case, the reference method generates pixels with zero value
	// because the transform from RGB->XYZ->RGB yields negative values
	delta[0] = 1;
	delta[1] = 5;
	params.invWhitePoint = 8.7875755e-011f;
	params.invGamma = 0.32298070f;
	params.multiplier = 1.4785858f;
	params.scale = 3.5021718f;
	params.isSRGB = false;
	pixel.r = 0.15293880f;
	pixel.g = 0.10678958f;
	pixel.b = 172735.61f;
	pixel.a = 0.82842076f;
	testBasic(pixel, params, delta);
}

#else

MTS_NAMESPACE_BEGIN

class TestTonemapperSSE : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(testNoSSE);
	MTS_END_TESTCASE()

	void testNoSSE() {
		Log(EError, "This test requries SSE support.");
	}
};

#endif /* MTS_SSE */


MTS_EXPORT_TESTCASE(TestTonemapperSSE,
	"Testcase for the SSE Tonemapper used by mtsgui")

MTS_NAMESPACE_END
