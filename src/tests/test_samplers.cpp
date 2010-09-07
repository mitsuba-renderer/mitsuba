#include <mitsuba/render/testcase.h>
#include <mitsuba/core/plugin.h>

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
				createObject(Sampler::m_theClass, Properties("halton")));

		/* MATLAB: p = haltonset(5); net(p,5) */
		Float comparison[] = {
							0,                 0,                 0,                 0,                 0,
			0.500000000000000, 0.333333333333333, 0.200000000000000, 0.142857142857143, 0.090909090909091,
			0.250000000000000, 0.666666666666667, 0.400000000000000, 0.285714285714286, 0.181818181818182,
			0.750000000000000, 0.111111111111111, 0.600000000000000, 0.428571428571429, 0.272727272727273,
			0.125000000000000, 0.444444444444444, 0.800000000000000, 0.571428571428571, 0.363636363636364
		};

		int pos = 0;
		sampler->generate();
		for (int i=0; i<5; ++i) {
			for (int j=0; j<5; ++j)
				Assert(std::abs(sampler->next1D() - comparison[pos++]) < 1e-7);
			sampler->advance();
		}
	}

	void test02_Hammersley() {
		Properties props("hammersley");
		props.setInteger("sampleCount", 5);

		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(Sampler::m_theClass, props));

		Float comparison[] = {
			0.0,						0,                 0,                 0,                 0,                 0,
			1.0/5.0,	0.500000000000000, 0.333333333333333, 0.200000000000000, 0.142857142857143, 0.090909090909091,
			2.0/5.0,	0.250000000000000, 0.666666666666667, 0.400000000000000, 0.285714285714286, 0.181818181818182,
			3.0/5.0,	0.750000000000000, 0.111111111111111, 0.600000000000000, 0.428571428571429, 0.272727272727273,
			4.0/5.0,	0.125000000000000, 0.444444444444444, 0.800000000000000, 0.571428571428571, 0.363636363636364
		};

		int pos = 0;
		sampler->generate();
		for (int i=0; i<5; ++i) {
			for (int j=0; j<6; ++j) 
				Assert(std::abs(sampler->next1D() - comparison[pos++]) < 1e-7);
			sampler->advance();
		}
	}

	void test03_radicalInverseIncr() {
		Float x = 0.0f;

		for (int i=0; i<20; ++i) {
			Assert(x == radicalInverse(2, i));
			x = radicalInverseIncremental(2, x);
		}
	}
};

MTS_EXPORT_TESTCASE(TestSamplers, "Testcase for Sampler implementations")

MTS_NAMESPACE_END
