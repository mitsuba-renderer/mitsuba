#include <mitsuba/mitsuba.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/photonmap.h>
#include <mitsuba/core/wavelet.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/grid.h>
#include "../medium/maxexp.h"
#include <fstream>

using namespace mitsuba;

void testHalton() {
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
			SAssert(std::abs(sampler->next1D() - comparison[pos++]) < 1e-7);
		sampler->advance();
	}
}

void testHammersley() {
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
			SAssert(std::abs(sampler->next1D() - comparison[pos++]) < 1e-7);
		sampler->advance();
	}
}

void testRadicalInverseIncr() {
	Float x = 0.0f;

	for (int i=0; i<20; ++i) {
		SAssert(x == radicalInverse(2, i));
		x = radicalInverseIncremental(2, x);
	}
}

void testSpotLuminaire() {
	Properties props("spot");
	props.setFloat("beamWidth", 20);
	props.setFloat("cutoffAngle", 25);
	ref<Luminaire> spot = static_cast<Luminaire *> (PluginManager::getInstance()->
		createObject(Luminaire::m_theClass, props));
	ref<Random> random = new Random();

	Spectrum power;
	EmissionRecord eRec;
	const int nSamples = 1000000;
	for (int i=0; i<nSamples; ++i) {
		Point2 areaSample(random->nextFloat(), random->nextFloat());
		Point2 dirSample(random->nextFloat(), random->nextFloat());
		spot->sampleEmission(eRec, areaSample, dirSample);
		power += eRec.P / (eRec.pdfArea*eRec.pdfDir);
	}
	power /= (Float) nSamples;
	cout << "Estimated power: " << power.toString() << endl;
	cout << "Actual power: " << spot->getPower().toString() << endl;
}

void testHG1() {
	Properties props("hg");
	props.setFloat("g", .8);

	ref<PhaseFunction> hg = static_cast<PhaseFunction *> (PluginManager::getInstance()->
		createObject(PhaseFunction::m_theClass, props));
	ref<Random> random = new Random();

	Float sum = 0;
	Vector dir1(0, 0, 1);
	const int nSamples = 100000;
	MediumSamplingRecord mRec;
	for (int i=0; i<nSamples; ++i) {
		Vector dir2(squareToSphere(Point2(random->nextFloat(), random->nextFloat())));
		sum += hg->f(mRec, dir2, dir2)[0];
	}
	sum *= (4 * M_PI) / nSamples;

	cout << sum << endl;
}

void testHG2() {
	const int res = 8;
	Float buf[res][res];
	Properties props("hg");
	props.setFloat("g", -.2);
	MediumSamplingRecord mRec;

	ref<PhaseFunction> hg = static_cast<PhaseFunction *> (PluginManager::getInstance()->
		createObject(PhaseFunction::m_theClass, props));
	ref<Random> random = new Random();

	const Vector dir1(0, 0, 1);
	memset(buf, 0, sizeof(Float)*res*res);
	Vector dir2;
	PhaseFunction::ESampledType sampledType;

	const int nSamples = 10000000;
	for (int i=0; i<nSamples; ++i) {
		Point2 sample(random->nextFloat(), random->nextFloat());
		hg->sample(mRec, dir1, dir2, sampledType, sample);
		Float pdf = hg->f(mRec, dir1, dir2)[0];
		Float dir2Phi = (std::atan2(dir2.y, dir2.x) + (Float) M_PI) / (2 * (Float) M_PI);
		Float dir2CosTheta = (dir2.z + 1) / (Float) 2;

		int pos1 = std::max(0, std::min((int) (dir2Phi * res), res-1));
		int pos2 = std::max(0, std::min((int) (dir2CosTheta* res), res-1));
		buf[pos1][pos2] += 1/pdf;

	}
	cout << "A=[";
	for (int y=0; y<res; ++y)  {
		for (int x=0; x<res; ++x) 
			cout << buf[y][x] / nSamples << " ";
		cout << ";" << endl;
	}
	cout << "];" << endl;
}

void testMATLAB() {
	int count = 256;
	int dim1 = 17, dim2 = 19;
	cout << "x=[";
	for (int i=1; i<=count; ++i) {
		cout << radicalInverse(dim1, i);
		if (i+1<=count)
			cout << " ";
	}
	cout << "];" << endl;
	cout << "y=[";
	for (int i=1; i<=count; ++i) {
		cout << radicalInverse(dim2, i);
		if (i+1<=count)
			cout << " ";
	}
	cout << "];" << endl;
}

Float sample(Float *sigma, int n, Float u, Float &prob) {
	SAssert(n >= 2);
	Float cdf[4], intervalStart[3];
	cdf[0] = 0;

	/* Requires the coefficients to be sorted in descending order */
	for (int i=0; i<n; ++i) {
		/* Integrate max(f_1(t), .., f_n(t)) on [0, \infty]*/
		Float lower = (i==0) ? -1 : -std::pow((sigma[i]/sigma[i-1]), 
					-sigma[i] / (sigma[i]-sigma[i-1]));
		Float upper = (i==n-1) ? 0 : -std::pow((sigma[i+1]/sigma[i]), 
					-sigma[i] / (sigma[i+1]-sigma[i]));
		cdf[i+1] = cdf[i] + (upper - lower);

		/* Store the interval covered by each f_i */
		intervalStart[i] = (i == 0) ? 0
			: std::log(sigma[i]/sigma[i-1]) / (sigma[i]-sigma[i-1]);
	}

	/* Turn into a discrete CDF and keep the normalization factor */
	Float normFactor = cdf[n];
	for (int i=0; i<=n; ++i) 
		cdf[i] /= normFactor;

	/* Find the f_i for this sample */
	Float *lowerBound = std::lower_bound(&cdf[0], &cdf[n+1], u);
	int index = std::max(0, (int) (lowerBound - &cdf[0]) - 1);
	SAssert(index >= 0 && index < n);
	
	/* Sample according to f_i */
	Float t = -std::log(std::exp(-intervalStart[index] * sigma[index]) 
		- normFactor * (u - cdf[index])) / sigma[index];

	/* Compute the probability of this sample */
	prob = sigma[index] * std::exp(-sigma[index] * t) / normFactor;

	/* CDF computation
		Float *lowerBound = std::lower_bound(&intervalStart[0], &intervalStart[n], t);
		int index = std::max(0, (int) (lowerBound - &intervalStart[0]) - 1);
		Float lower = (index==0) ? -1 : -std::pow((sigma[index]/sigma[index-1]), 
					-sigma[index] / (sigma[index]-sigma[index-1]));
		Float upper = -std::exp(-sigma[index] * t);
		Float integral = cdf[index] + (upper - lower) / normFactor;
	*/

	return t;
}

void testVariance() {
	ref<Random> random = new Random();
	Spectrum mean, meanSqr, variance;
	Spectrum sigmaT;
	sigmaT[0] = 0.7014;
	sigmaT[1] = 1.2225;
	sigmaT[2] = 1.9142;
	int nSamples = 100000;

	for (int chan=0; chan<3; ++chan) {
		cout << "Sampling with respect to channel " << chan << endl;
		mean = Spectrum(0.0f); meanSqr = Spectrum(0.0f);
		for (int i=0; i<nSamples; ++i) {
			Float desiredAttenuation = 1-random->nextFloat();

			Float t = -std::log(desiredAttenuation)/sigmaT[chan];
			Float prob = sigmaT[chan] * std::exp(-sigmaT[chan] * t);
			Spectrum value = (sigmaT * (-t)).exp() / prob;

			Spectrum delta = value - mean;
			mean += delta / (Float) (i+1);
			meanSqr += delta * (value - mean);
			variance = meanSqr / (Float) i;
		}
		cout << "Expectation     : " << (Spectrum(1)/sigmaT).toString() << endl;
		cout << "Mean            : " << mean.toString() << endl;
		cout << "Sample variance : " << variance.toString() << endl << endl;
	}

	cout << "Random sampling test" << endl;
	mean = Spectrum(0.0f); meanSqr = Spectrum(0.0f);
	for (int i=0; i<nSamples; ++i) {
		Float desiredAttenuation = 1-random->nextFloat();
		int chan = random->nextInteger(3);

		Float t = -std::log(desiredAttenuation)/sigmaT[chan];
		Spectrum prob = sigmaT * (-sigmaT * t).exp();
		Spectrum value = (sigmaT * (-t)).exp() / (prob[chan]);

		Spectrum delta = value - mean;
		mean += delta / (Float) (i+1);
		meanSqr += delta * (value - mean);
		variance = meanSqr / (Float) i;
	}
	cout << "Expectation     : " << (Spectrum(1)/sigmaT).toString() << endl;
	cout << "Mean            : " << mean.toString() << endl;
	cout << "Sample variance : " << variance.toString() << endl << endl;
	
	cout << "Single sample test" << endl;
	mean = Spectrum(0.0f); meanSqr = Spectrum(0.0f);
	for (int i=0; i<nSamples; ++i) {
		Float desiredAttenuation = 1-random->nextFloat();
		int chan = random->nextInteger(3);

		Float t = -std::log(desiredAttenuation)/sigmaT[chan];
		Spectrum prob = sigmaT * (-sigmaT * t).exp();
		Float weight = prob[chan] / (prob[0] + prob[1] + prob[2]);
		Spectrum value = (sigmaT * (-t)).exp() / (prob[chan] * 1/(Float) 3) * weight;

		Spectrum delta = value - mean;
		mean += delta / (Float) (i+1);
		meanSqr += delta * (value - mean);
		variance = meanSqr / (Float) i;
	}
	cout << "Expectation     : " << (Spectrum(1)/sigmaT).toString() << endl;
	cout << "Mean            : " << mean.toString() << endl;
	cout << "Sample variance : " << variance.toString() << endl << endl;


	cout << "New distribution test" << endl;
	mean = Spectrum(0.0f); meanSqr = Spectrum(0.0f);
	Float extinction[] = {1.9142, 1.2225, 0.7014};
	for (int i=0; i<nSamples; ++i) {
		Float prob;
		Float t = sample(extinction, 3, random->nextFloat(), prob);

		Spectrum value = (sigmaT * (-t)).exp() / prob;

		Spectrum delta = value - mean;
		mean += delta / (Float) (i+1);
		meanSqr += delta * (value - mean);
		variance = meanSqr / (Float) i;
	}

	cout << "Expectation     : " << (Spectrum(1)/sigmaT).toString() << endl;
	cout << "Mean            : " << mean.toString() << endl;
	cout << "Sample variance : " << variance.toString() << endl << endl;
}

void testPhotonMap() {
	int nPhotons = 3000000;
	ref<PhotonMap> map = new PhotonMap(nPhotons);
	Spectrum power;
	ref<Random> random = new Random();
	power.fromLinearRGB(1, 2, 3);
	power *= M_PI;

	bool storePhoton(const Point &pos, const Vector &dir, const Spectrum &power);
	for (int i=0; i<nPhotons; ++i) {
		Point2 disk = squareToDisk(Point2(random->nextFloat(), random->nextFloat()));
		map->storePhoton(Point(disk.x,disk.y,0), Normal(0,0,1), Vector(0, 0, -1), power, 1);
	}

	map->setScale(1/(Float) nPhotons);
	map->balance();

	for (int i=1; i<10; ++i) {
		Point2 disk = squareToDisk(Point2(random->nextFloat(), random->nextFloat())) * .2;

		cout << map->estimateIrradiance(Point(disk.x, disk.y, 0), Normal(0, 0, 1), .3, 10000).toString() << endl;
		cout << map->estimateIrradianceFiltered(Point(disk.x, disk.y, 0), Normal(0, 0, 1), .3, 10000).toString() << endl;
	}
}

void testMaxExp() {
	std::vector<Float> sigmaT;
	sigmaT.push_back(0.7014);
	sigmaT.push_back(1.2225);
	sigmaT.push_back(1.9142);

	MaxExpDist dist(sigmaT);

	for (int i=0; i<9; ++i) {
		Float U=.01f + i/10.0f, t, pdf, cdf;
		t = dist.sample(U, pdf);
		cdf = dist.cdf(t);

		cout << "Sampled U=" << U << " => t = " << t << ", pdf=" << pdf << ", cdf=" << cdf << endl;
	}
}

void testHeterogeneous() {
	ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
			createObject(Sampler::m_theClass, Properties("independent")));

	Properties props("heterogeneous");
	props.setString("filename", "scenes/cornell/smoke-density.65.vol");
	props.setSpectrum("sigmaS", Spectrum(2.62));
	props.setSpectrum("sigmaA", Spectrum(.05f));
	props.setString("strategy", "standard");
	props.setFloat("sizeMultiplier", 10);
	ref<Medium> standardMedium = static_cast<Medium *> (PluginManager::getInstance()->
			createObject(Medium::m_theClass, props));
	standardMedium->configure();

	props = Properties("heterogeneous");
	props.setString("filename", "scenes/cornell/smoke-density.65.vol");
	props.setString("strategy", "coleman");
	props.setSpectrum("sigmaS", Spectrum(2.62));
	props.setSpectrum("sigmaA", Spectrum(.05f));
	props.setFloat("sizeMultiplier", 10);
	ref<Medium> colemanMedium = static_cast<Medium *> (PluginManager::getInstance()->
			createObject(Medium::m_theClass, props));
	colemanMedium->configure();

	Ray ray(Point(0, 1, .4), Vector(1, 0, 0));

	MediumSamplingRecord mRec;

	std::ofstream distr("distr.m");
	int res = 300000, failures = 0;

	distr << "samplesStandard=[" << endl;
	for (int i=0; i<res; ++i) {
		if (standardMedium->sampleDistance(ray, std::numeric_limits<Float>::infinity(), mRec, sampler)) {
//			distr << mRec.t << " ";
			distr << mRec.sigmaS[0] * mRec.attenuation[0] / mRec.pdf << " " << endl;
		} else {
			++failures;
		}
	}
	distr << "];" << endl;
	cout << failures/(Float)res << endl;
	failures = 0;
	cout << failures/(Float)res << endl;
	distr << "samplesColeman=[" << endl;
	for (int i=0; i<res; ++i) {
		if (colemanMedium->sampleDistance(ray, std::numeric_limits<Float>::infinity(), mRec, sampler)) {
//			distr << mRec.t << " ";
			distr << mRec.sigmaS[0] * mRec.attenuation[0] / mRec.pdf << " " << endl;
		} else {
			++failures;
		}
	}
	distr << "];" << endl;
	distr << "hist(samplesStandard, 300);" << endl;
	distr << "title('Standard');" << endl;
	distr << "figure;" << endl;
	distr << "hist(samplesColeman, 300);" << endl;
	distr.close();
	cout << failures/(Float)res << endl;
}

void testWavelet() {
	ref<FileStream> stream = new FileStream("cat.png", FileStream::EReadOnly);
	ref<Bitmap> bitmap = new Bitmap(Bitmap::EPNG, stream);

	ref<Wavelet2D> wavelet = new Wavelet2D(bitmap);
	cout << "2D compression ratio: " << wavelet->compress(.015) << endl;
	ref<SparseWavelet2D> sw = wavelet->toSparseWavelet();
	ref<FileStream> stream3 = new FileStream("cat2.raw", FileStream::ETruncReadWrite);
	sw->serialize(stream3, NULL);
	stream3->close();
	ref<FileStream> stream4 = new FileStream("cat2.raw", FileStream::EReadOnly);
	ref<SparseWavelet2D> sw2 = new SparseWavelet2D(stream4, NULL);
	ref<Wavelet2D> wavelet2 = new Wavelet2D(sw2);
	ref<Bitmap> bitmap2 = new Bitmap(bitmap->getWidth(), bitmap->getHeight(), 8);

	wavelet2->decode(bitmap2);

	size_t nEntries = 512*512;
	Float mean = 0, sqrError = 0;
	for (size_t i=0; i<nEntries; ++i) {
		mean += bitmap->getData()[i] / (255.0f * nEntries);
		sqrError += std::pow(bitmap->getData()[i]/(255.0f)-bitmap2->getData()[i]/(255.0f), 2)/nEntries;
	}

	cout << "Mean: " << mean << ", error:" << std::sqrt(sqrError)/mean<< endl;

	ref<FileStream> stream2 = new FileStream("cat2.png", FileStream::ETruncReadWrite);
	bitmap2->save(Bitmap::EPNG, stream2);
}

void testLineIntegral() {
	ref<FileStream> stream = new FileStream("cat.png", FileStream::EReadOnly);
	ref<Bitmap> bitmap = new Bitmap(Bitmap::EPNG, stream);

	Point2 start(0.121, 0.566);
	Point2 end(0.815, 0.318);

	Vector2 dir = Vector2(end-start);

	for (int i=3; i<18; ++i) {
		size_t nsteps = (int) std::ldexp((Float) 1, i);
		uint8_t *data = bitmap->getData();
		double accum = 0;
		for (size_t i=0; i<nsteps; ++i) {
			Vector2 p = (start + dir * (i / (Float) (nsteps-1)));
			p *= 512;
			int x = (int) p.x, y = (int) p.y;

			accum += data[x+512*y] / 255.0f;
		}
		accum /= nsteps;
		accum *= dir.length();

		printf("%i steps=%f\n", (int) nsteps, accum);
	}
}

void testLineIntegralWavelet() {
	ref<FileStream> stream = new FileStream("cat.png", FileStream::EReadOnly);
	ref<Bitmap> bitmap = new Bitmap(Bitmap::EPNG, stream);
	ref<Wavelet2D> wavelet = new Wavelet2D(bitmap);
	ref<SparseWavelet2D> sparse = wavelet->toSparseWavelet();

	Point2 start(0.121, 0.566);
	Point2 end(0.815, 0.318);
	start*=512; end *= 512;

	cout << sparse->lineIntegral(start, end) << endl;
}

void testSHRotation() {
	ref<Random> random = new Random();
	int bands = 8;

	SHVector vec1(bands);
	for (int l=0; l<bands; ++l)
		for (int m=-l; m<=l; ++m)
			vec1(l, m) = random->nextFloat();

	Vector axis(squareToSphere(Point2(random->nextFloat(), random->nextFloat())));
	Transform trafo = Transform::rotate(axis, random->nextFloat()*360);
	Transform inv = trafo.inverse();
	SHRotation rot(vec1.getBands());

	SHVector::rotation(trafo, rot);
	SHVector vec2(bands);

	rot(vec1, vec2);

	for (int i=0; i<100; ++i) {
		Vector dir1(squareToSphere(Point2(random->nextFloat(), random->nextFloat()))), dir2;
		trafo(dir1, dir2);

		Float value1 = vec1.eval(dir2);
		Float value2 = vec2.eval(dir1);
		SAssert(std::abs(value1-value2) < Epsilon);
	}
	cout << "Passed." << endl;
}

struct ClampedCos {
	Vector axis;
	ClampedCos(Vector axis) : axis(axis) { }
	Float operator()(const Vector &w) const { return std::max((Float) 0, dot(w, axis)); }
};

void testSHSampler() {
	int bands = 25, numSamples = 100, depth = 12;

	Vector v = normalize(Vector(1, 2, 3));
	ref<Random> random = new Random();
	SHVector clampedCos = SHVector(bands);
	clampedCos.project(ClampedCos(v), numSamples);
	Float clampedCosError = clampedCos.l2Error(ClampedCos(v), numSamples);
	clampedCos.normalize();

	cout << "Projection error = " << clampedCosError << endl;
	cout << "Precomputing mip-maps" << endl;
	ref<SHSampler> sampler = new SHSampler(bands, depth);
	cout << "Done: "<< sampler->toString() << endl;
	Float accum = 0;
	int nsamples = 100;
	for (int i=0; i<=nsamples; ++i) {
		Point2 sample(random->nextFloat(), random->nextFloat());
		Float pdf1 = sampler->warp(clampedCos, sample);
		Float pdf2 = dot(v, sphericalDirection(sample.x, sample.y))/M_PI;
		Float relerr = std::abs(pdf1-pdf2)/pdf2;
		accum += relerr;
		SAssert(relerr < 0.04);
	}
	SAssert(accum / nsamples < 1);
}

struct GridFunctor {
	Float operator()(Float value, Float length) const {
		cout << "functor(value=" << value << ", length=" << length << ")" << endl;
		return 0;
	}
};

void testGrid() {
	Grid<Float> grid(Vector3i(3, 3, 1), AABB(Point(10, 0, 0), Point(11,1,1)));
	Ray ray(Point(10 + 1.0f/6.0f, .5, .5), Vector(1, 0, 0), 0, 2.0f/3.0f);

	grid(0, 1, 0) = 1;
	grid(1, 1, 0) = 2;
	grid(2, 1, 0) = 3;
	
	SAssert(std::abs(grid.lookup(Point(10 + 1.0/3.0f, .5, .5))-1.5) < 1e-6);

	GridFunctor functor;
	grid.rasterize(ray, functor);
	
	cout << "Expected: functor(value=1, length=" << 1.0f/6.0f << ")" << endl;
	cout << "Expected: functor(value=2, length=" << 1.0f/3.0f << ")" << endl;
	cout << "Expected: functor(value=3, length=" << 1.0f/6.0f << ")" << endl;

	SAssert(grid(0,1,0) == 0);
	SAssert(grid(1,1,0) == 0);
	SAssert(grid(2,1,0) == 0);
}

int main(int argc, char **argv) {
	Class::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	SHVector::staticInitialization();
	try {
		/*
		testHalton();
		testHammersley();
		testRadicalInverseIncr();
		testSpotLuminaire();
		testHG1();
		testHG2();
		testVariance();
		testPhotonMap();
		testMaxExp();
		testHeterogeneous();
		testLineIntegral();
		testLineIntegralWavelet();
		testWaveletBasic();
		testWavelet();
		testWavelet3D();
		testGrid();
		testSHRotation();
		testSHSampler();
		*/
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
		exit(-1);
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type! " << std::endl;
		exit(-1);
	}
	SHVector::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Class::staticShutdown();
	return 0;
}
