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

#include <mitsuba/render/scene.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

/*
 * When inverting an integral of the form f(t):=\int_0^t [...] dt
 * using composite Simpson's rule quadrature, the following constant
 * specified how many times the step size will be reduced when we
 * have almost reached the desired value.
 */
#define NUM_REFINES 10

/**
 * Heterogeneous medium class using trilinear interpolation, 
 * Simpson quadrature and one of several possible sampling
 * strategies. Data files have to be provided in an ASCII
 * format as follows:
 *
 * - The first three numbers determine the X,Y and Z resolution,
 *   each of which has to be larger than 2.
 * - The next six numbers determine the minimum X, Y and Z
 *   as well as the maximum X, Y and Z values of the enclosing
 *   axis-aligned bounding box.
 * - Afterwards, (xres*yres*zres) density samples follow in 
 *   XYZ order, (e.g. the second entry has coordinate x=2).
 */
class HeterogeneousMedium : public Medium {
public:
	/*
	 * We want to generate realizations of a random variable, which has a 
	 * cumulative distribution function given by
	 *
	 * F(t) := 1 - \int_0^t \sigma_t(s) \exp(\int_0^s \sigma_t(s') ds') ds
	 */
	enum ESamplingStrategy {
		/* Standard approach: generate a 'desired' accumulated density by 
		   sampling an exponentially distributed random variable. Afterwards,
		   try to find the distance, at which this much density has been 
		   accumulated. The composite Simpson's rule is used to integrate 
		   densities along the ray */
		EStandard,

		/* Naive variant for verification purposes: uniformly sample a distance
		   along the ray segment intersecting the volume */
		EUniform,

		/* Double integral approach - stupid and slow, but it also works.. */
		EDoubleIntegral,

		/* Rejection sampling approach by [Coleman et al., 1967]. Only for 
		   media with a wavelength-independent extinction coefficient 
		   (see below) */
		EColeman
	};

	HeterogeneousMedium(const Properties &props) 
		: Medium(props), m_data(NULL), m_empty(NULL) {
		std::string strategy = props.getString("strategy", "standard");

		if (strategy == "standard") {
			m_strategy = EStandard;
		} else if (strategy == "uniform") {
			m_strategy = EUniform;
		} else if (strategy == "double") {
			m_strategy = EDoubleIntegral;
		} else if (strategy == "coleman") {
			m_strategy = EColeman;
			for (size_t i=1; i<SPECTRUM_SAMPLES; ++i)
				if (m_sigmaT[i] != m_sigmaT[0])
					Log(EError, "The [Coleman et al.] sampling strategy is not supported "
						"for colored media - computing the attenuation "
						"of individual channels requires integrating volume densities, "
						"which this strategy was trying to avoid..");
		} else {
			Log(EError, "Specified an unknown sampling strategy");
		}

	
		std::string volData = props.getString("filename");
		volData = FileResolver::getInstance()->resolve(volData);

		/* Medium to world transformation - can't have nonuniform scales. Also note 
		   that a uniform scale factor of 100 will not reduce densities by that 
		   amount */
		m_mediumToWorld = props.getTransform("toWorld", Transform());

		Log(EInfo, "Loading volume data from \"%s\" ..", volData.c_str());
		std::ifstream is(volData.c_str());
		if (is.bad() || is.fail())
			Log(EError, "Invalid medium data file specified");

		/* While integrating density along a ray, approximately one sample
		   per voxel is taken - that number can be changed here */
		m_stepSizeFactor = props.getFloat("stepSizeFactor", 1.0f);

		Float proposedSigma = std::numeric_limits<Float>::infinity();
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			proposedSigma = std::min(proposedSigma, m_sigmaT[i]);
		
		/* Can be used to override the extinction coefficient used to sample distances 
			in the in-scatter line integral. By default, the smallest spectral sample of 
			<tt>sigmaA+sigmaT</tt> is used. */
		m_sigma = props.getFloat("sigma", proposedSigma);

		/* [Fedkiw et al.] volume density format */
		is >> m_res.x; 
		is >> m_res.y; 
		is >> m_res.z;
		is >> m_dataAABB.min.x; is >> m_dataAABB.min.y; is >> m_dataAABB.min.z;
		is >> m_dataAABB.max.x; is >> m_dataAABB.max.y; is >> m_dataAABB.max.z;

		m_aabb.reset();
		for (int i=0; i<8; ++i)
			m_aabb.expandBy(m_mediumToWorld(m_dataAABB.getCorner(i)));

		m_data = new Float[m_res.x * m_res.y * m_res.z];

		for (int i=0; i<m_res.x*m_res.y*m_res.z; ++i)
			is >> m_data[i];
	}

	/* Unserialize from a binary data stream */
	HeterogeneousMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_res = Vector3i(stream);
		m_mediumToWorld = Transform(stream);
		m_dataAABB = AABB(stream);
		m_strategy = (ESamplingStrategy) stream->readInt();
		m_sigma = stream->readFloat();
		m_stepSizeFactor = stream->readFloat();
		m_data = new Float[m_res.x * m_res.y * m_res.z];
		stream->readFloatArray(m_data, m_res.x*m_res.y*m_res.z);
		configure();
	}

	void configure() {
		Medium::configure();

		m_voxels.x = m_res.x-1; m_voxels.y = m_res.y-1; m_voxels.z = m_res.z-1;
		m_empty = new bool[m_voxels.x*m_voxels.y*m_voxels.z];
		m_extents = m_dataAABB.getExtents();

		for (int i=0; i<3; ++i)
			m_cellWidth[i] = m_extents[i] / (Float) m_voxels[i];

		size_t nPoints = m_res.x*m_res.y*m_res.z;
		m_maxDensity = 0;
		for (size_t i=0; i<nPoints; ++i) 
			m_maxDensity = std::max(m_maxDensity, m_data[i]);

		/* Take note of completely empty areas */
		int numEmpty = 0;
		for (int x1=0; x1<m_voxels.x; ++x1) {
			for (int y1=0; y1<m_voxels.y; ++y1) {
				for (int z1=0; z1<m_voxels.z; ++z1) {
					const int x2 = x1+1, y2 = y1+1, z2 = z1+1;
					const Float
						d000 = m_data[(z1*m_res.y + y1)*m_res.x + x1],
						d001 = m_data[(z1*m_res.y + y1)*m_res.x + x2],
						d010 = m_data[(z1*m_res.y + y2)*m_res.x + x1],
						d011 = m_data[(z1*m_res.y + y2)*m_res.x + x2],
						d100 = m_data[(z2*m_res.y + y1)*m_res.x + x1],
						d101 = m_data[(z2*m_res.y + y1)*m_res.x + x2],
						d110 = m_data[(z2*m_res.y + y2)*m_res.x + x1],
						d111 = m_data[(z2*m_res.y + y2)*m_res.x + x2];

					int pos = x1+m_voxels.x*(y1 + m_voxels.y * z1);
					if (d000 == 0 && d001 == 0 && d010 == 0 && d011 == 0 &&
						d100 == 0 && d101 == 0 && d110 == 0 && d111 == 0) {
						numEmpty++;
						m_empty[pos] = true;
					} else {
						m_empty[pos] = false;
					}
				}
			}
		}

		m_worldToMedium = m_mediumToWorld.inverse();
		Float scaleFactor = m_mediumToWorld(Vector(0,0,1)).length();
		m_stepSize = std::numeric_limits<Float>::infinity();
		for (int i=0; i<3; ++i)
			m_stepSize = std::min(m_stepSize, m_extents[i] / (Float) (m_res[i]-1));
		m_stepSize *= scaleFactor * m_stepSizeFactor;
		int voxelCount = m_voxels.x*m_voxels.y*m_voxels.z;
		Log(EInfo, "Resolution: %ix%ix%i, %s, %i voxels, %.1f%% empty", m_res.x, m_res.y, 
			m_res.z, m_dataAABB.toString().c_str(), voxelCount, 100*numEmpty/(Float) voxelCount);
	}

	virtual ~HeterogeneousMedium() {
		if (m_data)
			delete[] m_data;
		if (m_empty)	
			delete[] m_empty;
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);

		m_res.serialize(stream);
		m_mediumToWorld.serialize(stream);
		m_dataAABB.serialize(stream);
		stream->writeInt(m_strategy);
		stream->writeFloat(m_sigma);
		stream->writeFloat(m_stepSizeFactor);
		stream->writeFloatArray(m_data, m_res.x*m_res.y*m_res.z);
	}

	/* Evaluate the density using trilinear interpolation */
	inline Float eval(const Point &p) const {
		if (p.x < 0 || p.y < 0 || p.z < 0)
			return 0.0f;
		const int x1 = (int) p.x, y1 = (int) p.y, z1 = (int) p.z,
				  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x || y2 >= m_res.y || z2 >= m_res.z) {
			/* Do an integer bounds test (may seem redundant - this is
			   to avoid a segfault, should a NaN/Inf ever find its way here..) */
			return 0;
		}

		if (m_empty[x1 + m_voxels.x * (y1 + m_voxels.y * z1)]) 
			return 0.0f;

		const Float fx = p.x-x1, fy = p.y-y1, fz = p.z-z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f-fz;

		const Float
			d000 = m_data[(z1*m_res.y + y1)*m_res.x + x1],
			d001 = m_data[(z1*m_res.y + y1)*m_res.x + x2],
			d010 = m_data[(z1*m_res.y + y2)*m_res.x + x1],
			d011 = m_data[(z1*m_res.y + y2)*m_res.x + x2],
			d100 = m_data[(z2*m_res.y + y1)*m_res.x + x1],
			d101 = m_data[(z2*m_res.y + y1)*m_res.x + x2],
			d110 = m_data[(z2*m_res.y + y2)*m_res.x + x1],
			d111 = m_data[(z2*m_res.y + y2)*m_res.x + x2];

		return ((d000*_fx + d001*fx)*_fy +
				(d010*_fx + d011*fx)*fy)*_fz +
			   ((d100*_fx + d101*fx)*_fy +
				(d110*_fx + d111*fx)*fy)*fz;
	}

	Spectrum tau(const Ray &r) const {
		/* Ensure we're dealing with a normalized ray */
		Float rayLength = r.d.length();
		Ray ray;
		m_worldToMedium(r, ray);
		ray.setDirection(ray.d/rayLength);
		ray.maxt *= rayLength;
		ray.mint *= rayLength;

		/* Intersect against the AABB containing the medium.
	       Slightly offset to avoid zero samples at the endpoints */
		Float mint, maxt;
		if (!m_dataAABB.rayIntersect(ray, mint, maxt))
			return Spectrum(0.0f);

		mint = std::max(ray.mint, mint+Epsilon);
		maxt = std::min(ray.maxt, maxt-Epsilon);

		if (mint >= maxt) 
			return Spectrum(0.0f);

		/* Traverse directly in grid coordinates */
		ray.o = Point(ray(mint)-m_dataAABB.min);
		for (int i=0; i<3; ++i) {
			ray.o[i] /= m_cellWidth[i];
			ray.d[i] /= m_cellWidth[i];
		}
		maxt -= mint;

		/* Composite Simpson's rule */
		int nParts = (int) std::ceil(maxt/m_stepSize);
		nParts += nParts % 2;
		const Float stepSize = maxt/nParts;
		const Vector increment = ray.d * stepSize;

		Float m=4, accumulatedDensity = eval(ray.o) + eval(ray(maxt));
		ray.o += increment;
		for (int i=1; i<nParts; ++i) {
			accumulatedDensity += m*eval(ray.o);
			ray.o += increment;
			m = 6-m;
		}
		accumulatedDensity *= stepSize/3;

		return m_sigmaT * accumulatedDensity;
	}

	bool sampleDistance(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		if (m_strategy == EStandard)
			return importanceSampleStandard(r, maxDist, mRec, sampler);
		else if (m_strategy == EUniform)
			return uniformlySampleDistance(r, maxDist, mRec, sampler);
		else if (m_strategy == EDoubleIntegral)
			return importanceSampleDoubleIntegral(r, maxDist, mRec, sampler);
		else
			return importanceSampleColeman(r, maxDist, mRec, sampler);
	}

	bool importanceSampleStandard(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		/* Ensure we're dealing with a normalized ray */
		Float rayLength = r.d.length();
		Ray ray;
		m_worldToMedium(r, ray);
		ray.setDirection(ray.d/rayLength);
		maxDist *= rayLength;
		ray.mint *= rayLength;

		/* Intersect against the AABB containing the medium.
	       Slightly offset to avoid zero samples at the endpoints */
		Float mint, maxt;
		mRec.pdf = 1; mRec.miWeight = 1; mRec.attenuation = Spectrum(1.0f);
		if (!m_dataAABB.rayIntersect(ray, mint, maxt))
			return false;

		mint = std::max(ray.mint, mint+Epsilon);
		maxt = std::min(maxDist, maxt-Epsilon);

		if (mint >= maxt) 
			return false;

		/* Traverse directly in grid coordinates */
		ray.o = Point(ray(mint)-m_dataAABB.min);
		for (int i=0; i<3; ++i) {
			ray.o[i] /= m_cellWidth[i];
			ray.d[i] /= m_cellWidth[i];
		}
		maxt -= mint;
		Float sigmaT = m_sigma;

		int nParts = (int) std::ceil(maxt/m_stepSize);
		Float stepSize = maxt/nParts;
		Vector fullIncrement = ray.d * stepSize,
			   halfIncrement = fullIncrement * .5f;

		Float t=0, desiredDensity = -std::log(1-sampler->next1D())/sigmaT,
			  lastNode = eval(ray.o), accumulatedDensity = 0;

		bool success = false;
		int numRefines = 0;

		while (t <= maxt) {
			const Float node1 = lastNode,
						node2 = eval(ray.o + halfIncrement),
						node3 = eval(ray.o + fullIncrement);
			const Float newDensity = accumulatedDensity + stepSize / 6 * (node1+4*node2+node3);
			
			if (newDensity > desiredDensity) {
				if (++numRefines > NUM_REFINES) {
					lastNode = node3;
					success = true;
					break;
				}

				stepSize *= .5f;
				fullIncrement *= .5f;
				halfIncrement *= .5f;
				continue;
			}
		
			accumulatedDensity = newDensity;
			lastNode = node3;
			t += stepSize;
			ray.o += fullIncrement;
		}

		if (!success) {
			/* Could not achieve the desired density - hit the next surface instead */
			mRec.pdf = std::exp(-accumulatedDensity * sigmaT);
			mRec.attenuation = (m_sigmaT * (-accumulatedDensity)).exp();
			return false;
		}

		mRec.t = t+mint;
		mRec.p = r(mRec.t);
		mRec.pdf = std::exp(-accumulatedDensity * sigmaT) * lastNode * sigmaT;
		mRec.attenuation = (m_sigmaT * (-accumulatedDensity)).exp();
		mRec.sigmaA = m_sigmaA * lastNode;
		mRec.sigmaS = m_sigmaS * lastNode;
		mRec.albedo = m_albedo;
		mRec.medium = this;

		return true;
	}

	bool importanceSampleDoubleIntegral(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		/* Ensure we're dealing with a normalized ray */
		Float rayLength = r.d.length();
		Ray ray;
		m_worldToMedium(r, ray);
		ray.setDirection(ray.d/rayLength);
		maxDist *= rayLength;
		ray.mint *= rayLength;

		/* Intersect against the AABB containing the medium.
	       Slightly offset to avoid zero samples at the endpoints */
		Float mint, maxt;
		mRec.pdf = 1; mRec.miWeight = 1; mRec.attenuation = Spectrum(1.0f);
		if (!m_dataAABB.rayIntersect(ray, mint, maxt))
			return false;

		mint = std::max(ray.mint, mint+Epsilon);
		maxt = std::min(maxDist, maxt-Epsilon);

		if (mint >= maxt) 
			return false;

		/* Traverse directly in grid coordinates */
		ray.o = Point(ray(mint)-m_dataAABB.min);
		for (int i=0; i<3; ++i) {
			ray.o[i] /= m_cellWidth[i];
			ray.d[i] /= m_cellWidth[i];
		}
		maxt -= mint;
		Float sigmaT = m_sigma;
		Float innerIntegral = 0, outerIntegral = 0;

		int nParts = (int) std::ceil(maxt/m_stepSize);
		Float stepSize = maxt/nParts, halfStepSize = stepSize *.5f;
		Vector fullIncrement = ray.d * stepSize,
			   halfIncrement = fullIncrement * .5f,
			   quarterIncrement = halfIncrement * .5f;
		Float lastNode = eval(ray.o);
		Ray rayBackup(ray);

		/* First pass: normalization */
		for (int i=0; i<nParts; ++i) {
			/* Inner integral: integrate density \rho(t) */
			const Float node1 = lastNode,
						node2 = eval(ray.o + quarterIncrement),
						node3 = eval(ray.o + halfIncrement),
						node4 = eval(ray.o + halfIncrement + quarterIncrement),
						node5 = eval(ray.o + fullIncrement);

			/* Simpson's rule with twice the resolution */
			const Float innerIntegralMiddle = innerIntegral 
							+ halfStepSize / 6 * (node1 + 4*node2 + node3),
						innerIntegralRight = innerIntegralMiddle
							+ halfStepSize / 6 * (node3 + 4*node4 + node5);

			/* Outer integral: integrate \rho(s) exp(-\int_0^t \sigma_t \rho(s) ds) */
			const Float contributionLeft   = std::exp(-innerIntegral*sigmaT)       * node1,
						contributionMiddle = std::exp(-innerIntegralMiddle*sigmaT) * node3,
						contributionRight  = std::exp(-innerIntegralRight*sigmaT)  * node5;

			const Float outerIntegralRight = outerIntegral + 
				stepSize / 6 * (contributionLeft + 4*contributionMiddle + contributionRight);

			innerIntegral = innerIntegralRight;
			outerIntegral = outerIntegralRight;
			lastNode = node5;
			ray.o += fullIncrement;
		}

		Float t=0, xi = sampler->next1D();
		const Float totalAttenuation = std::exp(-innerIntegral * sigmaT);

		if (xi < totalAttenuation) {
			/* Sample the surface instead */
			mRec.pdf = totalAttenuation;
			mRec.attenuation = (m_sigmaT * (-innerIntegral)).exp();
			return false;
		} else {
			mRec.pdf = 1-totalAttenuation;
			xi = (xi-totalAttenuation)/mRec.pdf;
		}
		const Float normFactor = 1.0f / outerIntegral;

		/* Second pass: invert the outer integral */
		innerIntegral = 0; outerIntegral = 0;
		ray = rayBackup; lastNode = eval(ray.o);

		bool success = false;
		int numRefines = 0;

		while (t <= maxt) {
			/* Inner integral: integrate density \rho(t) */
			const Float node1 = lastNode,
						node2 = eval(ray.o + quarterIncrement),
						node3 = eval(ray.o + halfIncrement),
						node4 = eval(ray.o + halfIncrement + quarterIncrement),
						node5 = eval(ray.o + fullIncrement);
			
			/* Simpson's rule with twice the resolution */
			const Float innerIntegralMiddle = innerIntegral 
							+ halfStepSize / 6 * (node1+4*node2+node3),
						innerIntegralRight = innerIntegralMiddle
							+ halfStepSize / 6 * (node3+4*node4+node5);

			/* Outer integral: integrate \rho(s) exp(-\int_0^t \sigma_t \rho(s) ds) */
			const Float contributionLeft   = std::exp(-innerIntegral*sigmaT)       * node1 * normFactor,
						contributionMiddle = std::exp(-innerIntegralMiddle*sigmaT) * node3 * normFactor,
						contributionRight  = std::exp(-innerIntegralRight*sigmaT)  * node5 * normFactor;

			const Float outerIntegralRight = outerIntegral + 
				stepSize / 6 * (contributionLeft + 4*contributionMiddle + contributionRight);

			if (outerIntegralRight > xi) {
				if (++numRefines > NUM_REFINES) {
					outerIntegral = outerIntegralRight;
					mRec.pdf *= contributionRight;
					lastNode = node5;
					success = true;
					break;
				}

				stepSize *= .5f;
				halfStepSize *= .5f;
				fullIncrement *= .5f;
				halfIncrement *= .5f;
				quarterIncrement *= .5f;
				continue;
			}
		
			innerIntegral = innerIntegralRight;
			outerIntegral = outerIntegralRight;
			lastNode = node5;
			t += stepSize;
			ray.o += fullIncrement;
		}

		if (!success) {
			/* Left the integration domain due to round-off errors */
			mRec.pdf = 1;
			lastNode = 0;
		}

		mRec.t = t+mint;
		mRec.p = r(mRec.t);
		mRec.attenuation = (m_sigmaT * (-innerIntegral)).exp();
		mRec.sigmaA = m_sigmaA * lastNode;
		mRec.sigmaS = m_sigmaS * lastNode;
		mRec.albedo = m_albedo;
		mRec.medium = this;

		return true;
	}

	bool importanceSampleColeman(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		/* Ensure we're dealing with a normalized ray */
		Float rayLength = r.d.length();
		Ray ray;
		m_worldToMedium(r, ray);
		ray.setDirection(ray.d/rayLength);
		maxDist *= rayLength;
		ray.mint *= rayLength;

		/* Intersect against the AABB containing the medium.
	       Slightly offset to avoid zero samples at the endpoints */
		Float mint, maxt;
		mRec.pdf = 1; mRec.miWeight = 1; mRec.attenuation = Spectrum(1.0f);
		if (!m_dataAABB.rayIntersect(ray, mint, maxt))
			return false;

		mint = std::max(ray.mint, mint+Epsilon);
		maxt = std::min(maxDist, maxt-Epsilon);

		if (mint >= maxt) 
			return false;

		/* Traverse directly in grid coordinates */
		ray.o = Point(ray(mint)-m_dataAABB.min);
		for (int i=0; i<3; ++i) {
			ray.o[i] /= m_cellWidth[i];
			ray.d[i] /= m_cellWidth[i];
		}
		maxt -= mint;

		Float sigmaT = m_sigma;
		Float maxSigmaT = sigmaT * m_maxDensity;

		Float density, t = 0;
		while (true) {
			Float x = sampler->next1D(), y = sampler->next1D();
			t += -std::log(1-x) / maxSigmaT;

			if (t > maxt) {
				mRec.pdf = 1;
				mRec.attenuation = Spectrum(1.0f);
				return false;
			}
			density = eval(ray(t));

			if (density/m_maxDensity > y)
				break;
		}

		mRec.sigmaA = m_sigmaA * density;
		mRec.sigmaS = m_sigmaS * density;
		mRec.albedo = m_albedo;

		/* Intentionally not computing some of these,
		   since they will cancel out when computing
			   sigmaS * attenuation / pdf
		   (which is hopefully how this method is used :))
		*/
		mRec.attenuation = Spectrum(1);
		mRec.pdf = density * sigmaT;
		mRec.medium = this;
		mRec.t = t+mint;
		mRec.p = r(mRec.t);
		return true;
	}

	bool uniformlySampleDistance(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		mRec.pdf = 1; mRec.miWeight = 1; mRec.attenuation = Spectrum(1.0f);

		const Float probSurface = .5, probMedium = 1-probSurface;
		mRec.pdf = 1;

		Ray ray;
		m_worldToMedium(r, ray);

		/* Intersect against the AABB containing the medium.
		   Slightly offset to avoid zero samples at the endpoints */
		Float mint, maxt;
		if (!m_dataAABB.rayIntersect(ray, mint, maxt)) 
			return false;
		mint = std::max(ray.mint, mint+Epsilon);
		maxt = std::min(maxDist, maxt-Epsilon);
		if (mint >= maxt) 
			return false;

		Float U = sampler->next1D();
		if (U > probSurface) {
			U = (U - probSurface)/probMedium;
		} else {
			mRec.pdf = probSurface;
			mRec.attenuation = (-tau(Ray(r, r.mint, maxDist))).exp();
			return false;
		}

		mRec.t = U*(maxt-mint) + mint;
		mRec.p = r(mRec.t);

		Point gridPos(ray(mRec.t)-m_dataAABB.min);
		for (int i=0; i<3; ++i) 
			gridPos[i] /= m_cellWidth[i];

		Float fval = eval(gridPos);
		mRec.pdf = probMedium/(maxt-mint);
		mRec.sigmaA = m_sigmaA * fval;
		mRec.sigmaS = m_sigmaS * fval;
		mRec.albedo = m_albedo;
		mRec.medium = this;
		mRec.attenuation = (-tau(Ray(r, r.mint, mRec.t))).exp();

		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousMedium[" << endl
			<< "  sigmaA = " << m_sigmaA.toString() << "," << std::endl
			<< "  sigmaS = " << m_sigmaS.toString() << "," << std::endl
			<< "  sigmaT = " << m_sigmaT.toString() << "," << std::endl
			<< "  phase = " << indent(m_phaseFunction->toString()) << "," << std::endl
			<< "  res = " << m_res.x << "x" << m_res.y << "x" << m_res.z << "," << std::endl
			<< "  aabb = " << m_dataAABB.toString() << std::endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	ESamplingStrategy m_strategy;
	Float m_sigma;
	Vector3i m_res, m_voxels;
	Vector m_cellWidth;
	Transform m_worldToMedium;
	Transform m_mediumToWorld;
	Vector m_extents;
	Float *m_data;
	bool *m_empty;
	Float m_stepSize, m_stepSizeFactor, m_maxDensity;
	AABB m_dataAABB;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMedium, "Heterogeneous medium");
MTS_NAMESPACE_END
