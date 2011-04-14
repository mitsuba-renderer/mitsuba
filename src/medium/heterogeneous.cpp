/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief When the following line is uncommented, the medium implementation 
 * stops integrating densities when it is determined that the segment has a
 * throughput of less than 'Epsilon' (see \c mitsuba/core/constants.h)
 */
#define HETVOL_EARLY_EXIT 1

/// Generate a few statistics related to the implementation?
//#define HETVOL_STATISTICS 1

#if defined(HETVOL_STATISTICS)
static StatsCounter avgNewtonIterations("Heterogeneous volume", 
		"Avg. # of Newton-Bisection iterations", EAverage);
static StatsCounter avgRayMarchingStepsTransmission("Heterogeneous volume", 
		"Avg. # of ray marching steps (transmission)", EAverage);
static StatsCounter avgRayMarchingStepsSampling("Heterogeneous volume", 
		"Avg. # of ray marching steps (sampling)", EAverage);
static StatsCounter earlyExits("Heterogeneous volume", 
		"Number of early exits", EPercentage);
#endif

/**
 * Flexible heterogeneous medium implementation, which acquires its data from 
 * nested <tt>Volume</tt> instances. These can be constant, use a procedural 
 * function, or fetch data from disk, e.g. using a memory-mapped density grid.
 *
 * Instead of allowing separate volumes to be provided for the scattering
 * parameters sigma_s and sigma_t, this class instead takes the approach of 
 * enforcing a spectrally uniform sigma_t, which must be provided using a 
 * nested scalar-valued volume named 'density'.
 *
 * A nested spectrum-valued 'albedo' volume must also be provided, which is 
 * used to compute the parameter sigma_s using the expression
 * "sigma_s = density * albedo" (i.e. 'albedo' contains the single-scattering
 * albedo of the medium).
 *
 * Optionally, one can also provide an vector-valued 'orientation' volume,
 * which contains local particle orientations that will be passed to
 * scattering models such as Kajiya-Kay phase function.
 */
class HeterogeneousMedium : public Medium {
public:
	HeterogeneousMedium(const Properties &props) 
		: Medium(props) {
		m_stepSize = props.getFloat("stepSize", 0);
		if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
			Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
				"homogeneous media. Please use nested volume instances to supply "
				"these parameters");
	}

	/* Unserialize from a binary data stream */
	HeterogeneousMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_densities = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientations = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		configure();
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		manager->serialize(stream, m_densities.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientations.get());
		stream->writeFloat(m_stepSize);
	}

	void configure() {
		Medium::configure();
		if (m_densities.get() == NULL)
			Log(EError, "No densities specified!");
		if (m_albedo.get() == NULL)
			Log(EError, "No albedo specified!");
		m_aabb = m_densities->getAABB();

		if (m_stepSize == 0) {
			m_stepSize = std::min(
				m_densities->getStepSize(), m_albedo->getStepSize());
			if (m_orientations != NULL)
				m_stepSize = std::min(m_stepSize,
					m_orientations->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity()) 
				Log(EError, "Unable to infer a suitable step size, please specify one "
						"manually using the 'stepSize' parameter.");
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);

			if (name == "albedo") {
				Assert(volume->supportsSpectrumLookups());
				m_albedo = volume;
			} else if (name == "density") {
				Assert(volume->supportsFloatLookups());
				m_densities = volume;
			} else if (name == "orientations") {
				Assert(volume->supportsVectorLookups());
				m_orientations = volume;
			} else {
				Medium::addChild(name, child);
			}
		} else {
			Medium::addChild(name, child);
		}
	}

	/*
	 * This function uses Simpson quadrature to compute following 
	 * integral:
	 * 
	 *    \int_{ray.mint}^{ray.maxt} density(ray(x)) dx
	 * 
	 * The integration proceeds by splitting the function into
	 * approximately \c (ray.maxt-ray.mint)/m_stepSize segments,
	 * each of which are then approximated by a quadratic polynomial.
	 * The step size must be chosen so that this approximation is 
	 * valid given the behavior of the integrand.
	 *
	 * \param ray
	 *    Ray segment to be used for the integration
	 *
	 * \return
	 *    The integrated density
	 */
	Float integrateDensity(const Ray &ray) const {
		/* Determine the ray segment, along which the
		   density integration should take place */
		Float mint, maxt;
		if (!m_aabb.rayIntersect(ray, mint, maxt))
			return 0.0f;

		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt-mint, maxComp = 0;

		Point p = ray(mint), pLast = ray(maxt);

		for (int i=0; i<3; ++i) {
			maxComp = std::max(maxComp, std::abs(p[i]));
			maxComp = std::max(maxComp, std::abs(pLast[i]));
		}

		/* Ignore degenerate path segments */
		if (length < 1e-6f * maxComp) 
			return 0.0f;

		/* Compute a suitable step size */
		uint32_t nSteps = (uint32_t) std::ceil(length / m_stepSize);
		nSteps += nSteps % 2;
		const Float stepSize = length/nSteps;
		const Vector increment = ray.d * stepSize;

		#if defined(HETVOL_STATISTICS)
			avgRayMarchingStepsTransmission.incrementBase();
			earlyExits.incrementBase();
		#endif

		/* Perform lookups at the first and last node */
		Float integratedDensity = m_densities->lookupFloat(p) 
			+ m_densities->lookupFloat(pLast);

		#if defined(HETVOL_EARLY_EXIT)
			const Float stopAfterDensity = -std::log(Epsilon);
			const Float stopValue = stopAfterDensity*3.0f/(stepSize
					* m_densityMultiplier);
		#endif

		p += increment;

		Float m = 4;
		for (uint32_t i=1; i<nSteps; ++i) {
			integratedDensity += m * m_densities->lookupFloat(p);
			m = 6 - m;

			#if defined(HETVOL_STATISTICS)
				++avgRayMarchingStepsTransmission;
			#endif
			
			#if defined(HETVOL_EARLY_EXIT)
				if (integratedDensity > stopValue) {
					// Reached the threshold -- stop early
					#if defined(HETVOL_STATISTICS)
						++earlyExits;
					#endif
					return std::numeric_limits<Float>::infinity();
				}
			#endif

			Point next = p + increment;
			if (p == next) {
				Log(EWarn, "integrateDensity(): unable to make forward progress -- "
						"round-off error issues? The step size was %e, mint=%f, "
						"maxt=%f, nSteps=%i, ray=%s", stepSize, mint, maxt, nSteps, 
						ray.toString().c_str());
				break;
			}
			p = next;
		}

		return integratedDensity * m_densityMultiplier
			* stepSize * (1.0f / 3.0f);
	}

	/**
	 * This function uses composite Simpson quadrature to solve the 
	 * following integral equation for \a t:
	 * 
	 *    \int_{ray.mint}^t density(ray(x)) dx == desiredDensity
	 * 
	 * The integration proceeds by splitting the function into
	 * approximately \c (ray.maxt-ray.mint)/m_stepSize segments,
	 * each of which are then approximated by a quadratic polynomial.
	 * The step size must be chosen so that this approximation is 
	 * valid given the behavior of the integrand.
	 * 
	 * \param ray
	 *    Ray segment to be used for the integration
	 *
	 * \param desiredDensity
	 *    Right hand side of the above equation
	 *
	 * \param integratedDensity
	 *    Contains the final integrated density. Upon success, this value
	 *    should closely match \c desiredDensity. When the equation could
	 *    \a not be solved, the parameter contains the integrated density
	 *    from \c ray.mint to \c ray.maxt (which, in this case, must be 
	 *    less than \c desiredDensity).
	 *
	 * \param t
	 *    After calling this function, \c t will store the solution of the above
	 *    equation. When there is no solution, it will be set to zero.
	 *
	 * \param densityAtMinT
	 *    After calling this function, \c densityAtMinT will store the
	 *    underlying density function evaluated at \c ray(ray.mint).
	 *
	 * \param densityAtT
	 *    After calling this function, \c densityAtT will store the
	 *    underlying density function evaluated at \c ray(t). When
	 *    there is no solution, it will be set to zero.
	 *
	 * \return
	 *    When no solution can be found in [ray.mint, ray.maxt] the
	 *    function returns \c false.
	 */
	bool invertDensityIntegral(const Ray &ray, Float desiredDensity,
			Float &integratedDensity, Float &t, Float &densityAtMinT,
			Float &densityAtT) const {
		integratedDensity = densityAtMinT = densityAtT = 0.0f;

		/* Determine the ray segment, along which the
		   density integration should take place */
		Float mint, maxt;
		if (!m_aabb.rayIntersect(ray, mint, maxt))
			return false;
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt - mint, maxComp = 0;
		Point p = ray(mint), pLast = ray(maxt);

		for (int i=0; i<3; ++i) {
			maxComp = std::max(maxComp, std::abs(p[i]));
			maxComp = std::max(maxComp, std::abs(pLast[i]));
		}

		/* Ignore degenerate path segments */
		if (length < 1e-6f * maxComp) 
			return 0.0f;

		/* Compute a suitable step size */
		uint32_t nSteps = (uint32_t) std::ceil(length / m_stepSize);
		Float stepSize = length / nSteps,
			  multiplier = (1.0f / 6.0f) * stepSize
				  * m_densityMultiplier;
		Vector fullStep = ray.d * stepSize,
			   halfStep = fullStep * .5f;

		Float node1 = m_densities->lookupFloat(p);

		if (ray.mint == mint)
			densityAtMinT = node1 * m_densityMultiplier;
		else
			densityAtMinT = 0.0f;

		#if defined(HETVOL_STATISTICS)
			avgRayMarchingStepsSampling.incrementBase();
		#endif

		for (uint32_t i=0; i<nSteps; ++i) {
			Float node2 = m_densities->lookupFloat(p + halfStep),
				  node3 = m_densities->lookupFloat(p + fullStep),
				  newDensity = integratedDensity + multiplier * 
						(node1+node2*4+node3);
			#if defined(HETVOL_STATISTICS)
				++avgRayMarchingStepsSampling;
			#endif
			if (newDensity >= desiredDensity) {
				/* The integrated density of the last segment exceeds the desired
				   amount -- now use the Simpson quadrature expression and 
				   Newton-Bisection to find the precise location of the scattering
				   event. Note that no further density queries are performed after
				   this point; instead, the densities are modeled based on a 
				   quadratic polynomial that is fit to the last three lookups */

				Float a = 0, b = stepSize, x = a,
					  fx = integratedDensity - desiredDensity,
					  stepSizeSqr = stepSize * stepSize,
					  temp = m_densityMultiplier / stepSizeSqr;
				int it = 1;

				#if defined(HETVOL_STATISTICS)
					avgNewtonIterations.incrementBase();
				#endif
				while (true) {
					#if defined(HETVOL_STATISTICS)
						++avgNewtonIterations;
					#endif
					/* Lagrange polynomial from the Simpson quadrature */
					Float dfx = temp * (node1 * stepSizeSqr
						- (3*node1 - 4*node2 + node3)*stepSize*x
						+ 2*(node1 - 2*node2 + node3)*x*x);
					#if 0
						cout << "Iteration " << it << ":  a=" << a << ", b=" << b 
							<< ", x=" << x << ", fx=" << fx << ", dfx=" << dfx << endl;
					#endif

					x -= fx/dfx;

					if (EXPECT_NOT_TAKEN(x <= a || x >= b || dfx == 0)) 
						x = 0.5f * (b + a);

					/* Integrated version of the above Lagrange polynomial */
					Float intval = integratedDensity + temp * (1.0f / 6.0f) * (x *
						(6*node1*stepSizeSqr - 3*(3*node1 - 4*node2 + node3)*stepSize*x
						+ 4*(node1 - 2*node2 + node3)*x*x));
					fx = intval-desiredDensity;

					if (std::abs(fx) < 1e-6f) {
						t = mint + stepSize * i + x;
						integratedDensity = intval;
						densityAtT = temp * (node1 * stepSizeSqr
							- (3*node1 - 4*node2 + node3)*stepSize*x
							+ 2*(node1 - 2*node2 + node3)*x*x);
						return true;
					} else if (++it > 30) {
						Log(EWarn, "invertDensityIntegral(): stuck in Newton-Bisection -- "
							"round-off error issues? The step size was %e, fx=%f, dfx=%f, "
							"a=%f, b=%f", stepSize, fx, dfx, a, b);
						return false;
					}

					if (fx > 0)
						b = x;
					else
						a = x;
				}
			}

			Point next = p + fullStep;
			if (p == next) {
				Log(EWarn, "invertDensityIntegral(): unable to make forward progress -- "
						"round-off error issues? The step size was %e", stepSize);
				break;
			}
			integratedDensity = newDensity;
			node1 = node3;
			p = next;
		}

		return false;
	}

	Spectrum getTransmittance(const Ray &ray) const {
		return Spectrum(std::exp(-integrateDensity(ray)));
	}

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		Float desiredDensity = -std::log(1-sampler->next1D());
		Float integratedDensity, densityAtMinT, densityAtT;
		bool success = false;

		if (invertDensityIntegral(ray, desiredDensity, integratedDensity, 
				mRec.t, densityAtMinT, densityAtT)) {
			mRec.p = ray(mRec.t);
			success = true;
			Spectrum albedo = m_albedo->lookupSpectrum(mRec.p);
			mRec.sigmaS = albedo * densityAtT;
			mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
			mRec.albedo = albedo.max();
			mRec.orientation = m_orientations != NULL 
				? m_orientations->lookupVector(mRec.p) : Vector(0.0f);
		}

		Float expVal = std::exp(-integratedDensity);
		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * densityAtT;
		mRec.pdfSuccessRev = expVal * densityAtMinT;
		mRec.transmittance = Spectrum(expVal);

		return success;
	}

	void pdfDistance(const Ray &ray, MediumSamplingRecord &mRec) const {
		Float expVal = std::exp(-integrateDensity(ray));

		mRec.transmittance = Spectrum(expVal);
		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * 
			m_densities->lookupFloat(ray(ray.maxt)) * m_densityMultiplier;
		mRec.pdfSuccessRev = expVal * 
			m_densities->lookupFloat(ray(ray.mint)) * m_densityMultiplier;
	}

	bool isHomogeneous() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousMedium[" << endl
			<< "  albedo = " << indent(m_albedo.toString()) << "," << endl
			<< "  orientations = " << indent(m_orientations.toString()) << "," << endl
			<< "  densities = " << indent(m_densities.toString()) << "," << endl
			<< "  stepSize = " << m_stepSize << "," << endl
			<< "  densityMultiplier = " << m_densityMultiplier << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<VolumeDataSource> m_densities;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientations;
	Float m_stepSize;
	AABB m_aabb;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMedium, "Heterogeneous medium");
MTS_NAMESPACE_END
