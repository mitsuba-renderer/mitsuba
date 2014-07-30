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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/statistics.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief When the following line is uncommented, the medium implementation
 * stops integrating density when it is determined that the segment has a
 * throughput of less than 'Epsilon' (see \c mitsuba/core/constants.h)
 */
#define HETVOL_EARLY_EXIT 1

/// Generate a few statistics related to the implementation?
// #define HETVOL_STATISTICS 1

#if defined(HETVOL_STATISTICS)
static StatsCounter avgNewtonIterations("Heterogeneous volume",
		"Avg. # of Newton-Bisection iterations", EAverage);
static StatsCounter avgRayMarchingStepsTransmittance("Heterogeneous volume",
		"Avg. # of ray marching steps (transmittance)", EAverage);
static StatsCounter avgRayMarchingStepsSampling("Heterogeneous volume",
		"Avg. # of ray marching steps (sampling)", EAverage);
static StatsCounter earlyExits("Heterogeneous volume",
		"Number of early exits", EPercentage);
#endif

/*!\plugin{heterogeneous}{Heterogeneous participating medium}
 * \order{2}
 * \parameters{
 *     \parameter{method}{\String}{
 *         Specifies the sampling method that is used to generate
 *         scattering events within the medium.
 *         \begin{enumerate}[(i)]
 *             \item \code{simpson}: Sampling is done by inverting a
 *             deterministic quadrature rule based on composite
 *             Simpson integration over small ray segments. Benefits
 *             from the use of good sample generators (e.g. \pluginref{ldsampler}).
 *             \item \code{woodcock}: Generate samples using
 *             Woodcock tracking. This is usually faster and
 *             always unbiased, but has the disadvantages of not benefiting
 *             from good sample generators and not providing
 *             information that is required by bidirectional
 *             rendering techniques.
 *         \end{enumerate}
 *         Default: \texttt{woodcock}
 *     }
 *     \parameter{density}{\Volume}{
 *         Volumetric data source that supplies the medium densities
 *         (in inverse scene units)
 *     }
 *     \parameter{albedo}{\Volume}{
 *         Volumetric data source that supplies the
 *         single-scattering albedo
 *     }
 *     \parameter{orientation}{\Volume}{
 *         Optional: volumetric data source that supplies the
 *         local particle orientations throughout the medium
 *     }
 *     \parameter{scale}{\Float}{
 *         Optional scale factor that will be applied to the \code{density} parameter.
 *         Provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}
 *     }
 *     \parameter{\Unnamed}{\Phase}{
 *          A nested phase function that describes the directional
 *          scattering properties of the medium. When none is specified,
 *          the renderer will automatically use an instance of
 *          \pluginref{isotropic}.
 *     }
 * }
 *
 * \renderings{
 *     \medrendering{40}{medium_heterogeneous_density_40}
 *     \medrendering{200}{medium_heterogeneous_density_200}
 *     \medrendering{1000}{medium_heterogeneous_density_1000}
 *     \vspace{-2mm}
 *     \caption{Renderings of an index-matched medium using different scale factors (\lstref{hetvolume})}
 * }
 *
 * This plugin provides a flexible heterogeneous medium implementation, which
 * acquires its data from nested \code{volume} instances. These can be
 * constant, use a procedural function, or fetch data from disk, e.g. using a
 * memory-mapped density grid. See \secref{volumes} for details on volume data
 * sources.
 *
 * Instead of allowing separate volumes to be provided for the scattering
 * and absorption parameters \code{sigmaS} and \code{sigmaA} (as is done in
 * \pluginref{homogeneous}), this class instead takes the approach of
 * enforcing a spectrally uniform value of \code{sigmaT}, which must be
 * provided using a nested scalar-valued volume named \code{density}.
 *
 * Another nested spectrum-valued \code{albedo} volume must also be provided, which is
 * used to compute the scattering coefficient $\sigma_s$ using the expression
 * $\sigma_s = \code{scale} * \code{density} * \code{albedo}$ (i.e. 'albedo' contains the
 * single-scattering albedo of the medium.
 *
 * Optionally, one can also provide an vector-valued \code{orientation} volume,
 * which contains local particle orientation that will be passed to
 * scattering models that support this, such as a the Micro-flake or
 * Kajiya-Kay phase functions.
 *
 * \vspace{4mm}
 *
 * \begin{xml}[label=lst:hetvolume,caption=A simple heterogeneous medium backed by a grid volume]
 * <!-- Declare a heterogeneous participating medium named 'smoke' -->
 * <medium type="heterogeneous" id="smoke">
 *     <string name="method" value="simpson"/>
 *
 *     <!-- Acquire density values from an external data file -->
 *     <volume name="density" type="gridvolume">
 *         <string name="filename" value="frame_0150.vol"/>
 *     </volume>
 *
 *     <!-- The albedo is constant and set to 0.9 -->
 *     <volume name="albedo" type="constvolume">
 *         <spectrum name="value" value="0.9"/>
 *     </volume>
 *
 *     <!-- Use an isotropic phase function -->
 *     <phase type="isotropic"/>
 *
 *     <!-- Scale the density values as desired -->
 *     <float name="scale" value="200"/>
 *  </medium>
 *
 * <!-- Attach the index-matched medium to a shape in the scene -->
 * <shape type="obj">
 *     <!-- Load an OBJ file, which contains a mesh version
 *          of the axis-aligned box of the volume data file -->
 *     <string name="filename" value="bounds.obj"/>
 *
 *     <!-- Reference the medium by ID -->
 *     <ref name="interior" id="smoke"/>
 *
 *     <!-- If desired, this shape could also declare
 *          a BSDF to create an index-mismatched
 *          transition, e.g.
 *
 *     <bsdf type="dielectric"/>
 *     -->
 * </shape>
 * \end{xml}
 */
class HeterogeneousMedium : public Medium {
public:
	/// Possible integration modes
	enum EIntegrationMethod {
		/**
		 * \brief Use deterministic composite Simpson quadrature both
		 * to compute transmittances, and to sample scattering locations
		 */
		ESimpsonQuadrature = 0,

		/**
		 * \brief Use stochastic Woodcock tracking. This is potentially
		 * faster and more robust, but has the disadvantage of being
		 * incompatible with bidirectional rendering methods, which
		 * usually need to know the probability of a sample.
		 */
		EWoodcockTracking
	};

	HeterogeneousMedium(const Properties &props)
		: Medium(props) {
		m_stepSize = props.getFloat("stepSize", 0);
		m_scale = props.getFloat("scale", 1);
		if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
			Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
				"homogeneous media. Please use nested volume instances to supply "
				"these parameters");

		if (props.hasProperty("densityMultiplier"))
			Log(EError, "The 'densityMultiplier' parameter has been deprecated and is now called 'scale'.");

		std::string method = boost::to_lower_copy(props.getString("method", "woodcock"));
		if (method == "woodcock")
			m_method = EWoodcockTracking;
		else if (method == "simpson")
			m_method = ESimpsonQuadrature;
		else
			Log(EError, "Unsupported integration method \"%s\"!", method.c_str());
	}

	/* Unserialize from a binary data stream */
	HeterogeneousMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager) {
		m_method = (EIntegrationMethod) stream->readInt();
		m_scale = stream->readFloat();
		m_density = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientation = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		configure();
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeInt(m_method);
		stream->writeFloat(m_scale);
		manager->serialize(stream, m_density.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientation.get());
		stream->writeFloat(m_stepSize);
	}

	void configure() {
		Medium::configure();
		if (m_density.get() == NULL)
			Log(EError, "No density specified!");
		if (m_albedo.get() == NULL)
			Log(EError, "No albedo specified!");
		m_densityAABB = m_density->getAABB();
		m_anisotropicMedium =
			m_phaseFunction->needsDirectionallyVaryingCoefficients();

		/* Assumes that the density medium does not
		   contain values greater than one! */
		m_maxDensity = m_scale * m_density->getMaximumFloatValue();
		if (m_anisotropicMedium)
			m_maxDensity *= m_phaseFunction->sigmaDirMax();
		m_invMaxDensity = 1.0f/m_maxDensity;

		if (m_stepSize == 0) {
			m_stepSize = std::min(
				m_density->getStepSize(), m_albedo->getStepSize());
			if (m_orientation != NULL)
				m_stepSize = std::min(m_stepSize,
					m_orientation->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity())
				Log(EError, "Unable to infer a suitable step size for deterministic "
						"integration, please specify one manually using the 'stepSize' "
						"parameter.");
		}

		if (m_anisotropicMedium && m_orientation.get() == NULL)
			Log(EError, "Cannot use anisotropic phase function: "
				"did not specify a particle orientation field!");
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);

			if (name == "albedo") {
				Assert(volume->supportsSpectrumLookups());
				m_albedo = volume;
			} else if (name == "density") {
				Assert(volume->supportsFloatLookups());
				m_density = volume;
			} else if (name == "orientation") {
				Assert(volume->supportsVectorLookups());
				m_orientation = volume;
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
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return 0.0f;

		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt-mint, maxComp = 0;

		Point p = ray(mint), pLast = ray(maxt);

		/* Ignore degenerate path segments */
		for (int i=0; i<3; ++i)
			maxComp = std::max(std::max(maxComp,
				std::abs(p[i])), std::abs(pLast[i]));
		if (length < 1e-6f * maxComp)
			return 0.0f;

		/* Compute a suitable step size */
		uint32_t nSteps = (uint32_t) std::ceil(length / m_stepSize);
		nSteps += nSteps % 2;
		const Float stepSize = length/nSteps;
		const Vector increment = ray.d * stepSize;

		#if defined(HETVOL_STATISTICS)
			avgRayMarchingStepsTransmittance.incrementBase();
			earlyExits.incrementBase();
		#endif

		/* Perform lookups at the first and last node */
		Float integratedDensity = lookupDensity(p, ray.d)
			+ lookupDensity(pLast, ray.d);

		#if defined(HETVOL_EARLY_EXIT)
			const Float stopAfterDensity = -math::fastlog(Epsilon);
			const Float stopValue = stopAfterDensity*3.0f/(stepSize
					* m_scale);
		#endif

		p += increment;

		Float m = 4;
		for (uint32_t i=1; i<nSteps; ++i) {
			integratedDensity += m * lookupDensity(p, ray.d);
			m = 6 - m;

			#if defined(HETVOL_STATISTICS)
				++avgRayMarchingStepsTransmittance;
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

		return integratedDensity * m_scale
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
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return false;
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt - mint, maxComp = 0;
		Point p = ray(mint), pLast = ray(maxt);

		/* Ignore degenerate path segments */
		for (int i=0; i<3; ++i)
			maxComp = std::max(std::max(maxComp,
				std::abs(p[i])), std::abs(pLast[i]));
		if (length < 1e-6f * maxComp)
			return 0.0f;

		/* Compute a suitable step size (this routine samples the integrand
		   between steps, hence the factor of 2) */
		uint32_t nSteps = (uint32_t) std::ceil(length / (2*m_stepSize));
		Float stepSize = length / nSteps,
			  multiplier = (1.0f / 6.0f) * stepSize
				  * m_scale;
		Vector fullStep = ray.d * stepSize,
			   halfStep = fullStep * .5f;

		Float node1 = lookupDensity(p, ray.d);

		if (ray.mint == mint)
			densityAtMinT = node1 * m_scale;
		else
			densityAtMinT = 0.0f;

		#if defined(HETVOL_STATISTICS)
			avgRayMarchingStepsSampling.incrementBase();
		#endif

		for (uint32_t i=0; i<nSteps; ++i) {
			Float node2 = lookupDensity(p + halfStep, ray.d),
				  node3 = lookupDensity(p + fullStep, ray.d),
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
				   this point; instead, the density are modeled based on a
				   quadratic polynomial that is fit to the last three lookups */

				Float a = 0, b = stepSize, x = a,
					  fx = integratedDensity - desiredDensity,
					  stepSizeSqr = stepSize * stepSize,
					  temp = m_scale / stepSizeSqr;
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

	Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const {
		if (m_method == ESimpsonQuadrature || sampler == NULL) {
			return Spectrum(math::fastexp(-integrateDensity(ray)));
		} else {
			/* When Woodcock tracking is selected as the sampling method,
			   we can use this method to get a noisy (but unbiased) estimate
			   of the transmittance */
			Float mint, maxt;
			if (!m_densityAABB.rayIntersect(ray, mint, maxt))
				return Spectrum(1.0f);
			mint = std::max(mint, ray.mint);
			maxt = std::min(maxt, ray.maxt);

			#if defined(HETVOL_STATISTICS)
				avgRayMarchingStepsTransmittance.incrementBase();
			#endif
			int nSamples = 2; /// XXX make configurable
			Float result = 0;

			for (int i=0; i<nSamples; ++i) {
				Float t = mint;
				while (true) {
					t -= math::fastlog(1-sampler->next1D()) * m_invMaxDensity;
					if (t >= maxt) {
						result += 1;
						break;
					}

					Point p = ray(t);
					Float density = lookupDensity(p, ray.d) * m_scale;

					#if defined(HETVOL_STATISTICS)
						++avgRayMarchingStepsTransmittance;
					#endif

					if (density * m_invMaxDensity > sampler->next1D())
						break;
				}
			}
			return Spectrum(result/nSamples);
		}
	}

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		Float integratedDensity, densityAtMinT, densityAtT;
		bool success = false;

		if (m_method == ESimpsonQuadrature) {
			Float desiredDensity = -math::fastlog(1-sampler->next1D());
			if (invertDensityIntegral(ray, desiredDensity, integratedDensity,
					mRec.t, densityAtMinT, densityAtT)) {
				mRec.p = ray(mRec.t);
				success = true;
				Spectrum albedo = m_albedo->lookupSpectrum(mRec.p);
				mRec.sigmaS = albedo * densityAtT;
				mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
				mRec.orientation = m_orientation != NULL
					? m_orientation->lookupVector(mRec.p) : Vector(0.0f);
			}

			Float expVal = math::fastexp(-integratedDensity);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * densityAtT;
			mRec.pdfSuccessRev = expVal * densityAtMinT;
			mRec.transmittance = Spectrum(expVal);
			mRec.time = ray.time;
		} else {
			/* The following information is invalid when
			   using Woodcock-tracking */
			mRec.pdfFailure = 1.0f;
			mRec.pdfSuccess = 1.0f;
			mRec.pdfSuccessRev = 1.0f;
			mRec.transmittance = Spectrum(1.0f);
			mRec.time = ray.time;

			#if defined(HETVOL_STATISTICS)
				avgRayMarchingStepsSampling.incrementBase();
			#endif

			Float mint, maxt;
			if (!m_densityAABB.rayIntersect(ray, mint, maxt))
				return false;
			mint = std::max(mint, ray.mint);
			maxt = std::min(maxt, ray.maxt);

			Float t = mint, densityAtT = 0;
			while (true) {
				t -= math::fastlog(1-sampler->next1D()) * m_invMaxDensity;
				if (t >= maxt)
					break;

				Point p = ray(t);
				densityAtT = lookupDensity(p, ray.d) * m_scale;
				#if defined(HETVOL_STATISTICS)
					++avgRayMarchingStepsSampling;
				#endif
				if (densityAtT * m_invMaxDensity > sampler->next1D()) {
					mRec.t = t;
					mRec.p = p;
					Spectrum albedo = m_albedo->lookupSpectrum(p);
					mRec.sigmaS = albedo * densityAtT;
					mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
					mRec.transmittance = Spectrum(densityAtT != 0.0f ? 1.0f / densityAtT : 0);
					if (!std::isfinite(mRec.transmittance[0])) // prevent rare overflow warnings
						mRec.transmittance = Spectrum(0.0f);
					mRec.orientation = m_orientation != NULL
						? m_orientation->lookupVector(p) : Vector(0.0f);
					mRec.medium = this;
					success = true;
					break;
				}
			}
		}
		mRec.medium = this;

		return success && mRec.pdfSuccess > 0;
	}

	void eval(const Ray &ray, MediumSamplingRecord &mRec) const {
		if (m_method == ESimpsonQuadrature) {
			Float expVal = math::fastexp(-integrateDensity(ray));
			Float mintDensity = lookupDensity(ray(ray.mint), ray.d) * m_scale;
			Float maxtDensity = 0.0f;
			Spectrum maxtAlbedo(0.0f);
			if (ray.maxt < std::numeric_limits<Float>::infinity()) {
				Point p = ray(ray.maxt);
				maxtDensity = lookupDensity(p, ray.d) * m_scale;
				maxtAlbedo = m_albedo->lookupSpectrum(p);
			}
			mRec.transmittance = Spectrum(expVal);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * maxtDensity;
			mRec.pdfSuccessRev = expVal * mintDensity;
			mRec.sigmaS = maxtAlbedo * maxtDensity;
			mRec.sigmaA = Spectrum(maxtDensity) - mRec.sigmaS;
			mRec.time = ray.time;
			mRec.medium = this;
		} else {
			Log(EError, "eval(): unsupported integration method!");
		}
	}

	bool isHomogeneous() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousMedium[" << endl
			<< "  density = " << indent(m_density.toString()) << "," << endl
			<< "  albedo = " << indent(m_albedo.toString()) << "," << endl
			<< "  orientation = " << indent(m_orientation.toString()) << "," << endl
			<< "  stepSize = " << m_stepSize << "," << endl
			<< "  scale = " << m_scale << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	inline Float lookupDensity(const Point &p, const Vector &d) const {
		Float density = m_density->lookupFloat(p);
		if (m_anisotropicMedium && density != 0) {
			Vector orientation = m_orientation->lookupVector(p);
			if (!orientation.isZero())
				density *= m_phaseFunction->sigmaDir(dot(d, orientation));
			else
				return 0;
		}
		return density;
	}
protected:
	EIntegrationMethod m_method;
	ref<VolumeDataSource> m_density;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientation;
	Float m_scale;
	bool m_anisotropicMedium;
	Float m_stepSize;
	AABB m_densityAABB;
	Float m_maxDensity;
	Float m_invMaxDensity;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMedium, "Heterogeneous medium");
MTS_NAMESPACE_END
