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

#pragma once
#if !defined(__MITSUBA_CORE_SPECTRUM_H_)
#define __MITSUBA_CORE_SPECTRUM_H_

#include <mitsuba/mitsuba.h>

#if !defined(SPECTRUM_SAMPLES)
#error The desired number of spectral samples must be \
	specified in the configuration file!
#endif

#define SPECTRUM_MIN_WAVELENGTH   360
#define SPECTRUM_MAX_WAVELENGTH   830
#define SPECTRUM_RANGE                \
	(SPECTRUM_MAX_WAVELENGTH-SPECTRUM_MIN_WAVELENGTH)

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract continous spectral power distribution data type,
 * which supports evaluation at arbitrary wavelengths.
 *
 * Here, the term 'continous' doesn't necessarily mean that the
 * underlying spectrum is continous, but rather emphasizes the fact
 * that it is a function over the reals (as opposed to the discrete
 * spectrum, which only stores samples for a discrete set of wavelengths).
 *
 * \ingroup libpython
 * \ingroup libcore
 */
class MTS_EXPORT_CORE ContinuousSpectrum {
public:
	/**
	 * Evaluate the value of the spectral power distribution
	 * at the given wavelength.
	 *
	 * \param lambda  A wavelength in nanometers
	 */
	virtual Float eval(Float lambda) const = 0;

	/**
	 * \brief Integrate the spectral power distribution
	 * over a given interval and return the average value
	 *
	 * Unless overridden in a subclass, the integration is done
	 * using adaptive Gauss-Lobatto quadrature.
	 *
	 * \param lambdaMin
	 *     The lower interval bound in nanometers
	 *
	 * \param lambdaMax
	 *     The upper interval bound in nanometers
	 *
	 * \remark If \c lambdaMin >= \c lambdaMax, the
	 *     implementation will return zero.
	 */
	virtual Float average(Float lambdaMin, Float lambdaMax) const;

	/// \brief Return a string representation
	virtual std::string toString() const = 0;

	/// Virtual destructor
	virtual ~ContinuousSpectrum() { }
};

/**
 * \brief Spectral power distribution based on Planck's black body law
 *
 * Computes the spectral power distribution of a black body of the
 * specified temperature.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE BlackBodySpectrum : public ContinuousSpectrum {
public:
	/**
	 * \brief Construct a new black body spectrum given the emitter's
	 * temperature in Kelvin.
	 */
	inline BlackBodySpectrum(Float temperature) {
		m_temperature = temperature;
	}

	virtual ~BlackBodySpectrum() { }

	/** \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 *
	 * The units are Watts per unit surface area (m^-2)
	 * per unit wavelength (nm^-1) per steradian (sr^-1)
	 */
	virtual Float eval(Float lambda) const;

	/// Return a string representation
	std::string toString() const;
private:
	Float m_temperature;
};

/**
 * \brief Spectral distribution for rendering participating media
 * with Rayleigh scattering.
 *
 * This distribution captures the 1/lambda^4 wavelength dependence
 * of Rayleigh scattering. It can provide both the scattering and
 * extinction coefficient needed for simulating planetary
 * atmospheres with participating media.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE RayleighSpectrum : public ContinuousSpectrum {
public:
	enum EMode {
		/// Compute the scattering coefficient
		ESigmaS,
		/// Compute the extinction coefficient
		ESigmaT
	};

	/**
	 * \brief Create a Rayleigh spectrum instance
	 *
	 * \param mode        Specifies the requested type of spectrum
	 * \param eta         Refractive index of the medium (e.g. air)
	 * \param height      Height above sea level (in meters)
	 */
	RayleighSpectrum(EMode mode, Float eta = 1.000277f, Float height = 0);

	virtual ~RayleighSpectrum() { }

	/** \brief Evaluate the extinction/scattering coefficient for
	 * a specified wavelength.
	 *
	 * The returned value is in units of 1/meter.
	 */
	virtual Float eval(Float lambda) const;

	/// Return a string representation
	std::string toString() const;
private:
	Float m_precomp;
};

/**
 * \brief This spectral power distribution is defined as the
 * product of two other continuous spectra.
 */
class MTS_EXPORT_CORE ProductSpectrum : public ContinuousSpectrum {
public:
	/** \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 */
	ProductSpectrum(const ContinuousSpectrum &s1,
		const ContinuousSpectrum &s2) : m_spec1(s1),
	    m_spec2(s2) { }

	/** \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 */
	virtual Float eval(Float lambda) const;

	/// Virtual destructor
	virtual ~ProductSpectrum() { }

	/// Return a string representation
	std::string toString() const;
private:
	const ContinuousSpectrum &m_spec1;
	const ContinuousSpectrum &m_spec2;
};

/**
 * \brief Linearly interpolated spectral power distribution
 *
 * This class implements a linearly interpolated spectral
 * power distribution that is defined over a discrete set of
 * measurements at different wavelengths. Outside of the
 * specified range, the spectrum is assumed to be zero. Hence,
 * at least two entries are required to produce a nonzero
 * spectrum.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE InterpolatedSpectrum : public ContinuousSpectrum {
public:
	/**
	 * \brief Create a new interpolated spectrum with space
	 * for the specified number of samples
	 */
	InterpolatedSpectrum(size_t size = 0);

	/**
	 * \brief Create a interpolated spectrum instance from
	 * a float array
	 */
	InterpolatedSpectrum(const Float *wavelengths,
		const Float *values, size_t nEntries);

	/**
	 * \brief Read an interpolated spectrum from a simple
	 * ASCII format.
	 *
	 * Each line of the file should contain an entry of the form
	 * \verbatim
	 * <wavelength in nm> <value>
	 * \endverbatim
	 * Comments preceded by '#' are also valid.
	 */
	InterpolatedSpectrum(const fs::path &path);

	/**
	 * \brief Append an entry to the spectral power distribution.
	 *
	 * Entries must be added in order of increasing wavelength
	 */
	void append(Float lambda, Float value);

	/**
	 * \brief This function adds a zero entry before and after
	 * the stored wavelength range.
	 *
	 * This is useful when handling datasets that don't fall
	 * off to zero at the ends. The spacing of the added entries
	 * is determined by computing the average spacing of the
	 * existing samples.
	 */
	void zeroExtend();

	/// Clear all stored entries
	void clear();

	/**
	 * \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 */
	Float eval(Float lambda) const;

	/**
	 * \brief Integrate the spectral power distribution
	 * over a given interval and return the average value
	 *
	 * This method overrides the implementation in
	 * \ref ContinousSpectrum, since the integral can be
	 * analytically computed for linearly interpolated spectra.
	 *
	 * \param lambdaMin
	 *     The lower interval bound in nanometers
	 *
	 * \param lambdaMax
	 *     The upper interval bound in nanometers
	 *
	 * \remark If \c lambdaMin >= \c lambdaMax, the
	 *     implementation will return zero.
	 */
	Float average(Float lambdaMin, Float lambdaMax) const;

	/// \brief Return a string representation
	std::string toString() const;

	/// Virtual destructor
	virtual ~InterpolatedSpectrum() { }
protected:
	std::vector<Float> m_wavelengths, m_values;
};

/**
 * \brief Abstract spectral power distribution data type
 *
 * This class defines a vector-like data type that can be used for
 * computations involving radiance. A concrete instantiation for the
 * precision and spectral discretization chosen at compile time is
 * given by the \ref Spectrum data type.
 *
 * \ingroup libcore
 */
template <typename T, int N> struct TSpectrum {
public:
	typedef T          Scalar;

	/// Number of dimensions
	const static int dim = N;

	/// Create a new spectral power distribution, but don't initialize the contents
#if !defined(MTS_DEBUG_UNINITIALIZED)
	inline TSpectrum() { }
#else
	inline TSpectrum() {
		for (int i=0; i<N; i++)
			s[i] = std::numeric_limits<Scalar>::quiet_NaN();
	}
#endif

	/// Create a new spectral power distribution with all samples set to the given value
	explicit inline TSpectrum(Scalar v) {
		for (int i=0; i<N; i++)
			s[i] = v;
	}

	/// Copy a spectral power distribution
	explicit inline TSpectrum(Scalar spec[N]) {
		memcpy(s, spec, sizeof(Scalar)*N);
	}

	/// Unserialize a spectral power distribution from a binary data stream
	explicit inline TSpectrum(Stream *stream) {
		stream->readArray(s, N);
	}

	/// Initialize with a TSpectrum data type based on a alternate representation
	template <typename AltScalar> explicit TSpectrum(const TSpectrum<AltScalar, N> &v) {
		for (int i=0; i<N; ++i)
			s[i] = (Scalar) v[i];
	}

	/// Add two spectral power distributions
	inline TSpectrum operator+(const TSpectrum &spec) const {
		TSpectrum value = *this;
		for (int i=0; i<N; i++)
			value.s[i] += spec.s[i];
		return value;
	}

	/// Add a spectral power distribution to this instance
	inline TSpectrum& operator+=(const TSpectrum &spec) {
		for (int i=0; i<N; i++)
			s[i] += spec.s[i];
		return *this;
	}

	/// Subtract a spectral power distribution
	inline TSpectrum operator-(const TSpectrum &spec) const {
		TSpectrum value = *this;
		for (int i=0; i<N; i++)
			value.s[i] -= spec.s[i];
		return value;
	}

	/// Subtract a spectral power distribution from this instance
	inline TSpectrum& operator-=(const TSpectrum &spec) {
		for (int i=0; i<N; i++)
			s[i] -= spec.s[i];
		return *this;
	}

	/// Multiply by a scalar
	inline TSpectrum operator*(Scalar f) const {
		TSpectrum value = *this;
		for (int i=0; i<N; i++)
			value.s[i] *= f;
		return value;
	}

	/// Multiply by a scalar
	inline friend TSpectrum operator*(Scalar f, const TSpectrum &spec) {
		return spec * f;
	}

	/// Multiply by a scalar
	inline TSpectrum& operator*=(Scalar f) {
		for (int i=0; i<N; i++)
			s[i] *= f;
		return *this;
	}

	/// Perform a component-wise multiplication by another spectrum
	inline TSpectrum operator*(const TSpectrum &spec) const {
		TSpectrum value = *this;
		for (int i=0; i<N; i++)
			value.s[i] *= spec.s[i];
		return value;
	}

	/// Perform a component-wise multiplication by another spectrum
	inline TSpectrum& operator*=(const TSpectrum &spec) {
		for (int i=0; i<N; i++)
			s[i] *= spec.s[i];
		return *this;
	}

	/// Perform a component-wise division by another spectrum
	inline TSpectrum& operator/=(const TSpectrum &spec) {
		for (int i=0; i<N; i++)
			s[i] /= spec.s[i];
		return *this;
	}

	/// Perform a component-wise division by another spectrum
	inline TSpectrum operator/(const TSpectrum &spec) const {
		TSpectrum value = *this;
		for (int i=0; i<N; i++)
			value.s[i] /= spec.s[i];
		return value;
	}

	/// Divide by a scalar
	inline TSpectrum operator/(Scalar f) const {
		TSpectrum value = *this;
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "TSpectrum: Division by zero!");
#endif
		Scalar recip = 1.0f / f;
		for (int i=0; i<N; i++)
			value.s[i] *= recip;
		return value;
	}

	/// Equality test
	inline bool operator==(const TSpectrum &spec) const {
		for (int i=0; i<N; i++) {
			if (s[i] != spec.s[i])
				return false;
		}
		return true;
	}

	/// Inequality test
	inline bool operator!=(const TSpectrum &spec) const {
		return !operator==(spec);
	}

	/// Divide by a scalar
	inline friend TSpectrum operator/(Scalar f, TSpectrum &spec) {
		return TSpectrum(f) / spec;
	}

	/// Divide by a scalar
	inline TSpectrum& operator/=(Scalar f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "TTSpectrum: Division by zero!");
#endif
		Scalar recip = 1.0f / f;
		for (int i=0; i<N; i++)
			s[i] *= recip;
		return *this;
	}

	/// Check for NaNs
	inline bool isNaN() const {
		for (int i=0; i<N; i++)
			if (std::isnan(s[i]))
				return true;
		return false;
	}

	/// Returns whether the spectrum only contains valid (non-NaN, nonnegative) samples
	inline bool isValid() const {
		for (int i=0; i<N; i++)
			if (!std::isfinite(s[i]) || s[i] < 0.0f)
				return false;
		return true;
	}

	/// Multiply-accumulate operation, adds \a weight * \a spec
	inline void addWeighted(Scalar weight, const TSpectrum &spec) {
		for (int i=0; i<N; i++)
			s[i] += weight * spec.s[i];
	}

	/// Return the average over all wavelengths
	inline Scalar average() const {
		Scalar result = 0.0f;
		for (int i=0; i<N; i++)
			result += s[i];
		return result * (1.0f / N);
	}

	/// Component-wise absolute value
	inline TSpectrum abs() const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = std::abs(s[i]);
		return value;
	}

	/// Component-wise square root
	inline TSpectrum sqrt() const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = std::sqrt(s[i]);
		return value;
	}

	/// Component-wise square root
	inline TSpectrum safe_sqrt() const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = math::safe_sqrt(s[i]);
		return value;
	}

	/// Component-wise logarithm
	inline TSpectrum log() const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = math::fastlog(s[i]);
		return value;
	}

	/// Component-wise exponentation
	inline TSpectrum exp() const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = math::fastexp(s[i]);
		return value;
	}

	/// Component-wise power
	inline TSpectrum pow(Scalar f) const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = std::pow(s[i], f);
		return value;
	}

	/// Clamp negative values
	inline void clampNegative() {
		for (int i=0; i<N; i++)
			s[i] = std::max((Scalar) 0.0f, s[i]);
	}

	/// Return the highest-valued spectral sample
	inline Scalar max() const {
		Scalar result = s[0];
		for (int i=1; i<N; i++)
			result = std::max(result, s[i]);
		return result;
	}

	/// Return the lowest-valued spectral sample
	inline Scalar min() const {
		Scalar result = s[0];
		for (int i=1; i<N; i++)
			result = std::min(result, s[i]);
		return result;
	}

	/// Negate
	inline TSpectrum operator-() const {
		TSpectrum value;
		for (int i=0; i<N; i++)
			value.s[i] = -s[i];
		return value;
	}

	/// Indexing operator
	inline Scalar &operator[](int entry) {
		return s[entry];
	}

	/// Indexing operator
	inline Scalar operator[](int entry) const {
		return s[entry];
	}

	/// Check if this spectrum is zero at all wavelengths
	inline bool isZero() const {
		for (int i=0; i<N; i++) {
			if (s[i] != 0.0f)
				return false;
		}
		return true;
	}

	/// Serialize this spectrum to a stream
	inline void serialize(Stream *stream) const {
		stream->writeArray(s, N);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "[";
		for (int i=0; i<N; i++) {
			oss << s[i];
			if (i < N - 1)
				oss << ", ";
		}
		oss << "]";
		return oss.str();
	}

protected:
	Scalar s[N];
};


/** \brief RGB color data type
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE Color3 : public TSpectrum<Float, 3> {
public:
	typedef TSpectrum<Float, 3> Parent;

	/// Create a new color value, but don't initialize the contents
#if !defined(MTS_DEBUG_UNINITIALIZED)
	inline Color3() { }
#else
	inline Color3() {
		for (int i=0; i<3; i++)
			s[i] = std::numeric_limits<Scalar>::quiet_NaN();
	}
#endif

	/// Copy constructor
	inline Color3(const Parent &s) : Parent(s) { }

	/// Initialize to a constant value
	inline Color3(Float value) : Parent(value) { }

	/// Initialize to the given RGB value
	inline Color3(Float r, Float g, Float b) {
		s[0] = r; s[1] = g; s[2] = b;
	}

	/// Return the luminance (assuming the color value is expressed in linear sRGB)
	inline Float getLuminance() const {
		return s[0] * 0.212671f + s[1] * 0.715160f + s[2] * 0.072169f;
	}
};


/** \brief Discrete spectral power distribution based on a number
 * of wavelength bins over the 360-830 nm range.
 *
 * This class defines a vector-like data type that can be used for
 * computations involving radiance.
 *
 * When configured for spectral rendering (i.e. when the compile-time flag
 * \c SPECTRUM_SAMPLES is set to a value != 3), the implementation discretizes
 * the visible spectrum of light into a set of intervals, where the
 * distribution within each bin is modeled as being uniform.
 *
 * When SPECTRUM_SAMPLES == 3, the class reverts to a simple linear
 * RGB-based internal representation.
 *
 * The implementation of this class is based on PBRT.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE Spectrum : public TSpectrum<Float, SPECTRUM_SAMPLES> {
public:
	typedef TSpectrum<Float, SPECTRUM_SAMPLES> Parent;

	/**
	 * \brief When converting from RGB reflectance values to
	 * discretized color spectra, the following `intent' flag
	 * can be provided to improve the results of this highly
	 * under-constrained problem.
	 */
	enum EConversionIntent {
		/// Unitless reflectance data is converted
		EReflectance,

		/// Radiance-valued illumination data is converted
		EIlluminant
	};

	/// Create a new spectral power distribution, but don't initialize the contents
#if !defined(MTS_DEBUG_UNINITIALIZED)
	inline Spectrum() { }
#else
	inline Spectrum() {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] = std::numeric_limits<Scalar>::quiet_NaN();
	}
#endif

	/// Construct from a TSpectrum instance
	inline Spectrum(const Parent &s) : Parent(s) { }

	/// Initialize with a TSpectrum data type based on a alternate representation
	template <typename AltScalar> explicit Spectrum(const TSpectrum<AltScalar, SPECTRUM_SAMPLES> &v) {
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			s[i] = (Scalar) v[i];
	}

	/// Create a new spectral power distribution with all samples set to the given value
	explicit inline Spectrum(Float v) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] = v;
	}

	/// Copy a spectral power distribution
	explicit inline Spectrum(Float value[SPECTRUM_SAMPLES]) {
		memcpy(s, value, sizeof(Float)*SPECTRUM_SAMPLES);
	}

	/// Unserialize a spectral power distribution from a binary data stream
	explicit inline Spectrum(Stream *stream) : Parent(stream) { }

	/**
	 * \brief Evaluate the SPD for the given wavelength
	 * in nanometers.
	 */
	Float eval(Float lambda) const;

	/// \brief Return the wavelength range covered by a spectral bin
	static std::pair<Float, Float> getBinCoverage(size_t index);

	/// Return the luminance in candelas.
#if SPECTRUM_SAMPLES == 3
	inline Float getLuminance() const {
		return s[0] * 0.212671f + s[1] * 0.715160f + s[2] * 0.072169f;
	}
#else
	Float getLuminance() const;
#endif

	/**
	 * \brief Convert from a spectral power distribution to XYZ
	 * tristimulus values
	 *
	 * In the Python API, this function returns a 3-tuple
	 * with the result of the operation.
	 */
	void toXYZ(Float &x, Float &y, Float &z) const;

	/**
	 * \brief Convert XYZ tristimulus into a plausible spectral
	 * reflectance or spectral power distribution
	 *
	 * The \ref EConversionIntent parameter can be used to provide more
	 * information on how to solve this highly under-constrained problem.
	 * The default is \ref EReflectance.
	 */
	void fromXYZ(Float x, Float y, Float z,
			EConversionIntent intent = EReflectance);

	/**
	 * \brief Convert from a spectral power distribution to
	 * the perceptually uniform IPT color space by Ebner and Fairchild
	 *
	 * This is useful e.g. for computing color differences.
	 * \c I encodes intensity, \c P (protan) roughly encodes
	 * red-green color opponency, and \c T (tritan) encodes
	 * blue-red color opponency. For normalized input, the
	 * range of attainable values is given by
	 * \f$ I\in $[0,1], P,T\in [-1,1]\f$.
	 *
	 * In the Python API, this function returns a 3-tuple
	 * with the result of the operation.
	 */
	void toIPT(Float &I, Float &P, Float &T) const;

	/**
	 * \brief Convert a color value represented in the IPT
	 * space into a plausible spectral reflectance or
	 * spectral power distribution.
	 *
	 * The \ref EConversionIntent parameter can be used to provide more
	 * information on how to solve this highly under-constrained problem.
	 * The default is \ref EReflectance.
	 */
	void fromIPT(Float I, Float P, Float T,
			EConversionIntent intent = EReflectance);

#if SPECTRUM_SAMPLES == 3
	/**
	 * \brief Convert to linear RGB
	 *
	 * In the Python API, this function returns a 3-tuple
	 * with the result of the operation.
	 */
	inline void toLinearRGB(Float &r, Float &g, Float &b) const {
		/* Nothing to do -- the renderer is in RGB mode */
		r = s[0]; g = s[1]; b = s[2];
	}

	/// Convert from linear RGB
	inline void fromLinearRGB(Float r, Float g, Float b,
			EConversionIntent intent = EReflectance /* unused */) {
		/* Nothing to do -- the renderer is in RGB mode */
		s[0] = r; s[1] = g; s[2] = b;
	}
#else
	/**
	 * \brief Convert to linear RGB
	 *
	 * In the Python API, this function returns a 3-tuple
	 * with the result of the operation.
	 */
	void toLinearRGB(Float &r, Float &g, Float &b) const;

	/**
	 * \brief Convert linear RGB colors into a plausible
	 * spectral power distribution
	 *
	 * The \ref EConversionIntent parameter can be used to provide more
	 * information on how to solve this highly under-constrained problem.
	 * The default is \ref EReflectance.
	 */
	void fromLinearRGB(Float r, Float g, Float b,
			EConversionIntent intent = EReflectance);
#endif

	/**
	 * \brief Convert to sRGB
	 *
	 * In the Python API, this function returns a 3-tuple
	 * with the result of the operation.
	 */
	void toSRGB(Float &r, Float &g, Float &b) const;

	/**
	 * \brief Convert sRGB color values into a plausible spectral
	 * power distribution
	 *
	 * Note that compared to \ref fromLinearRGB, no \c intent parameter
	 * is available. For sRGB colors, it is assumed that the intent is
	 * always \ref EReflectance.
	 */
	void fromSRGB(Float r, Float g, Float b);

	/**
	 * \brief Convert linear RGBE colors into a plausible
	 * spectral power distribution
	 *
	 * Based on code by Bruce Walter and Greg ward.
	 *
	 * The \ref EConversionIntent parameter can be used to provide more
	 * information on how to solve this highly under-constrained problem.
	 * For RGBE values, the default is \ref EIlluminant.
	 */
	void fromRGBE(const uint8_t rgbe[4], EConversionIntent intent = EIlluminant);

	/// Linear RGBE conversion based on Bruce Walter's and Greg Ward's code
	void toRGBE(uint8_t rgbe[4]) const;

	/// Initialize with spectral values from a smooth spectrum representation
	void fromContinuousSpectrum(const ContinuousSpectrum &smooth);

	/// Equality test
	inline bool operator==(const Spectrum &val) const {
		for (int i=0; i<SPECTRUM_SAMPLES; i++) {
			if (s[i] != val.s[i])
				return false;
		}
		return true;
	}

	/// Inequality test
	inline bool operator!=(const Spectrum &val) const {
		return !operator==(val);
	}

	/// Return a string representation
	std::string toString() const;

	/**
	 * \brief Return a spectral color distribution of the
	 * D65 white point (with unit luminance)
	 */
	inline static const Spectrum &getD65() { return CIE_D65; }

	/**
	 * \brief Static initialization (should be called once during the
	 * application's initialization phase)
	 *
	 * This function is responsible for choosing the wavelengths
	 * that will be used during rendering. It also pre-integrates
	 * the CIE matching curves so that sampled spectra can
	 * efficiently be converted to XYZ tristimulus values.
	 * Finally, it sets up pre-integrated color spectra for conversions
	 * from linear RGB to plausible spectral color distributions.
	 */
	static void staticInitialization();
	static void staticShutdown();
protected:
	#if SPECTRUM_SAMPLES != 3
	/// Configured wavelengths bins in nanometers
	static Float m_wavelengths[SPECTRUM_SAMPLES+1];

	/// @{ \name Pre-integrated CIE 1931 XYZ color matching functions.
	static Spectrum CIE_X;
	static Spectrum CIE_Y;
	static Spectrum CIE_Z;
	static Float CIE_normalization;
	/// @}

	/**
	 * @{ \name Pre-integrated Smits-style RGB to Spectrum
	 * conversion spectra, data by Karl vom Berge
	 */
	static Spectrum rgbRefl2SpecWhite;
	static Spectrum rgbRefl2SpecCyan;
	static Spectrum rgbRefl2SpecMagenta;
	static Spectrum rgbRefl2SpecYellow;
	static Spectrum rgbRefl2SpecRed;
	static Spectrum rgbRefl2SpecGreen;
	static Spectrum rgbRefl2SpecBlue;
	static Spectrum rgbIllum2SpecWhite;
	static Spectrum rgbIllum2SpecCyan;
	static Spectrum rgbIllum2SpecMagenta;
	static Spectrum rgbIllum2SpecYellow;
	static Spectrum rgbIllum2SpecRed;
	static Spectrum rgbIllum2SpecGreen;
	static Spectrum rgbIllum2SpecBlue;
	/// @}
	#endif

	/// Pre-integrated D65 illuminant
	static Spectrum CIE_D65;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SPECTRUM_H_ */
