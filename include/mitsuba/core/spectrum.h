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

#if !defined(__SPECTRUM_H)
#define __SPECTRUM_H

#include <mitsuba/mitsuba.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

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
struct MTS_EXPORT_CORE Spectrum {
public:
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
			s[i] = std::numeric_limits<double>::quiet_NaN();
	}
#endif

	/// Create a new spectral power distribution with all samples set to the given value
	explicit inline Spectrum(Float v) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] = v;
	}

	/// Copy a spectral power distribution
	explicit inline Spectrum(Float spd[SPECTRUM_SAMPLES]) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] = spd[i];
	}

	/// Unserialize a spectral power distribution from a binary data stream
	explicit inline Spectrum(Stream *stream) {
		stream->readFloatArray(s, SPECTRUM_SAMPLES);
	}

	/// Add two spectral power distributions
	inline Spectrum operator+(const Spectrum &spd) const {
		Spectrum value = *this;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] += spd.s[i];
		return value;
	}

	/// Add a spectral power distribution to this instance
	inline Spectrum& operator+=(const Spectrum &spd) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] += spd.s[i];
		return *this;
	}

	/// Subtract a spectral power distribution
	inline Spectrum operator-(const Spectrum &spd) const {
		Spectrum value = *this;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] -= spd.s[i];
		return value;
	}

	/// Subtract a spectral power distribution from this instance
	inline Spectrum& operator-=(const Spectrum &spd) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] -= spd.s[i];
		return *this;
	}

	/// Multiply by a scalar
	inline Spectrum operator*(Float f) const {
		Spectrum value = *this;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] *= f;
		return value;
	}

	/// Multiply by a scalar
	inline friend Spectrum operator*(Float f, Spectrum &spd) {
		return spd * f;
	}

	/// Multiply by a scalar
	inline Spectrum& operator*=(Float f) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] *= f;
		return *this;
	}

	/// Perform a component-wise multiplication by another spectrum
	inline Spectrum operator*(const Spectrum &spd) const {
		Spectrum value = *this;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] *= spd.s[i];
		return value;
	}

	/// Perform a component-wise multiplication by another spectrum
	inline Spectrum& operator*=(const Spectrum &spd) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] *= spd.s[i];
		return *this;
	}

	/// Perform a component-wise division by another spectrum
	inline Spectrum& operator/=(const Spectrum &spd) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] /= spd.s[i];
		return *this;
	}
	
	/// Perform a component-wise division by another spectrum
	inline Spectrum operator/(Spectrum spd) const {
		Spectrum value = *this;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] /= spd.s[i];
		return value;
	}

	/// Divide by a scalar
	inline Spectrum operator/(Float f) const {
		Spectrum value = *this;
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Spectrum: Division by zero!");
#endif
		Float recip = 1.0f / f;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] *= recip;
		return value;
	}

	/// Equality test
	inline bool operator==(Spectrum spd) const {
		for (int i=0; i<SPECTRUM_SAMPLES; i++) {
			if (s[i] != spd.s[i])
				return false;
		}
		return true;
	}

	/// Inequality test
	inline bool operator!=(Spectrum spd) const {
		return !operator==(spd);
	}

	/// Divide by a scalar
	inline friend Spectrum operator/(Float f, Spectrum &spd) {
		return spd / f;
	}

	/// Divide by a scalar
	inline Spectrum& operator/=(Float f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Spectrum: Division by zero!");
#endif
		Float recip = 1.0f / f;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] *= recip;
		return *this;
	}

	/// Check for NaNs
	inline bool isNaN() const {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			if (mts_isnan(s[i]))
				return true;
		return false;
	}

	/// Returns whether the spectrum only contains valid (non-NaN, nonnegative) samples
	inline bool isValid() const {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			if (mts_isnan(s[i]) || s[i] < 0.0f)
				return false;
		return true;
	}

	/// Multiply-accumulate operation, adds \a weight * \a spd
	inline void addWeighted(Float weight, const Spectrum &spd) {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] += weight * spd.s[i];
	}

	/// Return the average over all wavelengths
	inline Float average() const {
		Float result = 0.0f;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			result += s[i];
		return result * (1.0f / SPECTRUM_SAMPLES);
	}

	/// Component-wise square root
	inline Spectrum sqrt() const {
		Spectrum value;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] = std::sqrt(s[i]);
		return value;
	}

	/// Component-wise exponentation
	inline Spectrum exp() const {
		Spectrum value;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] = std::fastexp(s[i]);
		return value;
	}

	/// Component-wise power
	inline Spectrum pow(Float f) const {
		Spectrum value;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] = std::pow(s[i], f);
		return value;
	}

	/// Clamp negative values
	inline void clampNegative() {
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			s[i] = std::max((Float) 0.0f, s[i]);
	}

	/// Return the highest-valued spectral sample
	inline Float max() const {
		Float result = s[0];
		for (int i=1; i<SPECTRUM_SAMPLES; i++)
			result = std::max(result, s[i]);
		return result;
	}

	/// Return the lowest-valued spectral sample
	inline Float min() const {
		Float result = s[0];
		for (int i=1; i<SPECTRUM_SAMPLES; i++)
			result = std::min(result, s[i]);
		return result;
	}

	/// Negate
	inline Spectrum operator-() const {
		Spectrum value;
		for (int i=0; i<SPECTRUM_SAMPLES; i++)
			value.s[i] = -s[i];
		return value;
	}

	/// Indexing operator
	inline Float &operator[](int entry) {
		return s[entry];
	}

	/// Indexing operator
	inline Float operator[](int entry) const {
		return s[entry];
	}

	/// Check if this spectrum is zero at all wavelengths
	inline bool isZero() const {
		for (int i=0; i<SPECTRUM_SAMPLES; i++) {
			if (s[i] != 0.0f)
				return false;
		}	
		return true;
	}

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

	/// Convert from a spectral power distribution to XYZ tristimulus values
	void toXYZ(Float &x, Float &y, Float &z) const;

	/**
	 * \brief Convert XYZ tristimulus into a plausible spectral
	 * power distribution
	 *
	 * The \ref EConversionIntent parameter can be used to provide more 
	 * information on how to solve this highly under-constrained problem.
	 * The default is \ref EReflectance.
	 */
	void fromXYZ(Float x, Float y, Float z, 
			EConversionIntent intent = EReflectance);

#if SPECTRUM_SAMPLES == 3
	/// Convert to linear RGB
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
	/// Convert to linear RGB
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

	/// Convert to sRGB
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

	/// Serialize this spectrum to a stream
	inline void serialize(Stream *stream) const {
		stream->writeFloatArray(s, SPECTRUM_SAMPLES);
	}

	/// Return a string representation
	std::string toString() const;

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
	Float s[SPECTRUM_SAMPLES];

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
};

MTS_NAMESPACE_END

#endif /* __SPECTRUM_H */
