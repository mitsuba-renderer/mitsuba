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

#if !defined(__SPECTRUM_H)
#define __SPECTRUM_H

#include <mitsuba/mitsuba.h>

#define SPECTRUM_MIN_WAVELENGTH   400
#define SPECTRUM_MAX_WAVELENGTH   700
#define SPECTRUM_RANGE            (SPECTRUM_MAX_WAVELENGTH-SPECTRUM_MIN_WAVELENGTH+1)
#define SPECTRUM_SAMPLES          3

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract smooth spectral power distribution data type,
 * which supports evaluation at arbitrary wavelengths
 */
class MTS_EXPORT_CORE SmoothSpectrum {
public:
	/**
	 * Evaluate the value of the spectral power distribution
	 * at the given wavelength.
	 */
	virtual Float eval(Float lambda) const = 0;
	
	virtual ~SmoothSpectrum() { }
};

/**
 * \brief Spectral power distribution based on Planck's black body law
 *
 * Computes the spectral power distribution of a black body of the 
 * specified temperature.
 */
class MTS_EXPORT_CORE BlackBodySpectrum : public SmoothSpectrum {
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
	 */
	virtual Float eval(Float lambda) const;
private:
	Float m_temperature;
};

/**
 * \brief Linearly interpolated spectral power distribution
 */
class MTS_EXPORT_CORE InterpolatedSpectrum : public SmoothSpectrum {
public:
	/**
	 * \brief Create a new interpolated spectrum with space 
	 * for the specified number of samples
	 */
	inline InterpolatedSpectrum(size_t size = 0) {
		m_wavelength.reserve(size);
		m_value.reserve(size);
	}

	/**
	 * \brief Append an entry to the spectral power distribution.
	 *
	 * Entries must be added in order of increasing wavelength
	 */
	void appendSample(Float lambda, Float value);

	/**
	 * \brief Return the value of the spectral power distribution
	 * at the given wavelength.
	 */
	virtual Float eval(Float lambda) const;
	
	virtual ~InterpolatedSpectrum() { }
private:
	std::vector<Float> m_wavelength, m_value;
};

/** \brief Discrete spectral power distribution 
 * based on a (usually small) number of wavelength samples. 
 *
 * When SPECTRUM_SAMPLES is set to 3 (the default), this class 
 * falls back to linear RGB as its internal representation.
 */
struct MTS_EXPORT_CORE Spectrum {
public:
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
		if (f == 0) {
			SLog(EWarn, "Spectrum: Division by zero!");
		//	exit(-1);
		}
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
		if (f == 0) {
			SLog(EWarn, "Spectrum: Division by zero!");
			//exit(-1);
		}
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
		return result / SPECTRUM_SAMPLES;
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
			value.s[i] = std::exp(s[i]);
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
	 * \brief Evaluate the SPD at an arbitrary wavelength
	 * (uses interpolation)
	 */
	Float eval(Float lambda) const;

	/// Return the luminance in candelas.
#if SPECTRUM_SAMPLES == 3
	inline Float getLuminance() const {
		return s[0] * 0.212671f + s[1] * 0.715160f + s[2] * 0.072169f;
	}
#else
	Float getLuminance() const;
#endif

	/// Convert from a spectral power distribution to XYZ colors
	void toXYZ(Float &x, Float &y, Float &z) const;

	/// Convert from XYZ to a spectral power distribution
	void fromXYZ(Float x, Float y, Float z);

#if SPECTRUM_SAMPLES == 3
	/// Convert to linear RGB
	inline void toLinearRGB(Float &r, Float &g, Float &b) const {
		r = s[0];
		g = s[1];
		b = s[2];
	}

	/// Convert from linear RGB
	inline void fromLinearRGB(Float r, Float g, Float b) {
		s[0] = r;
		s[1] = g;
		s[2] = b;
	}
#else
	/// Convert to linear RGB
	void toLinearRGB(Float &r, Float &g, Float &b) const;

	/// Convert from linear RGB
	void fromLinearRGB(Float r, Float g, Float b);	
#endif

	/// Convert to sRGB
	void toSRGB(Float &r, Float &g, Float &b) const;

	/// Convert from sRGB
	void fromSRGB(Float r, Float g, Float b);

	/// Linear RGBE conversion based on Bruce Walter's and Greg Ward's code
	void fromRGBE(const uint8_t rgbe[4]);

	/// Linear RGBE conversion based on Bruce Walter's and Greg Ward's code 
	void toRGBE(uint8_t rgbe[4]) const;

	/// Initialize with spectral values from a smooth spectrum representation
	void fromSmoothSpectrum(const SmoothSpectrum *smooth);

	/// Serialize this spectrum to a stream
	inline void serialize(Stream *stream) const {
		stream->writeFloatArray(s, SPECTRUM_SAMPLES);
	}

	/// Return the wavelength corresponding to an index
	inline static Float getWavelength(int index) {
		SAssert(index < SPECTRUM_SAMPLES);
		return m_wavelengths[index];
	}

	/// Return a string representation
	std::string toString() const;

	/** 
	 * Static initialization (should be called once during the
	 * application's initialization phase
	 */
	static void staticInitialization();
	static void staticShutdown();
protected:
	Float s[SPECTRUM_SAMPLES];

	/// Configured wavelengths in nanometers
	static Float m_wavelengths[SPECTRUM_SAMPLES];

	/// Normalization factor for XYZ<->RGB conversion
	static Float m_normalization;

	/// Inverse of \ref m_normalization
	static Float m_invNormalization;

	/**
	 * @{ \name CIE 1931 XYZ color matching functions. 
	 * From http://www.cvrl.org/database/data/cmfs/ciexyz31_1.txt
	 */
	static const int   CIE_start = 360;
	static const int   CIE_end   = 830;
	static const int   CIE_count = CIE_end - CIE_start + 1;
	static const Float CIE_normalization;
	static const Float CIE_X[CIE_count];
	static const Float CIE_Y[CIE_count];
	static const Float CIE_Z[CIE_count];
	/// @}
};

MTS_NAMESPACE_END

#endif /* __SPECTRUM_H */
