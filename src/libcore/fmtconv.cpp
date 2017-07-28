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

// The default mpl::vector limit of 20 is not enough
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

#include <mitsuba/core/bitmap.h>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/transform.hpp>

MTS_NAMESPACE_BEGIN

namespace mpl = boost::mpl;

/****************************************************************************/
/*  This file contains an implementation of the FormatConverter interface   */
/*  in mitsuba/render/bitmap.h. Ideally we'd like to have optimized code    */
/*  for each possible combination of source & target pixel and component    */
/*  formats. The switch() and Boost MPL craziness below does exactly this:  */
/*  it produces code for each possible pair                                 */
/****************************************************************************/
/*  Note that some cases in this file really scream for SSE optimization    */
/*  (in case they ever become a bottleneck, that is..)                      */
/****************************************************************************/

namespace detail {
    /* Mapping from a C++ type to Bitmap::EComponentFormat */

    template <typename T> struct get_pixelformat {
        enum { value = Bitmap::EInvalid };
    };
    template <> struct get_pixelformat<uint8_t> {
        enum { value = Bitmap::EUInt8 };
    };
    template <> struct get_pixelformat<uint16_t> {
        enum { value = Bitmap::EUInt16 };
    };
    template <> struct get_pixelformat<uint32_t> {
        enum { value = Bitmap::EUInt32 };
    };
    template <> struct get_pixelformat<half> {
        enum { value = Bitmap::EFloat16 };
    };
    template <> struct get_pixelformat<float> {
        enum { value = Bitmap::EFloat32 };
    };
    template <> struct get_pixelformat<double> {
        enum { value = Bitmap::EFloat64 };
    };

    // Safely convert from size_t to other types, avoiding downcasting warning
    template <typename T, typename S> inline T safe_cast(S a) {
        return static_cast<T>(a);
    }
    template <> inline half safe_cast(size_t a) {
        return static_cast<half>(static_cast<float>(a));
    }
    template <> inline half safe_cast(double a) {
        return static_cast<half>(static_cast<float>(a));
    }
}

template <typename T> struct FormatConverterImpl : public FormatConverter {
    typedef typename T::first  SourceFormat;
    typedef typename T::second DestFormat;

    template <typename FormatType>
    struct format_traits
    {
        // Is the type a floating point-based format
        static const bool is_float =
            boost::is_same<FormatType, double>::value ||
            boost::is_same<FormatType, float>::value ||
            boost::is_same<FormatType, half>::value;

        // Does the type have a "compact" representation? (i.e. uint8_t/uint16_t)
        static const bool is_compact =
            boost::is_same<FormatType, uint8_t>::value ||
            boost::is_same<FormatType, uint16_t>::value;
    };

    virtual Conversion getConversion() const {
        return std::make_pair(
            (Bitmap::EComponentFormat) detail::get_pixelformat<SourceFormat>::value,
            (Bitmap::EComponentFormat) detail::get_pixelformat<DestFormat>::value
        );
    }

    virtual void convert(
            Bitmap::EPixelFormat sourceFormat, Float sourceGamma, const void *_source,
            Bitmap::EPixelFormat destFormat, Float destGamma, void *_dest,
            size_t count, Float multiplier, Spectrum::EConversionIntent intent, int channelCount) const {

        #if 0
            std::ostringstream oss;
            oss << "FormatConverter::convert([" << sourceFormat << ", "
                << (Bitmap::EComponentFormat) detail::get_pixelformat<SourceFormat>::value
                << "] -> [" << destFormat << ", "
                << (Bitmap::EComponentFormat) detail::get_pixelformat<DestFormat>::value
                << "], gamma = " << sourceGamma << " -> " << destGamma
                << ", count = " << count << ")";
            SLog(EInfo, "%s", oss.str().c_str());
        #endif

        /* Revert to memcpy when the underlying data needs no transformation */
        if ((int) detail::get_pixelformat<SourceFormat>::value == (int) detail::get_pixelformat<DestFormat>::value &&
            sourceFormat == destFormat && sourceGamma == destGamma && multiplier == 1.0) {
            switch (sourceFormat) {
                case Bitmap::ELuminance:            channelCount = 1; break;
                case Bitmap::ELuminanceAlpha:       channelCount = 2; break;
                case Bitmap::ERGB:
                case Bitmap::EXYZ:                  channelCount = 3; break;
                case Bitmap::EXYZA:
                case Bitmap::ERGBA:                 channelCount = 4; break;
                case Bitmap::ESpectrum:             channelCount = SPECTRUM_SAMPLES; break;
                case Bitmap::ESpectrumAlpha:        channelCount = SPECTRUM_SAMPLES + 1; break;
                case Bitmap::ESpectrumAlphaWeight:  channelCount = SPECTRUM_SAMPLES + 2; break;
                case Bitmap::EMultiChannel:         break;
                default:
                    SLog(EError, "Unsupported source/target pixel format!");
                    return;
            }
            memcpy(_dest, _source, sizeof(SourceFormat) * channelCount * count);
            return;
        }

        const SourceFormat *source = reinterpret_cast<const SourceFormat *>(_source);
        DestFormat *dest = reinterpret_cast<DestFormat *>(_dest);
        const Float invDestGamma = 1.0f / destGamma;
        const size_t maxValue = (size_t) std::numeric_limits<SourceFormat>::max();

        DestFormat *precomp = NULL;
        if (format_traits<SourceFormat>::is_compact && count > maxValue) {
            /* When transforming an uint8_t or uint16_t-based image, it
               will generally be much cheaper to precompute a table for this ahead of time */
            precomp = (DestFormat *) alloca(sizeof(DestFormat) * (maxValue+1));
            for (size_t i=0; i <= maxValue; ++i)
                precomp[i] = convertScalar<DestFormat>(detail::safe_cast<SourceFormat>(i), sourceGamma, NULL, multiplier, invDestGamma);
        }

        const DestFormat one = convertScalar<DestFormat>(1.0f);

        Spectrum spec;

        switch (sourceFormat) {
            case Bitmap::ELuminance: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i)
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = value; *dest++ = value; *dest++ = value;
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = value; *dest++ = value; *dest++ = value; *dest++ = one;
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                Float value = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(value * 0.950456f, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value * 1.08875f, 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                Float value = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(value * 0.950456f, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value * 1.08875f, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = value;
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = value;
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = value;
                                *dest++ = one; *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::ELuminanceAlpha: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = value; *dest++ = value; *dest++ = value;
                                source++;
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = value; *dest++ = value; *dest++ = value;
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                Float value = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(value * 0.950456f, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value * 1.08875f, 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                Float value = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(value * 0.950456f, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(value * 1.08875f, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = value;
                                source++;
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = value;
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                DestFormat value = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = value;
                                *dest++ = convertScalar<DestFormat>(*source++);
                                *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::ERGB: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Float luminance = RGBToLuminance(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Float luminance = RGBToLuminance(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Color3 xyz = RGBToXYZ(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(xyz[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[2], 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Color3 xyz = RGBToXYZ(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(xyz[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[2], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromLinearRGB(r, g, b, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromLinearRGB(r, g, b, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromLinearRGB(r, g, b, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one; *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::ERGBA: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Float luminance = RGBToLuminance(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Float luminance = RGBToLuminance(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Color3 xyz = RGBToXYZ(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(xyz[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[2], 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                Color3 xyz = RGBToXYZ(Color3(r, g, b));
                                *dest++ = convertScalar<DestFormat>(xyz[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(xyz[2], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromLinearRGB(r, g, b, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromLinearRGB(r, g, b, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                Float r = convertScalar<Float>(*source++, sourceGamma);
                                Float g = convertScalar<Float>(*source++, sourceGamma);
                                Float b = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromLinearRGB(r, g, b, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                                *dest++ = convertScalar<DestFormat>(1.0f);
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::EXYZ: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                Float luminance = convertScalar<Float>(source[1], sourceGamma);
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                source += 3;
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float luminance = convertScalar<Float>(source[1], sourceGamma);
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                                source += 3;
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                Color3 rgb = XYZToRGB(Color3(x, y, z));
                                *dest++ = convertScalar<DestFormat>(rgb[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[2], 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                Color3 rgb = XYZToRGB(Color3(x, y, z));
                                *dest++ = convertScalar<DestFormat>(rgb[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[2], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromXYZ(x, y, z, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromXYZ(x, y, z, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromXYZ(x, y, z, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one; *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::EXYZA: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                Float luminance = convertScalar<Float>(source[1], sourceGamma);
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                source += 4;
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float luminance = convertScalar<Float>(source[1], sourceGamma);
                                *dest++ = convertScalar<DestFormat>(luminance, 1.0f, NULL, multiplier, invDestGamma);
                                source += 3;
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                Color3 rgb = XYZToRGB(Color3(x, y, z));
                                *dest++ = convertScalar<DestFormat>(rgb[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[2], 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                Color3 rgb = XYZToRGB(Color3(x, y, z));
                                *dest++ = convertScalar<DestFormat>(rgb[0], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[1], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(rgb[2], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                source++;
                                spec.fromXYZ(x, y, z, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromXYZ(x, y, z, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                Float x = convertScalar<Float>(*source++, sourceGamma);
                                Float y = convertScalar<Float>(*source++, sourceGamma);
                                Float z = convertScalar<Float>(*source++, sourceGamma);
                                spec.fromXYZ(x, y, z, intent);
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                                *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;


            case Bitmap::ESpectrum: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(spec.getLuminance(), 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(spec.getLuminance(), 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float r, g, b;
                                spec.toLinearRGB(r, g, b);
                                *dest++ = convertScalar<DestFormat>(r, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(g, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(b, 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float r, g, b;
                                spec.toLinearRGB(r, g, b);
                                *dest++ = convertScalar<DestFormat>(r, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(g, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(b, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float x, y, z;
                                spec.toXYZ(x, y, z);
                                *dest++ = convertScalar<DestFormat>(x, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(y, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(z, 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float x, y, z;
                                spec.toXYZ(x, y, z);
                                *dest++ = convertScalar<DestFormat>(x, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(y, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(z, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0, n = count*SPECTRUM_SAMPLES; i<n; ++i)
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = one;
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = one; *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::ESpectrumAlpha: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(spec.getLuminance(), 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                *dest++ = convertScalar<DestFormat>(spec.getLuminance(), 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float r, g, b;
                                spec.toLinearRGB(r, g, b);
                                *dest++ = convertScalar<DestFormat>(r, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(g, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(b, 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float r, g, b;
                                spec.toLinearRGB(r, g, b);
                                *dest++ = convertScalar<DestFormat>(r, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(g, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(b, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float x, y, z;
                                spec.toXYZ(x, y, z);
                                *dest++ = convertScalar<DestFormat>(x, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(y, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(z, 1.0f, NULL, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float x, y, z;
                                spec.toXYZ(x, y, z);
                                *dest++ = convertScalar<DestFormat>(x, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(y, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(z, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                source++;
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                                *dest++ = one;
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

            case Bitmap::ESpectrumAlphaWeight: {
                    switch (destFormat) {
                        case Bitmap::ELuminance:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                source++;
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                *dest++ = convertScalar<DestFormat>(spec.getLuminance()*invWeight, 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ELuminanceAlpha:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float alpha = convertScalar<Float>(*source++);
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                *dest++ = convertScalar<DestFormat>(spec.getLuminance()*invWeight, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(alpha * invWeight);
                            }
                            break;

                        case Bitmap::ERGB:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                source++;
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                Float r, g, b;
                                Spectrum(spec * invWeight).toLinearRGB(r, g, b);
                                *dest++ = convertScalar<DestFormat>(r, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(g, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(b, 1.0f, NULL, multiplier, invDestGamma);
                            }
                            break;

                        case Bitmap::ERGBA:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float alpha = convertScalar<Float>(*source++);
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                Float r, g, b;
                                Spectrum(spec * invWeight).toLinearRGB(r, g, b);
                                *dest++ = convertScalar<DestFormat>(r, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(g, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(b, 1.0f, NULL, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(alpha * invWeight);
                            }
                            break;

                        case Bitmap::EXYZ:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                source++;
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                Float x, y, z;
                                Spectrum(spec * invWeight * multiplier).toXYZ(x, y, z);
                                *dest++ = convertScalar<DestFormat>(x, 1.0f, NULL, 1.0f, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(y, 1.0f, NULL, 1.0f, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(z, 1.0f, NULL, 1.0f, invDestGamma);
                            }
                            break;

                        case Bitmap::EXYZA:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float alpha = convertScalar<Float>(*source++);
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                Float x, y, z;
                                Spectrum(spec * invWeight * multiplier).toXYZ(x, y, z);
                                *dest++ = convertScalar<DestFormat>(x, 1.0f, NULL, 1.0f, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(y, 1.0f, NULL, 1.0f, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(z, 1.0f, NULL, 1.0f, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(alpha * invWeight);
                            }
                            break;

                        case Bitmap::ESpectrum:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                ++source;
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier*invWeight, invDestGamma);
                            }
                            break;

                        case Bitmap::ESpectrumAlpha:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    spec[j] = convertScalar<Float>(*source++, sourceGamma);
                                Float alpha = convertScalar<Float>(*source++);
                                Float weight = convertScalar<Float>(*source++), invWeight = (weight != 0) ? 1 / weight : weight;
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(spec[j], 1.0f, NULL, multiplier*invWeight, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(alpha * invWeight);
                            }
                            break;

                        case Bitmap::ESpectrumAlphaWeight:
                            for (size_t i=0; i<count; ++i) {
                                for (int j=0; j<SPECTRUM_SAMPLES; ++j)
                                    *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                                *dest++ = convertScalar<DestFormat>(*source++);
                                *dest++ = convertScalar<DestFormat>(*source++);
                            }
                            break;

                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;

                case Bitmap::EMultiChannel: {
                    switch (destFormat) {
                        case Bitmap::EMultiChannel:
                            for (size_t i=0; i<count*channelCount; ++i)
                                *dest++ = convertScalar<DestFormat>(*source++, sourceGamma, precomp, multiplier, invDestGamma);
                            break;
                        default:
                            SLog(EError, "Unsupported destination pixel format!");
                    }
                }
                break;


            default:
                SLog(EError, "Unsupported source pixel format!");
        }
    }

private:
    static Float undoGamma(Float value, Float gamma) {
        if (gamma == -1) {
            if (value <= (Float) 0.04045)
                return value * (Float) (1.0 / 12.92);
            else
                return std::pow((Float) ((value + (Float) 0.055) * (Float) (1.0 / 1.055)), (Float) 2.4);
        } else {
            return std::pow(value, gamma);
        }
    }

    static Float applyGamma(Float value, Float invGamma) {
        if (invGamma == -1) {
            return (value <= (Float) 0.0031308) ? ((Float) 12.92 * value)
                : ((Float) 1.055 * std::pow(value, (Float) (1.0/2.4)) - (Float) 0.055);
        } else {
            return std::pow(value, invGamma);
        }
    }

    /// Convert ITU-R Rec. BT.709 linear RGB to a luminance value
    inline static Float RGBToLuminance(const Color3 &value) {
        return value[0] * (Float) 0.212671 + value[1] * (Float) 0.715160 + value[2] * (Float) 0.072169;
    }

    /// Convert ITU-R Rec. BT.709 linear RGB to XYZ tristimulus values
    inline static Color3 RGBToXYZ(const Color3 &value) {
        Color3 result;
        result[0] = value[0] * (Float) 0.412453 + value[1] * (Float) 0.357580 + value[2] * (Float) 0.180423;
        result[1] = value[0] * (Float) 0.212671 + value[1] * (Float) 0.715160 + value[2] * (Float) 0.072169;
        result[2] = value[0] * (Float) 0.019334 + value[1] * (Float) 0.119193 + value[2] * (Float) 0.950227;
        return result;
    }

    /// Convert from XYZ tristimulus values to ITU-R Rec. BT.709 linear RGB
    inline static Color3 XYZToRGB(const Color3 &value) {
        Color3 result;
        result[0] = value[0] * (Float)  3.240479 + value[1] * (Float) -1.537150 + value[2] * (Float) -0.498535;
        result[1] = value[0] * (Float) -0.969256 + value[1] * (Float)  1.875991 + value[2] * (Float)  0.041556;
        result[2] = value[0] * (Float)  0.055648 + value[1] * (Float) -0.204043 + value[2] * (Float)  1.057311;
        return result;
    }

    template <typename DestFmt, typename SourceFmt>
    inline static DestFmt convertScalar(SourceFmt source, Float sourceGamma = 1.0f,
            DestFmt *precomp = NULL, Float multiplier = 1.0f, Float invDestGamma = 1.0f) {
        if (format_traits<SourceFmt>::is_compact && precomp)
            return precomp[(size_t) source];

        Float value = static_cast<Float>(source);

        if (!format_traits<SourceFmt>::is_float)
            value *= static_cast<Float>(1.0f/std::numeric_limits<SourceFmt>::max());

        if (sourceGamma != 1)
            value = undoGamma(value, sourceGamma);

        value *= multiplier;

        if (invDestGamma != 1)
            value = applyGamma(value, invDestGamma);

        if (format_traits<DestFmt>::is_float)
            return detail::safe_cast<DestFmt> (value);
        else /* Round to nearest value and clamp to representable range */
            return detail::safe_cast<DestFmt> (std::min(static_cast<Float>(std::numeric_limits<DestFmt>::max()),
                std::max((Float) 0, value * (Float) std::numeric_limits<DestFmt>::max() + (Float) 0.5f)));
    }
};

/* ================================================
    The following Boost MPL magic is responsible
    for generating code that efficiently converts
    from any image type to any other image type,
    while applying a specified Gamma transformation
 * ================================================ */

/* List of supported conversion types */
typedef mpl::vector<
    uint8_t, uint16_t, uint32_t, half, float, double> SupportedTypes;

/* Meta-function that generates all type-pairs with T
   as the first element */
template <typename T, typename Initial> struct ConversionFn
    : mpl::fold<SupportedTypes, Initial,
        mpl::push_back<mpl::_1, mpl::pair<T, mpl::_2> > > { };

/* Generate all possible type pairs */
typedef mpl::fold<
    SupportedTypes, mpl::vector<>,
        mpl::lambda<ConversionFn<mpl::_2, mpl::_1> > >::type
    SupportedConversions;

/* Instantiate the converter class with all possible arguments */
typedef mpl::transform<
    SupportedConversions, mpl::lambda<FormatConverterImpl<mpl::_1> > >::type
    ConverterImplementations;

/// Helper class, which registers a FormatConverter at program startup
struct RegisterConverter {
public:
    inline RegisterConverter(FormatConverter::ConverterMap &map) : m_converters(map) { }

    template <class T> void operator() (T) const {
        FormatConverter *converter = new T();
        FormatConverter::Conversion cType = converter->getConversion();
        m_converters[cType] = converter;
    }
private:
    FormatConverter::ConverterMap &m_converters;
};

FormatConverter::ConverterMap FormatConverter::m_converters;

void FormatConverter::staticInitialization() {
    mpl::for_each<ConverterImplementations>(RegisterConverter(m_converters));
}

void FormatConverter::staticShutdown() {
    for (ConverterMap::iterator it = m_converters.begin();
            it != m_converters.end(); ++it)
        delete it->second;
    m_converters.clear();
}

const FormatConverter *FormatConverter::getInstance(Conversion types) {
    if (m_converters.find(types) == m_converters.end()) {
        std::ostringstream oss;
        oss << "Unable to find a format converter from '" << types.first << "' to '" << types.second << "'!";
        SLog(EError, "%s", oss.str().c_str());
    }

    return m_converters[types];
}

MTS_NAMESPACE_END
