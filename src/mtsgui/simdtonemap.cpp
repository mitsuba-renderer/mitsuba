/*============================================================================
  HDRITools - High Dynamic Range Image Tools
  Copyright 2008-2012 Program of Computer Graphics, Cornell University

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
 -----------------------------------------------------------------------------
 Primary author:
     Edgar Velazquez-Armendariz <cs#cornell#edu - eva5>
============================================================================*/

#include "simdtonemap.h"

#if defined(MTS_SSE)

#include <mitsuba/core/sse.h>
#include <mitsuba/core/ssemath.h>
#include <mitsuba/core/ssevector.h>
#include <mitsuba/core/bitmap.h>


namespace
{

typedef mitsuba::math::SSEVector4f V4f;
typedef mitsuba::math::SSEVector4i V4i;

// Method for scaling the luminance
enum ELuminanceMethod {
    EExposure,
    EReinhard02
};

// Non-linear display transform
enum EDisplayMethod {
    ESRGB,
    EGamma
};

// Group of 4 AoS RGBA32F pixels: R,G,B,A, R,G,B,A, R,G,B,A, R,G,B,A
struct PixelRGBA32FGroup
{
    V4f p[4];

    inline static PixelRGBA32FGroup zero() {
        PixelRGBA32FGroup group;
        group.p[0] = group.p[1] = group.p[2] = group.p[3] = V4f::zero();
        return group;
    }
};

union PixelRGBA8
{
    struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
        uint8_t a;
    };
    uint32_t abgr;
};

union PixelRGBA8Group
{
    PixelRGBA8 p[4];
    __m128i packed;

    inline static PixelRGBA8Group zero() {
        PixelRGBA8Group p;
        p.packed = V4i::zero();
        return p;
    }
};


// Inline version of fastlog_ps for extra performance
inline mitsuba::math::SSEVector4f am_log(const mitsuba::math::SSEVector4f& x) {
    typedef mitsuba::math::SSEVector4f v4f;
    typedef mitsuba::math::SSEVector4i v4i;
    using mitsuba::math::castAsFloat;
    using mitsuba::math::clamp_ps;
    using mitsuba::math::fastpow_ps;

    // Constants
    const v4f min_normal(castAsFloat(v4i::constant<0x00800000>()));
    const v4f inv_mantissa_mask(castAsFloat(v4i::constant<~0x7f800000>()));
    const v4f const_1(1.0f);
    const v4i const_127(v4i::constant<127>());

    const v4f log_p0(-7.89580278884799154124e-1f);
    const v4f log_p1( 1.63866645699558079767e1f);
    const v4f log_p2(-6.41409952958715622951e1f);

    const v4f log_q0(-3.56722798256324312549e1f);
    const v4f log_q1( 3.12093766372244180303e2f);
    const v4f log_q2(-7.69691943550460008604e2f);

    const v4f log_c0( 0.693147180559945f);


    // Kill invalid values (ignoring NaN and Inf)
    const v4f x0 = max(x, min_normal);

    // Kill the exponent and combine with the exponent of 1.0f to get the
    // actual embedded mantissa as a valid floating point value:
    // a value in the range [1.0, 2.0)
    const v4f mantissa = (x0 & inv_mantissa_mask) | const_1;

    const v4f v_min1  = mantissa - const_1;
    const v4f v_plus1 = mantissa + const_1;

    // Extract the original exponent and undo the bias
    const v4i biasedExponent = srl(castAsInt(x0), 23);
    const v4f origExponent = toFloat(biasedExponent - const_127);

    v4f vFrac = v_min1 * rcp(v_plus1); // Is it worth it to use rcp_nr?
    vFrac += vFrac;
    const v4f vFracSqr = vFrac * vFrac;

    // Evaluate the polynomial
    const v4f polyP = ((((log_p0 * vFracSqr) + log_p1) * vFracSqr)
                                             + log_p2) * vFracSqr;
    const v4f polyQ =  (((log_q0 * vFracSqr) + log_q1) * vFracSqr) + log_q2;

    const v4f poly = polyP * rcp(polyQ); // Use rcp_nr?
    const v4f logApprox = poly * vFrac;

    // Scale by log(2) to get the natural logarithm of the exponent part
    const v4f logExpPart = origExponent * log_c0;

    // Combine the different parts
    const v4f result = logApprox + vFrac + logExpPart;

    return result;
}


/* Fast sRGB curve calculation using a rational approximation for the
non-linear term. Assumes the input in in the range [0,1]

Mathematica input:
-------------------------------------------------------------------------------
<< FunctionApproximations`
MiniMaxApproximation[-(11/200) + (211 x^(5/12))/
  200, {x, {0.0031308, 1}, 4, 3}]
CForm[HornerForm[%[[2, 1]]]]
-------------------------------------------------------------------------------
Result:
(-0.016036752726326525 + x*(23.24653363361083 +
       x*(1832.6027368173256 + x*(10602.877994753313 + 2764.157524016198*x))))/
 (1. + x*(255.72605859770067 + x*(5415.6856291461045 + 9542.777488625074*x)))

Maximum relative error: -0.000504731
Abscissas where the relative error is a local maximum:
   0.0031308, 0.00392026, 0.00737952, 0.0181293, 0.0499015,
   0.139431,  0.361749,   0.751249,   1.0
*/
inline V4f sRGB_Remez43(const V4f& x) throw()
{
    // Constants
    const V4f P0(    -0.016036752726326525f );
    const V4f P1(    23.24653363361083f );
    const V4f P2(  1832.6027368173256f );
    const V4f P3( 10602.877994753313f );
    const V4f P4(  2764.157524016198f );

    const V4f Q0(     1.f );
    const V4f Q1(   255.72605859770067f );
    const V4f Q2( 5415.6856291461045f  );
    const V4f Q3( 9542.777488625074f  );

    const V4f CUTOFF(0.0031308f);
    const V4f LINEAR_FACTOR(12.92f);

    // Rational approximation using the Remez algorithm in Mathematica to
    //   1.055 * pow(x, 1.0/2.4) - 0.055
    const V4f num = (P0 + x*(P1 + x*(P2 + x*(P3 + x*P4))));
    const V4f den = (Q0 + x*(Q1 + x*(Q2 + x*Q3)));
    const V4f p   = num * rcp(den);

    const V4f result = select(x < CUTOFF, x*LINEAR_FACTOR, p);
    return result;
}


// Applies the global Reinhard-2002 TMO. The parameters are calculated
// separately by a different process.
//
// The canonical approach is
//   a. Transform sRGB to xyY
//   b. Apply the TMO to Y
//   c. Transform x,y,TMO(Y) back to sRGB
//
// However, having only Y and assuming that TMO(Y) == k*Y, then the
// result of all the transformation is just k*[r,g,b]
// Thus:
//         (key/avgLogLum) * (1 + (key/avgLogLum)/pow(Lwhite,2) * Y)
//    k == ---------------------------------------------------------
//                        1 + (key/avgLogLum)*Y
//
//    k == (P * (R + Q*(P*Y)) / (R + P*Y)
//    P == key / avgLogLum   --> scale
//    Q == 1 / pow(Lwhite,2) --> invWp2
//    R == 1
//
inline void reinhard(V4f& r, V4f& g, V4f& b, const V4f& multiplier,
    const V4f& scale, const V4f& invWp2)
{
    const V4f ONE(1.0f);
    const V4f LVec0(0.357580f);
    const V4f LVec1(0.715160f);
    const V4f LVec2(0.119193f);

    // Get the luminance
    const V4f Y = multiplier * (LVec0*r + LVec1*g + LVec2*b);

    // Compute the scale factor
    const V4f Lp = scale * Y;
    V4f k = (scale * (ONE + invWp2*Lp)) * rcp_nr(ONE + Lp);
    k *= multiplier;

    // And apply
    r *= k;
    g *= k;
    b *= k;
}



template <ELuminanceMethod luminanceMethod, EDisplayMethod displayMethod>
void tonemap(const PixelRGBA32FGroup* const begin,
    const PixelRGBA32FGroup* const end, PixelRGBA8Group* const target,
    const TonemapCPU::Params& params)
{
    // Use 24 KiB for temporary block storage
    const ptrdiff_t BLOCK_SIZE = ((24*1024)/4) / sizeof(V4f);

    // SoA temporary storage
    V4f buffer[BLOCK_SIZE * 4];

    // Constants
    const V4f invGamma(params.invGamma);
    const V4f invWhitePoint(params.invWhitePoint);
    const V4f invWp2(params.invWhitePoint * params.invWhitePoint);
    const V4f multiplier(params.multiplier);
    const V4f scale(params.scale);

    const V4f const_1f(1.0f);
    const V4f const_255f(255.0f);

    // Actual output iterator
    PixelRGBA8Group* out = target;

    for (const PixelRGBA32FGroup* it = begin; it != end;) {
         const ptrdiff_t numIter = std::min((end-it), BLOCK_SIZE);
         const ptrdiff_t numColorIter = 3 * numIter;
         const ptrdiff_t rOffset = 0;
         const ptrdiff_t gOffset = rOffset + numIter;
         const ptrdiff_t bOffset = gOffset + numIter;
         const ptrdiff_t aOffset = bOffset + numIter;

        // Convert the input data from AoS to SoA, advancing the input iterator
        for (ptrdiff_t i = 0; i != numIter; ++i, ++it) {
            buffer[rOffset+i] = it->p[0];
            buffer[gOffset+i] = it->p[1];
            buffer[bOffset+i] = it->p[2];
            buffer[aOffset+i] = it->p[3];
            transpose(buffer[rOffset+i], buffer[gOffset+i],
                buffer[bOffset+i], buffer[aOffset+i]);
            buffer[aOffset+i] = clamp_ps(buffer[aOffset+i],
                V4f::zero(), const_1f);
        }

        // Scale luminance using only the color components
        if (luminanceMethod == EReinhard02) {
            for (ptrdiff_t i = 0; i != numIter; ++i) {
                reinhard(buffer[rOffset+i],buffer[gOffset+i],buffer[bOffset+i],
                    multiplier, scale, invWp2);
            }
        } else if (luminanceMethod == EExposure) {
            for (ptrdiff_t i = 0; i != numColorIter; ++i) {
                buffer[i] *= invWhitePoint;
            }
        }

        // Clamp and nonlinear display transform of the color components
        if (displayMethod == ESRGB) {
            for (ptrdiff_t i = 0; i != numColorIter; ++i) {
                const V4f s = clamp_ps(buffer[i], V4f::zero(), const_1f);
                buffer[i] = sRGB_Remez43(s);
            }
        } else if (displayMethod == EGamma) {
            for (ptrdiff_t i = 0; i != numColorIter; ++i) {
                const V4f s = clamp_ps(buffer[i], V4f::zero(), const_1f);
                buffer[i] = fastpow_ps(s, invGamma);
            }
        }

        // Quantize and build the pixel, advancing the output iterator
        for (ptrdiff_t i = 0; i != numIter; ++i, ++out) {
            V4i rQ = roundToInt(const_255f * buffer[rOffset+i]);
            V4i gQ = roundToInt(const_255f * buffer[gOffset+i]);
            V4i bQ = roundToInt(const_255f * buffer[bOffset+i]);
            V4i aQ = roundToInt(const_255f * buffer[aOffset+i]);

            V4i gShift = sll(gQ,  8);
            V4i bShift = sll(bQ, 16);
            V4i aShift = sll(aQ, 24);
            V4i pixel = aShift | bShift | gShift | rQ;
            stream(&(out->packed), pixel);
        }
    }
}



// Tonemap dispatcher
template <ELuminanceMethod luminanceMethod, EDisplayMethod displayMethod>
bool tonemap(const mitsuba::Bitmap* source, mitsuba::Bitmap* target,
    const TonemapCPU::Params& params)
{
    using namespace mitsuba;

    if (source->getSize() != target->getSize()) {
        SLog(EWarn, "TonemapCPU: images size missmatch");
        return false;
    }

    // Check the format
    if (source->getPixelFormat() != Bitmap::ERGBA ||
        target->getPixelFormat() != Bitmap::ERGBA) {
        SLog(EWarn, "TonemapCPU: the images are not in RGBA format");
        return false;
    }
    else if (source->getComponentFormat() != Bitmap::EFloat32) {
        SLog(EWarn, "TonemapCPU: the source component format is not Float32");
        return false;
    }
    else if (target->getComponentFormat() != Bitmap::EUInt8) {
        SLog(EWarn, "TonemapCPU: the target component format is not UInt8");
        return false;
    }

    // Raw pointers to the data, checking for alignment
    const float *sourceData = source->getFloat32Data();
    uint8_t *targetData = static_cast<uint8_t*>(target->getData());
    if (reinterpret_cast<uintptr_t>(sourceData) % 16 != 0) {
        SLog(EWarn, "TonemapCPU: the source data is not 16-byte aligned");
        return false;
    }
    else if (reinterpret_cast<uintptr_t>(targetData) % 16 != 0) {
        SLog(EWarn, "TonemapCPU: the target data is not 16-byte aligned");
        return false;
    }

    // So far so good, now reinterpret the data as groups of 4 pixels. We may
    // have an incomplete final group which needs separate handling
    const size_t pixelCount = source->getPixelCount();
    const PixelRGBA32FGroup* begin =
        reinterpret_cast<const PixelRGBA32FGroup*>(sourceData);
    const PixelRGBA32FGroup* end   = begin + (pixelCount / 4);
    PixelRGBA8Group* dest = reinterpret_cast<PixelRGBA8Group*>(targetData);

    // Main processing
    tonemap<luminanceMethod, displayMethod> (begin, end, dest, params);

    if (pixelCount % 4 != 0) {
        // Individual final group
        PixelRGBA32FGroup last = PixelRGBA32FGroup::zero();
        const size_t offset = pixelCount & ~0x3;
        const V4f* pixels = reinterpret_cast<const V4f*>(sourceData) + offset;
        for (int i = 0; i < static_cast<int>(pixelCount % 4); ++i) {
            last.p[i] = pixels[i];
        }

        PixelRGBA8Group lastTarget = PixelRGBA8Group::zero();
        tonemap<luminanceMethod, displayMethod> (&last, (&last)+1,
            &lastTarget, params);

        // Copy the final pixels to the target
        PixelRGBA8* dest = reinterpret_cast<PixelRGBA8*>(targetData) + offset;
        for (int i = 0; i < static_cast<int>(pixelCount % 4); ++i) {
            dest[i].abgr = lastTarget.p[i].abgr;
        }
    }

    return true;
};



struct LuminanceResult
{
    V4f sumLogLuminance;
    V4f maxLuminance;
};


// Compute the luminance mapping NaNs, negatives and the Mitsuba logo to 0.0
inline V4f computeLuminance(const V4f& r, const V4f& g, const V4f& b,
    const V4f& multiplier)
{
    const V4f L0(0.212671f);
    const V4f L1(0.715160f);
    const V4f L2(0.072169f);
    const V4f const_1024f(1024.0f);

    const V4f Y = multiplier * (L0*r + L1*g + L2*b);
    const V4f invalidMask = (Y < V4f::zero()) | isnan(Y) | (Y == const_1024f);

    // 0.0 is just zeros, so the mask works
    const V4f result = andnot(invalidMask, Y);
    return result;
}

LuminanceResult luminance(const PixelRGBA32FGroup* const begin,
    const PixelRGBA32FGroup* const end, float inMultiplier)
{
    // Use 24 KiB for temporary block storage
    const ptrdiff_t BLOCK_SIZE = (24*1024) / sizeof(V4f);
    V4f buffer[BLOCK_SIZE];

    // Constants
    const V4f multiplier(inMultiplier);
    const V4f const_1e_3f(1e-3f);

    // Totals
    V4f sumLogLuminance = V4f::zero();
    V4f maxLuminance(-1.0f);        // All luminances will be at least zero

    for (const PixelRGBA32FGroup* it = begin; it != end;) {
        const ptrdiff_t numIter = std::min((end - it), BLOCK_SIZE);

        // Compute the luminance and advance the input iterator
        for (ptrdiff_t i = 0; i != numIter; ++i, ++it) {
            V4f r    = it->p[0];
            V4f g    = it->p[1];
            V4f b    = it->p[2];
            V4f junk = it->p[3];
            transpose(r, g, b, junk);
            const V4f luminance = computeLuminance(r, g, b, multiplier);
            maxLuminance = max(luminance, maxLuminance);
            buffer[i] = luminance;
        }

        // Calculate and accumulate the log-luminance
        for (ptrdiff_t i = 0; i != numIter; ++i) {
            const V4f logLuminance = am_log(buffer[i] + const_1e_3f);
            sumLogLuminance += logLuminance;
        }
    }

    LuminanceResult result;
    result.sumLogLuminance = sumLogLuminance;
    result.maxLuminance    = maxLuminance;
    return result;
}



// Luminance dispatcher
bool luminance(const mitsuba::Bitmap* source, const float multiplier,
    float& outMaxLuminance, float &outAvgLogLuminance)
{
    using namespace mitsuba;

    // Check the format
    if (source->getPixelFormat() != Bitmap::ERGBA) {
        SLog(EWarn, "TonemapCPU: the image is not in RGBA format");
        return false;
    }
    else if (source->getComponentFormat() != Bitmap::EFloat32) {
        SLog(EWarn, "TonemapCPU: the image component format is not Float32");
        return false;
    }

    // Raw pointers to the data, checking for alignment
    const float *sourceData = source->getFloat32Data();
    if (reinterpret_cast<uintptr_t>(sourceData) % 16 != 0) {
        SLog(EWarn, "TonemapCPU: the source data is not 16-byte aligned");
        return false;
    }

    // So far so good, now reinterpret the data as groups of 4 pixels. We may
    // have an incomplete final group which needs separate handling
    const size_t pixelCount = source->getPixelCount();
    const PixelRGBA32FGroup* begin =
        reinterpret_cast<const PixelRGBA32FGroup*>(sourceData);
    const PixelRGBA32FGroup* end   = begin + (pixelCount / 4);

    // Main processing
    LuminanceResult result = luminance(begin, end, multiplier);

    if (pixelCount % 4 != 0) {
        // Individual final group
        PixelRGBA32FGroup last = PixelRGBA32FGroup::zero();
        const size_t offset = pixelCount & ~0x3;
        const V4f* pixels = reinterpret_cast<const V4f*>(sourceData) + offset;
        for (int i = 0; i < static_cast<int>(pixelCount % 4); ++i) {
            last.p[i] = pixels[i];
        }
        LuminanceResult rTail = luminance(&last, (&last)+1, multiplier);

        // Remove the invalid results
        V4i tailMask = V4i::zero();
        switch (pixelCount % 4) {
        case 1:
            tailMask = V4i::constant<0, 0, 0, -1>();
            break;
        case 2:
            tailMask = V4i::constant<0, 0, -1, -1>();
            break;
        case 3:
            tailMask = V4i::constant<0, -1, -1, -1>();
            break;
        }
        // maxLuminance contains zeros in the invalid positions which is OK
        rTail.sumLogLuminance &= castAsFloat(tailMask);

        result.sumLogLuminance += rTail.sumLogLuminance;
        result.maxLuminance = max(rTail.maxLuminance, result.maxLuminance);
    }

    using mitsuba::math::hmax_ps;
    using mitsuba::math::hsum_ps;
    using mitsuba::math::fastexp;

    outMaxLuminance = hmax_ps(result.maxLuminance);
    const float sumLogLuminance = hsum_ps(result.sumLogLuminance);
    outAvgLogLuminance =  fastexp(sumLogLuminance / pixelCount);
    return true;
};

} // namespace



bool TonemapCPU::gammaTonemap(const mitsuba::Bitmap* source,
    mitsuba::Bitmap* target) const
{
    if (m_params.isSRGB) {
        return tonemap<EExposure, ESRGB> (source, target, m_params);
    } else {
        return tonemap<EExposure, EGamma>(source, target, m_params);
    }
}


bool TonemapCPU::reinhardTonemap(const mitsuba::Bitmap* source,
    mitsuba::Bitmap* target) const
{
    if (m_params.isSRGB) {
        return tonemap<EReinhard02, ESRGB> (source, target, m_params);
    } else {
        return tonemap<EReinhard02, EGamma>(source, target, m_params);
    }
}


bool TonemapCPU::setLuminanceInfo(const mitsuba::Bitmap* source,
    mitsuba::Float multiplier)
{
    return luminance(source, static_cast<float>(multiplier),
        m_params.maxLum, m_params.avgLogLum);
}


using mitsuba::Object;
MTS_IMPLEMENT_CLASS_I(TonemapCPU, false/*abstract*/, Object)

#endif
