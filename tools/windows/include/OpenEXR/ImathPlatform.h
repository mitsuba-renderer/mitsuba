///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2002, Industrial Light & Magic, a division of Lucas
// Digital Ltd. LLC
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// *       Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// *       Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
// *       Neither the name of Industrial Light & Magic nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////



#ifndef INCLUDED_IMATHPLATFORM_H
#define INCLUDED_IMATHPLATFORM_H

//----------------------------------------------------------------------------
//
//	ImathPlatform.h
//
//	This file contains functions and constants which aren't 
//	provided by the system libraries, compilers, or includes on 
//	certain platforms.
//----------------------------------------------------------------------------

#if defined(_MSC_VER) 
/* [i_a] MSVC again: see the comment from Microsoft's math.h below (MSVC2005).

   Other compilers may also lack M_PI / M_PI_2 which both are used in the
   OpenEXR (test) code.

   Microsoft says:

	Define _USE_MATH_DEFINES before including math.h to expose these macro
	definitions for common math constants.  These are placed under an #ifdef
	since these commonly-defined names are not part of the C/C++ standards.

   End of quote.

   The other defines have been added for completeness sake (they exist on
   BSD and Linux at least and other code may expect these).

   Microsoft doesn't define M_2PI ever, other compilers may lack some of these
   too, hence the sequence as it is: load math.h, then see what's lacking still.
*/
#define _USE_MATH_DEFINES 1
#endif

#include <math.h>

// #if defined _WIN32 || defined _WIN64 -- [i_a] some other platforms (MingW?) don't have this either, it seems.
    #ifndef M_PI
        #define M_PI    3.14159265358979323846 /* pi */
    #endif

	#ifndef M_PI_2
        #define M_PI_2  1.57079632679489661923 /* pi/2 */
    #endif

/* Some extra useful constants.  */

#ifndef M_E
#define M_E            2.7182818284590452354   /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E        1.4426950408889634074   /* log_2 e */
#endif

#ifndef M_LOG10E
#define M_LOG10E       0.43429448190325182765  /* log_10 e */
#endif

#ifndef M_LN2
#define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif

#ifndef M_LN10
#define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif

#ifndef M_PI_4
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI         0.31830988618379067154  /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI         0.63661977236758134308  /* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif

#ifndef M_2PI
#define M_2PI           6.283185307179586232    /* 2*pi  */
#endif

// #endif

//-----------------------------------------------------------------------------
//
//    Some, but not all, C++ compilers support the C99 restrict
//    keyword or some variant of it, for example, __restrict.
//
//-----------------------------------------------------------------------------

#if defined __GNUC__

    //
    // supports __restrict
    //

    #define IMATH_RESTRICT __restrict

#elif defined (__INTEL_COMPILER) || \
      defined(__ICL) || \
      defined(__ICC) || \
      defined(__ECC)

    //
    // supports restrict
    //

    #define IMATH_RESTRICT restrict

#elif defined __sgi

    //
    // supports restrict
    //

    #define IMATH_RESTRICT restrict

#elif defined _MSC_VER

    //
    // supports __restrict
    //

    #define IMATH_RESTRICT __restrict

#else

    //
    // restrict / __restrict not supported
    //

    #define IMATH_RESTRICT

#endif

#endif // INCLUDED_IMATHPLATFORM_H
