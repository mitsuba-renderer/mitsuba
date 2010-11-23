/* s_nextafterf.c -- float version of s_nextafter.c.
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */

/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

#include <mitsuba/mitsuba.h>

/* Use strict IEEE 754 floating point computations 
   for the following code */
#pragma float_control(precise, on)

typedef union {
  float value;
  uint32_t word;
} ieee_float_shape_type;

#define GET_FLOAT_WORD(i,d) \
do { \
  ieee_float_shape_type gf_u; \
  gf_u.value = (d); \
  (i) = gf_u.word; \
} while (0)

#define SET_FLOAT_WORD(d,i) \
do { \
  ieee_float_shape_type sf_u; \
  sf_u.word = (i); \
  (d) = sf_u.value; \
} while (0)

float nextafterf(float x, float y) {
	int32_t hx, hy, ix, iy;

	GET_FLOAT_WORD(hx, x);
	GET_FLOAT_WORD(hy, y);
	ix = hx & 0x7fffffff;		/* |x| */
	iy = hy & 0x7fffffff;		/* |y| */

	/* x is nan or y is nan? */
	if ((ix > 0x7f800000) || (iy > 0x7f800000))
		return x + y;

	if (x == y)
		return y;

	if (ix == 0) { /* x == 0? */
		SET_FLOAT_WORD(x, (hy & 0x80000000) | 1);
		return x;
	}

	if (hx >= 0) { /* x > 0 */
		if (hx > hy) { /* x > y: x -= ulp */
			hx -= 1;
		} else { /* x < y: x += ulp */
			hx += 1;
		}
	} else { /* x < 0 */
		if (hy >= 0 || hx > hy) { /* x < y: x -= ulp */
			hx -= 1;
		} else { /* x > y: x += ulp */
			hx += 1;
		}
	}
	hy = hx & 0x7f800000;
	if (hy >= 0x7f800000) {
		x = x + x; /* overflow */
		return x; /* overflow */
	}
	SET_FLOAT_WORD(x, hx);
	return x;
}
