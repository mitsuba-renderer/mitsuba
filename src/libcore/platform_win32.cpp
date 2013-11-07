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

#if defined(_MSC_VER) && _MSC_VER < 1800

/* Use strict IEEE 754 floating point computations
   for the following code */
#pragma float_control(precise, on)

typedef union {
  float value;
  uint32_t word;
} ieee_float_shape_type;

typedef union {
  double value;
  struct {
    uint32_t lsw;
    uint32_t msw;
  } parts;
} ieee_double_shape_type;

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


#define EXTRACT_WORDS(ix0,ix1,d) \
do { \
  ieee_double_shape_type ew_u; \
  ew_u.value = (d); \
  (ix0) = ew_u.parts.msw; \
  (ix1) = ew_u.parts.lsw; \
} while (0)


#define INSERT_WORDS(d,ix0,ix1) \
do { \
  ieee_double_shape_type iw_u; \
  iw_u.parts.msw = (ix0); \
  iw_u.parts.lsw = (ix1); \
  (d) = iw_u.value; \
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

double nextafter(double x, double y) {
	int32_t hx,hy,ix,iy;
	uint32_t lx,ly;

	EXTRACT_WORDS(hx,lx,x);
	EXTRACT_WORDS(hy,ly,y);

	ix = hx & 0x7fffffff; /* |x| */
	iy = hy & 0x7fffffff; /* |y| */

	if (((ix>=0x7ff00000) && ((ix-0x7ff00000)|lx) != 0) ||   /* x is nan */
	    ((iy>=0x7ff00000) && ((iy-0x7ff00000)|ly) != 0))     /* y is nan */
	   return x+y;

	if (x==y)
		return x; /* x=y, return x */

	if ((ix|lx) == 0) { /* x == 0 */
	    INSERT_WORDS(x,hy & 0x80000000,1);	/* return +-minsubnormal */
	    y = x*x;
	    if (y==x)
			return y;
		else
			return x; /* raise underflow flag */
	}

	if (hx>=0) { /* x > 0 */
	    if (hx>hy || ((hx==hy) && (lx>ly))) { /* x > y, x -= ulp */
			if (lx==0)
				hx -= 1;
			lx -= 1;
	    } else {				/* x < y, x += ulp */
			lx += 1;
			if (lx==0)
				hx += 1;
	    }
	} else { /* x < 0 */
	    if (hy>=0 || hx>hy || ((hx==hy) && (lx>ly))) { /* x < y, x -= ulp */
			if (lx==0)
				hx -= 1;
			lx -= 1;
	    } else { /* x > y, x += ulp */
			lx += 1;
			if (lx==0)
				hx += 1;
		}
	}

	hy = hx & 0x7ff00000;

	if (hy >= 0x7ff00000)
		return x+x;	/* overflow  */

	if (hy < 0x00100000) { /* underflow */
	    y = x*x;
	    if (y!=x) {
			/* raise underflow flag */
	        INSERT_WORDS(y,hx,lx);
			return y;
	    }
	}

	INSERT_WORDS(x,hx,lx);

	return x;
}
#endif
