#if !defined(__CONSTANTS_H)
#define __CONSTANTS_H

/* Choice of precision */
#ifdef DOUBLE_PRECISION
#define Float double
#define Epsilon 1e-6
#else
#ifndef SINGLE_PRECISION
#define SINGLE_PRECISION
#endif
#define Float float
#define Epsilon 1e-4f
#endif

/// relative eps. for shadow rays
#define ShadowEpsilon 1e-3 

#ifdef M_PI
#undef M_PI
#endif

#ifdef INFINITY
#undef INFINITY
#endif

#ifdef SINGLE_PRECISION
#define M_PI         3.14159265358979323846f
#define INV_PI       0.31830988618379067154f
#define INV_TWOPI    0.15915494309189533577f
#define SQRT_TWO     1.41421356237309504880f
#define INV_SQRT_TWO 0.70710678118654752440f
#else
#define M_PI         3.14159265358979323846
#define INV_PI       0.31830988618379067154
#define INV_TWOPI    0.15915494309189533577
#define SQRT_TWO     1.41421356237309504880
#define INV_SQRT_TWO 0.70710678118654752440
#endif

#endif /* __CONSTANTS_H */
