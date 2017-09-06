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

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

namespace math {

Float erfinv(Float x) {
    // Based on "Approximating the erfinv function" by Mark Giles
    Float w = -math::fastlog(((Float) 1 - x)*((Float) 1 + x));
    Float p;
    if (w < (Float) 5) {
        w = w - (Float) 2.5;
        p = (Float) 2.81022636e-08;
        p = (Float) 3.43273939e-07 + p*w;
        p = (Float) -3.5233877e-06 + p*w;
        p = (Float) -4.39150654e-06 + p*w;
        p = (Float) 0.00021858087 + p*w;
        p = (Float) -0.00125372503 + p*w;
        p = (Float) -0.00417768164 + p*w;
        p = (Float) 0.246640727 + p*w;
        p = (Float) 1.50140941 + p*w;
    } else {
        w = std::sqrt(w) - (Float) 3;
        p = (Float) -0.000200214257;
        p = (Float) 0.000100950558 + p*w;
        p = (Float) 0.00134934322 + p*w;
        p = (Float) -0.00367342844 + p*w;
        p = (Float) 0.00573950773 + p*w;
        p = (Float) -0.0076224613 + p*w;
        p = (Float) 0.00943887047 + p*w;
        p = (Float) 1.00167406 + p*w;
        p = (Float) 2.83297682 + p*w;
    }
    return p*x;
}

Float erf(Float x) {
    Float a1 = (Float)  0.254829592;
    Float a2 = (Float) -0.284496736;
    Float a3 = (Float)  1.421413741;
    Float a4 = (Float) -1.453152027;
    Float a5 = (Float)  1.061405429;
    Float p  = (Float)  0.3275911;

    // Save the sign of x
    Float sign = math::signum(x);
    x = std::abs(x);

    // A&S formula 7.1.26
    Float t = (Float) 1.0 / ((Float) 1.0 + p*x);
    Float y = (Float) 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math::fastexp(-x*x);

    return sign*y;
}

float hypot2(float a, float b) {
    float r;
    if (std::abs(a) > std::abs(b)) {
        r = b / a;
        r = std::abs(a) * std::sqrt(1.0f + r*r);
    } else if (b != 0.0f) {
        r = a / b;
        r = std::abs(b) * std::sqrt(1.0f + r*r);
    } else {
        r = 0.0f;
    }
    return r;
}


double hypot2(double a, double b) {
    double r;
    if (std::abs(a) > std::abs(b)) {
        r = b / a;
        r = std::abs(a) * std::sqrt(1.0 + r*r);
    } else if (b != 0.0) {
        r = a / b;
        r = std::abs(b) * std::sqrt(1.0 + r*r);
    } else {
        r = 0.0;
    }
    return r;
}

float log2(float value) {
    const float invLn2 = 1.0f / std::log(2.0f);
    return fastlog(value) * invLn2;
}

double log2(double value) {
    const double invLn2 = 1.0 / std::log(2.0);
    return fastlog(value) * invLn2;
}

int log2i(uint32_t value) {
    int r = 0;
    while ((value >> r) != 0)
        r++;
    return r-1;
}

int log2i(uint64_t value) {
    int r = 0;
    while ((value >> r) != 0)
        r++;
    return r-1;
}

/* Fast rounding & power-of-two test algorithms from PBRT */
uint32_t roundToPowerOfTwo(uint32_t i) {
    i--;
    i |= i >> 1; i |= i >> 2;
    i |= i >> 4; i |= i >> 8;
    i |= i >> 16;
    return i+1;
}

uint64_t roundToPowerOfTwo(uint64_t i) {
    i--;
    i |= i >> 1;  i |= i >> 2;
    i |= i >> 4;  i |= i >> 8;
    i |= i >> 16; i |= i >> 32;
    return i+1;
}

};

MTS_NAMESPACE_END
