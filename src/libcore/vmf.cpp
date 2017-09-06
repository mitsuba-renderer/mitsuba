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

#include <mitsuba/core/vmf.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/brent.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

Float VonMisesFisherDistr::eval(Float cosTheta) const {
    if (m_kappa == 0.0f)
        return INV_FOURPI;
#if 0
    return math::fastexp(cosTheta * m_kappa)
        * m_kappa / (4 * M_PI * std::sinh(m_kappa));
#else
    /* Numerically stable version */
    return math::fastexp(m_kappa * std::min((Float)0, cosTheta - 1))
        * m_kappa / (2 * M_PI * (1-math::fastexp(-2*m_kappa)));
#endif
}

Vector VonMisesFisherDistr::sample(const Point2 &sample) const {
    if (m_kappa == 0)
        return warp::squareToUniformSphere(sample);

#if 0
    Float cosTheta = math::fastlog(math::fastexp(-m_kappa) + 2 *
                        sample.x * std::sinh(m_kappa)) / m_kappa;
#else
    /* Numerically stable version */
    Float cosTheta = 1 + (math::fastlog(sample.x +
        math::fastexp(-2 * m_kappa) * (1 - sample.x))) / m_kappa;
#endif

    Float sinTheta = math::safe_sqrt(1-cosTheta*cosTheta),
          sinPhi, cosPhi;

    math::sincos(2*M_PI * sample.y, &sinPhi, &cosPhi);

    return Vector(cosPhi * sinTheta,
        sinPhi * sinTheta, cosTheta);
}

Float VonMisesFisherDistr::getMeanCosine() const {
    if (m_kappa == 0)
        return 0;
    Float coth = m_kappa > 6 ? 1 : ((std::exp(2*m_kappa)+1)/(std::exp(2*m_kappa)-1));
    return coth-1/m_kappa;
}

static Float A3(Float kappa) {
    return 1/ std::tanh(kappa) - 1 / kappa;
}

std::string VonMisesFisherDistr::toString() const {
    std::ostringstream oss;
    oss << "VonMisesFisherDistr[kappa=" << m_kappa << "]";
    return oss.str();
}

static Float dA3(Float kappa) {
    Float csch = 2.0f /
        (math::fastexp(kappa)-math::fastexp(-kappa));
    return 1/(kappa*kappa) - csch*csch;
}

static Float A3inv(Float y, Float guess) {
    Float x = guess;
    int it = 1;

    while (true) {
        Float residual = A3(x)-y,
              deriv = dA3(x);
        x -= residual/deriv;

        if (++it > 20) {
            SLog(EWarn, "VanMisesFisherDistr::convolve(): Newton's method "
                " did not converge!");
            return guess;
        }

        if (std::abs(residual) < 1e-5f)
            break;
    }
    return x;
}

Float VonMisesFisherDistr::convolve(Float kappa1, Float kappa2) {
    return A3inv(A3(kappa1) * A3(kappa2), std::min(kappa1, kappa2));
}

Float VonMisesFisherDistr::forPeakValue(Float x) {
    if (x < INV_FOURPI) {
        return 0.0f;
    } else if (x > 0.795) {
        return 2 * M_PI * x;
    } else {
        return std::max((Float) 0.0f,
            (168.479f * x * x + 16.4585f * x - 2.39942f) /
            (-1.12718f * x * x + 29.1433f * x + 1.0f));
    }
}

static Float meanCosineFunctor(Float kappa, Float g) {
    return VonMisesFisherDistr(kappa).getMeanCosine()-g;
}

Float VonMisesFisherDistr::forMeanLength(Float l) {
    return (3*l - l*l*l) / (1-l*l);
}

Float VonMisesFisherDistr::forMeanCosine(Float g) {
    if (g == 0)
        return 0;
    else if (g < 0)
        SLog(EError, "Error: vMF distribution cannot be created for g<0.");

    BrentSolver brentSolver(100, 1e-6f);
    BrentSolver::Result result = brentSolver.solve(
        boost::bind(&meanCosineFunctor, _1, g), 0, 1000);
    SAssert(result.success);
    return result.x;
}

MTS_NAMESPACE_END
