/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    This particular file is based on code by Piti Irawan, which is
    redistributed with permission.

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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/noise.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/qmc.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/bitmap.h>
#include "irawan.h"
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{irawan}{Irawan \& Marschner woven cloth BRDF}
 *
 * \parameters{
 *     \parameter{filename}{\String}{Path to a weave pattern description}
 *     \parameter{repeatU, repeatV}{\Float}{Specifies the number
 *         of weave pattern repetitions over a $[0,1]^2$ region of the UV
 *         parameterization}
 *     \parameter{(\emph{Additional parameters})}{\Spectrum\Or\Float}{
 *         Weave pattern files may define their own custom parameters; this is
 *         useful for instance to support changing the color of a weave
 *         without having to create a new file every time.
 *         These parameters must be specified directly to the plugin
 *         so that they can be appropriately resolved when the pattern file is loaded.
 *     }
 * }
 *
 * This plugin implements the Irawan \& Marschner BRDF,
 * a realistic model for rendering woven materials.
 * This spatially-varying reflectance model uses an explicit description
 * of the underlying weave pattern to create fine-scale texture and
 * realistic reflections across a wide range of different weave types.
 * To use the model, you must provide a special weave pattern
 * file---for an example of what these look like, see the
 * examples scenes available on the Mitsuba website.
 *
 * A detailed explanation of the model is beyond the scope of this manual.
 * For reference, it is described in detail in the PhD thesis of
 * Piti Irawan (``The Appearance of Woven Cloth'' \cite{IrawanThesis}).
 * The code in Mitsuba a modified port of a previous Java implementation
 * by Piti, which has been extended with a simple domain-specific weave
 * pattern description language. \vspace{8mm}
 *
 * \renderings{
 *     \unframedmedrendering{Silk charmeuse}{bsdf_irawan_charmeuse}
 *     \unframedmedrendering{Cotton denim}{bsdf_irawan_denim}
 *     \unframedmedrendering{Wool gabardine}{bsdf_irawan_gabardine}
 * }
 * \renderings{
 *     \setcounter{subfigure}{3}
 *     \unframedmedrendering{Polyester lining cloth}{bsdf_irawan_polyester}
 *     \unframedmedrendering{Silk shantung}{bsdf_irawan_shantung}
 *     \unframedmedrendering{Cotton twill}{bsdf_irawan_twill}
 * }
 */
class IrawanClothBRDF : public BSDF {
public:
    IrawanClothBRDF(const Properties &props)
        : BSDF(props), m_specularNormalization(0) {

        FileResolver *fResolver = Thread::getThread()->getFileResolver();
        fs::path path = fResolver->resolve(props.getString("filename"));
        if (!fs::exists(path))
            Log(EError, "Weave pattern file \"%s\" could not be found!",
                path.string().c_str());
        fs::ifstream in(path);
        typedef spirit::istream_iterator iterator_type;
        iterator_type end, begin(in);

        WeavePatternGrammar<iterator_type> g(props);
        SkipGrammar<iterator_type> sg;

        bool success = phrase_parse(begin, end, g, sg, m_pattern);
        if (!success)
            Log(EError, "Unable to parse the weave pattern file \"%s\"!",
                path.string().c_str());

        /* Some sanity checks */
        SAssert(m_pattern.pattern.size() ==
                m_pattern.tileWidth * m_pattern.tileHeight);
        for (size_t i=0; i<m_pattern.pattern.size(); ++i)
            SAssert(m_pattern.pattern[i] > 0 &&
                    m_pattern.pattern[i] <= m_pattern.yarns.size());

        /* U and V tile count */
        m_repeatU = props.getFloat("repeatU");
        m_repeatV = props.getFloat("repeatV");

        if (props.hasProperty("ksMultiplier") || props.hasProperty("kdMultiplier"))
            Log(EError, "The 'ksMultiplier' and 'kdMultiplier' parameters were "
                "replaced by a normalization scheme. Please remove them and "
                "appropriately set the 'kd' and 'ks'-values used in your model.");
    }

    IrawanClothBRDF(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_pattern = WeavePattern(stream);
        m_repeatU = stream->readFloat();
        m_repeatV = stream->readFloat();
        m_specularNormalization = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        m_pattern.serialize(stream);
        stream->writeFloat(m_repeatU);
        stream->writeFloat(m_repeatV);
        stream->writeFloat(m_specularNormalization);
    }

    void configure() {
        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide
            | EAnisotropic | ESpatiallyVarying);
        m_components.push_back(EDiffuseReflection | EFrontSide
            | ESpatiallyVarying);

        /* Estimate the average reflectance under diffuse
           illumination and use it to normalize the specular
           component */
        ref<Random> random = new Random();
        size_t nSamples = 10000;

        if (m_specularNormalization == 0) {
            Intersection its;
            BSDFSamplingRecord bRec(its, NULL, ERadiance);
            Spectrum result(0.0f);
            m_initialization = true;
            for (size_t i=0; i<nSamples; ++i) {
                bRec.wi = warp::squareToCosineHemisphere(Point2(random->nextFloat(), random->nextFloat()));
                bRec.wo = warp::squareToCosineHemisphere(Point2(random->nextFloat(), random->nextFloat()));
                its.uv = Point2(random->nextFloat(), random->nextFloat());

                result += eval(bRec, ESolidAngle) / Frame::cosTheta(bRec.wo);
            }
            m_initialization = false;

            if (result.max() == 0)
                m_specularNormalization = 0;
            else
                m_specularNormalization = nSamples / (result.max() * M_PI);
        }

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        Point2 uv = Point2(its.uv.x * m_repeatU,
            (1 - its.uv.y) * m_repeatV);
        Point2 xy(uv.x * m_pattern.tileWidth, uv.y * m_pattern.tileHeight);
        Point2i lookup(
            math::modulo((int) xy.x, m_pattern.tileWidth),
            math::modulo((int) xy.y, m_pattern.tileHeight));
        int yarnID = m_pattern.pattern[lookup.x + lookup.y * m_pattern.tileWidth] - 1;
        const Yarn &yarn = m_pattern.yarns.at(yarnID);

        return yarn.kd;
    }


    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            (!hasDiffuse && !hasSpecular) ||
            measure != ESolidAngle)
            return Spectrum(0.0f);

        Point2 uv = Point2(bRec.its.uv.x * m_repeatU,
            (1 - bRec.its.uv.y) * m_repeatV);
        Point2 xy(uv.x * m_pattern.tileWidth, uv.y * m_pattern.tileHeight);

        Point2i lookup(
            math::modulo((int) xy.x, m_pattern.tileWidth),
            math::modulo((int) xy.y, m_pattern.tileHeight));

        int yarnID = m_pattern.pattern[lookup.x + lookup.y * m_pattern.tileWidth] - 1;

        const Yarn &yarn = m_pattern.yarns.at(yarnID);
        // store center of the yarn segment
        Point2 center
            (((int) xy.x / m_pattern.tileWidth) * m_pattern.tileWidth
                + yarn.centerU * m_pattern.tileWidth,
             ((int) xy.y / m_pattern.tileHeight) * m_pattern.tileHeight
                + (1 - yarn.centerV) * m_pattern.tileHeight);

        // transform x and y to new coordinate system with (0,0) at the
        // center of the yarn segment
        xy.x =    xy.x - center.x;
        xy.y = - (xy.y - center.y);

        int type = yarn.type;
        Float w = yarn.width;
        Float l = yarn.length;

        // Get incident and exitant directions.
        Vector om_i = bRec.wi;
        Vector om_r = bRec.wo;

        Float psi = yarn.psi;
        Float umax = yarn.umax;
        Float kappa = yarn.kappa;

        Float dUmaxOverDWarp, dUmaxOverDWeft;
        if (type == Yarn::EWarp) {
            dUmaxOverDWarp = m_pattern.dWarpUmaxOverDWarp;
            dUmaxOverDWeft = m_pattern.dWarpUmaxOverDWeft;
        } else { // type == EWeft
            dUmaxOverDWarp = m_pattern.dWeftUmaxOverDWarp;
            dUmaxOverDWeft = m_pattern.dWeftUmaxOverDWeft;
            // Rotate xy, incident, and exitant directions pi/2 radian about z-axis
            Float tmp = xy.x;
            xy.x = -xy.y;
            xy.y = tmp;
            tmp = om_i.x;
            om_i.x = -om_i.y;
            om_i.y = tmp;
            tmp = om_r.x;
            om_r.x = -om_r.y;
            om_r.y = tmp;
        }

        // Correlated (Perlin) noise.
        Float random1 = 1.0f;
        Float random2 = 1.0f;

        /* Number of TEA iterations (the more, the better the
           quality of the pseudorandom floats) */
        const int teaIterations = 8;

        if (m_pattern.period > 0.0f) {
            // generate 1 seed per yarn segment
            Point2u pos(center);

            random1 = Noise::perlinNoise(Point(
                (center.x * (m_pattern.tileHeight * m_repeatV
                    + sampleTEAFloat(pos.x, 2*pos.y, teaIterations)) + center.y) / m_pattern.period, 0, 0));
            random2 = Noise::perlinNoise(Point(
                (center.y * (m_pattern.tileWidth * m_repeatU
                    + sampleTEAFloat(pos.x, 2*pos.y+1, teaIterations)) + center.x) / m_pattern.period, 0, 0));
            umax = umax + random1 * dUmaxOverDWarp + random2 * dUmaxOverDWeft;
        }

        // Compute u and v.
        // See Chapter 6.
        Float u = xy.y / (l / 2.0f) * umax;
        Float v = xy.x * M_PI / w;

        // Compute specular contribution.
        Spectrum result(0.0f);
        if (hasSpecular) {
            Float integrand;
            if (psi != 0.0f)
                integrand = evalStapleIntegrand(u, v, om_i, om_r, m_pattern.alpha,
                        m_pattern.beta, psi, umax, kappa, w, l);
            else
                integrand = evalFilamentIntegrand(u, v, om_i, om_r, m_pattern.alpha,
                        m_pattern.beta, m_pattern.ss, umax, kappa, w, l);

            // Initialize random number generator based on texture location.
            Float intensityVariation = 1.0f;
            if (m_pattern.fineness > 0.0f) {
                // Compute random variation and scale specular component.
                // Generate fineness^2 seeds per 1 unit of texture.
                uint32_t index1 = (uint32_t) ((center.x + xy.x) * m_pattern.fineness);
                uint32_t index2 = (uint32_t) ((center.y + xy.y) * m_pattern.fineness);

                Float xi = sampleTEAFloat(index1, index2, teaIterations);
                intensityVariation = std::min(-math::fastlog(xi), (Float) 10.0f);
            }

            if (!m_initialization)
                result = yarn.ks * (intensityVariation * integrand * m_specularNormalization);
            else
                result = Spectrum(intensityVariation * integrand);

            if (type == Yarn::EWarp)
                result *= (m_pattern.warpArea + m_pattern.weftArea) / m_pattern.warpArea;
            else
                result *= (m_pattern.warpArea + m_pattern.weftArea) / m_pattern.weftArea;
        }

        if (hasDiffuse && !m_initialization)
            result += yarn.kd * INV_PI;

        return result * Frame::cosTheta(bRec.wo);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            (!hasDiffuse && !hasSpecular) ||
            measure != ESolidAngle)
            return 0.0f;

        return warp::squareToCosineHemispherePdf(bRec.wo);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (Frame::cosTheta(bRec.wi) <= 0 ||
            (!hasDiffuse && !hasSpecular))
            return Spectrum(0.0f);

        /* Lacking a better sampling method, generate cosine-weighted directions */
        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;
        return eval(bRec, ESolidAngle) * M_PI / Frame::cosTheta(bRec.wo);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (Frame::cosTheta(bRec.wi) <= 0 ||
            (!hasDiffuse && !hasSpecular))
            return Spectrum(0.0f);

        /* Lacking a better sampling method, generate cosine-weighted directions */
        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;
        pdf = warp::squareToCosineHemispherePdf(bRec.wo);
        return eval(bRec, ESolidAngle) / pdf;
    }

    /** parameters:
     *  u    to be compared to u(v) in texturing
     *  v    for filament, we compute u(v)
     *  om_i  incident direction
     *  om_r  exitant direction
     *  ss  filament smoothing parameter
     *  fiber properties
     *  alpha uniform scattering
     *  beta  forward scattering
     *  yarn geometry
     *  psi   fiber twist angle; because this is filament, psi = pi/2
     *  umax  maximum inclination angle
     *  kappa spine curvature parameter
     *  weave pattern
     *  w    width of segment rectangle
     *  l    length of segment rectangle
     */
    Float evalFilamentIntegrand(Float u, Float v, const Vector &om_i,
            const Vector &om_r, Float alpha, Float beta, Float ss,
            Float umax, Float kappa, Float w, Float l) const {
        // 0 <= ss < 1.0
        if (ss < 0.0f || ss >= 1.0f)
            return 0.0f;

        // w * sin(umax) < l
        if (w * std::sin(umax) >= l)
            return 0.0f;

        // -1 < kappa < inf
        if (kappa < -1.0f)
            return 0.0f;

        // h is the half vector
        Vector h = normalize(om_r + om_i);

        // u_of_v is location of specular reflection.
        Float u_of_v = std::atan(h.y / h.z);

        // Check if u_of_v within the range of valid u values
        if (std::abs(u_of_v) < umax) {
            // n is normal to the yarn surface
            // t is tangent of the fibers.
            Normal n = normalize(Normal(std::sin(v), std::sin(u_of_v) * std::cos(v),
                    std::cos(u_of_v) * std::cos(v)));
            Vector t = normalize(Vector(0.0f, std::cos(u_of_v), -std::sin(u_of_v)));

            // R is radius of curvature.
            Float R = radiusOfCurvature(std::min(std::abs(u_of_v),
                (1-ss)*umax), (1-ss)*umax, kappa, w, l);

            // G is geometry factor.
            Float a = 0.5f * w;
            Vector om_i_plus_om_r = om_i + om_r,
                   t_cross_h = cross(t, h);
            Float Gu = a * (R + a * std::cos(v)) /
                (om_i_plus_om_r.length() * std::abs(t_cross_h.x));

            // fc is phase function
            Float fc = alpha + vonMises(-dot(om_i, om_r), beta);

            // A is attenuation function without smoothing.
            // As is attenuation function with smoothing.
            Float A = seeliger(dot(n, om_i), dot(n, om_r), 0, 1);
            Float As;
            if (ss == 0.0f)
                As = A;
            else
                As = A * (1.0f - math::smoothStep((Float) 0, (Float) 1, (std::abs(u_of_v)
                    - (1.0f - ss) * umax) / (ss * umax)));

            // fs is scattering function.
            Float fs = Gu * fc * As;

            // Domain transform.
            fs = fs * M_PI * l;

            // Highlight has constant width delta_y on screen.
            Float delta_y = l * m_pattern.hWidth;

            // Clamp y_of_v between -(l - delta_y)/2 and (l - delta_y)/2.
            Float y_of_v = u_of_v * 0.5f * l / umax;
            if (y_of_v > 0.5f * (l - delta_y))
                y_of_v = 0.5f * (l - delta_y);
            else if (y_of_v < 0.5f * (delta_y - l))
                y_of_v = 0.5f * (delta_y - l);

            // Check if |y(u(v)) - y(u)| < delta_y/2.
            if (std::abs(y_of_v - u * 0.5f * l / umax) < 0.5f * delta_y)
                return fs / delta_y;
        }
        return 0.0f;
    }

    /** parameters:
     *  u    for staple, we compute v(u)
     *  v    to be compared to v(u) in texturing
     *  om_i  incident direction
     *  om_r  exitant direction
     *  fiber properties
     *  alpha uniform scattering
     *  beta  forward scattering
     *  yarn geometry
     *  psi   fiber twist angle
     *  umax  maximum inclination angle
     *  kappa spine curvature parameter
     *  weave pattern
     *  w    width of segment rectangle
     *  l    length of segment rectangle
     */
    Float evalStapleIntegrand(Float u, Float v, const Vector &om_i,
            const Vector &om_r, Float alpha, Float beta, Float psi,
            Float umax, Float kappa, Float w, Float l) const {
        // w * sin(umax) < l
        if (w * std::sin(umax) >= l)
            return 0.0f;

        // -1 < kappa < inf
        if (kappa < -1.0f)
            return 0.0f;

        // h is the half vector
        Vector h = normalize(om_i + om_r);

        // v_of_u is location of specular reflection.
        Float D = (h.y*std::cos(u) - h.z*std::sin(u))
            / (std::sqrt(h.x * h.x + std::pow(h.y * std::sin(u) + h.z * std::cos(u), (Float) 2.0f)) * std::tan(psi));
        Float v_of_u = std::atan2(-h.y * std::sin(u) - h.z * std::cos(u), h.x) + math::safe_acos(D);

        // Check if v_of_u within the range of valid v values
        if (std::abs(D) < 1.0f && std::abs(v_of_u) < M_PI / 2.0f) {
            // n is normal to the yarn surface.
            // t is tangent of the fibers.

            Normal n = normalize(Normal(std::sin(v_of_u), std::sin(u)
                    * std::cos(v_of_u), std::cos(u) * std::cos(v_of_u)));

            /*Vector t = normalize(Vector(-std::cos(v_of_u) * std::sin(psi),
                    std::cos(u) * std::cos(psi) + std::sin(u) * std::sin(v_of_u) * std::sin(psi),
                    -std::sin(u) * std::cos(psi) + std::cos(u) * std::sin(v_of_u) * std::sin(psi))); */

            // R is radius of curvature.
            Float R = radiusOfCurvature(std::abs(u), umax, kappa, w, l);

            // G is geometry factor.
            Float a = 0.5f * w;
            Vector om_i_plus_om_r(om_i + om_r);
            Float Gv = a * (R + a * std::cos(v_of_u))
                / (om_i_plus_om_r.length() * dot(n, h) * std::abs(std::sin(psi)));

            // fc is phase function.
            Float fc = alpha + vonMises(-dot(om_i, om_r), beta);

            // A is attenuation function without smoothing.
            Float A = seeliger(dot(n, om_i), dot(n, om_r), 0, 1);

            // fs is scattering function.
            Float fs = Gv * fc * A;

            // Domain transform.
            fs = fs * 2.0f * w * umax;

            // Highlight has constant width delta_x on screen.
            Float delta_x = w * m_pattern.hWidth;

            // Clamp x_of_u between (w - delta_x)/2 and -(w - delta_x)/2.
            Float x_of_u = v_of_u * w / M_PI;
            if (x_of_u > 0.5f * (w - delta_x))
                x_of_u = 0.5f * (w - delta_x);
            else if (x_of_u < 0.5f * (delta_x - w))
                x_of_u = 0.5f * (delta_x - w);

            // Check if |x(v(u)) - x(v)| < delta_x/2.
            if (std::abs(x_of_u - v * w / M_PI) < 0.5f * delta_x)
                return fs / delta_x;
        }
        return 0.0f;
    }

    Float radiusOfCurvature(Float u, Float umax, Float kappa, Float w, Float l) const {
        // rhat determines whether the spine is a segment
        // of an ellipse, a parabole, or a hyperbola.
        // See Section 5.3.
        Float rhat = 1.0f + kappa * (1.0f + 1.0f / std::tan(umax));

        Float a = 0.5f * w, R = 0;
        if (rhat == 1.0f) { // circle; see Subsection 5.3.1.
            R = (0.5f * l - a * std::sin(umax)) / std::sin(umax);
        } else if (rhat > 0.0f) {
            Float tmax = std::atan(rhat * std::tan(umax));
            Float bhat = (0.5f * l - a * std::sin(umax)) / std::sin(tmax);
            Float ahat = bhat / rhat;
            Float t = std::atan(rhat * std::tan(u));
            R = std::pow(bhat * bhat * std::cos(t) * std::cos(t)
              + ahat * ahat * std::sin(t) * std::sin(t),(Float) 1.5f) / (ahat * bhat);
        } else if (rhat < 0.0f) { // hyperbola; see Subsection 5.3.3.
            Float tmax = -atanh(rhat * std::tan(umax));
            Float bhat = (0.5f * l - a * std::sin(umax)) / std::sinh(tmax);
            Float ahat = bhat / rhat;
            Float t = -atanh(rhat * std::tan(u));
            R = -std::pow(bhat * bhat * std::cosh(t) * std::cosh(t)
                + ahat * ahat * std::sinh(t) * std::sinh(t), (Float) 1.5f) / (ahat * bhat);
        } else { // rhat == 0  // parabola; see Subsection 5.3.2.
            Float tmax = std::tan(umax);
            Float ahat = (0.5f * l - a * std::sin(umax)) / (2 * tmax);
            Float t = std::tan(u);
            R = 2 * ahat * std::pow(1 + t * t, (Float) 1.5f);
        }
        return R;
    }

    inline Float atanh(Float arg) const {
        return math::fastlog((1.0f + arg) / (1.0f - arg)) / 2.0f;
    }

    // von Mises Distribution
    Float vonMises(Float cos_x, Float b) const {
        // assumes a = 0, b > 0 is a concentration parameter.

        Float I0, absB = std::abs(b);
        if (std::abs(b) <= 3.75f) {
            Float t = absB / 3.75f;
            t = t * t;
            I0 = 1.0f + t*(3.5156229f + t*(3.0899424f + t*(1.2067492f
                    + t*(0.2659732f + t*(0.0360768f + t*0.0045813f)))));
        } else {
            Float t = 3.75f / absB;
            I0 = math::fastexp(absB) / std::sqrt(absB) * (0.39894228f + t*(0.01328592f
                + t*(0.00225319f + t*(-0.00157565f + t*(0.00916281f + t*(-0.02057706f
                + t*(0.02635537f + t*(-0.01647633f + t*0.00392377f))))))));
        }

        return math::fastexp(b * cos_x) / (2 * M_PI * I0);
    }

    /// Attenuation term
    Float seeliger(Float cos_th1, Float cos_th2, Float sg_a, Float sg_s) const {
        Float al = sg_s / (sg_a + sg_s); // albedo
        Float c1 = std::max((Float) 0, cos_th1);
        Float c2 = std::max((Float) 0, cos_th2);
        if (c1 == 0.0f || c2 == 0.0f)
            return 0.0f;
        return al / (4.0f * M_PI) * c1 * c2 / (c1 + c2);
    }

    Float getRoughness(const Intersection &its, int component) const {
        /* For lack of a better value, treat this material as diffuse
           in Manifold Exploration */
        return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "IrawanClothBRDF[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  weavePattern = " << indent(m_pattern.toString()) << "," << endl
            << "  repeatU = " << m_repeatU << "," << endl
            << "  repeatV = " << m_repeatV << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    WeavePattern m_pattern;
    Float m_repeatU, m_repeatV;
    Float m_specularNormalization;
    bool m_initialization;
};

// ================ Hardware shader implementation ================

/**
 * In place of a real shader, let's just show the average
 * diffuse albedo for now..
 */
class IrawanShader : public Shader {
public:
    IrawanShader(Renderer *renderer, Spectrum albedo)
        : Shader(renderer, EBSDFShader), m_albedo(albedo) {
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_albedo", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_albedo);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_albedo;" << endl
            << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << evalName << "_albedo * inv_pi * cosTheta(wo);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_albedo;
};

Shader *IrawanClothBRDF::createShader(Renderer *renderer) const {
    Spectrum albedo(0.0f);
    for (size_t i=0; i<m_pattern.yarns.size(); ++i)
        albedo += m_pattern.yarns[i].kd;
    albedo /= (Float) m_pattern.yarns.size();
    return new IrawanShader(renderer, albedo);
}


MTS_IMPLEMENT_CLASS(IrawanShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(IrawanClothBRDF, false, BSDF)
MTS_EXPORT_PLUGIN(IrawanClothBRDF, "Irawan & Marschner woven cloth BRDF")
MTS_NAMESPACE_END
