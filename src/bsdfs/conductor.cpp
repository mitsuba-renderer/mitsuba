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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/hw/basicshader.h>
#include <boost/algorithm/string.hpp>
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{conductor}{Smooth conductor}
 * \order{6}
 * \icon{bsdf_conductor}
 * \parameters{
 *     \parameter{material}{\String}{Name of a material preset, see
 *           \tblref{conductor-iors}.\!\default{\texttt{Cu} / copper}}
 *     \parameter{eta, k}{\Spectrum}{Real and imaginary components of the material's index of
 *             refraction \default{based on the value of \texttt{material}}}
 *     \parameter{extEta}{\Float\Or\String}{
 *           Real-valued index of refraction of the surrounding dielectric,
 *           or a material name of a dielectric \default{\code{air}}
 *     }
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 * }
 * \renderings{
 *     \rendering{Measured copper material (the default), rendered using 30
 *     spectral samples between 360 and 830$nm$}
 *         {bsdf_conductor_copper.jpg}
 *     \rendering{Measured gold material (\lstref{conductor-gold})}
 *         {bsdf_conductor_gold.jpg}
 * }
 *
 * This plugin implements a perfectly smooth interface to a conducting material,
 * such as a metal. For a similar model that instead describes a rough surface
 * microstructure, take a look at the separately available
 * \pluginref{roughconductor} plugin.

 * In contrast to dielectric materials, conductors do not transmit
 * any light. Their index of refraction is complex-valued and tends to undergo
 * considerable changes throughout the visible color spectrum.
 *
 * To facilitate the tedious task of specifying spectrally-varying index of
 * refraction information, Mitsuba ships with a set of measured data for
 * several materials, where visible-spectrum information was publicly
 * available\footnote{
 *   These index of refraction values are identical to the data distributed
 *   with PBRT. They are originally from the Luxpop database
 *   (\url{www.luxpop.com}) and are based on data by Palik et al.
 *   \cite{Palik1998Handbook} and measurements of atomic scattering factors
 *   made by the Center For X-Ray Optics (CXRO) at Berkeley and the
 *   Lawrence Livermore National Laboratory (LLNL).
 * }.
 *
 * Note that \tblref{conductor-iors} also includes several popular optical
 * coatings, which are not actually conductors. These materials can also
 * be used with this plugin, though note that the plugin will ignore any
 * refraction component that the actual material might have had.
 * There is also a special material profile named \code{none}, which disables
 * the computation of Fresnel reflectances and produces an idealized
 * 100% reflecting mirror.
 *
 * When using this plugin, you should ideally compile Mitsuba with support for
 * spectral rendering to get the most accurate results. While it also works
 * in RGB mode, the computations will be more approximate in nature.
 * Also note that this material is one-sided---that is, observed from the
 * back side, it will be completely black. If this is undesirable,
 * consider using the \pluginref{twosided} BRDF adapter plugin.\vspace{4mm}
 *
 * \begin{xml}[caption=A material configuration for a smooth conductor with
 *    measured gold data, label=lst:conductor-gold]
 * <shape type="...">
 *     <bsdf type="conductor">
 *         <string name="material" value="Au"/>
 *     </bsdf>
 * <shape>
 * \end{xml}
 * \vspace{5mm}
 * It is also possible to load spectrally varying index of refraction data from
 * two external files containing the real and imaginary components,
 * respectively (see \secref{format-spectra} for details on the file
 * format):
 * \begin{xml}[caption=Rendering a smooth conductor with custom data]
 * <shape type="...">
 *     <bsdf type="conductor">
 *         <spectrum name="eta" filename="conductorIOR.eta.spd"/>
 *         <spectrum name="k" filename="conductorIOR.k.spd"/>
 *     </bsdf>
 * <shape>
 * \end{xml}
 * \vspace{1.5cm}
 * \begin{table}[hb!]
 * \centering
 * \scriptsize
 * \begin{tabular}{>{\ttfamily}llp{1mm}>{\ttfamily}ll}
 * \toprule
 * \rmfamily\textbf{Preset(s)} & \textbf{Description} &&
 * \rmfamily\textbf{Preset(s)} & \textbf{Description}\\
 * \cmidrule{1-2} \cmidrule{4-5}
 * a-C                   & Amorphous carbon &&            Na\_palik             & Sodium \\
 * Ag                    & Silver &&                      Nb, Nb\_palik         & Niobium \\
 * Al                    & Aluminium &&                   Ni\_palik             & Nickel \\
 * AlAs, AlAs\_palik     & Cubic aluminium arsenide &&    Rh, Rh\_palik         & Rhodium \\
 * AlSb, AlSb\_palik     & Cubic aluminium antimonide &&  Se, Se\_palik         & Selenium \\
 * Au                    & Gold &&                        SiC, SiC\_palik       & Hexagonal silicon carbide \\
 * Be, Be\_palik         & Polycrystalline beryllium &&   SnTe, SnTe\_palik     & Tin telluride\\
 * Cr                    & Chromium &&                    Ta, Ta\_palik         & Tantalum \\
 * CsI, CsI\_palik       & Cubic caesium iodide &&        Te, Te\_palik         & Trigonal tellurium \\
 * Cu, Cu\_palik         & Copper &&                      ThF4, ThF4\_palik     & Polycryst. thorium (IV) fluoride \\
 * Cu2O, Cu2O\_palik     & Copper (I) oxide &&            TiC, TiC\_palik       & Polycrystalline titanium carbide \\
 * CuO, CuO\_palik       & Copper (II) oxide &&           TiN, TiN\_palik       & Titanium nitride \\
 * d-C, d-C\_palik       & Cubic diamond &&               TiO2, TiO2\_palik     & Tetragonal titan. dioxide\\
 * Hg, Hg\_palik         & Mercury &&                     VC, VC\_palik         & Vanadium carbide \\
 * HgTe, HgTe\_palik     & Mercury telluride &&           V\_palik              & Vanadium \\
 * Ir, Ir\_palik         & Iridium &&                     VN, VN\_palik         & Vanadium nitride \\
 * K, K\_palik           & Polycrystalline potassium &&   W                     & Tungsten \\
 * Li, Li\_palik         & Lithium &&                     \\
 * MgO, MgO\_palik       & Magnesium oxide &&             \\
 * Mo, Mo\_palik         & Molybdenum &&                  none & No mat. profile ($\to$ 100% reflecting mirror)\\
 * \bottomrule
 * \end{tabular}
 * \caption{
 *     \label{tbl:conductor-iors}
 *      This table lists all supported materials that can be passed into the
 *      \pluginref{conductor} and \pluginref{roughconductor} plugins. Note that
 *      some of them are not actually conductors---this is not a problem,
 *      they can be used regardless (though only the reflection component and
 *      no transmission will be simulated). In most cases, there are
 *      multiple entries for each material, which represent measurements by
 *      different authors.
 * }
 * \end{table}
 */
class SmoothConductor : public BSDF {
public:
    SmoothConductor(const Properties &props) : BSDF(props) {
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));

        std::string materialName = props.getString("material", "Cu");

        Spectrum intEta, intK;
        if (boost::to_lower_copy(materialName) == "none") {
            intEta = Spectrum(0.0f);
            intK = Spectrum(1.0f);
        } else {
            intEta.fromContinuousSpectrum(InterpolatedSpectrum(
                fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
            intK.fromContinuousSpectrum(InterpolatedSpectrum(
                fResolver->resolve("data/ior/" + materialName + ".k.spd")));
        }

        Float extEta = lookupIOR(props, "extEta", "air");

        m_eta = props.getSpectrum("eta", intEta) / extEta;
        m_k   = props.getSpectrum("k", intK) / extEta;
    }

    SmoothConductor(Stream *stream, InstanceManager *manager)
            : BSDF(stream, manager) {
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_eta = Spectrum(stream);
        m_k = Spectrum(stream);

        configure();
    }

    void configure() {
        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
            m_specularReflectance, "specularReflectance", 1.0f);

        m_usesRayDifferentials =
            m_specularReflectance->usesRayDifferentials();

        m_components.clear();
        m_components.push_back(EDeltaReflection | EFrontSide
            | (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));

        BSDF::configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);
        manager->serialize(stream, m_specularReflectance.get());
        m_eta.serialize(stream);
        m_k.serialize(stream);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularReflectance") {
            m_specularReflectance = static_cast<Texture *>(child);
            m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
        } else {
            BSDF::addChild(name, child);
        }
    }

    /// Reflection in local coordinates
    inline Vector reflect(const Vector &wi) const {
        return Vector(-wi.x, -wi.y, wi.z);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);

        /* Verify that the provided direction pair matches an ideal
           specular reflection; tolerate some roundoff errors */
        if (!sampleReflection || measure != EDiscrete ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            std::abs(dot(reflect(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
            return Spectrum(0.0f);

        return m_specularReflectance->eval(bRec.its) *
            fresnelConductorExact(Frame::cosTheta(bRec.wi), m_eta, m_k);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);

        /* Verify that the provided direction pair matches an ideal
           specular reflection; tolerate some roundoff errors */
        if (!sampleReflection || measure != EDiscrete ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            std::abs(dot(reflect(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
            return 0.0f;

        return 1.0f;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);

        if (!sampleReflection || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.sampledComponent = 0;
        bRec.sampledType = EDeltaReflection;
        bRec.wo = reflect(bRec.wi);
        bRec.eta = 1.0f;

        return m_specularReflectance->eval(bRec.its) *
            fresnelConductorExact(Frame::cosTheta(bRec.wi), m_eta, m_k);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);

        if (!sampleReflection || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.sampledComponent = 0;
        bRec.sampledType = EDeltaReflection;
        bRec.wo = reflect(bRec.wi);
        bRec.eta = 1.0f;
        pdf = 1;

        return m_specularReflectance->eval(bRec.its) *
            fresnelConductorExact(Frame::cosTheta(bRec.wi), m_eta, m_k);
    }

    Float getRoughness(const Intersection &its, int component) const {
        return 0.0f;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SmoothConductor[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  eta = " << m_eta.toString() << "," << endl
            << "  k = " << m_k.toString() << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_specularReflectance;
    Spectrum m_eta;
    Spectrum m_k;
};

/* Smooth conductor shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least
   something that suggests the presence of a specularly-reflecting
   conductor.

   The code below is an isotropic version of the shader in
   roughconductor.cpp, with \alpha fixed to 0.4f
*/
class SmoothConductorShader : public Shader {
public:
    SmoothConductorShader(Renderer *renderer, const Texture *specularReflectance,
            const Spectrum &eta, const Spectrum &k) : Shader(renderer, EBSDFShader),
            m_specularReflectance(specularReflectance) {
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());

        /* Compute the reflectance at perpendicular incidence */
        m_R0 = fresnelConductorExact(1.0f, eta, k);

        m_alpha = 0.4f;
    }

    bool isComplete() const {
        return m_specularReflectanceShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_specularReflectanceShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_specularReflectance.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_R0);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_R0;" << endl
            << endl
            << "float " << evalName << "_D(vec3 m, float alpha) {" << endl
            << "    alpha = 2 / (alpha * alpha) - 2;" << endl
            << "    return (alpha + 2) * 0.15915 * pow(cosTheta(m), alpha);" << endl
            << "}" << endl
            << endl
            << "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
            << "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
            << "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
            << "        return 0.0;" << endl
            << "    float nDotM = cosTheta(m);" << endl
            << "    return min(1.0, min(" << endl
            << "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
            << "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_schlick(float ct) {" << endl
            << "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
            << "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
            << "        return vec3(0.0);" << endl
            << "   vec3 H = normalize(wi + wo);" << endl
            << "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
            << "   float D = " << evalName << "_D(H, " << m_alpha << ")" << ";" << endl
            << "   float G = " << evalName << "_G(H, wi, wo);" << endl
            << "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);"<< endl
            << "}" << endl;
    }
    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_specularReflectance;
    ref<Shader> m_specularReflectanceShader;
    Spectrum m_R0;
    Float m_alpha;
};

Shader *SmoothConductor::createShader(Renderer *renderer) const {
    return new SmoothConductorShader(renderer,
        m_specularReflectance.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(SmoothConductorShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothConductor, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothConductor, "Smooth conductor");
MTS_NAMESPACE_END
