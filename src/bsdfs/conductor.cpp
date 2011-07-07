/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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
#include <mitsuba/render/texture.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{conductor}{Smooth conductor}
 * \order{5}
 * \parameters{
 *     \parameter{preset}{\String}{Name of a material preset, see 
 *           \tblref{conductor-iors}.\!\default{\texttt{Cu} / copper}}
 *     \parameter{eta}{\Spectrum}{Real part of the material's index 
 *           of refraction \default{based on the value of \texttt{preset}}}
 *     \parameter{k}{\Spectrum}{Imaginary part of the material's index of 
 *             refraction, also known as absorption coefficient.
 *             \default{based on the value of \texttt{preset}}}
 *     \lastparameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *        Optional factor used to modulate the reflectance component
 *        \default{1.0}}
 * }
 * \renderings{
 *     \rendering{Measured copper material (the default)}
 *         {bsdf_conductor_copper.jpg}
 *     \rendering{Measured gold material (\lstref{conductor-gold})}
 *         {bsdf_conductor_gold.jpg}
 * }
 
 * This plugin implements a perfectly smooth interface to a conducting material, 
 * such as a metal. For a similar model that instead describes a rough surface 
 * microstructure, take a look at the seperately available 
 * \pluginref{roughconductor} plugin.

 * In contrast to dielectric materials, conductors do not transmit 
 * any light. Their index of refraction is complex-valued and tends to undergo
 * considerable changes throughout the visible color spectrum. 
 * 
 * To faciliate the tedious task of specifying spectrally-varying index of 
 * refraction information, Mitsuba ships with a set of measured data for a 
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
 * The table also contains a few birefingent materials, which are split into
 * separate measurements correponding to their two indices of 
 * refraction (named ``ordinary'' and ``extraordinary ray'').
 *
 * When using this plugin, you should ideally compile Mitsuba with support for 
 * spectral renderings to get the most accurate results. While it also works 
 * in RGB mode, the computations will be much more approximate in this case.
 *
 * \begin{xml}[caption=Material configuration for a smooth conductor with 
 *    measured gold data, label=lst:conductor-gold]
 * <shape type="...">
 *     <bsdf type="conductor">
 *         <string name="preset" value="Au"/>
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
 * AlSb, AlSb\_palik     & Cubic aluminium antimonide &&  Se, Se\_palik         & Selenium (ord. ray) \\
 * Au                    & Gold &&                        Se-e, Se-e\_palik     & Selenium (extr. ray) \\
 * Be, Be\_palik         & Polycrystalline beryllium &&   SiC, SiC\_palik       & Hexagonal silicon carbide \\
 * Cr                    & Chromium &&                    SnTe, SnTe\_palik     & Tin telluride\\
 * CsI, CsI\_palik       & Cubic caesium iodide &&        Ta, Ta\_palik         & Tantalum \\
 * Cu, Cu\_palik         & Copper &&                      Te, Te\_palik         & Trigonal tellurium (ord. ray) \\
 * Cu2O, Cu2O\_palik     & Copper (I) oxide &&            Te-e, Te-e\_palik     & Trigonal tellurium (extr. ray) \\
 * CuO, CuO\_palik       & Copper (II) oxide &&           ThF4, ThF4\_palik     & Polycryst. thorium (IV) fluoride \\
 * d-C, d-C\_palik       & Cubic diamond &&               TiC, TiC\_palik       & Polycrystalline titanium carbide \\
 * Hg, Hg\_palik         & Mercury &&                     TiN, TiN\_palik       & Titanium nitride \\
 * HgTe, HgTe\_palik     & Mercury telluride &&           TiO2, TiO2\_palik     & Tetragonal titan. dioxide (ord. ray) \\
 * Ir, Ir\_palik         & Iridium &&                     TiO2-e, TiO2-e\_palik & Tetragonal titan. dioxide (extr. ray) \\
 * K, K\_palik           & Polycrystalline potassium &&   VC, VC\_palik         & Vanadium carbide \\
 * Li, Li\_palik         & Lithium &&                     V\_palik              & Vanadium \\
 * MgO, MgO\_palik       & Magnesium oxide &&             VN, VN\_palik         & Vanadium nitride \\
 * Mo, Mo\_palik         & Molybdenum &&                  W                     & Tungsten\\
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

		std::string preset = props.getString("preset", "Cu");
		Spectrum presetEta, presetK;
		presetEta.fromContinuousSpectrum(InterpolatedSpectrum(
			fResolver->resolve("data/ior/" + preset + ".eta.spd")));
		presetK.fromContinuousSpectrum(InterpolatedSpectrum(
			fResolver->resolve("data/ior/" + preset + ".k.spd")));

		m_eta = props.getSpectrum("eta", presetEta);
		m_k = props.getSpectrum("k", presetK);

		m_components.push_back(EDeltaReflection | EFrontSide);
		m_usesRayDifferentials = false;
	}

	SmoothConductor(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = Spectrum(stream);
		m_k = Spectrum(stream);
		m_components.push_back(EDeltaReflection | EFrontSide);
		m_usesRayDifferentials = 
			m_specularReflectance->usesRayDifferentials();
	}

	virtual ~SmoothConductor() { }

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

	void configure() {
		BSDF::configure();
		/* Verify the input parameter and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);

		/* Verify that the provided direction pair matches an ideal
		   specular reflection; tolerate some roundoff errors */
		if (!sampleReflection || measure != EDiscrete ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			std::abs(1 - dot(reflect(bRec.wi), bRec.wo)) > Epsilon)
			return Spectrum(0.0f);

		return m_specularReflectance->getValue(bRec.its) *
			fresnelConductor(Frame::cosTheta(bRec.wi), m_eta, m_k);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);

		/* Verify that the provided direction pair matches an ideal
		   specular reflection; tolerate some roundoff errors */
		if (!sampleReflection || measure != EDiscrete ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			std::abs(1 - dot(reflect(bRec.wi), bRec.wo)) > Epsilon)
			return 0.0f;

		return 1.0f;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		
		if (!sampleReflection || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.sampledComponent = 0;
		bRec.sampledType = EDeltaReflection;
		bRec.wo = reflect(bRec.wi);

		return m_specularReflectance->getValue(bRec.its) *
			fresnelConductor(Frame::cosTheta(bRec.wi), m_eta, m_k);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		
		if (!sampleReflection || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.sampledComponent = 0;
		bRec.sampledType = EDeltaReflection;
		bRec.wo = reflect(bRec.wi);
		pdf = 1;

		return m_specularReflectance->getValue(bRec.its) *
			fresnelConductor(Frame::cosTheta(bRec.wi), m_eta, m_k);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothConductor[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  eta = " << m_eta.toString() << "," << endl
			<< "  k = " << m_k.toString() << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_specularReflectance;
	Spectrum m_eta;
	Spectrum m_k;
};

MTS_IMPLEMENT_CLASS_S(SmoothConductor, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothConductor, "Smooth conductor");
MTS_NAMESPACE_END
