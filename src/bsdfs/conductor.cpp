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
#include <mitsuba/render/consttexture.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*! \plugin{conductor}{Smooth conductor}
 *
 * \begin{table}
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
 *     \label{tbl:dielectric-iors}
 *      This table lists all supported material names that can be passed into the
 *      \pluginref{conductor} plugin. In most cases, there are two separate
 *       measurements of the same material made using different approaches.
 * }
 * \end{table}
 **
 */
class SmoothConductor : public BSDF {
public:
	SmoothConductor(const Properties &props) : BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		std::string preset = props.getString("preset", "Cu");

		//m_eta = props.getSpectrum("eta", presetEta);
//		m_k = props.getSpectrum("k", presetK);

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

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		return Spectrum(0.0f);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		return Spectrum(0.0f);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		return 0.0f;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothConductor[" << endl
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
