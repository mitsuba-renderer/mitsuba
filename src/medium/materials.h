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

#if !defined(__MATERIAL_DATA_H)
#define __MATERIAL_DATA_H

#include <mitsuba/mitsuba.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

struct MaterialEntry {
	const char *name;
	Float sigmaS[3];
	Float sigmaA[3];
	Float eta;
};

/* Fitted data from "A Practical Model for Subsurface scattering" (Jensen et al.) */
static MaterialEntry materialData[] = {
	{ "apple",       { 2.29,  2.39,  1.97 }, { 0.0030, 0.0034,  0.046  },  1.3 },
	{ "chicken1",    { 0.15,  0.21,  0.38 }, { 0.0015, 0.077,   0.19   },  1.3 },
	{ "chicken2",    { 0.19,  0.25,  0.32 }, { 0.0018, 0.088,   0.20   },  1.3 },
	{ "cream",       { 7.38,  5.47,  3.15 }, { 0.0002, 0.0028,  0.0163 },  1.3 },
	{ "ketchup",     { 0.18,  0.07,  0.03 }, { 0.061,  0.97,    1.45   },  1.3 },
	{ "marble",      { 2.19,  2.62,  3.00 }, { 0.0021, 0.0041,  0.0071 },  1.5 },
	{ "potato",      { 0.68,  0.70,  0.55 }, { 0.0024, 0.0090,  0.12   },  1.3 },
	{ "skimmilk",    { 0.70,  1.22,  1.90 }, { 0.0014, 0.0025,  0.0142 },  1.3 },
	{ "skin1",       { 0.74,  0.88,  1.01 }, { 0.032,  0.17,    0.48   },  1.3 },
	{ "skin2",       { 1.09,  1.59,  1.79 }, { 0.013,  0.070,   0.145  },  1.3 },
	{ "spectralon",  { 11.6,  20.4,  14.9 }, { 0.00,   0.00,    0.00   },  1.3 },
	{ "wholemilk",   { 2.55,  3.21,  3.77 }, { 0.0011, 0.0024,  0.014  },  1.3 },
	{ NULL,          { 0.00,  0.00,  0.00 }, { 0.00,   0.00,    0.00   },  0.0 }
};

static void lookupMaterial(const Properties &props, Spectrum &sigmaS, Spectrum &sigmaA, Float *eta = NULL, bool requireValues = true) {
	bool hasSigmaAS = props.hasProperty("sigmaS") || props.hasProperty("sigmaA"),
		hasSigmaTAlbedo = props.hasProperty("sigmaT") || props.hasProperty("albedo"),
		manual = hasSigmaAS || hasSigmaTAlbedo,
		hasPreset = props.hasProperty("material");

	if (manual && hasPreset)
		SLog(EError, "Please specify either a preset material or "
			"scattering coefficients (you provided both!)");

	SAssertEx(!props.hasProperty("sizeMultiplier"),
		"The sizeMultiplier property was deprecated!");

	Float densityMultiplier = props.getFloat("densityMultiplier", 1.0f);

	if (manual) {
		if (hasSigmaAS && hasSigmaTAlbedo) {
			SLog(EError, "You can either specify sigmaS & sigmaA *or* "
				"sigmaT & albedo, but no other combinations!");
		}
		if (hasSigmaAS) {
			if (requireValues) {
				sigmaS = props.getSpectrum("sigmaS") * densityMultiplier;
				sigmaA = props.getSpectrum("sigmaA") * densityMultiplier;
			} else {
				sigmaS = props.getSpectrum("sigmaS", Spectrum(0.0f)) * densityMultiplier;
				sigmaA = props.getSpectrum("sigmaA", Spectrum(0.0f)) * densityMultiplier;
			}
		} else {
			Spectrum albedo, sigmaT;
			if (requireValues) {
				albedo = props.getSpectrum("albedo");
				sigmaT = props.getSpectrum("sigmaT") * densityMultiplier;
			} else {
				albedo = props.getSpectrum("albedo", Spectrum(0.0f));
				sigmaT = props.getSpectrum("sigmaT", Spectrum(0.0f)) * densityMultiplier;
			}
			sigmaS = albedo * sigmaT;
			sigmaA = sigmaT - sigmaS;
		}
		if (eta)
			*eta = props.getFloat("eta", 1.3f);
	} else {
		std::string material = 
			boost::to_lower_copy(props.getString("material", "skin1"));
		densityMultiplier *= 100;

		MaterialEntry *matEntry = materialData;
		while (matEntry->name) {
			if (material == matEntry->name) {
				sigmaS.fromLinearRGB(
					matEntry->sigmaS[0],
					matEntry->sigmaS[1],
					matEntry->sigmaS[2]);
				sigmaA.fromLinearRGB(
					matEntry->sigmaA[0],
					matEntry->sigmaA[1],
					matEntry->sigmaA[2]);
				sigmaS *= densityMultiplier;
				sigmaA *= densityMultiplier;
				if (eta)
					*eta = matEntry->eta;
				return;
			}
			++matEntry;
		}

		std::ostringstream oss;
		oss << "Unable to find a material preset for \"" << material
			<< "\"! Valid choices are:";

		/* Unable to find the material preset by name -- print an error
		   message that lists all possible options */
		for (matEntry = materialData; matEntry->name != NULL; ++matEntry) {
			oss << matEntry->name;
			if ((matEntry+1)->name)
				oss << ", ";
		}

		SLog(EError, "%s", oss.str().c_str());
	}
}

MTS_NAMESPACE_END

#endif /* __MATERIAL_DATA_H */
