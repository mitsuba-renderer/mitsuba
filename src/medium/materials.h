/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include "../bsdfs/ior.h"

MTS_NAMESPACE_BEGIN

struct MaterialEntry {
	const char *name;
	Float sigmaS[3];
	Float sigmaA[3];
	Float eta;
};

/* Fitted data from "A Practical Model for Subsurface scattering" (Jensen et al.) */
static MaterialEntry materialData[] = {
	{ "apple",       { 2.29f, 2.39f, 1.97f }, { 0.0030f, 0.0034f, 0.046f  }, 1.3f },
	{ "chicken1",    { 0.15f, 0.21f, 0.38f }, { 0.0015f, 0.077f,  0.19f   }, 1.3f },
	{ "chicken2",    { 0.19f, 0.25f, 0.32f }, { 0.0018f, 0.088f,  0.20f   }, 1.3f },
	{ "cream",       { 7.38f, 5.47f, 3.15f }, { 0.0002f, 0.0028f, 0.0163f }, 1.3f },
	{ "ketchup",     { 0.18f, 0.07f, 0.03f }, { 0.061f,  0.97f,   1.45f   }, 1.3f },
	{ "marble",      { 2.19f, 2.62f, 3.00f }, { 0.0021f, 0.0041f, 0.0071f }, 1.5f },
	{ "potato",      { 0.68f, 0.70f, 0.55f }, { 0.0024f, 0.0090f, 0.12f   }, 1.3f },
	{ "skimmilk",    { 0.70f, 1.22f, 1.90f }, { 0.0014f, 0.0025f, 0.0142f }, 1.3f },
	{ "skin1",       { 0.74f, 0.88f, 1.01f }, { 0.032f,  0.17f,   0.48f   }, 1.3f },
	{ "skin2",       { 1.09f, 1.59f, 1.79f }, { 0.013f,  0.070f,  0.145f  }, 1.3f },
	{ "spectralon",  { 11.6f, 20.4f, 14.9f }, { 0.00f,   0.00f,   0.00f   }, 1.3f },
	{ "wholemilk",   { 2.55f, 3.21f, 3.77f }, { 0.0011f, 0.0024f, 0.014f  }, 1.3f },
	{ NULL,          { 0.00f, 0.00f, 0.00f }, { 0.00f,   0.00f,   0.00f   }, 0.0f }
};

static void lookupMaterial(const Properties &props, Spectrum &sigmaS, Spectrum &sigmaA, Float *eta = NULL) {
	bool hasSigmaAS = props.hasProperty("sigmaS") || props.hasProperty("sigmaA"),
		hasSigmaTAlbedo = props.hasProperty("sigmaT") || props.hasProperty("albedo"),
		hasIOR = props.hasProperty("intIOR") || props.hasProperty("extIOR"),
		manual = hasSigmaAS || hasSigmaTAlbedo,
		hasPreset = props.hasProperty("material");

	if (manual && hasPreset)
		SLog(EError, "Please specify either a preset material or "
			"scattering coefficients (you provided both!)");

	if (props.hasProperty("densityMultiplier") || props.hasProperty("sizeMultiplier"))
		SLog(EError, "The 'densityMultiplier' parameter has been deprecated and is now called 'scale'.");

	if (hasSigmaAS && hasSigmaTAlbedo)
		SLog(EError, "You can either specify sigmaS & sigmaA *or* "
			"sigmaT & albedo, but no other combinations!");

	std::string material =
		boost::to_lower_copy(props.getString("material", "skin1"));

	/* Start with a preset */
	bool found = false;
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
			sigmaS *= 100;
			sigmaA *= 100;
			if (eta)
				*eta = matEntry->eta;
			found = true;
			break;
		}
		++matEntry;
	}

	if (!found) {
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

	if (hasSigmaAS) {
		sigmaS = props.getSpectrum("sigmaS", sigmaS);
		sigmaA = props.getSpectrum("sigmaA", sigmaA);
	} else if (hasSigmaTAlbedo) {
		Spectrum albedo, sigmaT;
		sigmaT = props.getSpectrum("sigmaT", sigmaA + sigmaS);
		albedo = props.getSpectrum("albedo", sigmaS / (sigmaS + sigmaA));
		sigmaS = albedo * sigmaT;
		sigmaA = sigmaT - sigmaS;
	}

	if (eta && hasIOR) {
		/* Specifies the internal index of refraction at the interface */
		Float intIOR = lookupIOR(props, "intIOR", "bk7");

		/* Specifies the external index of refraction at the interface */
		Float extIOR = lookupIOR(props, "extIOR", "air");

		if (intIOR < 0 || extIOR < 0)
			SLog(EError, "The interior and exterior indices of "
				"refraction must be positive!");

		*eta = intIOR / extIOR;
	}

	Float scale = props.getFloat("scale", 1.0f);
	sigmaS *= scale;
	sigmaA *= scale;

	std::ostringstream oss;
	oss << "Medium parameters: sigmaS=" << sigmaS.toString() << ", sigmaA=" << sigmaA.toString();
	SLog(EDebug, "%s", oss.str().c_str());
}

MTS_NAMESPACE_END

#endif /* __MATERIAL_DATA_H */
