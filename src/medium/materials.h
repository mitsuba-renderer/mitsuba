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
	Float g[3];
	Float eta;
};

static MaterialEntry materialData[] = {
	/* Fitted data from "A Practical Model for Subsurface scattering" (Jensen et al.). No anisotropy data available. */
	{ "Apple",                      { 2.29f, 2.39f, 1.97f }, { 0.0030f, 0.0034f, 0.046f  }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Chicken1",                   { 0.15f, 0.21f, 0.38f }, { 0.0015f, 0.077f,  0.19f   }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Chicken2",                   { 0.19f, 0.25f, 0.32f }, { 0.0018f, 0.088f,  0.20f   }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Cream",                      { 7.38f, 5.47f, 3.15f }, { 0.0002f, 0.0028f, 0.0163f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Ketchup",                    { 0.18f, 0.07f, 0.03f }, { 0.061f,  0.97f,   1.45f   }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Marble",                     { 2.19f, 2.62f, 3.00f }, { 0.0021f, 0.0041f, 0.0071f }, { 0.0f, 0.0f, 0.0f }, 1.5f },
	{ "Potato",                     { 0.68f, 0.70f, 0.55f }, { 0.0024f, 0.0090f, 0.12f   }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Skimmilk",                   { 0.70f, 1.22f, 1.90f }, { 0.0014f, 0.0025f, 0.0142f }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Skin1",                      { 0.74f, 0.88f, 1.01f }, { 0.032f,  0.17f,   0.48f   }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Skin2",                      { 1.09f, 1.59f, 1.79f }, { 0.013f,  0.070f,  0.145f  }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Spectralon",                 { 11.6f, 20.4f, 14.9f }, { 0.00f,   0.00f,   0.00f   }, { 0.0f, 0.0f, 0.0f }, 1.3f },
	{ "Wholemilk",                  { 2.55f, 3.21f, 3.77f }, { 0.0011f, 0.0024f, 0.014f  }, { 0.0f, 0.0f, 0.0f }, 1.3f },

	/* From "Acquiring Scattering Properties of Participating Media by Dilution"
	   by Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen (SIGGRAPH 2006) */
	{ "Lowfat Milk",                { 0.9124, 1.0744, 1.2492 }, { 0.0002, 0.0004, 0.0008 }, { 0.9320, 0.9020, 0.8590 }, 1.33f },
	{ "Reduced Milk",               { 1.0748, 1.2209, 1.3931 }, { 0.0002, 0.0004, 0.0010 }, { 0.8190, 0.7970, 0.7460 }, 1.33f },
	{ "Regular Milk",               { 1.1873, 1.3293, 1.4589 }, { 0.0001, 0.0003, 0.0013 }, { 0.7500, 0.7140, 0.6810 }, 1.33f },
	{ "Espresso",                   { 0.2707, 0.2828, 0.2970 }, { 0.1669, 0.2287, 0.3078 }, { 0.9070, 0.8960, 0.8800 }, 1.33f },
	{ "Mint Mocha Coffee",          { 0.0916, 0.1081, 0.1460 }, { 0.0984, 0.1519, 0.2040 }, { 0.9100, 0.9070, 0.9140 }, 1.33f },
	{ "Lowfat Soy Milk",            { 0.1418, 0.1620, 0.2715 }, { 0.0001, 0.0005, 0.0025 }, { 0.8500, 0.8530, 0.8420 }, 1.33f },
	{ "Regular Soy Milk",           { 0.2433, 0.2714, 0.4563 }, { 0.0001, 0.0005, 0.0034 }, { 0.8730, 0.8580, 0.8320 }, 1.33f },
	{ "Lowfat Chocolate Milk",      { 0.4277, 0.4998, 0.5723 }, { 0.0005, 0.0016, 0.0068 }, { 0.9340, 0.9270, 0.9160 }, 1.33f },
	{ "Regular Chocolate Milk",     { 0.7352, 0.9142, 1.0588 }, { 0.0007, 0.0030, 0.0100 }, { 0.8620, 0.8380, 0.8060 }, 1.33f },
	{ "Coke",                       { 0.0177, 0.0208, 0.0000 }, { 0.6966, 1.1480, 1.7169 }, { 0.9650, 0.9720, 0.9685 }, 1.33f },
	{ "Pepsi",                      { 0.0058, 0.0141, 0.0000 }, { 0.6375, 0.9849, 1.4420 }, { 0.9260, 0.9790, 0.9525 }, 1.33f },
	{ "Sprite",                     { 0.0069, 0.0089, 0.0089 }, { 0.1230, 0.1194, 0.1306 }, { 0.9430, 0.9530, 0.9520 }, 1.33f },
	{ "Gatorade",                   { 0.2392, 0.2927, 0.3745 }, { 0.1617, 0.1258, 0.0579 }, { 0.9330, 0.9330, 0.9350 }, 1.33f },
	{ "Chardonnay",                 { 0.0030, 0.0047, 0.0069 }, { 0.1547, 0.1701, 0.3443 }, { 0.9140, 0.9580, 0.9750 }, 1.33f },
	{ "White Zinfandel",            { 0.0031, 0.0048, 0.0066 }, { 0.1732, 0.2322, 0.2847 }, { 0.9190, 0.9430, 0.9720 }, 1.33f },
	{ "Merlot",                     { 0.0053, 0.0000, 0.0000 }, { 0.7586, 1.6429, 1.9196 }, { 0.9740, 0.9740, 0.9740 }, 1.33f },
	{ "Budweiser Beer",             { 0.0037, 0.0069, 0.0074 }, { 0.1449, 0.3141, 0.7286 }, { 0.9170, 0.9560, 0.9820 }, 1.33f },
	{ "Coors Light Beer",           { 0.0027, 0.0055, 0.0000 }, { 0.0268, 0.0608, 0.1521 }, { 0.9180, 0.9660, 0.9420 }, 1.33f },
	{ "Clorox",                     { 0.1425, 0.1723, 0.1928 }, { 0.0175, 0.0777, 0.1372 }, { 0.9120, 0.9050, 0.8920 }, 1.33f },
	{ "Apple Juice",                { 0.0201, 0.0243, 0.0323 }, { 0.1014, 0.1858, 0.4084 }, { 0.9470, 0.9490, 0.9450 }, 1.33f },
	{ "Cranberry Juice",            { 0.0128, 0.0155, 0.0196 }, { 0.2572, 0.6145, 0.8104 }, { 0.9470, 0.9510, 0.9740 }, 1.33f },
	{ "Grape Juice",                { 0.0072, 0.0000, 0.0000 }, { 0.5428, 1.2500, 1.5300 }, { 0.9610, 0.9610, 0.9610 }, 1.33f },
	{ "Ruby Grapefruit Juice",      { 0.1617, 0.1606, 0.1669 }, { 0.0896, 0.1911, 0.2636 }, { 0.9290, 0.9290, 0.9310 }, 1.33f },
	{ "White Grapefruit Juice",     { 0.3513, 0.3669, 0.5237 }, { 0.0096, 0.0131, 0.0395 }, { 0.5480, 0.5450, 0.5650 }, 1.33f },
	{ "Shampoo",                    { 0.0104, 0.0114, 0.0147 }, { 0.0184, 0.0596, 0.0805 }, { 0.9100, 0.9050, 0.9200 }, 1.33f },
	{ "Strawberry Shampoo",         { 0.0028, 0.0032, 0.0033 }, { 0.0189, 0.0756, 0.0989 }, { 0.9270, 0.9350, 0.9940 }, 1.33f },
	{ "Head & Shoulders Shampoo",   { 0.2791, 0.2890, 0.3086 }, { 0.0883, 0.1637, 0.2125 }, { 0.9110, 0.8960, 0.8840 }, 1.33f },
	{ "Lemon Tea Powder",           { 0.0798, 0.0898, 0.1073 }, { 0.2602, 0.4902, 0.7727 }, { 0.9460, 0.9460, 0.9490 }, 1.33f },
	{ "Orange Juice Powder",        { 0.1928, 0.2132, 0.2259 }, { 0.1449, 0.3441, 0.7863 }, { 0.9190, 0.9180, 0.9220 }, 1.33f },
	{ "Pink Lemonade Powder",       { 0.1235, 0.1334, 0.1305 }, { 0.1165, 0.2366, 0.3195 }, { 0.9020, 0.9020, 0.9040 }, 1.33f },
	{ "Cappuccino Powder",          { 0.0654, 0.0882, 0.1568 }, { 0.1920, 0.2654, 0.3272 }, { 0.8490, 0.8430, 0.9260 }, 1.33f },
	{ "Salt Powder",                { 0.2485, 0.2822, 0.3216 }, { 0.5115, 0.5863, 0.6147 }, { 0.8020, 0.7930, 0.8210 }, 1.33f },
	{ "Sugar Powder",               { 0.0145, 0.0162, 0.0202 }, { 0.0650, 0.1597, 0.2578 }, { 0.9210, 0.9190, 0.9310 }, 1.33f },
	{ "Suisse Mocha",               { 0.3223, 0.3583, 0.4148 }, { 0.1875, 0.2893, 0.3796 }, { 0.9070, 0.8940, 0.8880 }, 1.33f },

	{ NULL,                         { 0.0f,  0.0f,   0.0f   }, { 0.0f,   0.0f,   0.0f   }, { 0.0f,   0.0f,  0.0f   }, 0.0f  }
};

static void lookupMaterial(const Properties &props, Spectrum &sigmaS, Spectrum &sigmaA, Spectrum &g, Float *eta = NULL) {
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
		boost::to_lower_copy(props.getString("material", "Skin1"));

	/* Start with a preset */
	bool found = false;
	MaterialEntry *matEntry = materialData;
	while (matEntry->name) {
		if (material == boost::to_lower_copy(std::string(matEntry->name))) {
			sigmaS.fromLinearRGB(
				matEntry->sigmaS[0],
				matEntry->sigmaS[1],
				matEntry->sigmaS[2]);
			sigmaA.fromLinearRGB(
				matEntry->sigmaA[0],
				matEntry->sigmaA[1],
				matEntry->sigmaA[2]);
			g.fromLinearRGB(
				matEntry->g[0],
				matEntry->g[1],
				matEntry->g[2]);
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

	if (props.hasProperty("g")) {
		if (props.getType("g") == Properties::ESpectrum)
			g = props.getSpectrum("g");
		else
			g = Spectrum(props.getFloat("g"));
	}

	if (g.min() <= -1 || g.max() >= 1)
		SLog(EError, "The anisotropy parameter 'g' must be in the range (-1, 1)!");

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
	oss << "Medium parameters: sigmaS=" << sigmaS.toString() << ", sigmaA=" << sigmaA.toString() << ", g=" << g.average();
	SLog(EDebug, "%s", oss.str().c_str());
}

MTS_NAMESPACE_END

#endif /* __MATERIAL_DATA_H */
