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
	{ "Lowfat Milk",                { 0.9124f, 1.0744f, 1.2492f }, { 0.0002f, 0.0004f, 0.0008f }, { 0.9320f, 0.9020f, 0.8590f }, 1.33f },
	{ "Reduced Milk",               { 1.0748f, 1.2209f, 1.3931f }, { 0.0002f, 0.0004f, 0.0010f }, { 0.8190f, 0.7970f, 0.7460f }, 1.33f },
	{ "Regular Milk",               { 1.1873f, 1.3293f, 1.4589f }, { 0.0001f, 0.0003f, 0.0013f }, { 0.7500f, 0.7140f, 0.6810f }, 1.33f },
	{ "Espresso",                   { 0.2707f, 0.2828f, 0.2970f }, { 0.1669f, 0.2287f, 0.3078f }, { 0.9070f, 0.8960f, 0.8800f }, 1.33f },
	{ "Mint Mocha Coffee",          { 0.0916f, 0.1081f, 0.1460f }, { 0.0984f, 0.1519f, 0.2040f }, { 0.9100f, 0.9070f, 0.9140f }, 1.33f },
	{ "Lowfat Soy Milk",            { 0.1418f, 0.1620f, 0.2715f }, { 0.0001f, 0.0005f, 0.0025f }, { 0.8500f, 0.8530f, 0.8420f }, 1.33f },
	{ "Regular Soy Milk",           { 0.2433f, 0.2714f, 0.4563f }, { 0.0001f, 0.0005f, 0.0034f }, { 0.8730f, 0.8580f, 0.8320f }, 1.33f },
	{ "Lowfat Chocolate Milk",      { 0.4277f, 0.4998f, 0.5723f }, { 0.0005f, 0.0016f, 0.0068f }, { 0.9340f, 0.9270f, 0.9160f }, 1.33f },
	{ "Regular Chocolate Milk",     { 0.7352f, 0.9142f, 1.0588f }, { 0.0007f, 0.0030f, 0.0100f }, { 0.8620f, 0.8380f, 0.8060f }, 1.33f },
	{ "Coke",                       { 0.0177f, 0.0208f, 0.0000f }, { 0.6966f, 1.1480f, 1.7169f }, { 0.9650f, 0.9720f, 0.9685f }, 1.33f },
	{ "Pepsi",                      { 0.0058f, 0.0141f, 0.0000f }, { 0.6375f, 0.9849f, 1.4420f }, { 0.9260f, 0.9790f, 0.9525f }, 1.33f },
	{ "Sprite",                     { 0.0069f, 0.0089f, 0.0089f }, { 0.1230f, 0.1194f, 0.1306f }, { 0.9430f, 0.9530f, 0.9520f }, 1.33f },
	{ "Gatorade",                   { 0.2392f, 0.2927f, 0.3745f }, { 0.1617f, 0.1258f, 0.0579f }, { 0.9330f, 0.9330f, 0.9350f }, 1.33f },
	{ "Chardonnay",                 { 0.0030f, 0.0047f, 0.0069f }, { 0.1547f, 0.1701f, 0.3443f }, { 0.9140f, 0.9580f, 0.9750f }, 1.33f },
	{ "White Zinfandel",            { 0.0031f, 0.0048f, 0.0066f }, { 0.1732f, 0.2322f, 0.2847f }, { 0.9190f, 0.9430f, 0.9720f }, 1.33f },
	{ "Merlot",                     { 0.0053f, 0.0000f, 0.0000f }, { 0.7586f, 1.6429f, 1.9196f }, { 0.9740f, 0.9740f, 0.9740f }, 1.33f },
	{ "Budweiser Beer",             { 0.0037f, 0.0069f, 0.0074f }, { 0.1449f, 0.3141f, 0.7286f }, { 0.9170f, 0.9560f, 0.9820f }, 1.33f },
	{ "Coors Light Beer",           { 0.0027f, 0.0055f, 0.0000f }, { 0.0268f, 0.0608f, 0.1521f }, { 0.9180f, 0.9660f, 0.9420f }, 1.33f },
	{ "Clorox",                     { 0.1425f, 0.1723f, 0.1928f }, { 0.0175f, 0.0777f, 0.1372f }, { 0.9120f, 0.9050f, 0.8920f }, 1.33f },
	{ "Apple Juice",                { 0.0201f, 0.0243f, 0.0323f }, { 0.1014f, 0.1858f, 0.4084f }, { 0.9470f, 0.9490f, 0.9450f }, 1.33f },
	{ "Cranberry Juice",            { 0.0128f, 0.0155f, 0.0196f }, { 0.2572f, 0.6145f, 0.8104f }, { 0.9470f, 0.9510f, 0.9740f }, 1.33f },
	{ "Grape Juice",                { 0.0072f, 0.0000f, 0.0000f }, { 0.5428f, 1.2500f, 1.5300f }, { 0.9610f, 0.9610f, 0.9610f }, 1.33f },
	{ "Ruby Grapefruit Juice",      { 0.1617f, 0.1606f, 0.1669f }, { 0.0896f, 0.1911f, 0.2636f }, { 0.9290f, 0.9290f, 0.9310f }, 1.33f },
	{ "White Grapefruit Juice",     { 0.3513f, 0.3669f, 0.5237f }, { 0.0096f, 0.0131f, 0.0395f }, { 0.5480f, 0.5450f, 0.5650f }, 1.33f },
	{ "Shampoo",                    { 0.0104f, 0.0114f, 0.0147f }, { 0.0184f, 0.0596f, 0.0805f }, { 0.9100f, 0.9050f, 0.9200f }, 1.33f },
	{ "Strawberry Shampoo",         { 0.0028f, 0.0032f, 0.0033f }, { 0.0189f, 0.0756f, 0.0989f }, { 0.9270f, 0.9350f, 0.9940f }, 1.33f },
	{ "Head & Shoulders Shampoo",   { 0.2791f, 0.2890f, 0.3086f }, { 0.0883f, 0.1637f, 0.2125f }, { 0.9110f, 0.8960f, 0.8840f }, 1.33f },
	{ "Lemon Tea Powder",           { 0.0798f, 0.0898f, 0.1073f }, { 0.2602f, 0.4902f, 0.7727f }, { 0.9460f, 0.9460f, 0.9490f }, 1.33f },
	{ "Orange Juice Powder",        { 0.1928f, 0.2132f, 0.2259f }, { 0.1449f, 0.3441f, 0.7863f }, { 0.9190f, 0.9180f, 0.9220f }, 1.33f },
	{ "Pink Lemonade Powder",       { 0.1235f, 0.1334f, 0.1305f }, { 0.1165f, 0.2366f, 0.3195f }, { 0.9020f, 0.9020f, 0.9040f }, 1.33f },
	{ "Cappuccino Powder",          { 0.0654f, 0.0882f, 0.1568f }, { 0.1920f, 0.2654f, 0.3272f }, { 0.8490f, 0.8430f, 0.9260f }, 1.33f },
	{ "Salt Powder",                { 0.2485f, 0.2822f, 0.3216f }, { 0.5115f, 0.5863f, 0.6147f }, { 0.8020f, 0.7930f, 0.8210f }, 1.33f },
	{ "Sugar Powder",               { 0.0145f, 0.0162f, 0.0202f }, { 0.0650f, 0.1597f, 0.2578f }, { 0.9210f, 0.9190f, 0.9310f }, 1.33f },
	{ "Suisse Mocha",               { 0.3223f, 0.3583f, 0.4148f }, { 0.1875f, 0.2893f, 0.3796f }, { 0.9070f, 0.8940f, 0.8880f }, 1.33f },

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
