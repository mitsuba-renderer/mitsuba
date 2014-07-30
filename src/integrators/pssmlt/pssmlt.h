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

#if !defined(__PSSMLT_H)
#define __PSSMLT_H

#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/core/bitmap.h>

/// Use Kelemen-style mutations in random number space?
#define KELEMEN_STYLE_MUTATIONS 1

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters used by
 * the MLT rendering implementation
 */
struct PSSMLTConfiguration {
	PathSampler::ETechnique technique;
	int maxDepth;
	bool directSampling;
	int rrDepth;
	bool separateDirect;
	Float luminance;
	Float pLarge;
	int workUnits;
	int directSamples;
	int luminanceSamples;
	size_t nMutations;
	bool kelemenStyleWeights;
	Float mutationSizeLow;
	Float mutationSizeHigh;
	bool twoStage;
	bool firstStage;
	int firstStageSizeReduction;
	size_t timeout;
	ref<Bitmap> importanceMap;

	inline PSSMLTConfiguration() { }

	void dump() const {
		SLog(EDebug, "PSSMLT configuration:");
		SLog(EDebug, "   Maximum path depth          : %i", maxDepth);
		SLog(EDebug, "   Bidirectional path tracing  : %s",
			(technique == PathSampler::EBidirectional) ? "yes" : "no");
		SLog(EDebug, "   Two-stage MLT               : %s",
			twoStage ? "yes" : "no");
		if (twoStage)
			SLog(EDebug, "   First-stage size reduction  : %i",
				firstStageSizeReduction);
		SLog(EDebug, "   Direct illum. samples       : %i", directSamples);
		SLog(EDebug, "   Separate direct illum.      : %s",
			separateDirect ? "yes" : "no");
		SLog(EDebug, "   Direct sampling strategies  : %s",
			directSampling ? "yes" : "no");
		SLog(EDebug, "   Russian roulette depth      : %i", rrDepth);
		SLog(EDebug, "   Large step probability      : %f", pLarge);
		SLog(EDebug, "   Kelemen et al. weights      : %s",
			kelemenStyleWeights ? "yes" : "no");
		SLog(EDebug, "   Mutation size               : [%f, %f]",
			mutationSizeLow, mutationSizeHigh);
		SLog(EDebug, "   Overall MLT image luminance : %f (%i samples)",
			luminance, luminanceSamples);
		SLog(EDebug, "   Total number of work units  : %i", workUnits);
		SLog(EDebug, "   Mutations per work unit     : " SIZE_T_FMT, nMutations);
		if (timeout)
			SLog(EDebug, "   Timeout                     : " SIZE_T_FMT,  timeout);
	}

	inline PSSMLTConfiguration(Stream *stream) {
		technique = (PathSampler::ETechnique) stream->readUInt();
		maxDepth = stream->readInt();
		directSampling = stream->readBool();
		rrDepth = stream->readInt();
		separateDirect = stream->readBool();
		luminance = stream->readFloat();
		pLarge = stream->readFloat();
		workUnits = stream->readInt();
		directSamples = stream->readInt();
		luminanceSamples = stream->readInt();
		nMutations = stream->readSize();
		kelemenStyleWeights = stream->readBool();
		mutationSizeLow = stream->readFloat();
		mutationSizeHigh = stream->readFloat();
		twoStage = stream->readBool();
		firstStage = stream->readBool();
		firstStageSizeReduction = stream->readInt();
		Vector2i size(stream);
		if (size != Vector2i(0)) {
			importanceMap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat, size);
			stream->readFloatArray(importanceMap->getFloatData(),
				(size_t) size.x * (size_t) size.y);
		}
		timeout = stream->readSize();
	}

	inline void serialize(Stream *stream) const {
		stream->writeUInt((uint32_t) technique);
		stream->writeInt(maxDepth);
		stream->writeBool(directSampling);
		stream->writeInt(rrDepth);
		stream->writeBool(separateDirect);
		stream->writeFloat(luminance);
		stream->writeFloat(pLarge);
		stream->writeInt(workUnits);
		stream->writeInt(directSamples);
		stream->writeInt(luminanceSamples);
		stream->writeSize(nMutations);
		stream->writeBool(kelemenStyleWeights);
		stream->writeFloat(mutationSizeLow);
		stream->writeFloat(mutationSizeHigh);
		stream->writeBool(twoStage);
		stream->writeBool(firstStage);
		stream->writeInt(firstStageSizeReduction);
		if (importanceMap.get()) {
			importanceMap->getSize().serialize(stream);
			stream->writeFloatArray(importanceMap->getFloatData(),
				(size_t) importanceMap->getWidth() * (size_t) importanceMap->getHeight());
		} else {
			Vector2i(0, 0).serialize(stream);
		}
		stream->writeSize(timeout);
	}
};

MTS_NAMESPACE_END

#endif /* __PSSMLT_H */
