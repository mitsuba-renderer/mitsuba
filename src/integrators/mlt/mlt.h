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

#if !defined(__MLT_H)
#define __MLT_H

#include <mitsuba/core/bitmap.h>
#include <mitsuba/bidir/manifold.h>
#include <mitsuba/bidir/mut_manifold.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters used
 * by the MLT rendering implementation
 */
struct MLTConfiguration {
    int maxDepth;
    bool separateDirect;
    bool bidirectionalMutation;
    bool causticPerturbation;
    bool lensPerturbation;
    bool multiChainPerturbation;
    bool manifoldPerturbation;
    Float luminance;
    Float probFactor;
    int workUnits;
    int directSamples;
    int luminanceSamples;
    size_t nMutations;
    bool twoStage;
    bool firstStage;
    int firstStageSizeReduction;
    ref<Bitmap> importanceMap;
    size_t timeout;

    inline MLTConfiguration() { }

    void dump() const {
        std::ostringstream oss;
        if (bidirectionalMutation)
            oss << "bidir ";
        if (causticPerturbation)
            oss << "caustic ";
        if (lensPerturbation)
            oss << "lens ";
        if (multiChainPerturbation)
            oss << "multiChain ";
        if (manifoldPerturbation)
            oss << "manifold ";
        SLog(EDebug, "Veach-MLT configuration:");
        SLog(EDebug, "   Maximum path length         : %i", maxDepth);
        SLog(EDebug, "   Active mutators             : %s", oss.str().c_str());
        SLog(EDebug, "   Two-stage MLT               : %s",
            twoStage ? "yes" : "no");
        if (twoStage)
            SLog(EDebug, "   First-stage size reduction  : %i",
                firstStageSizeReduction);
        SLog(EDebug, "   Separate direct illum.      : %s",
            separateDirect ? formatString("%i samples", directSamples).c_str() : "no");

        SLog(EDebug, "   Overall MLT image luminance : %f (%i samples)",
            luminance, luminanceSamples);
        SLog(EDebug, "   Total number of work units  : %i", workUnits);
        SLog(EDebug, "   Mutations per work unit     : " SIZE_T_FMT, nMutations);
        SLog(EDebug, "   Universal perturb. factor   : %f", probFactor);
        SLog(EDebug, "   Manifold max iterations     : %i", MTS_MANIFOLD_MAX_ITERATIONS);
        SLog(EDebug, "   Quantiles                   : %f (surfaces), %f (media)",
                MTS_MANIFOLD_QUANTILE_SURFACE, MTS_MANIFOLD_QUANTILE_MEDIUM);
        if (timeout)
            SLog(EDebug, "   Timeout                     : " SIZE_T_FMT,  timeout);
    }

    inline MLTConfiguration(Stream *stream) {
        maxDepth = stream->readInt();
        separateDirect = stream->readBool();
        bidirectionalMutation = stream->readBool();
        causticPerturbation = stream->readBool();
        lensPerturbation = stream->readBool();
        multiChainPerturbation = stream->readBool();
        manifoldPerturbation = stream->readBool();
        luminance = stream->readFloat();
        probFactor = stream->readFloat();
        workUnits = stream->readInt();
        directSamples = stream->readInt();
        luminanceSamples = stream->readInt();
        nMutations = stream->readSize();
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
        stream->writeInt(maxDepth);
        stream->writeBool(separateDirect);
        stream->writeBool(bidirectionalMutation);
        stream->writeBool(causticPerturbation);
        stream->writeBool(lensPerturbation);
        stream->writeBool(multiChainPerturbation);
        stream->writeBool(manifoldPerturbation);
        stream->writeFloat(luminance);
        stream->writeFloat(probFactor);
        stream->writeInt(workUnits);
        stream->writeInt(directSamples);
        stream->writeInt(luminanceSamples);
        stream->writeSize(nMutations);
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

#endif /* __MLT_H */
