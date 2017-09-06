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

#if !defined(__ERPT_H)
#define __ERPT_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/bidir/manifold.h>
#include <mitsuba/bidir/mut_manifold.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters used
 * by the ERPT rendering implementation
 */
struct ERPTConfiguration {
    int maxDepth;
    bool separateDirect;
    int directSamples;
    bool bidirectionalMutation;
    bool causticPerturbation;
    bool lensPerturbation;
    bool multiChainPerturbation;
    bool manifoldPerturbation;
    Float probFactor;
    size_t chainLength;
    Float numChains;
    int blockSize;
    Float luminance;
    size_t luminanceSamples;
    Float avgAngleChangeSurface;
    Float avgAngleChangeMedium;
    int maxChains;

    inline ERPTConfiguration() { }

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
        SLog(EDebug, "ERPT configuration:");
        SLog(EDebug, "   Maximum path length         : %i", maxDepth);
        SLog(EDebug, "   Chain length                : " SIZE_T_FMT, chainLength);
        SLog(EDebug, "   Average number of chains    : %f", numChains);
        SLog(EDebug, "   Separate direct illum.      : %s",
            separateDirect ? formatString("%i samples", directSamples).c_str() : "no");
        SLog(EDebug, "   Active mutators             : %s", oss.str().c_str());
        SLog(EDebug, "   Block size                  : %i", blockSize);
        SLog(EDebug, "   Overall sample luminance    : %f (%i samples)",
            luminance, luminanceSamples);
        SLog(EDebug, "   Universal perturb. factor   : %f", probFactor);
        SLog(EDebug, "   Manifold max iterations     : %i", MTS_MANIFOLD_MAX_ITERATIONS);
        SLog(EDebug, "   Quantiles                   : %f (surfaces), %f (media)",
                MTS_MANIFOLD_QUANTILE_SURFACE, MTS_MANIFOLD_QUANTILE_MEDIUM);
    }

    inline ERPTConfiguration(Stream *stream) {
        maxDepth = stream->readInt();
        separateDirect = stream->readBool();
        directSamples = stream->readInt();
        bidirectionalMutation = stream->readBool();
        causticPerturbation = stream->readBool();
        lensPerturbation = stream->readBool();
        multiChainPerturbation = stream->readBool();
        manifoldPerturbation = stream->readBool();
        probFactor = stream->readFloat();
        chainLength = stream->readSize();
        numChains = stream->readFloat();
        blockSize = stream->readInt();
        luminance = stream->readFloat();
        luminanceSamples = stream->readSize();
        avgAngleChangeSurface = stream->readFloat();
        avgAngleChangeMedium = stream->readFloat();
        maxChains = stream->readInt();
    }

    inline void serialize(Stream *stream) const {
        stream->writeInt(maxDepth);
        stream->writeBool(separateDirect);
        stream->writeInt(directSamples);
        stream->writeBool(bidirectionalMutation);
        stream->writeBool(causticPerturbation);
        stream->writeBool(lensPerturbation);
        stream->writeBool(multiChainPerturbation);
        stream->writeBool(manifoldPerturbation);
        stream->writeFloat(probFactor);
        stream->writeSize(chainLength);
        stream->writeFloat(numChains);
        stream->writeInt(blockSize);
        stream->writeFloat(luminance);
        stream->writeSize(luminanceSamples);
        stream->writeFloat(avgAngleChangeSurface);
        stream->writeFloat(avgAngleChangeMedium);
        stream->writeInt(maxChains);
    }
};

MTS_NAMESPACE_END

#endif /* __ERPT_H */
