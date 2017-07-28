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

#include <mitsuba/core/statistics.h>
#include <mitsuba/bidir/mut_mchain.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("Multi-chain perturbation",
        "Acceptance rate", EPercentage);
static StatsCounter statsGenerated("Multi-chain perturbation",
        "Successful generation rate", EPercentage);

MultiChainPerturbation::MultiChainPerturbation(const Scene *scene, Sampler *sampler,
        MemoryPool &pool, Float minJump, Float coveredArea) :
    m_scene(scene), m_sampler(sampler), m_pool(pool) {

    if (!scene->getSensor()->getClass()->derivesFrom(MTS_CLASS(PerspectiveCamera)))
        Log(EError, "The multi-chain perturbation requires a perspective camera.");

    Vector2i sizeInPixels = scene->getFilm()->getCropSize();
    m_filmRes = Vector2((Float) sizeInPixels.x, (Float) sizeInPixels.y);
    m_imagePlaneArea = m_filmRes.x * m_filmRes.y;

    /* Pixel jump range (in pixels) [Veach, p.354] */
    m_r1 = minJump;
    m_r2 = std::sqrt(coveredArea * m_filmRes.x * m_filmRes.y / M_PI);
    m_theta1 = degToRad(0.0001f);
    m_theta2 = degToRad(0.1f);
    m_logRatio = -math::fastlog(m_r2/m_r1);
    m_thetaLogRatio = -math::fastlog(m_theta2/m_theta1);
}

MultiChainPerturbation::~MultiChainPerturbation() { }

Mutator::EMutationType MultiChainPerturbation::getType() const {
    return EMultiChainPerturbation;
}

Float MultiChainPerturbation::suitability(const Path &path) const {
    int k = path.length(), m = k - 1, l = m-1, nChains = 1;

    while (l-1 >= 0 && (!path.vertex(l)->isConnectable()
            || !path.vertex(l-1)->isConnectable())) {
        if (path.vertex(l)->isConnectable())
            ++nChains;
        --l;
    }

    return (l-1 >= 0) && nChains >= 2;
}

bool MultiChainPerturbation::sampleMutation(
        Path &source, Path &proposal, MutationRecord &muRec, const MutationRecord& sourceMuRec) {
    int k = source.length(), m = k - 1, l = m-1, nChains = 1;

    while (l-1 >= 0 && (!source.vertex(l)->isConnectable()
            || !source.vertex(l-1)->isConnectable())) {
        if (source.vertex(l)->isConnectable())
            ++nChains;
        --l;
    }
    --l;
    if (l < 0 || nChains < 2)
        return false;

    muRec = MutationRecord(EMultiChainPerturbation, l, m, m-l,
        source.getPrefixSuffixWeight(l, m));
    statsAccepted.incrementBase();
    statsGenerated.incrementBase();

    /* Generate a screen-space offset */
    Float r = m_r2 * math::fastexp(m_logRatio * m_sampler->next1D());
    Float phi = m_sampler->next1D() * 2 * M_PI;
    Vector2 offset(r*std::cos(phi), r*std::sin(phi));

    Point2 proposalSamplePosition = source.getSamplePosition() + offset;

    /* Immediately reject if we went off the image plane */
    if (proposalSamplePosition.x <= 0 || proposalSamplePosition.x >= m_filmRes.x
     || proposalSamplePosition.y <= 0 || proposalSamplePosition.y >= m_filmRes.y)
        return false;

    const PerspectiveCamera *sensor = static_cast<const PerspectiveCamera *>(m_scene->getSensor());

    Ray ray;
    if (sensor->sampleRay(ray, proposalSamplePosition, Point2(0.5f), 0.0f).isZero())
        return false;

    Float focusDistance = sensor->getFocusDistance() /
        absDot(sensor->getWorldTransform(0)(Vector(0,0,1)), ray.d);

    /* Correct direction based on the current aperture sample.
       This is necessary to support thin lens cameras */
    Vector d = normalize(ray(focusDistance) - source.vertex(m)->getPosition());

    /* Allocate memory for the proposed path */
    proposal.clear();
    proposal.append(source, 0, l);
    if (l > 0)
        proposal.append(source.edge(l-1));
    proposal.append(source.vertex(l)->clone(m_pool));
    for (int i=l; i<m-1; ++i) {
        proposal.append(m_pool.allocEdge());
        proposal.append(m_pool.allocVertex());
    }
    proposal.append(m_pool.allocEdge());
    proposal.append(source.vertex(m)->clone(m_pool));
    proposal.append(source.edge(m));
    proposal.append(source.vertex(k));

    BDAssert(proposal.vertexCount() == source.vertexCount());
    BDAssert(proposal.edgeCount() == source.edgeCount());

    Float dist = source.edge(m-1)->length
        + perturbMediumDistance(m_sampler, source.vertex(m-1));

    /* Sample a perturbation and propagate it through specular interactions */
    if (!proposal.vertex(m)->perturbDirection(m_scene,
            proposal.vertex(k), proposal.edge(m),
            proposal.edge(m-1), proposal.vertex(m-1), d, dist,
            source.vertex(m-1)->getType(), ERadiance)) {
        proposal.release(l, m+1, m_pool);
        return false;
    }

    /* If necessary, propagate the perturbation through a sequence of
       ideally specular interactions */
    for (int i=m-1; i>l+1; --i) {
        if (!source.vertex(i)->isConnectable()) {
            Float dist = source.edge(i-1)->length +
                perturbMediumDistance(m_sampler, source.vertex(i-1));

            if (!proposal.vertex(i)->propagatePerturbation(m_scene,
                    proposal.vertex(i+1), proposal.edge(i),
                    proposal.edge(i-1), proposal.vertex(i-1),
                    source.vertex(i)->getComponentType(), dist,
                    source.vertex(i-1)->getType(), ERadiance)) {
                proposal.release(l, m+1, m_pool);
                return false;
            }
        } else {
            Vector oldD = source.vertex(i-1)->getPosition()
                - source.vertex(i)->getPosition();
            Float dist = oldD.length();
            oldD /= dist;

            Float theta = m_theta2 * math::fastexp(m_thetaLogRatio * m_sampler->next1D());
            Float phi = 2 * M_PI * m_sampler->next1D();
            Vector newD = Frame(oldD).toWorld(sphericalDirection(theta, phi));

            dist += perturbMediumDistance(m_sampler, source.vertex(i-1));

            if (!proposal.vertex(i)->perturbDirection(m_scene,
                    proposal.vertex(i+1), proposal.edge(i),
                    proposal.edge(i-1), proposal.vertex(i-1),
                    newD, dist, source.vertex(i-1)->getType(),
                    ERadiance)) {
                proposal.release(l, m+1, m_pool);
                return false;
            }
        }
    }

    if (!PathVertex::connect(m_scene,
            l > 0 ? proposal.vertex(l-1) : NULL,
            l > 0 ? proposal.edge(l-1) : NULL,
            proposal.vertex(l),
            proposal.edge(l),
            proposal.vertex(l+1),
            proposal.edge(l+1),
            proposal.vertex(l+2))) {
        proposal.release(l, m+1, m_pool);
        return false;
    }

    proposal.vertex(k-1)->updateSamplePosition(
        proposal.vertex(k-2));

    BDAssert(proposal.matchesConfiguration(source));

    ++statsGenerated;

    return true;
}

Float MultiChainPerturbation::Q(const Path &source, const Path &proposal,
        const MutationRecord &muRec) const {
    int l = muRec.l, m = muRec.m;

    Spectrum weight = muRec.weight *
        proposal.edge(l)->evalCached(proposal.vertex(l), proposal.vertex(l+1),
                PathEdge::EEverything);

    for (int i=m; i>l+1; --i) {
        const PathVertex *v0 = proposal.vertex(i-1),
              *v1 = proposal.vertex(i);
        const PathEdge *edge = proposal.edge(i-1);

        weight *= edge->evalCached(v0, v1,
                PathEdge::ETransmittance | (i != m ? PathEdge::EValueCosineRad : 0));

        if (v0->isMediumInteraction())
            weight /= pdfMediumPerturbation(source.vertex(i-1),
                    source.edge(i-1), edge);
    }

    return 1.0f / weight.getLuminance();
}

void MultiChainPerturbation::accept(const MutationRecord &) {
    ++statsAccepted;
}

MTS_IMPLEMENT_CLASS(MultiChainPerturbation, false, Mutator)
MTS_NAMESPACE_END

