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
#include <mitsuba/bidir/mut_caustic.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("Caustic perturbation",
        "Acceptance rate", EPercentage);
static StatsCounter statsGenerated("Caustic perturbation",
        "Successful generation rate", EPercentage);

CausticPerturbation::CausticPerturbation(const Scene *scene, Sampler *sampler,
        MemoryPool &pool, Float minJump, Float coveredArea) :
    m_scene(scene), m_sampler(sampler), m_pool(pool) {

    if (!scene->getSensor()->getClass()->derivesFrom(MTS_CLASS(PerspectiveCamera)))
        Log(EError, "The caustic perturbation requires a perspective camera.");

    const PerspectiveCamera *camera = static_cast<const PerspectiveCamera *>(scene->getSensor());
    Vector2i filmSize = camera->getFilm()->getSize(),
             cropSize = camera->getFilm()->getCropSize();

    /* Simple heuristic for choosing a jump size: assumes that each
       pixel on the camera subtends the same area on the sphere */
    Float degPerPixel = std::min(
                camera->getXFov() / filmSize.x,
                camera->getYFov() / filmSize.y),
          radPerPixel = degPerPixel * M_PI / 180.0f;

    Float r1 = minJump,
          r2 = std::sqrt(coveredArea * cropSize.x*cropSize.y / M_PI); /* [Veach, p. 354] */

    /* These represent the *desired* angle change range as seen from the camera */
    m_theta1 = radPerPixel * r1;
    m_theta2 = radPerPixel * r2;
    m_logRatio = -math::fastlog(m_theta2 / m_theta1);
}

CausticPerturbation::~CausticPerturbation() { }

Mutator::EMutationType CausticPerturbation::getType() const {
    return ECausticPerturbation;
}

Float CausticPerturbation::suitability(const Path &path) const {
    int k = path.length(), m = k - 1, l = m - 1;

    if (k < 4 || !path.vertex(l)->isConnectable())
        return false;

    --l;
    while (l >= 0 && !path.vertex(l)->isConnectable())
        --l;

    return l >= 1 ? 1.0f : 0.0f;
}

bool CausticPerturbation::sampleMutation(
        Path &source, Path &proposal, MutationRecord &muRec, const MutationRecord& sourceMuRec) {
    int k = source.length(), m = k - 1, l = m - 1;

    if (k < 4 || !source.vertex(l)->isConnectable())
        return false;
    --l;

    while (l >= 0 && !source.vertex(l)->isConnectable())
        --l;

    if (l < 1)
        return false;

    muRec = MutationRecord(ECausticPerturbation, l, m, m-l,
        source.getPrefixSuffixWeight(l, m));
    statsAccepted.incrementBase();
    statsGenerated.incrementBase();

    /* Heuristic perturbation size computation (Veach, p.354) */
    Float lengthE = source.edge(m-1)->length;
    Float lengthL = 0;
    for (int i=l; i<m-1; ++i)
        lengthL += source.edge(i)->length;
    Float factor = lengthE/lengthL,
        theta1 = m_theta1 * factor,
        theta2 = m_theta2 * factor;

    Vector woSource = normalize(source.vertex(l+1)->getPosition()
            - source.vertex(l)->getPosition());
    Float phi = m_sampler->next1D() * 2 * M_PI;
    Float theta = theta2 * math::fastexp(m_logRatio * m_sampler->next1D());
    Vector wo = Frame(woSource).toWorld(sphericalDirection(theta, phi));

    /* Allocate memory for the proposed path */
    proposal.clear();
    proposal.append(source, 0, l+1);
    proposal.append(m_pool.allocEdge());
    for (int i=l+1; i<m; ++i) {
        proposal.append(m_pool.allocVertex());
        proposal.append(m_pool.allocEdge());
    }
    proposal.append(source, m, k+1);
    proposal.vertex(l) = proposal.vertex(l)->clone(m_pool);
    proposal.vertex(m) = proposal.vertex(m)->clone(m_pool);
    BDAssert(proposal.vertexCount() == source.vertexCount());
    BDAssert(proposal.edgeCount() == source.edgeCount());

    Float dist = source.edge(l)->length +
        perturbMediumDistance(m_sampler, source.vertex(l+1));

    /* Sample a perturbation and propagate it through specular interactions */
    if (!proposal.vertex(l)->perturbDirection(m_scene,
            proposal.vertex(l-1), proposal.edge(l-1),
            proposal.edge(l), proposal.vertex(l+1), wo, dist,
            source.vertex(l+1)->getType(), EImportance)) {
        proposal.release(l, m+1, m_pool);
        return false;
    }

    Vector woProposal = normalize(proposal.vertex(l+1)->getPosition()
            - source.vertex(l)->getPosition());
    theta = unitAngle(woSource, woProposal);
    if (theta >= theta2 || theta <= theta1) {
        proposal.release(l, m+1, m_pool);
        return false;
    }

    /* If necessary, propagate the perturbation through a sequence of
       ideally specular interactions */
    for (int i=l+1; i<m-1; ++i) {
        Float dist = source.edge(i)->length +
            perturbMediumDistance(m_sampler, source.vertex(i+1));

        if (!proposal.vertex(i)->propagatePerturbation(m_scene,
                proposal.vertex(i-1), proposal.edge(i-1),
                proposal.edge(i), proposal.vertex(i+1),
                source.vertex(i)->getComponentType(), dist,
                source.vertex(i+1)->getType(), EImportance)) {
            proposal.release(l, m+1, m_pool);
            return false;
        }
    }

    if (!PathVertex::connect(m_scene,
            proposal.vertex(m-2),
            proposal.edge(m-2),
            proposal.vertex(m-1),
            proposal.edge(m-1),
            proposal.vertex(m),
            proposal.edge(m),
            proposal.vertex(m+1))) {
        proposal.release(l, m+1, m_pool);
        return false;
    }

    proposal.vertex(k-1)->updateSamplePosition(
        proposal.vertex(k-2));

    ++statsGenerated;
    return true;
}

Float CausticPerturbation::Q(const Path &source, const Path &proposal,
        const MutationRecord &muRec) const {
    int m = muRec.m, l = muRec.l;

    /* Heuristic perturbation size computation (Veach, p.354) */
    Float lengthE = source.edge(m-1)->length;
    Float lengthL = 0;
    for (int i=l; i<m-1; ++i)
        lengthL += source.edge(i)->length;
    Float factor = lengthE/lengthL,
        theta1 = m_theta1 * factor,
        theta2 = m_theta2 * factor;

    Vector d1 = normalize(source.vertex(l+1)->getPosition()   - source.vertex(l)->getPosition());
    Vector d2 = normalize(proposal.vertex(l+1)->getPosition() - source.vertex(l)->getPosition());
    Float theta = unitAngle(d1, d2);
    if (theta >= theta2 || theta <= theta1)
        return 0.0f;

    Float solidAngleDensity = 1.0f / (2*M_PI * -m_logRatio * std::sin(theta) * theta);

    Spectrum weight = muRec.weight * proposal.edge(m-1)->evalCached(
        proposal.vertex(m-1), proposal.vertex(m), PathEdge::EEverything);

    for (int i=l; i<m-1; ++i) {
        const PathVertex *v0 = proposal.vertex(i),
              *v1 = proposal.vertex(i+1);
        const PathEdge *edge = proposal.edge(i);

        weight *= edge->evalCached(v0, v1,
            PathEdge::ETransmittance | PathEdge::EValueCosineImp);

        if (v1->isMediumInteraction())
            weight /= pdfMediumPerturbation(source.vertex(i+1),
                    source.edge(i), edge);
    }

    const Float lumWeight = weight.getLuminance();
    if (lumWeight <= RCPOVERFLOW)
        return 0.f;

    return solidAngleDensity / lumWeight;
}

void CausticPerturbation::accept(const MutationRecord &) {
    ++statsAccepted;
}

MTS_IMPLEMENT_CLASS(CausticPerturbation, false, Mutator)
MTS_NAMESPACE_END
