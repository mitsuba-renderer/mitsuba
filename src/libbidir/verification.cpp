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

#include <mitsuba/bidir/path.h>

#define MTS_BD_MAXERR 1e-4f

MTS_NAMESPACE_BEGIN

static bool validateValue(const std::string name, const Spectrum &value, const Spectrum &cached, std::ostream &os) {
    bool valid = true;

    for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
        Float err = std::abs(cached[i]-value[i]);
        Float mag = std::max(std::abs(cached[i]), std::abs(value[i]));

        if (std::isnan(err))
            valid = false;
        else if (mag < MTS_BD_MAXERR && err > MTS_BD_MAXERR) // absolute error threshold
            valid = false;
        else if (mag > MTS_BD_MAXERR && err/mag > MTS_BD_MAXERR) // relative error threshold
            valid = false;
    }

    if (!valid)
        os << "  " << name << " mismatch: cached=" << cached.toString() << ", computed=" << value.toString() << endl;

    return valid;
}

static bool validateValue(const std::string name, Float value, Float cached, std::ostream &os) {
    bool valid = true;

    Float err = std::abs(cached-value);
    Float mag = std::max(std::abs(cached), std::abs(value));

    if (std::isnan(err))
        valid = false;
    else if (mag < MTS_BD_MAXERR && err > MTS_BD_MAXERR) // absolute error threshold
        valid = false;
    else if (mag > MTS_BD_MAXERR && err/mag > MTS_BD_MAXERR) // relative error threshold
        valid = false;

    if (!valid)
        os << "  " << name << " mismatch: cached=" << cached << ", computed=" << value << endl;

    return valid;
}

static bool validateValue(const std::string name, Vector value, Vector cached, std::ostream &os) {
    bool valid = true;

    Float err = (cached-value).length();
    Float mag = std::max(cached.length(), value.length());

    if (std::isnan(err))
        valid = false;
    else if (mag < MTS_BD_MAXERR && err > MTS_BD_MAXERR) // absolute error threshold
        valid = false;
    else if (mag > MTS_BD_MAXERR && err/mag > MTS_BD_MAXERR) // relative error threshold
        valid = false;

    if (!valid)
        os << "  " << name << " mismatch: cached=" << cached.toString() << ", computed=" << value.toString() << endl;

    return valid;
}

bool PathVertex::verify(const Scene *scene, const PathVertex *pred, const PathVertex *succ,
        ETransportMode mode, std::ostream &os) const {
    if (mode == ERadiance)
        std::swap(pred, succ);

    EMeasure measure = (EMeasure) this->measure;

    Float pdfL = evalPdf(scene, pred, succ, EImportance, measure),
          pdfE = evalPdf(scene, succ, pred, ERadiance, measure);

    bool valid = true;
    valid &= validateValue("pdf[ERadiance]",   pdfE, pdf[ERadiance], os);
    valid &= validateValue("pdf[EImportance]", pdfL, pdf[EImportance], os);

    Spectrum weightL = eval(scene, pred, succ, EImportance, measure),
             weightE = eval(scene, succ, pred, ERadiance, measure);

    if (!isSupernode()) {
        pdfL = evalPdf(scene, pred, succ, EImportance, measure == EArea ? ESolidAngle : measure);
        pdfE = evalPdf(scene, succ, pred, ERadiance, measure == EArea ? ESolidAngle : measure);

        if (isOnSurface() && measure != EDiscrete) {
            if (!succ->isSupernode())
                weightL *= absDot(getShadingNormal(), normalize(succ->getPosition() - getPosition()));
            if (!pred->isSupernode())
                weightE *= absDot(getShadingNormal(), normalize(pred->getPosition() - getPosition()));
        }
    }

    weightL = pdfL != 0 ? (weightL / pdfL) : Spectrum(0.0f);
    weightE = pdfE != 0 ? (weightE / pdfE) : Spectrum(0.0f);

    valid &= validateValue("weight[ERadiance]",   weightE, weight[ERadiance], os);
    valid &= validateValue("weight[EImportance]", weightL, weight[EImportance], os);

    bool isDegenerate;
    switch (type) {
        case EEmitterSupernode:
            isDegenerate = scene->hasDegenerateEmitters();
            break;
        case ESensorSupernode:
            isDegenerate = scene->hasDegenerateSensor();
            break;
        case EEmitterSample:
        case ESensorSample: {
                const PositionSamplingRecord &pRec = getPositionSamplingRecord();
                isDegenerate = static_cast<const AbstractEmitter *>(
                    pRec.object)->getType() & AbstractEmitter::EDeltaDirection;
            }
            break;
        case ESurfaceInteraction: {
                const Shape *shape = getIntersection().shape;
                const BSDF *bsdf = shape->getBSDF();
                isDegenerate = !(bsdf->hasComponent(BSDF::ESmooth) || shape->isEmitter() || shape->isSensor());
            }
            break;
        case EMediumInteraction:
            isDegenerate = false;
            break;
        default:
            SLog(EError, "PathVertex::verify(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return false;
    }

    if (degenerate != isDegenerate) {
        os << "  degeneracy mismatch: cached=" << degenerate << ", computed=" << isDegenerate << endl;
        valid = false;
    }

    return valid;
}

bool PathEdge::verify(const Scene *scene, const PathVertex *pred,
        const PathVertex *succ, ETransportMode mode, std::ostream &os) const {
    if (mode == ERadiance)
        std::swap(pred, succ);

    Spectrum weightL = evalTransmittance(pred, succ),
             weightE = evalTransmittance(succ, pred);
    Float    pdfL    = evalPdf(pred, succ),
             pdfE    = evalPdf(succ, pred);

    weightL = pdfL != 0 ? (weightL / pdfL) : Spectrum(0.0f);
    weightE = pdfE != 0 ? (weightE / pdfE) : Spectrum(0.0f);

    bool valid = true;

    if (!pred->isSupernode() && !succ->isSupernode()) {
        Vector refDireciton = succ->getPosition()-pred->getPosition();
        Float refLength = refDireciton.length();
        refDireciton /= refLength;
        valid &= validateValue("length", refLength, length, os);
        valid &= validateValue("d", refDireciton, d, os);
    }


    valid &= validateValue("weight[ERadiance]",   weightE, weight[ERadiance], os);
    valid &= validateValue("weight[EImportance]", weightL, weight[EImportance], os);
    valid &= validateValue("pdf[ERadiance]",   pdfE, pdf[ERadiance], os);
    valid &= validateValue("pdf[EImportance]", pdfL, pdf[EImportance], os);
    return valid;
}

bool Path::verify(const Scene *scene, ETransportMode mode, std::ostream &os) const {
    std::ostringstream oss;
    bool valid = true;

    for (size_t i=0; i<m_vertices.size(); ++i) {
        const PathVertex *pred = i   > 0 ? m_vertices[i-1] : NULL;
        const PathVertex *succ = i+1 < m_vertices.size() ? m_vertices[i+1] : NULL;
        oss << "Vertex " << i << ":" << endl;
        if (!m_vertices[i]->verify(scene, pred, succ, mode, oss))
            valid = false;
        if (i > 0 && i < m_edges.size())
        {
            oss << "Edge " << i << ":" << endl;
            if (!m_edges[i-1]->verify(scene, pred, m_vertices[i], mode, oss))
                valid = false;
        }
    }

    if (!valid) {
        os << "Detected an inconsistency in the path " << endl;
        os << toString() << endl;
        os << "Inconsistency list:" << endl;
        os << oss.str() << endl;
    }
    return valid;
}

MTS_NAMESPACE_END
