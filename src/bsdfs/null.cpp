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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

class Null : public BSDF {
public:
    Null(const Properties &props)
        : BSDF(props) { }

    Null(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);
    }

    void configure() {
        m_components.clear();
        m_components.push_back(ENull | EFrontSide | EBackSide);
        m_usesRayDifferentials = false;
        BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        return Spectrum(((bRec.typeMask & ENull) && measure == EDiscrete) ? 1.0f : 0.0f);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        return ((bRec.typeMask & ENull) && measure == EDiscrete) ? 1.0f : 0.0f;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
        if (bRec.typeMask & ENull) {
            bRec.wo = -bRec.wi;
            bRec.sampledComponent = 0;
            bRec.sampledType = ENull;
            bRec.eta = 1.0f;
            return Spectrum(1.0f);
        } else {
            return Spectrum(0.0f);
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
        if (bRec.typeMask & ENull) {
            bRec.wo = -bRec.wi;
            bRec.sampledComponent = 0;
            bRec.sampledType = ENull;
            bRec.eta = 1.0f;
            pdf = 1;
            return Spectrum(1.0f);
        } else {
            return Spectrum(0.0f);
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
        return 0.0f;
    }

    std::string toString() const {
        return "Null[]";
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
};

// ================ Hardware shader implementation ================

/* Null shader-- render as a 'black box' */
class NullShader : public Shader {
public:
    NullShader(Renderer *renderer) :
        Shader(renderer, EBSDFShader) {
        m_flags = ETransparent;
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return vec3(0.0);" << endl
            << "}" << endl;
        oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return vec3(0.0);" << endl
            << "}" << endl;
    }
    MTS_DECLARE_CLASS()
};

Shader *Null::createShader(Renderer *renderer) const {
    return new NullShader(renderer);
}

MTS_IMPLEMENT_CLASS(NullShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Null, false, BSDF)
MTS_EXPORT_PLUGIN(Null, "Null BSDF");
MTS_NAMESPACE_END
