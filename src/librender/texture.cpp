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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/mipmap.h>

MTS_NAMESPACE_BEGIN

namespace stats {
    StatsCounter mipStorage("Texture system", "Cumulative MIP map memory allocations", EByteCount);
    StatsCounter clampedAnisotropy("Texture system", "Lookups with clamped anisotropy", EPercentage);
    StatsCounter avgEWASamples("Texture system", "Average EWA samples / lookup", EAverage);
    StatsCounter filteredLookups("Texture system", "Filtered texture lookups", EPercentage);

}

Texture::Texture(const Properties &props)
 : ConfigurableObject(props) {
}

Texture::Texture(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
}

Vector3i Texture::getResolution() const {
    return Vector3i(0);
}

Spectrum Texture::eval(const Intersection &its, bool filter) const { NotImplementedError("eval"); }
Spectrum Texture::getAverage() const { NotImplementedError("getAverage"); }
Spectrum Texture::getMinimum() const { NotImplementedError("getMinimum"); }
Spectrum Texture::getMaximum() const { NotImplementedError("getMaximum"); }
bool Texture::isConstant() const { NotImplementedError("isConstant"); }
bool Texture::isMonochromatic() const { NotImplementedError("isMonochromatic"); }
bool Texture::usesRayDifferentials() const { NotImplementedError("usesRayDifferentials"); }
ref<Bitmap> Texture::getBitmap(const Vector2i &) const { NotImplementedError("getBitmap"); }

ref<Texture> Texture::expand() {
    return this;
}

void Texture::evalGradient(const Intersection &_its, Spectrum *gradient) const {
    const Float eps = Epsilon;
    Intersection its(_its);

    Spectrum value = eval(its, false);

    its.p = _its.p + its.dpdu * eps;
    its.uv = _its.uv + Point2(eps, 0);
    Spectrum valueU = eval(its, false);

    its.p = _its.p + its.dpdv * eps;
    its.uv = _its.uv + Point2(0, eps);
    Spectrum valueV = eval(its, false);

    gradient[0] = (valueU - value)*(1/eps);
    gradient[1] = (valueV - value)*(1/eps);
}

Texture::~Texture() { }

void Texture::serialize(Stream *stream, InstanceManager *manager) const {
    ConfigurableObject::serialize(stream, manager);
}

Texture2D::Texture2D(const Properties &props) : Texture(props) {
    if (props.getString("coordinates", "uv") == "uv") {
        m_uvOffset = Point2(
            props.getFloat("uoffset", 0.0f),
            props.getFloat("voffset", 0.0f)
        );
        Float uvscale = props.getFloat("uvscale", 1.0f);
        m_uvScale = Vector2(
            props.getFloat("uscale", uvscale),
            props.getFloat("vscale", uvscale)
        );
    } else {
        Log(EError, "Only UV coordinates are supported at the moment!");
    }
}

Texture2D::Texture2D(Stream *stream, InstanceManager *manager)
 : Texture(stream, manager) {
    m_uvOffset = Point2(stream);
    m_uvScale = Vector2(stream);
}

Texture2D::~Texture2D() {
}

void Texture2D::serialize(Stream *stream, InstanceManager *manager) const {
    Texture::serialize(stream, manager);
    m_uvOffset.serialize(stream);
    m_uvScale.serialize(stream);
}

Spectrum Texture2D::eval(const Intersection &its, bool filter) const {
    Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;
    if (its.hasUVPartials && filter) {
        return eval(uv,
            Vector2(its.dudx * m_uvScale.x, its.dvdx * m_uvScale.y),
            Vector2(its.dudy * m_uvScale.x, its.dvdy * m_uvScale.y));
    } else {
        return eval(uv);
    }
}

void Texture2D::evalGradient(const Intersection &its, Spectrum *gradient) const {
    Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;

    evalGradient(uv, gradient);

    gradient[0] *= m_uvScale.x;
    gradient[1] *= m_uvScale.y;
}

void Texture2D::evalGradient(const Point2 &uv, Spectrum *gradient) const {
    const Float eps = Epsilon;

    Spectrum value = eval(uv);
    Spectrum valueU = eval(uv + Vector2(eps, 0));
    Spectrum valueV = eval(uv + Vector2(0, eps));

    gradient[0] = (valueU - value)*(1/eps);
    gradient[1] = (valueV - value)*(1/eps);
}

ref<Bitmap> Texture2D::getBitmap(const Vector2i &sizeHint) const {
    Vector2i res(sizeHint);
    if (res.x <= 0 || res.y <= 0)
        res = Vector2i(32);

    Float invX = 1.0f / res.x, invY = 1.0f / res.y;

    ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, res);
    Spectrum *target = (Spectrum *) bitmap->getFloatData();
    for (int y=0; y<res.y; ++y)
        for (int x=0; x<res.x; ++x)
            *target++ = eval(Point2((x + 0.5f) * invX, (y + 0.5f) * invY));
    return bitmap;
}

MTS_IMPLEMENT_CLASS(Texture, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(Texture2D, true, Texture)
MTS_NAMESPACE_END
