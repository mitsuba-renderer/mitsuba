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

#pragma once
#if !defined(__MITSUBA_HW_BASICSHADER_H_)
#define __MITSUBA_HW_BASICSHADER_H_

#include <mitsuba/render/texture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/* ============================================================ */
/*    Some vary basic texture and shader definitions            */
/*    These classes are in libhw instead of librender, since    */
/*    they link to some functions in libhw.                     */
/* ============================================================ */

/**
 * \brief Constant spectrum-valued texture
 *
 * Includes a \ref Shader implementation for hardware rendering
 * \ingroup libhw
 */
class MTS_EXPORT_HW ConstantSpectrumTexture : public Texture {
public:
    inline ConstantSpectrumTexture(const Spectrum &value)
        : Texture(Properties()), m_value(value) {
    }

    ConstantSpectrumTexture(Stream *stream, InstanceManager *manager);

    inline Spectrum eval(const Intersection &its, bool /* unused */) const {
        return m_value;
    }

    inline Spectrum getAverage() const {
        return m_value;
    }

    inline Spectrum getMaximum() const {
        return m_value;
    }

    inline Spectrum getMinimum() const {
        return m_value;
    }

    inline bool isConstant() const {
        return true;
    }

    inline std::string toString() const {
        return m_value.toString();
    }

    inline bool usesRayDifferentials() const {
        return false;
    }

    inline bool isMonochromatic() const {
        return m_value == Spectrum(m_value[0]);
    }

    Shader *createShader(Renderer *renderer) const;

    ref<Bitmap> getBitmap(const Vector2i &resolutionHint) const;

    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    Spectrum m_value;
};

/**
 * \brief Constant float-valued texture
 *
 * Includes a \ref Shader implementation for hardware rendering
 * \ingroup libhw
 */
class MTS_EXPORT_HW ConstantFloatTexture : public Texture {
public:
    inline ConstantFloatTexture(const Float &value)
        : Texture(Properties()), m_value(value) {
    }

    ConstantFloatTexture(Stream *stream, InstanceManager *manager);

    inline Spectrum eval(const Intersection &its, bool /* unused */) const {
        return Spectrum(m_value);
    }

    inline Spectrum getAverage() const {
        return Spectrum(m_value);
    }

    inline Spectrum getMaximum() const {
        return Spectrum(m_value);
    }

    inline Spectrum getMinimum() const {
        return Spectrum(m_value);
    }

    inline bool isConstant() const {
        return true;
    }

    inline std::string toString() const {
        std::ostringstream oss;
        oss << m_value;
        return oss.str();
    }

    inline bool usesRayDifferentials() const {
        return false;
    }

    inline bool isMonochromatic() const {
        return true;
    }

    Shader *createShader(Renderer *renderer) const;

    ref<Bitmap> getBitmap(const Vector2i &resolutionHint) const;

    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    Float m_value;
};

/**
 * \brief Componentwise addition of two textures
 *
 * Includes a \ref Shader implementation for hardware rendering
 * \ingroup libhw
 */
class MTS_EXPORT_HW SpectrumAdditionTexture : public Texture {
public:
    inline SpectrumAdditionTexture(const Texture *a, const Texture *b)
        : Texture(Properties()), m_a(a), m_b(b) { }

    SpectrumAdditionTexture(Stream *stream, InstanceManager *manager);

    inline Spectrum eval(const Intersection &its, bool /* unused */) const {
        return m_a->eval(its) + m_b->eval(its);
    }

    inline Spectrum getAverage() const {
        return m_a->getAverage() + m_b->getAverage();
    }

    inline Spectrum getMaximum() const {
        // Return a conservative estimate
        return m_a->getMaximum() + m_b->getMaximum();
    }

    inline Spectrum getMinimum() const {
        // Return a conservative estimate
        return m_a->getMinimum() + m_b->getMinimum();
    }

    inline bool isConstant() const {
        return m_a->isConstant() && m_b->isConstant();
    }

    inline std::string toString() const {
        std::ostringstream oss;
        oss << "SpectrumAdditionTexture[" << endl
            << "  a = " << indent(m_a->toString()) << "," << endl
            << "  b = " << indent(m_a->toString()) << endl
            << "]";
        return oss.str();
    }

    inline bool usesRayDifferentials() const {
        return m_a->usesRayDifferentials() || m_b->usesRayDifferentials();
    }

    inline bool isMonochromatic() const {
        return m_a->isMonochromatic() && m_b->isMonochromatic();
    }

    Shader *createShader(Renderer *renderer) const;

    ref<Bitmap> getBitmap(const Vector2i &resolutionHint) const;

    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    ref<const Texture> m_a, m_b;
};

/**
 * \brief Componentwise subtraction of two textures
 *
 * Includes a \ref Shader implementation for hardware rendering
 * \ingroup libhw
 */
class MTS_EXPORT_HW SpectrumSubtractionTexture : public Texture {
public:
    inline SpectrumSubtractionTexture(const Texture *a, const Texture *b)
        : Texture(Properties()), m_a(a), m_b(b) { }

    SpectrumSubtractionTexture(Stream *stream, InstanceManager *manager);

    inline Spectrum eval(const Intersection &its, bool /* unused */) const {
        return m_a->eval(its) - m_b->eval(its);
    }

    inline Spectrum getAverage() const {
        return m_a->getAverage() - m_b->getAverage();
    }

    inline Spectrum getMaximum() const {
        // Return a conservative estimate
        return m_a->getMaximum() - m_b->getMinimum();
    }

    inline Spectrum getMinimum() const {
        // Return a conservative estimate
        return m_a->getMinimum() - m_b->getMaximum();
    }

    inline bool isConstant() const {
        return m_a->isConstant() && m_b->isConstant();
    }

    inline std::string toString() const {
        std::ostringstream oss;
        oss << "SpectrumSubtractionTexture[" << endl
            << "  a = " << indent(m_a->toString()) << "," << endl
            << "  b = " << indent(m_b->toString()) << endl
            << "]";
        return oss.str();
    }

    inline bool usesRayDifferentials() const {
        return m_a->usesRayDifferentials() || m_b->usesRayDifferentials();
    }

    inline bool isMonochromatic() const {
        return m_a->isMonochromatic() && m_b->isMonochromatic();
    }

    Shader *createShader(Renderer *renderer) const;

    ref<Bitmap> getBitmap(const Vector2i &resolutionHint) const;

    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    ref<const Texture> m_a, m_b;
};

/**
 * \brief Componentwise product of two textures
 *
 * Includes a \ref Shader implementation for hardware rendering
 * \ingroup libhw
 */
class MTS_EXPORT_HW SpectrumProductTexture : public Texture {
public:
    inline SpectrumProductTexture(const Texture *a, const Texture *b)
        : Texture(Properties()), m_a(a), m_b(b) { }

    SpectrumProductTexture(Stream *stream, InstanceManager *manager);

    inline Spectrum eval(const Intersection &its, bool /* unused */) const {
        return m_a->eval(its) * m_b->eval(its);
    }

    inline Spectrum getAverage() const {
        SLog(EError, "SpectrumProductTexture::getAverage() -- information unavailable!");
        return Spectrum(0.0f);
    }

    inline Spectrum getMaximum() const {
        // Return a conservative estimate
        return m_a->getMaximum() * m_b->getMaximum();
    }

    inline Spectrum getMinimum() const {
        // Return a conservative estimate
        return m_a->getMinimum() * m_b->getMinimum();
    }

    inline bool isConstant() const {
        return m_a->isConstant() && m_b->isConstant();
    }

    inline std::string toString() const {
        std::ostringstream oss;
        oss << "SpectrumProductTexture[" << endl
            << "  a = " << indent(m_a->toString()) << "," << endl
            << "  b = " << indent(m_b->toString()) << endl
            << "]";
        return oss.str();
    }

    inline bool usesRayDifferentials() const {
        return m_a->usesRayDifferentials() || m_b->usesRayDifferentials();
    }

    inline bool isMonochromatic() const {
        return m_a->isMonochromatic() && m_b->isMonochromatic();
    }

    Shader *createShader(Renderer *renderer) const;

    void serialize(Stream *stream, InstanceManager *manager) const;

    ref<Bitmap> getBitmap(const Vector2i &resolutionHint) const;

    MTS_DECLARE_CLASS()
protected:
    ref<const Texture> m_a, m_b;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_BASICSHADER_H_ */
