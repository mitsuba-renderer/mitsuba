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

#include <mitsuba/hw/gltexture.h>

MTS_NAMESPACE_BEGIN

GPUTexture::GPUTexture(const std::string &name, Bitmap *bitmap)
 : m_name(name) {
    m_filterType = EMipMapLinear;
    m_wrapTypeU = m_wrapTypeV = EClampToEdge;
    m_mipmapped = true;
    m_maxAnisotropy = 0.0f;
    m_fbType = ENone;
    m_samples = 1;
    m_depthMode = ECompare;
    m_borderColor = Color3(static_cast<Float>(0));
    m_size = Point3i(0);

    if (bitmap != NULL) {
        setBitmap(0, bitmap);
    } else {
        m_type = ETexture2D;
        m_pixelFormat = ERGB;
        m_componentFormat = EUInt8;
    }
}

GPUTexture::~GPUTexture() {
    for (size_t i=0; i<m_bitmaps.size(); i++) {
        if (m_bitmaps[i] != NULL)
            m_bitmaps[i]->decRef();
    }
}

void GPUTexture::initAndRelease() {
    init();
    release();
}

void GPUTexture::release() {
    for (size_t i=0; i<m_bitmaps.size(); i++) {
        if (m_bitmaps[i] == NULL)
            continue;
        m_bitmaps[i]->decRef();
        m_bitmaps[i] = NULL;
    }
}

void GPUTexture::setBitmap(unsigned int slot, Bitmap *bitmap) {
    while (slot >= m_bitmaps.size())
        m_bitmaps.push_back(NULL);

    if (slot == 0 && bitmap != NULL) {
        m_size = Point3i(
            bitmap->getWidth(),
            bitmap->getHeight(), 1);

        if (bitmap->getWidth() == 1 || bitmap->getHeight() == 1)
            m_type = ETexture1D;
        else
            m_type = ETexture2D;

        switch (bitmap->getPixelFormat()) {
            case Bitmap::ELuminance: m_pixelFormat = ELuminance; break;
            case Bitmap::ELuminanceAlpha: m_pixelFormat = ELuminanceAlpha; break;
            case Bitmap::ERGB: m_pixelFormat = ERGB; break;
            case Bitmap::ERGBA: m_pixelFormat = ERGBA; break;
#if SPECTRUM_SAMPLES == 3
            case Bitmap::ESpectrum: m_pixelFormat = ERGB; break;
            case Bitmap::ESpectrumAlpha: m_pixelFormat = ERGBA; break;
#endif
            default:
                Log(EError, "Unsupported pixel format %i!",
                    (int) bitmap->getPixelFormat());
        }

        switch (bitmap->getComponentFormat()) {
            case Bitmap::EUInt8: m_componentFormat = EUInt8; break;
            case Bitmap::EUInt16: m_componentFormat = EUInt16; break;
            case Bitmap::EUInt32: m_componentFormat = EUInt32; break;
            case Bitmap::EFloat16: m_componentFormat = EFloat16; break;
            case Bitmap::EFloat32: m_componentFormat = EFloat32; break;
            case Bitmap::EFloat64: m_componentFormat = EFloat64; break;
            default:
                Log(EError, "Unsupported component format %i!",
                    (int) bitmap->getComponentFormat());
        }
    }

    if (m_bitmaps[slot] != NULL)
        m_bitmaps[slot]->decRef();

    m_bitmaps[slot] = bitmap;

    if (bitmap != NULL)
        bitmap->incRef();
}

void GPUTexture::setFrameBufferType(EFrameBufferType type) {
    m_fbType = type;

    switch (m_fbType) {
        case EColorBuffer:
        case EColorAndDepthBuffer:
            break;
        case EDepthBuffer:
            m_pixelFormat = EDepth;
            m_mipmapped = false;
            m_wrapTypeU = m_wrapTypeV = EClamp;
            m_filterType = ELinear;
            break;
        default:
            Log(EError, "Invalid buffer type!");
    }
}

Bitmap *GPUTexture::getBitmap(unsigned int slot) {
    if (slot >= m_bitmaps.size())
        return NULL;
    return m_bitmaps[slot];
}

const Bitmap *GPUTexture::getBitmap(unsigned int slot) const {
    if (slot >= m_bitmaps.size())
        return NULL;
    return m_bitmaps[slot];
}

namespace detail {
    static const char *toString(GPUTexture::ETextureType type) {
        switch (type) {
            case GPUTexture::ETexture1D: return "texture1D";
            case GPUTexture::ETexture2D: return "texture2D";
            case GPUTexture::ETexture3D: return "texture3D";
            case GPUTexture::ETextureCubeMap: return "textureCubeMap";
            default: SLog(EError, "Invalid texture type"); return NULL;
        }
    }

    static const char *toString(GPUTexture::EPixelFormat format) {
        switch (format) {
            case GPUTexture::EDepth: return "depth";
            case GPUTexture::ELuminance: return "luminance";
            case GPUTexture::ELuminanceAlpha: return "luminance-alpha";
            case GPUTexture::ERGB: return "rgb";
            case GPUTexture::ERGBA: return "rgba";
            default: SLog(EError, "Invalid pixel format"); return NULL;
        }
    }

    static const char *toString(GPUTexture::EComponentFormat format) {
        switch (format) {
            case GPUTexture::EUInt8: return "uint8";
            case GPUTexture::EUInt16: return "uint16";
            case GPUTexture::EUInt32: return "uint32";
            case GPUTexture::EFloat16: return "float16";
            case GPUTexture::EFloat32: return "float32";
            case GPUTexture::EFloat64: return "float64";
            default: SLog(EError, "Invalid component format"); return NULL;
        }
    }

    static const char *toString(GPUTexture::EFilterType filter) {
        switch (filter) {
            case GPUTexture::ENearest: return "nearest";
            case GPUTexture::ELinear: return "linear";
            case GPUTexture::EMipMapNearest: return "mipMapNearest";
            case GPUTexture::EMipMapLinear: return "mipMapLinear";
            default: SLog(EError, "Invalid texture filter type"); return NULL;
        }
    }

    static const char *toString(GPUTexture::EWrapType wrap) {
        switch (wrap) {
            case GPUTexture::EClamp: return "clamp";
            case GPUTexture::EClampToEdge: return "clampToEdge";
            case GPUTexture::EClampToBorder: return "clampToBorder";
            case GPUTexture::ERepeat: return "repeat";
            case GPUTexture::EMirror: return "mirroredRepeat";
            default: SLog(EError, "Invalid texture wrap type"); return NULL;
        }
    }

    static const char *toString(GPUTexture::EFrameBufferType fbType) {
        switch (fbType) {
            case GPUTexture::ENone: return "none";
            case GPUTexture::EDepthBuffer: return "depthBuffer";
            case GPUTexture::EColorBuffer: return "colorBuffer";
            case GPUTexture::EColorAndDepthBuffer: return "colorAndDepthBuffer";
            default: SLog(EError, "Invalid framebuffer type"); return NULL;
        }
    }

    static const char *toString(GPUTexture::EDepthMode depthMode) {
        switch (depthMode) {
            case GPUTexture::ENormal: return "normal";
            case GPUTexture::ECompare: return "compare";
            default: SLog(EError, "Invalid depth read mode"); return NULL;
        }
    }
};

std::string GPUTexture::toString() const {
    std::ostringstream oss;
    oss << "GPUTexture[" << endl
        << "  name = '" << m_name << "'," << endl
        << "  type = " << detail::toString(m_type) << "," << endl
        << "  pixelFormat = " << detail::toString(m_pixelFormat) << "," << endl
        << "  componentFormat = " << detail::toString(m_componentFormat) << "," << endl
        << "  fbType = " << detail::toString(m_fbType) << "," << endl
        << "  size = " << m_size.toString() << "," << endl
        << "  filterType = " << detail::toString(m_filterType) << "," << endl
        << "  wrapType = [" << detail::toString(m_wrapTypeU)
        << ", " << detail::toString(m_wrapTypeV) << "]," << endl;

    if (m_samples > 1)
        oss << "  samples = " << m_samples << "," << endl;
    if (m_fbType == EDepthBuffer)
        oss << "  depthMode = " << detail::toString(m_depthMode) << "," << endl;

    oss << "  mipmapped = " << m_mipmapped << "," << endl
        << "  samples = " << m_samples << "," << endl
        << "  maxAnisotropy = " << m_maxAnisotropy;
    if (m_bitmaps.size() > 0) {
        oss << "," << endl;
        oss << "  bitmaps = {" << endl;
        for (size_t i=0; i<m_bitmaps.size(); i++) {
            oss << "    " << i << " => ";
            if (m_bitmaps[i] == NULL)
                oss << "null";
            else
                oss << indent(m_bitmaps[i]->toString(), 2);
            if (i != m_bitmaps.size() - 1)
                oss << ",";
            oss << endl;
        }
        oss << "  }" << endl;
    } else {
        oss << endl;
    }
    oss << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS(GPUTexture, true, Object)
MTS_NAMESPACE_END
