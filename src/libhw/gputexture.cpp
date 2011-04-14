/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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
	m_wrapType = EClampToEdge;
	m_mipmapped = true;
	m_maxAnisotropy = 0.0f;
	m_fbType = ENone;
	m_samples = 1;
	m_depthMode = ECompare;

	if (bitmap != NULL) {
		setBitmap(0, bitmap);
	} else {
		m_type = ETexture2D;
		m_format = ER8G8B8;
	}
}

GPUTexture::~GPUTexture() {
	for (unsigned int i=0; i<m_bitmaps.size(); i++) {
		if (m_bitmaps[i] != NULL)
			m_bitmaps[i]->decRef();
	}
}

void GPUTexture::disassociate() {
	m_textureUnits.cleanup();
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

		switch (bitmap->getBitsPerPixel()) {
			case 8: m_format = EL8; break;
			case 16: m_format = EL8A8; break;
			case 24: m_format = ER8G8B8; break;
			case 32: m_format = ER8G8B8A8; break;
			case 96: m_format = EFloat32RGB; break;
			case 128: m_format = EFloat32RGBA; break;
			default:
				Log(EError, "Invalid bpp for texture creation!");
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
			m_format = EDepth;
			m_mipmapped = false;
			m_wrapType = EClamp;
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

static const char *toString(GPUTexture::ETextureType type) {
    switch (type) {
		case GPUTexture::ETexture1D: return "texture1D";
		case GPUTexture::ETexture2D: return "texture2D";
		case GPUTexture::ETexture3D: return "texture3D";
		case GPUTexture::ETextureCubeMap: return "textureCubeMap";
        default: SLog(EError, "Invalid texture type"); return NULL;
    }
}

static const char *toString(GPUTexture::ETextureFormat format) {
    switch (format) {
		case GPUTexture::EFloat32L: return "Float32L";
		case GPUTexture::EFloat16L: return "Float16L";
		case GPUTexture::EFloat32RGBA: return "Float32RGBA";
		case GPUTexture::EFloat16RGBA: return "Float16RGBA";
		case GPUTexture::EFloat32RGB: return "Float32RGB";
		case GPUTexture::EFloat16RGB: return "Float16RGB";
		case GPUTexture::ER8G8B8: return "R8G8B8";
		case GPUTexture::ER8G8B8A8: return "R8G8B8A8";
		case GPUTexture::EL8: return "L8";
		case GPUTexture::EL8A8: return "L8A8";
		case GPUTexture::EDepth: return "depth";
        default: SLog(EError, "Invalid texture format"); return NULL;
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
		case GPUTexture::EMirroredRepeat: return "mirroredRepeat";
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

std::string GPUTexture::toString() const {
	std::ostringstream oss;
	oss << "GPUTexture[" << endl
		<< "  name = '" << m_name << "'," << endl
		<< "  type = " << mitsuba::toString(m_type) << "," << endl
		<< "  fbType = " << mitsuba::toString(m_fbType) << "," << endl
		<< "  format = " << mitsuba::toString(m_format) << "," << endl
		<< "  size = " << m_size.toString() << "," << endl
		<< "  filterType = " << mitsuba::toString(m_filterType) << "," << endl
		<< "  wrapType = " << mitsuba::toString(m_wrapType) << "," << endl;

	if (m_samples > 1)
		oss << "  samples = " << m_samples << "," << endl;
	if (m_fbType == EDepthBuffer)
		oss << "  depthMode = " << mitsuba::toString(m_depthMode) << "," << endl;

	oss << "  mipmapped = " << m_mipmapped << "," << endl
		<< "  maxAnisotropy = " << m_maxAnisotropy;
	if (m_bitmaps.size() > 0) {
		oss << "," << endl;
		oss << "  bitmaps = {" << endl;
		for (unsigned int i=0; i<m_bitmaps.size(); i++) {
			oss << "    " << i << " => ";
			if (m_bitmaps[i] == NULL)
				oss << "null";
			else
				oss << indent(m_bitmaps[i]->toString());
			if (i != m_bitmaps.size() - 1)
				oss << ",";
			oss << endl;
		}
		oss	<< "  }" << endl;
	} else {
		oss << endl;
	}
	oss << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(GPUTexture, true, Object)
MTS_NAMESPACE_END
