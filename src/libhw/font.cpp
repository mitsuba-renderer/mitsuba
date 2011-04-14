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

#include <mitsuba/hw/font.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/bitmap.h>
#include "veramono14_dsc.h"
#include "veramono14_png.h"
#include "vera14_dsc.h"
#include "vera14_png.h"

MTS_NAMESPACE_BEGIN

Font::Font(EFont font) {
	uint32_t png_size = 0, dsc_size = 0;
	uint8_t *png_ptr = NULL, *dsc_ptr = NULL;

	switch (font) {
		case EBitstreamVera14:
			m_name = "Bitstream Vera 14";
			dsc_ptr = vera14_dsc;
			dsc_size = vera14_dsc_size;
			png_ptr = vera14_png;
			png_size = vera14_png_size;
			break;
		case EBitstreamVeraMono14:
			m_name = "Bitstream Vera Mono 14";
			dsc_ptr = veramono14_dsc;
			dsc_size = veramono14_dsc_size;
			png_ptr = veramono14_png;
			png_size = veramono14_png_size;
			break;
		default:
			Log(EError, "Font is not available!");
	}

	ref<MemoryStream> pngStream = new MemoryStream(png_ptr, png_size);
	ref<MemoryStream> dscStream = new MemoryStream(dsc_ptr, dsc_size);
	dscStream->setByteOrder(Stream::ENetworkByteOrder);
	m_maxVerticalBearing = 0;

	m_bitmap = new Bitmap(Bitmap::EPNG, pngStream);
	for (int i=0; i<256; ++i) {
		Glyph &g = m_glyphs[i];
		
		g.tx.x = dscStream->readSingle();
		g.tx.y = dscStream->readSingle();
		g.ts.x = dscStream->readSingle();
		g.ts.y = dscStream->readSingle();
		g.size = Vector2i(dscStream);
		g.horizontalBearing = dscStream->readInt();
		g.verticalBearing = dscStream->readInt();
		g.horizontalAdvance = dscStream->readInt();
		m_maxVerticalBearing = std::max(m_maxVerticalBearing, g.verticalBearing);
	}
	dscStream->read(m_kerningMatrix, 256*256);
	Assert(dscStream->getPos() == dscStream->getSize());
}

void Font::init(Renderer *renderer) {
	m_texture = renderer->createGPUTexture(m_name, m_bitmap);
	m_texture->setFilterType(GPUTexture::ENearest);
	m_texture->setMipMapped(false);
	m_texture->init();
}

void Font::cleanup() {
	m_texture->cleanup();
}

Font::~Font() {
}

MTS_IMPLEMENT_CLASS(Font, false, Object)
MTS_NAMESPACE_END
