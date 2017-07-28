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

#include <mitsuba/hw/font.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/zstream.h>
#include "data/veramono14_dsc.h"
#include "data/veramono14_png.h"
#include "data/vera14_dsc.h"
#include "data/vera14_png.h"

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

    ref<Stream> pngStream = new MemoryStream(png_ptr, png_size);
    ref<Stream> dscStream = new ZStream(
        new MemoryStream(dsc_ptr, dsc_size), ZStream::EGZipStream);
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
}

void Font::convert(Bitmap::EPixelFormat pixelFormat, Bitmap::EComponentFormat componentFormat, Float gamma) {
    m_bitmap = m_bitmap->convert(pixelFormat, componentFormat, gamma);
}

void Font::drawText(Bitmap *dest, Point2i pos, const std::string &text) const {
    int initial = pos.x;

    for (size_t i=0; i<text.length(); i++) {
        char character = text[i];
        if (character == '\r')
            continue;
        if (character == '\n') {
            pos.x = initial;
            pos.y += (int) (getMaxVerticalBearing()*4.0/3.0);
            continue;
        }

        const Font::Glyph &glyph = getGlyph(character);

        Point2i targetOffset = pos + Vector2i(
            glyph.horizontalBearing,
            getMaxVerticalBearing() - glyph.verticalBearing - 1
        );

        Point2i sourceOffset(
            (int) (glyph.tx.x * m_bitmap->getWidth()),
            (int) (glyph.tx.y * m_bitmap->getHeight()));

        dest->accumulate(m_bitmap.get(), sourceOffset, targetOffset, glyph.size);

        pos.x += glyph.horizontalAdvance;

        if (i+1 < text.length())
            pos.x += getKerning(character, text[i+1]);
    }
}

Vector2i Font::getSize(const std::string &text) const {
    Vector2i size(0, getMaxVerticalBearing());
    int pos = 0;

    for (size_t i=0; i<text.length(); i++) {
        char character = text[i];
        if (character == '\r')
            continue;
        if (character == '\n') {
            size.y += (int) (getMaxVerticalBearing()*(4.0 / 3.0));
            size.x = std::max(size.x, pos);
            pos = 0;
            continue;
        }

        const Font::Glyph &glyph = getGlyph(character);

        pos += glyph.horizontalAdvance;

        if (i+1 < text.length())
            pos += getKerning(character, text[i+1]);
    }
    size.x = std::max(size.x, pos);

    return size;
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
