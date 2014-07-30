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
#if !defined(__MITSUBA_HW_FONT_H_)
#define __MITSUBA_HW_FONT_H_

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Utility class used to render text inside OpenGL programs
 * using pre-rasterized TrueType fonts stored as textures.
 *
 * A FreeType2-based generation tool is located in the directory
 * 'tools/linux/fontgen'. Only Latin-1 is supported at the moment.
 * \ingroup libhw
 */
class MTS_EXPORT_HW Font : public Object {
public:
	/// Glyph metrics data structure
	struct Glyph {
		/// Position on the font texture
		Point2 tx;

		/// Size on the font texture
		Vector2 ts;

		/// Glyph size in pixels
		Vector2i size;

		/** \brief Horizontal bearing of this glyph
		 *
		 * (# of pixels between the pen position before
		 * having drawn the glyph and the left bounding
		 * box edge)
		 */
		int32_t horizontalBearing;

		/** \brief Vertical bearing of this glyph
		 *
		 * (Vertical distance between the baseline and the
		 * top of the glyph bounding box)
		 */
		int32_t verticalBearing;

		/** \brief Horizontal advance value of this glyph
		 *
		 * (# of pixels the pen must be advanced
		 * after rendering a glyph)
		 */
		int32_t horizontalAdvance;
	};

	/// List of supplied fonts
	enum EFont {
		EBitstreamVera14,
		EBitstreamVeraMono14
	};

	/// Allocate memory for a certain font
	Font(EFont font);

	/// Draw text to the specified bitmap
	void drawText(Bitmap *dest, Point2i pos, const std::string &text) const;

	/// Compute the size covered by the given string when rendered using this font
	Vector2i getSize(const std::string &text) const;

	/// Upload the font to the GPU
	void init(Renderer *renderer);

	/// Free the GPU memory
	void cleanup();

	/// Convert the underlying bitmap to a different pixel format
	void convert(Bitmap::EPixelFormat pixelFormat,
		Bitmap::EComponentFormat componentFormat, Float gamma);

	/// Return the name of this font
	inline const std::string &getName() const { return m_name; }

	/// Return an entry from the kerning table
	inline int8_t getKerning(char i, char o) const {
		return m_kerningMatrix[(uint8_t) i + (uint8_t) o*256];
	}

	/// Return the associated texture
	inline GPUTexture *getTexture() { return m_texture; }

	/// Return the associated texture (const version)
	inline const GPUTexture *getTexture() const { return m_texture.get(); }

	/// Return the glyph data structure for a specified character
	inline const Glyph &getGlyph(char c) const { return m_glyphs[(uint8_t) c]; }

	/// Return the max. vertical bearing
	inline int getMaxVerticalBearing() const { return m_maxVerticalBearing; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Font();
private:
	std::string m_name;
	ref<GPUTexture> m_texture;
	ref<Bitmap> m_bitmap;
	Glyph m_glyphs[256];
	int8_t m_kerningMatrix[256*256];
	int m_maxVerticalBearing;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_FONT_H_ */
