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

#if !defined(__MIPMAP_H)
#define __MIPMAP_H

#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN

#define MIPMAP_LUTSIZE 128

/** \brief Isotropic/anisotropic EWA mip-map texture map class based on PBRT 
 * \ingroup librender
 */
class MTS_EXPORT_RENDER MIPMap : public Object {
public:
	enum EWrapMode {
		ERepeat = 0,
		EBlack,
		EWhite,
		EClamp
	};

	enum EFilterType {
		/// Elliptically weighted average
		EEWA,
		/// Trilinear filtering
		ETrilinear,
		/// No filtering
		ENone
	};

	/**
	 * Construct a new mip-map from the given texture. Does not
	 * need to have a power-of-two size.
	 */
	MIPMap(int width, int height, Spectrum *pixels, 
		EFilterType filterType = EEWA, EWrapMode wrapMode = ERepeat,
		Float maxAnisotropy = 8.0f);

	/// Construct a mip map from a HDR bitmap
	static ref<MIPMap> fromBitmap(Bitmap *bitmap, 
		EFilterType filterType = EEWA, EWrapMode wrapMode = ERepeat,
		Float maxAnisotropy = 8.0f, 
		Spectrum::EConversionIntent intent = Spectrum::EReflectance);

	/// Do a mip-map lookup at the appropriate level
	Spectrum getValue(Float u, Float v,
		Float dudx, Float dudy, Float dvdx, Float dvdy) const;

	/// Return the number of mip-map levels
	inline int getLevels() const { return m_levels; }

	/// Bilinear interpolation using a triangle filter
	Spectrum triangle(int level, Float x, Float y) const;

	/// Return the width of the represented texture
	inline int getWidth() const { return m_width; }

	/// Return the height of the represented texture
	inline int getHeight() const { return m_height; }

	/// Return a pointer to internal image representation at full resolution
	inline const Spectrum *getImageData() const { return m_pyramid[0]; }
	
	/// Return a pointer to internal image representation at the specified resolution
	inline const Spectrum *getImageData(int level) const { return m_pyramid[level]; }

	/// Return the resolution of the specified level
	inline const Vector2i getLevelResolution(int level) const {
		return Vector2i(m_levelWidth[level], m_levelHeight[level]);
	}

	/// Get the component-wise maximum at the zero level
	inline const Spectrum &getMinimum() const { return m_minimum; }

	/// Get the component-wise minimum
	inline const Spectrum &getMaximum() const { return m_maximum; }

	/// Get the component-wise average
	inline const Spectrum &getAverage() const { return m_average; }

	/// Return a bitmap representation of the full-resolution image
	Bitmap *getBitmap() const;

	/// Return a bitmap representation of the full-resolution image (8 bit/color)
	Bitmap *getLDRBitmap() const;

	MTS_DECLARE_CLASS()
protected:
	/// \cond
	struct ResampleWeight {
		int firstTexel;
		Float weight[4];
	};
	/// \endcond

	/// Calculate weights for up-sampling a texture
	ResampleWeight *resampleWeights(int oldRes, int newRes) const;

	/// Look up a texel at the given hierarchy level
	Spectrum getTexel(int level, int x, int y) const;

	/**
	 * Calculate the elliptically weighted average of a sample with
     * differential uv information
	 */
	Spectrum EWA(Float u, Float v, Float dudx, Float dudy, Float dvdx, 
		Float dvdy, int level) const;

	/* Virtual destructor */
	virtual ~MIPMap();
private:
	int m_width, m_height;
	int m_levels;
	int *m_levelWidth;
	int *m_levelHeight;
	Spectrum **m_pyramid;
	EFilterType m_filterType;
	EWrapMode m_wrapMode;
	Float *m_weightLut;
	Float m_maxAnisotropy;
	Spectrum m_minimum;
	Spectrum m_maximum;
	Spectrum m_average;
};

MTS_NAMESPACE_END

#endif /* __MIPMAP_H */
