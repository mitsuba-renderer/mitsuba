#if !defined(__MIPMAP_H)
#define __MIPMAP_H

#include <mitsuba/render/records.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN

#define MIPMAP_LUTSIZE 128

/** \brief Isotropic/anisotropic EWA mip-map texture map class based on PBRT 
 */
class MTS_EXPORT_RENDER MIPMap : public Object {
public:
	enum EWrapMode {
		ERepeat,
		EBlack,
		EClamp
	};

	/**
	 * Construct a new mip-map from the given texture. Does not
	 * need to have a power-of-two size.
	 */
	MIPMap(int width, int height, Spectrum *pixels, 
		bool isotropic = false, EWrapMode wrapMode = ERepeat,
		Float maxAnisotropy = 8.0f);

	/// Construct a mip map from a HDR bitmap
	static ref<MIPMap> fromBitmap(Bitmap *bitmap);

	/**
	 * Do a mip-map lookup at the appropriate level
	 */
	Spectrum getValue(const Intersection &its) const;

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

	/// Get the component-wise maximum at the zero level
	Spectrum getMaximum() const;

	/// Return a bitmap representation of the full-resolution image
	Bitmap *getBitmap() const;

	/// Return a bitmap representation of the full-resolution image (8 bit/color)
	Bitmap *getLDRBitmap() const;
protected:
	struct ResampleWeight {
		int firstTexel;
		Float weight[4];
	};

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

	MTS_DECLARE_CLASS()
private:
	int m_width, m_height;
	int m_levels;
	int *m_levelWidth;
	int *m_levelHeight;
	Spectrum **m_pyramid;
	bool m_isotropic;
	EWrapMode m_wrapMode;
	Float *m_weightLut;
	Float m_maxAnisotropy;
};

MTS_NAMESPACE_END

#endif /* __MIPMAP_H */
