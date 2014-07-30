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
#if !defined(__MITSUBA_HW_GLTEXTURE_H_)
#define __MITSUBA_HW_GLTEXTURE_H_

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL-based GPUTexture implementation
 * \ingroup libhw
 */
class MTS_EXPORT_HW GLTexture : public GPUTexture {
public:
	/// Create a new GLTexture with the given name and bitmap
	GLTexture(const std::string &name, Bitmap *bitmap);

	/// Upload the texture
	void init();

	/// Refresh (re-upload) the texture
	void refresh();

	/**
	 * \brief Refresh (re-upload) a subregion of the texture
	 *
	 * Note: this is only implemented for 2D textures
	 */
	void refresh(const Point2i &offset, const Vector2i &size);

	/// Free the texture from GPU memory
	void cleanup();

	/**
	 * \brief Bind the texture and enable texturing
	 *
	 * \param textureUnit
	 *     Specifies the unit to which this texture should be bound
	 * \param textureIndex
	 *     When this texture has multiple sub-textures (e.g.
	 *     a color and depth map in the case of a
	 *     \ref EColorAndDepthBuffer texture), this parameter
	 *     specifies the one to be bound
	 */
	void bind(int textureUnit = 0, int textureIndex = 0) const;

	/// Download the texture (only for render target textures)
	void download(Bitmap *bitmap = NULL);

	/// Unbind the texture and disable texturing
	void unbind() const;

	/// Activate the render target
	void activateTarget();

	/// Activate a certain face of a cube map as render target
	void activateSide(int side);

	/// Restrict rendering to a sub-region of the texture
	void setTargetRegion(const Point2i &offset, const Vector2i &size);

	/// Deactivate the render target
	void releaseTarget();

	/**
	 * \brief Blit a render buffer into another render buffer
	 *
	 * \param target
	 *     Specifies the target render buffer (or NULL for the framebuffer)
	 * \param what
	 *     A bitwise-OR of the components in \ref EFrameBufferType to copy
	 */
	void blit(GPUTexture *target, int what) const;

	/**
	 * \brief Blit a render buffer into another render buffer
	 *
	 * \param target
	 *     Specifies the target render buffer (or NULL for the framebuffer)
	 * \param what
	 *     A bitwise-OR of the components in \ref EFrameBufferType to copy
	 * \param sourceOffset
	 *     Offset in the source render buffer
	 * \param sourceOffset
	 *     Size of the region to be copied from the source render buffer
	 * \param destOffset
	 *     Offset in the destination render buffer
	 * \param destOffset
	 *     Size of the region to be copied into the dest destination buffer
	 */
	void blit(GPUTexture *target, int what, const Point2i &sourceOffset,
		const Vector2i &sourceSize, const Point2i &destOffset,
		const Vector2i &destSize) const;

	/// Clear (assuming that this is a render buffer)
	void clear();

	/// Assuming that this is a 2D RGB framebuffer, read a single pixel from the GPU
	Color3 getPixel(int x, int y) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLTexture();
protected:
	/// Look up relevant constants
	void lookupGLConstants();

	/// Configure texture filtering
	void configureTexture();

	uint32_t m_id;
	uint32_t m_glType;
	uint32_t m_format;
	uint32_t m_internalFormat;
	uint32_t m_dataFormat;
	uint32_t m_fboId, m_depthId;
	mutable bool m_needsUpdate;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_GLTEXTURE_H_ */
