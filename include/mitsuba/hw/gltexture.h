/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__GLTEXTURE_H)
#define __GLTEXTURE_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL-based GPUTexture implementation
 */
class MTS_EXPORT_HW GLTexture : public GPUTexture {
public:
	/// Create a new GLTexture with the given name and bitmap
	GLTexture(const std::string &name, Bitmap *bitmap);

	/// Upload the texture
	void init();

	/// Refresh (re-upload) the texture
	void refresh();

	/// Free the texture from GPU memory
	void cleanup();
	
	/// Bind the texture and enable texturing
	void bind(int textureUnit = 0) const;

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

	/// Blit a float render buffer into another render buffer
	void blit(GPUTexture *texture) const;

	/// Clear (assuming that this is a render buffer)
	void clear();

	/// Assuming that this is a 2D RGB framebuffer, read a single pixel from the GPU
	Spectrum getPixel(int x, int y) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLTexture();
protected:
	/// Look up relevant constants
	void lookupGLConstants();

	/// Configure texture filtering
	void configureTexture();

	GLuint m_id;
	GLuint m_glType;
	GLuint m_format;
	GLuint m_internalFormat;
	GLuint m_dataFormat;
	/* For render targets */
	GLuint m_fboId, m_depthId;
	GLuint m_extraId;
	mutable bool m_needsUpdate;
};

MTS_NAMESPACE_END

#endif /* __GLTEXTURE_H */
