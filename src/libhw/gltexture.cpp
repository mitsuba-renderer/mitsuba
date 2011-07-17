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

#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/gltexture.h>

#ifndef GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT
#define GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT 0x8CD8
#endif

MTS_NAMESPACE_BEGIN

GLTexture::GLTexture(const std::string &name, Bitmap *bitmap)
 : GPUTexture(name, bitmap), m_id(0), m_needsUpdate(true) {
}

void GLTexture::init() {
	if (m_fbType == ENone)
		Log(ETrace, "Uploading a texture : %s", toString().c_str());
	else
		Log(ETrace, "Creating a framebuffer : %s", toString().c_str());

	lookupGLConstants();

	/* Generate an identifier */
	glGenTextures(1, &m_id);

	/* Bind to the texture */
	glBindTexture(m_glType, m_id);

	/* Set the texture filtering / wrapping modes 
	   (don't do this for multisample textures)*/
	if (!((m_fbType & EColorBuffer) && m_samples > 1))
		configureTexture(); /* Multisample textures don't have parameters */

	if (m_samples > 1) {
		int samples;
		glGetIntegerv(GL_MAX_SAMPLES_EXT, &samples);
		if (m_samples > samples) {
			Log(EWarn, "Attempted to create a multisample framebuffer "
				"with an unsupported # of samples (requested=%i, supported=%i)", 
				m_samples, samples);
			m_samples = samples;
		}
	}


	if (m_fbType == ENone) {
		Assert(m_samples == 1);
		refresh();
	} else {
		/* Create the FBO and bind it */
		glGenFramebuffersEXT(1, &m_fboId);
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fboId);

		AssertEx(glIsFramebufferEXT(m_fboId), "Creating an FBO failed");
		bool depthAsTexture = m_fbType & EDepthBuffer;

		switch (m_fbType) {
			case EColorAndDepthBuffer: 
			case EColorBuffer: {
					if (m_type == ETexture2D) {
						if (!depthAsTexture) {
							glGenRenderbuffersEXT(1, &m_depthId);
							glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, m_depthId);
							if (m_samples == 1)
								glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, 
									GL_DEPTH_COMPONENT, m_size.x, m_size.y);
							else
								glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, 
									m_samples, GL_DEPTH_COMPONENT, m_size.x, m_size.y);
							glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, 
								GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_depthId);
						} else {
							glGenTextures(1, &m_depthId);
							glBindTexture(m_glType, m_depthId);
							configureTexture();
							glTexParameteri(m_glType, GL_TEXTURE_COMPARE_MODE, GL_NONE);
							glTexParameteri(m_glType, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
							if (m_samples == 1)
								glTexImage2D(m_glType, 0, GL_DEPTH_COMPONENT24, m_size.x, m_size.y,
									0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
							else
								glTexImage2DMultisample(m_glType, 
									m_samples, GL_DEPTH_COMPONENT24, m_size.x, m_size.y, GL_FALSE);
							glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, 
								GL_DEPTH_ATTACHMENT, m_glType, m_depthId, 0);
							glBindTexture(m_glType, m_id);
						}

						if (m_samples == 1) 
							glTexImage2D(m_glType, 0, m_internalFormat, m_size.x, m_size.y,
								0, m_format, m_dataFormat, NULL);
						else
							glTexImage2DMultisample(m_glType, 
								m_samples, m_internalFormat, m_size.x, m_size.y, GL_FALSE);

						if (isMipMapped())
							glGenerateMipmapEXT(m_glType);

						glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, 
							GL_COLOR_ATTACHMENT0_EXT, m_glType, m_id, 0);
					} else if (m_type == ETextureCubeMap) {
						Assert(m_size.x == m_size.y && isPowerOfTwo(m_size.x));
						Assert(m_fbType == EColorBuffer);
						Assert(m_samples == 1);

						for (int i=0; i<6; i++) 
							glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, m_internalFormat, 
								m_size.x, m_size.y, 0, m_format, m_dataFormat, NULL);

						if (isMipMapped())
							glGenerateMipmapEXT(m_glType);

						/* Generate an identifier */
						glGenTextures(1, &m_depthId);
						glBindTexture(m_glType, m_depthId);
						glTexParameteri(m_glType, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
						glTexParameteri(m_glType, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

						for (int i=0; i<6; i++) 
							glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_DEPTH_COMPONENT24, 
								m_size.x, m_size.y, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);

						if (glewIsSupported("GL_EXT_geometry_shader4"))
							activateSide(-1);
						else
							activateSide(0);
					} else {
						Log(EError, "Unsupported texture type!");
					}
				}
				break;
			case EDepthBuffer:
				Assert(m_samples == 1);
				if (m_depthMode == ECompare) {
					glTexParameteri(m_glType, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
					glTexParameteri(m_glType, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
				}

				if (m_type == ETexture2D) {
					/* Allocate the texture memory */
					glTexImage2D(m_glType, 0, m_internalFormat, 
						m_size.x, m_size.y, 0, GL_DEPTH_COMPONENT, 
						GL_UNSIGNED_BYTE, NULL);

					/* Attach the texture as a depth target */
					glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, 
						GL_DEPTH_ATTACHMENT_EXT, m_glType, m_id, 0);
				} else if (m_type == ETextureCubeMap) {
					Assert(m_size.x == m_size.y && isPowerOfTwo(m_size.x));
					for (int i=0; i<6; i++)
						glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, m_internalFormat, 
							m_size.x, m_size.y, 0, m_format, m_dataFormat, NULL);

					if (glewIsSupported("GL_EXT_geometry_shader4"))
						activateSide(-1);
					else
						activateSide(0);
				} else {
					Log(EError, "Unsupported texture type!");
				}

				glDrawBuffer(GL_NONE);
				glReadBuffer(GL_NONE);
				break;
			default:
				Log(EError, "Invalid render buffer type!");
		}
	
		GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
		switch (status) {
			case GL_FRAMEBUFFER_COMPLETE_EXT:
				break;
			case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
				Log(EError, "FBO Error: Incomplete attachment!");
			case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
				Log(EError, "FBO Error: Unsupported framebuffer format!");
			case GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT:
				Log(EError, "FBO Error: Incomplete framebuffer - duplicate attachment!");
			case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
				Log(EError, "FBO Error: Incomplete framebuffer - missing attachment!");
			case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
				Log(EError, "FBO Error: Incomplete framebuffer - invalid dimensions!");
			case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
				Log(EError, "FBO Error: Incomplete framebuffer - no draw buffer!");
			case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
				Log(EError, "FBO Error: Incomplete framebuffer - invalid formats!");
			case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
				Log(EError, "FBO Error: Incomplete framebuffer - no readbuffer!");
			case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE_EXT:
				Log(EError, "FBO Error: Incomplete multisample framebuffer!");
			case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS:
				Log(EError, "FBO Error: Incomplete layer targets!");
			default:
				Log(EError, "FBO Error: Unknown error status (0x%x)!", status);
		}
	
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, GL_NONE);
	}
}

void GLTexture::refresh() {
	Bitmap *bitmap = getBitmap();

	/* Bind to the texture */
	glBindTexture(m_glType, m_id);
	
	if (m_type == ETexture1D) {
		Assert((isPowerOfTwo(m_size.x) && m_size.y == 1)
			|| (isPowerOfTwo(m_size.y) && m_size.x == 1));

		if (isMipMapped()) {
			/* Let GLU generate mipmaps for us */
			gluBuild1DMipmaps(m_glType, m_internalFormat, m_size.x == 1 ? m_size.y : m_size.x,
				m_format, m_dataFormat, bitmap->getData());
		} else {
			glTexImage1D(m_glType, 0, m_internalFormat, m_size.x == 1 ? m_size.y : m_size.x,
				0, m_format, m_dataFormat, bitmap->getData());
		}
	} else if (m_type == ETexture2D) {
		//Assert(m_size.x == m_size.y);

		/* Anisotropic texture filtering */
		float anisotropy = (float) getMaxAnisotropy();
		if (anisotropy > 1.0f) {
			if (isMipMapped() && m_filterType == EMipMapLinear) {
				/* Only use anisotropy when it makes sense - otherwise some
				   GL implementations will enforce mipmapping */
				glTexParameterf(m_glType, GL_TEXTURE_MAX_ANISOTROPY_EXT, 
					anisotropy);
			}
		}

		glTexImage2D(m_glType, 0, m_internalFormat, m_size.x, m_size.y,
			0, m_format, m_dataFormat, bitmap->getData());

		if (isMipMapped())
			glGenerateMipmapEXT(m_glType);
	} else if (m_type == ETextureCubeMap) {
		Assert(bitmap != NULL);
		Assert(bitmap->getWidth() == bitmap->getHeight());
		Assert(isPowerOfTwo(bitmap->getWidth()));
		if (isMipMapped())
			glTexParameteri(m_glType, GL_GENERATE_MIPMAP, GL_TRUE);

		for (int i=0; i<6; i++) {
			Bitmap *bitmap = getBitmap(i);

			GLuint pos;
			switch (i) {
				case ECubeMapPositiveX: pos = GL_TEXTURE_CUBE_MAP_POSITIVE_X; break;
				case ECubeMapNegativeX: pos = GL_TEXTURE_CUBE_MAP_NEGATIVE_X; break;
				case ECubeMapPositiveY: pos = GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT; break;
				case ECubeMapNegativeY: pos = GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT; break;
				case ECubeMapPositiveZ: pos = GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT; break;
				case ECubeMapNegativeZ: pos = GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT; break;
				default: Log(EError, "Unknown cube map index"); return;
			};

			glTexImage2D(pos, 0, m_internalFormat, bitmap->getWidth(), bitmap->getHeight(),
				0, m_format, m_dataFormat, bitmap->getData());
		}
	} else {
		Log(EError, "Unknown texture type!");
	}
}

void GLTexture::lookupGLConstants() {
	m_dataFormat = GL_UNSIGNED_BYTE;

	/* Convert the texture type */
	switch (m_type) {
		case ETexture1D: m_glType = GL_TEXTURE_1D; break;
		case ETexture2D: 
			if (m_samples == 1)
				m_glType = GL_TEXTURE_2D;
			else
				m_glType = GL_TEXTURE_2D_MULTISAMPLE;
			break;
		case ETexture3D: m_glType = GL_TEXTURE_3D; break;
		case ETextureCubeMap: m_glType = GL_TEXTURE_CUBE_MAP_EXT; break;
		default: Log(EError, "Invalid texture type specified"); return;
	}

	/* Convert the texture format */
	switch (getFormat()) {
		case ER8G8B8:
			m_format = GL_RGB;
			m_internalFormat = GL_RGB;
			break;
		case ER8G8B8A8:
			m_format = GL_RGBA;
			m_internalFormat = GL_RGBA;
			break;
		case EL8:
			m_format = GL_LUMINANCE;
			m_internalFormat = GL_LUMINANCE8;
			break;
		case EL8A8:
			m_format = GL_LUMINANCE_ALPHA;
			m_internalFormat = GL_LUMINANCE8_ALPHA8;
			break;
		case EDepth:
			m_format = GL_DEPTH_COMPONENT;
			m_internalFormat = GL_DEPTH_COMPONENT32;
			break;
		case EFloat16L:
			m_format = GL_LUMINANCE;
			m_internalFormat = GL_RGBA16F_ARB;
			m_dataFormat = GL_FLOAT;
			break;
		case EFloat16RGB:
			m_format = GL_RGB;
			m_internalFormat = GL_RGBA16F_ARB;
			m_dataFormat = GL_FLOAT;
			break;
		case EFloat16RGBA:
			m_format = GL_RGBA;
			m_internalFormat = GL_RGBA16F_ARB;
			m_dataFormat = GL_FLOAT;
			break;
		case EFloat32L:
			m_format = GL_LUMINANCE;
			m_internalFormat = GL_RGBA32F_ARB;
			m_dataFormat = GL_FLOAT;
			break;
		case EFloat32RGB:
			m_format = GL_RGB;
			m_internalFormat = GL_RGBA32F_ARB;
			m_dataFormat = GL_FLOAT;
			break;
		case EFloat32RGBA:
			m_format = GL_RGBA;
			m_internalFormat = GL_RGBA32F_ARB;
			m_dataFormat = GL_FLOAT;
			break;
		default: Log(EError, "Invalid texture format specified"); return;
	}

}


void GLTexture::configureTexture() {
	int wrap, mag_filter, min_filter;

	/* Convert the texture coordinate wrapping type */
	switch (getWrapType()) {
		case EClamp: wrap = GL_CLAMP; break;
		case EClampToEdge: wrap = GL_CLAMP_TO_EDGE; break;
		case EClampToBorder: wrap = GL_CLAMP_TO_BORDER; break;
		case ERepeat: wrap = GL_REPEAT; break;
		case EMirroredRepeat: wrap = GL_MIRRORED_REPEAT_ARB; break;
		default: Log(EError, "Invalid texture wrap type specified"); return;
	}

	/* Convert the texture filter type */
	switch (m_filterType) {
		case ENearest:
			mag_filter = GL_NEAREST;
			min_filter = isMipMapped() ? GL_NEAREST_MIPMAP_NEAREST : GL_NEAREST;
			break;
		case ELinear:
			mag_filter = GL_LINEAR;
			min_filter = isMipMapped() ? GL_NEAREST_MIPMAP_LINEAR : GL_LINEAR;
			break;
		case EMipMapNearest:
			if (!isMipMapped()) {
				mag_filter = GL_LINEAR;
				min_filter = GL_LINEAR;
				break;
			}

			mag_filter = GL_LINEAR;
			min_filter = GL_LINEAR_MIPMAP_NEAREST;
			break;
		case EMipMapLinear:
			if (!isMipMapped()) {
				mag_filter = GL_LINEAR;
				min_filter = GL_LINEAR;
				break;
			}

			mag_filter = GL_LINEAR;
			min_filter = GL_LINEAR_MIPMAP_LINEAR;
			break;
		default:
			Log(EError, "Invalid filter type specified");
			return;
	}

	/* Set the filter type */
	glTexParameteri(m_glType, GL_TEXTURE_MAG_FILTER, mag_filter);
	glTexParameteri(m_glType, GL_TEXTURE_MIN_FILTER, min_filter);

	/* Set the texcoord wrapping type */
	if (m_type == ETexture1D) {
		glTexParameteri(m_glType, GL_TEXTURE_WRAP_S, wrap);
	} else if (m_type == ETexture2D) {
		glTexParameteri(m_glType, GL_TEXTURE_WRAP_S, wrap);
		glTexParameteri(m_glType, GL_TEXTURE_WRAP_T, wrap);
	} else if (m_type == ETextureCubeMap) {
		glTexParameteri(m_glType, GL_TEXTURE_WRAP_S, wrap);
		glTexParameteri(m_glType, GL_TEXTURE_WRAP_T, wrap);
		glTexParameteri(m_glType, GL_TEXTURE_WRAP_R, wrap);
	}
}

void GLTexture::download(Bitmap *bitmap) {
	Assert(m_type == ETexture2D); // Only type supported so far

	if (bitmap == NULL)
		bitmap = getBitmap();

	Assert(bitmap != NULL);
	unsigned char *temp, *pixels = bitmap->getData();

	activateTarget();
	GLenum format = 0, dataFormat = GL_UNSIGNED_BYTE;
	switch (bitmap->getBitsPerPixel()) {
		case 8: 
			if (m_fbType == EDepthBuffer)
				format = GL_DEPTH_COMPONENT;
			else
				format = GL_LUMINANCE;
			break;
		case 24: format = GL_RGB; break;
		case 32: format = GL_RGBA; break;
		case 96: format = GL_RGB; dataFormat = GL_FLOAT; break;
		case 128: format = GL_RGBA; dataFormat = GL_FLOAT; break;
	}

	glReadPixels(0, 0, bitmap->getWidth(), bitmap->getHeight(), format, 
		dataFormat, pixels);
	releaseTarget();

	/* OpenGL has associates (0, 0) with the lower left position
	   and the resulting bitmap is thus vertically flipped. */
	int rowSize = bitmap->getWidth() * bitmap->getBitsPerPixel() / 8,
	    halfHeight= bitmap->getHeight()/2;
	temp = new unsigned char[rowSize];
	for (int i=0, j=bitmap->getHeight()-1; i<halfHeight; ++i) {
		memcpy(temp, pixels + i * rowSize, rowSize);
		memcpy(pixels + i * rowSize, pixels + j * rowSize, rowSize);
		memcpy(pixels + j * rowSize, temp, rowSize);
		j--;
	}
	delete[] temp;
}

Spectrum GLTexture::getPixel(int x, int y) const {
	Assert(m_fbType == EColorBuffer);
	float pixels[3];
	Spectrum result;

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fboId);
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, m_size.x, m_size.y);
	glReadPixels(x, y, 1, 1, GL_RGB, GL_FLOAT, &pixels);
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, GL_NONE);

	result.fromLinearRGB(pixels[0], pixels[1], pixels[2]);
	return result;
}


void GLTexture::activateTarget() {
	Assert(m_fbType != ENone);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fboId);
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, m_size.x, m_size.y);
}

void GLTexture::activateSide(int side) {
	if (side == -1) {
		if (m_fbType == EColorBuffer) {
			glFramebufferTextureEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, m_id, 0);
			glFramebufferTextureEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, m_depthId, 0);
		} else if (m_fbType == EDepthBuffer) {
			glFramebufferTextureEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, m_id, 0);
		} else {
			Log(EError, "Unsupported framebuffer type!");
		}
	} else {
		if (m_fbType == EColorBuffer) {
			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
					GL_TEXTURE_CUBE_MAP_POSITIVE_X + side, m_id, 0);
			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
					GL_TEXTURE_CUBE_MAP_POSITIVE_X + side, m_depthId, 0);
		} else if (m_fbType == EDepthBuffer) {
			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
					GL_TEXTURE_CUBE_MAP_POSITIVE_X + side, m_id, 0);
		} else {
			Log(EError, "Unsupported framebuffer type!");
		}
	}
}
	
void GLTexture::setTargetRegion(const Point2i &offset, const Vector2i &size) {
	glViewport(offset.x, offset.y, size.x, size.y);
}

void GLTexture::releaseTarget() {
	Assert(m_fbType != ENone);
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, GL_NONE);
	if (isMipMapped())
		m_needsUpdate = true;
}

void GLTexture::bind(int textureUnit, int textureIndex) const {
	/* Bind to the texture */
	if (GLEW_VERSION_1_3) {
		m_textureUnits.get().insert(textureUnit);
		glActiveTexture(GL_TEXTURE0 + textureUnit);
	} else {
		if (textureUnit != 0)
			Log(EError, "Multitexturing is not supported");
	}

	glEnable(m_glType);
	if (textureIndex == 1 && m_fbType == EColorAndDepthBuffer) {
		glBindTexture(m_glType, m_depthId);
	} else {
		glBindTexture(m_glType, m_id);
	}

	if (isMipMapped() && m_needsUpdate) {
		glGenerateMipmapEXT(m_glType);
		m_needsUpdate = false;
	}
}

void GLTexture::unbind() const {
	if (GLEW_VERSION_1_3) {
		std::set<int> &textureUnits = m_textureUnits.get();
		for (std::set<int>::iterator it = textureUnits.begin();
				it != textureUnits.end(); ++it) {
			glActiveTexture(GL_TEXTURE0 + *it);
			glDisable(m_glType);
		}
		textureUnits.clear();
	} else {
		glDisable(m_glType);
	}
}

void GLTexture::cleanup() {
	if (m_id == 0)
		return;
	if (m_fbType != ENone) {
		Log(ETrace, "Freeing framebuffer \"%s\"", m_name.c_str());
		if (m_fbType == EColorAndDepthBuffer || (m_fbType == EColorBuffer && m_type == ETextureCubeMap)) {
			glDeleteTextures(1, &m_depthId);
		} else if (m_fbType == EColorBuffer) {
			glDeleteRenderbuffersEXT(1, &m_depthId);
		}

		glDeleteFramebuffersEXT(1, &m_fboId);
	} else {
		Log(ETrace, "Freeing texture \"%s\"", m_name.c_str());
	}
	glDeleteTextures(1, &m_id);
	m_id = 0;
}

void GLTexture::blit(GPUTexture *target, int what) const {
	GLTexture *dest = static_cast<GLTexture *>(target);
	Assert(m_fbType != ENone && (dest == NULL || dest->m_fbType != ENone));

	if (!GLEW_EXT_framebuffer_blit)
		Log(EError, "Your OpenGL driver does not support fast "
			"framebuffer blitting!");

	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, m_fboId);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 
		(dest == NULL) ? GL_NONE : dest->m_fboId);

	int flags = 0;
	if (what & EDepthBuffer)
		flags |= GL_DEPTH_BUFFER_BIT;
	if (what & EColorBuffer)
		flags |= GL_COLOR_BUFFER_BIT;

	if (!dest) {
		glBlitFramebufferEXT(0, 0, m_size.x, m_size.y, 0, 0, 
			m_size.x, m_size.y, flags, GL_NEAREST);
	} else {
		glBlitFramebufferEXT(0, 0, m_size.x, m_size.y, 0, 0, 
			dest->m_size.x, dest->m_size.y, flags, 
			(m_size == dest->m_size) ? GL_NEAREST : GL_LINEAR);
	}
	
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, GL_NONE);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, GL_NONE);
}

void GLTexture::blit(GPUTexture *target, int what, 
		const Point2i &sourceOffset, const Vector2i &sourceSize,
		const Point2i &destOffset, const Vector2i &destSize) const {
	GLTexture *dest = static_cast<GLTexture *>(target);
	Assert(m_fbType != ENone && (dest == NULL || dest->m_fbType != ENone));

	if (!GLEW_EXT_framebuffer_blit)
		Log(EError, "Your OpenGL driver does not support fast "
			"framebuffer blitting!");

	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, m_fboId);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 
		(dest == NULL) ? GL_NONE : dest->m_fboId);

	int flags = 0;
	if (what & EDepthBuffer)
		flags |= GL_DEPTH_BUFFER_BIT;
	if (what & EColorBuffer)
		flags |= GL_COLOR_BUFFER_BIT;

	glBlitFramebufferEXT(sourceOffset.x, sourceOffset.y, 
		sourceOffset.x + sourceSize.x, sourceOffset.x + sourceSize.y,
		destOffset.x, destOffset.y, destOffset.x + destSize.x, 
		destOffset.y + destSize.y, flags, 
		(sourceSize == destSize) ? GL_NEAREST : GL_LINEAR);

	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, GL_NONE);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, GL_NONE);
}

void GLTexture::clear() {
	Assert(m_fbType != ENone);
	glClear(GL_DEPTH_BUFFER_BIT 
		| ((m_fbType & EColorBuffer) ? GL_COLOR_BUFFER_BIT : 0));
}
	
GLTexture::~GLTexture() {
	cleanup();
}

MTS_IMPLEMENT_CLASS(GLTexture, false, Texture)
MTS_NAMESPACE_END
