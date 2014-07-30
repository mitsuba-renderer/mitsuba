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
#if !defined(__MITSUBA_RENDER_SHADER_H_)
#define __MITSUBA_RENDER_SHADER_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class Shader;
class Renderer;
class GPUProgram;

/**
 * \brief Abstract hardware resource.
 *
 * Implementations provides support for functionality that are also able to
 * live on the GPU. By default, the method 'createShader' just returns \c NULL,
 * which means that the BSDF/Light source/Texture/.. has not yet been ported
 * to the GPU-based renderer.
 */
class MTS_EXPORT_RENDER HWResource {
public:
	virtual Shader *createShader(Renderer *renderer) const;

	virtual ~HWResource() { }
};

/**
 * \brief %Shader base class for use with a VPL-style renderer.
 *
 * Subclasses can implement one of various things, such as a BSDF,
 * a light source, or a texture.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Shader : public Object {
public:
	enum EShaderType {
		EBSDFShader = 0,
		ETextureShader,
		EEmitterShader
	};

	enum EFlags {
		ETransparent = 0x01
	};

	// Return the type of shader represented by this instance
	inline EShaderType getType() const { return m_type; }

	/// Return a list of flags
	inline uint32_t getFlags() const { return m_flags; }

	/**
	 * \brief For transparent objects, this function returns
	 * the alpha blending weight
	 */
	virtual Float getAlpha() const;

	// List other shaders, on which this instance depends
	virtual void putDependencies(std::vector<Shader *> &deps);

	/**
	 * \brief Is this shader complete?
	 *
	 * This is mainly useful to check whether all dependencies
	 * could be constructed successfully. The default
	 * implementation returns \c true.
	 */
	virtual bool isComplete() const;

	/**
	 * \brief Generate a string version of this shader's evaluation
	 * routine.
	 *
	 * The appended string should assign the name \c evalName to this
	 * function. The function names of depedencies (as specified by
	 * \ref putDependencies), are supplied in the parameter \c depNames
	 * in identical order.
	 */
	virtual void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const = 0;

	/**
	 * \brief This function can optionally be implemented to resolve named
	 * program parameters to numerical IDs for increased performance.
	 *
	 * The int array returned here will later be passed to \ref bind().
	 * The default implementation does nothing.
	 */
	virtual void resolve(const GPUProgram *program, const std::string &evalName,
		std::vector<int> &parameterIDs) const;

	/**
	 * \brief Configure the the associated GPU program.
	 *
	 * This function is typically used to bind textures and
	 * to set program pararameters.
	 */
	virtual void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const;

	/// Release any bound resources.
	virtual void unbind() const;

	/// Release all resources
	virtual void cleanup(Renderer *renderer);

	MTS_DECLARE_CLASS()
protected:
	Shader(Renderer *renderer, EShaderType type);

	/// Virtual destructor
	virtual ~Shader();
protected:
	EShaderType m_type;
	uint32_t m_flags;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_SHADER_H_ */
