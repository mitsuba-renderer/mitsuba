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

#if !defined(__SHADER_H)
#define __SHADER_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class Shader;
class Renderer;
class GPUProgram;

/**
 * Abstract hardware resource -- provides support for functionality,
 * which can alternatively also live on the GPU. By default, the
 * method 'createShader' just returns NULL, which means that the
 * BSDF/Light source/Texture/.. cannot be used with a GPU-based renderer.
 */
class MTS_EXPORT_RENDER HWResource {
public:
	virtual Shader *createShader(Renderer *renderer) const;
};

/**
 * Shader base class for use with a VPL-style renderer. Could implement
 * one of various things, such as a BSDF, a light source, or a texture.
 */
class MTS_EXPORT_RENDER Shader : public Object {
public:
	enum EShaderType {
		EBSDFShader = 0,
		ETextureShader,
		ELuminaireShader
	};

	/**
	 * Return the type of shader represented by this instance
	 */
	inline EShaderType getType() const { return m_type; }

	/**
	 * List other shaders, on which this instance depends
	 */
	virtual void putDependencies(std::vector<Shader *> &deps);

	/**
	 * Is this shader complete? This is mainly useful to
	 * check whether all dependencies could be constructed
	 * successfully. The default implementation returns true
	 */
	virtual bool isComplete() const;

	/**
	 * Generate a string version of this shader's evaluation
	 * routine. The appended string should assign the name
	 * 'evalName' to this function. The function names of 
	 * depedencies (as specified by 'putDependencies'), are
	 * supplied in the parameter 'depNames' in identical order.
	 */
	virtual void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const = 0;

	/**
	 * This function can optionally be implemented to resolve named 
	 * program parameters to numerical IDs for increased performance.
	 * The int array returned here will later be passed to bind().
	 * The default implementation does nothing.
	 */
	virtual void resolve(const GPUProgram *program, const std::string &evalName,
		std::vector<int> &parameterIDs) const;

	/**
	 * Configure the the associated GPU program. This
	 * function is typically used to bind textures and
	 * to set program pararameters.
	 */
	virtual void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const;

	/**
	 * Release any bound resources.
	 */
	virtual void unbind() const;

	/**
	 * Release all resources
	 */
	virtual void cleanup(Renderer *renderer);

	MTS_DECLARE_CLASS()
protected:
	Shader(Renderer *renderer, EShaderType type);

	/// Virtual destructor
	virtual ~Shader();
private:
	EShaderType m_type;
};

MTS_NAMESPACE_END

#endif /* __SHADER_H */
