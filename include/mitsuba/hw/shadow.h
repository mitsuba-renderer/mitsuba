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
#if !defined(__MITSUBA_HW_SHADOW_H_)
#define __MITSUBA_HW_SHADOW_H_

#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/** \brief Utility class for creating different kinds of shadow maps (cube,
 * directional, and paraboloid shadow maps) using hardware rasterization
 *
 * \ingroup libhw
 */
class MTS_EXPORT_HW ShadowMapGenerator : public Object {
public:
	enum EShadowMapType {
		/// Directional (orthographic) shadow map
		EDirectional = 0,

		/// Nonlinear paraboloid shadow map (captures one hemisphere)
		EParaboloid,

		/// Omnidirectional shadow map, 6 passes
		ECube,

		/// Omnidirectional shadow map, 1 pass with geometry shader
		ECubeSinglePass,

		/// Hemispherical cube shadow map, 5 passes
		EHemicube,

		/// Hemispherical cube shadow map, 1 pass with geometry shader
		EHemicubeSinglePass,

		/// Unused
		ETypeCount
	};

	/// Create a new shadow map generator
	ShadowMapGenerator(Renderer *renderer);

	/// Register the associated GPU programs
	void init();

	/**
	 * \brief Allocate a texture that is suitable for storing
	 * the requested type of shadow map and return it
	 */
	ref<GPUTexture> allocate(Renderer *renderer,
		EShadowMapType type, int resolution);

 	/**
	 * \brief Render a shadow map using the desired technique
	 *
	 * \param renderer
	 *    Underlying \ref Renderer instance
	 * \param shadowMap
	 *    Target texture, which was previously created using \ref allocate()
	 * \param type
	 *    Desired shadow map type
	 * \param trafo
	 *    View transformation of the source
	 */
	void render(Renderer *renderer, GPUTexture *shadowMap, EShadowMapType type,
		const Transform &trafo, Float minDepth, Float maxDepth,
		const std::vector<Renderer::TransformedGPUGeometry> &geo);

	/**
	 * \brief Convenience function for computing a transformation that
	 * is suitable for creating a directional shadow map
	 *
	 * \param aabb
	 *     Axis-aligned bounding box of the world-space content that
	 *     should be covered
	 *
	 * \param d
	 *     Viewing direction of the shadow map
	 */
	Transform directionalFindGoodFrame(const AABB &aabb, const Vector &d) const;

	/// Release allocated resources
	void cleanup();

	/// Return the number of resident shaders
	size_t getShaderCount() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~ShadowMapGenerator() { }
protected:
	ref<GPUProgram> m_program[ETypeCount];
	bool m_cubeDepthMapsSupported;

	/* Cached parameter indices (ECube) */
	int m_cubeTransform;
	int m_cubeProjDir;

	/* Cached parameter indices (ECubeSinglePass) */
	int m_cubeSinglePassTransform[6];
	int m_cubeSinglePassProjDir[6];

	/* Cached parameter indices (EHemicubeSinglePass) */
	int m_hemicubeSinglePassTransform[5];
	int m_hemicubeSinglePassProjDir[5];

	/* Cached parameter indices (EParaboloid) */
	int m_paraboloidMinDepth;
	int m_paraboloidInvDepthRange;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_SHADOW_H_ */
