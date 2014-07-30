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
#if !defined(__MITSUBA_RENDER_SUBSURFACE_H_)
#define __MITSUBA_RENDER_SUBSURFACE_H_

#include <mitsuba/core/netobject.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract subsurface scattering models
 *
 * Can be attached to an arbitrary shape to compute exitant
 * radiance due to internal scattering. How that is done is
 * completely up to the implementation. It might for instance
 * recursively trace rays or perform lookups into a precomputed
 * point cloud radiance representation.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Subsurface : public NetworkedObject {
public:
	/**
	 * \brief Possibly perform a pre-process task.
	 *
	 * The last three parameters are resource IDs of the associated scene,
	 * camera and sample generator, which have been made available to all
	 * local and remote workers.
	 */
	virtual bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) = 0;

	/// Selectively activate/deactivate the subsurface integrator
	inline void setActive(bool active) { m_active = active; }

	/// Return whether or not the subsurface integrator is currently active
	inline bool isActive() const { return m_active; }

	/// Cancel any running pre-process tasks
	virtual void cancel();

	/// Return the list of shapes associated with this subsurface integrator
	inline const std::vector<Shape *> getShapes() const { return m_shapes; }

	/// Get the exitant radiance for a point on the surface
	virtual Spectrum Lo(const Scene *scene, Sampler *sampler,
		const Intersection &its, const Vector &d, int depth = 0) const = 0;

	/// Serialize this subsurface integrator to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Set the parent node of the subsurface integrator
	void setParent(ConfigurableObject *parent);

	MTS_DECLARE_CLASS()
protected:
	/// Create a new subsurface scattering class
	Subsurface(const Properties &props);

	/// Unserialize a subsurface integrator from a binary data stream
	Subsurface(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Subsurface();
protected:
	std::vector<Shape *> m_shapes;
	bool m_active;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_SUBSURFACE_H_ */
