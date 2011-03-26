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

#if !defined(__GIZMO_H)
#define __GIZMO_H

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/**
 * This class is responsible for rendering `gizmos', which are
 * intuitive controls for navigating through 3D space.
 */
class MTS_EXPORT_HW Gizmo : public Object {
public:
	Gizmo();

	/// Initialize the gizmo for the specified bounding sphere
	void init(const BSphere &bsphere);

	/// Start a new drag operation
	void startDrag(const Ray &ray);

	//// Perform a drag to the specified ray
	void dragTo(const Ray &ray, const Camera *camera);

	/// Check whether the gizmo is currently active
	inline bool isActive() const { return m_active; }
	
	/// Check whether the gizmo is currently used in a drag operation
	inline bool isDragging() const { return m_drag; }

	/// Check whether a certain ray intersects the gizmo
	inline void rayIntersect(const Ray &ray, Float &t) const {
		Float farT;
		if (!m_active || !m_bsphere.rayIntersect(ray, t, farT)) {
			t = std::numeric_limits<Float>::infinity();
		} else {
			if (t <= 0)
				t = farT;
			if (t <= 0)
				t = std::numeric_limits<Float>::infinity();
		}
	}

	/// Draw the gizmo
	void draw(Renderer *renderer, const Camera *camera);

	/// Return the transformation associated with the gizmo
	Transform getTransform() const;

	/// Reset the gizmo
	inline void reset() { m_active = m_drag = false; }
	
	/// Stop dragging
	inline void stopDrag() { m_drag = false; }

	/// Return the bounding sphere associated with the gizmo
	inline const BSphere &getBSphere() const { return m_bsphere; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Gizmo();
private:
	BSphere m_bsphere;
	Point m_dragStart, m_dragEnd;
	bool m_active, m_drag;
};

MTS_NAMESPACE_END

#endif /* __GIZMO_H */
