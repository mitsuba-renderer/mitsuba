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

#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/emitter.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief "Fake" shape that groups sub-shapes into a
 * separate KD-tree.
 *
 * This shape doesn't actually generate any intersectable
 * geometry on its own. Instead, the "instance" plugin must
 * be used to create references to the geometry stored inside it.
 */
class ShapeGroup : public Shape {
public:
	/// Create a new shape group
	ShapeGroup(const Properties &props);

	/// Unserialize from a binary data stream
	ShapeGroup(Stream *stream, InstanceManager *manager);

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Build the internal KD-tree
	void configure();

	/// Add a child object
	void addChild(const std::string &name, ConfigurableObject *child);
	/// Add an unnamed child
	inline void addChild(ConfigurableObject *child) { addChild("", child); }

	/// Return whether or not the shape is a compound object
	bool isCompound() const;

	/// Returns an invalid AABB
	AABB getAABB() const;

	/// Returns the surface area
	Float getSurfaceArea() const;

	/// Return a pointer to the internal KD-tree
	inline const ShapeKDTree *getKDTree() const { return m_kdtree.get(); }

	/// Return the primitive count of the nested shapes
	size_t getPrimitiveCount() const;

	/// Return the effective primitive count of this shape (always zero)
	size_t getEffectivePrimitiveCount() const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
private:
	ref<ShapeKDTree> m_kdtree;
};

MTS_NAMESPACE_END
