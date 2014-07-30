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
#if !defined(__MITSUBA_BIDIR_MANIFOLD_H_)
#define __MITSUBA_BIDIR_MANIFOLD_H_

#include <mitsuba/bidir/vertex.h>

#define MTS_MANIFOLD_DEBUG          0
#define MTS_MANIFOLD_EPSILON        Epsilon
#define MTS_MANIFOLD_MAX_ITERATIONS 20

MTS_NAMESPACE_BEGIN

/**
 * \brief Utility class for perturbing paths located on a specular manifold.
 * \author Wenzel Jakob
 */
class MTS_EXPORT_BIDIR SpecularManifold : public Object {
public:
	/// Construct an unitialized specular manifold data structure
	SpecularManifold(const Scene *scene, int maxIterations = -1);

	/**
	 * \brief Initialize the specular manifold with the specified
	 * path segment
	 */
	bool init(const Path &path, int start, int end);

	/**
	 * \brief Update the provided path segment based on the stored
	 * specular manifold configuration
	 */
	bool update(Path &path, int start, int end);

	/// Attempt to move the movable endpoint vertex to position \c target
	bool move(const Point &target, const Normal &normal);

	/**
	 * \brief Compute the generalized geometric term between 'a' and 'b'
	 */
	Float G(const Path &path, int a, int b);

	/**
	 * \brief Compute a product of standard and generalized geometric
	 * terms between 'a' and 'b' depending on whether vertices are
	 * specular or non-specular.
	 */
	Float multiG(const Path &path, int a, int b);

	Float det(const Path &path, int a, int b, int c);

	/// Return the number of iterations used by \ref move()
	inline int getIterationCount() const { return m_iterations; }

	/// Return the position of a vertex
	inline const Point &getPosition(int i) { return m_vertices[i].p; }

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~SpecularManifold() { }

private:
	enum EType {
		EPinnedPosition = 0,
		EPinnedDirection,
		EReflection,
		ERefraction,
		EMedium,
		EMovable
	};

	/// Describes a single interaction on the path
	struct SimpleVertex {
		bool degenerate : 1;
		EType type : 31;

		/* Position and partials */
		Point p;
		Vector dpdu, dpdv;

		/* Normal and partials */
		Normal n, gn;
		Vector dndu, dndv;

		/* "Microfacet" normal */
		Normal m;

		/* Further information about the vertex */
		Float eta;
		const Object *object;

		/* Scratch space for matrix assembly */
		Matrix2x2 a, b, c, u;

		/* Manifold tangent space projected onto this vertex */
		Matrix2x2 Tp;

		/// Initialize certain fields to zero by default
		inline SimpleVertex(EType type, const Point &p) :
			degenerate(false), type(type), p(p), dpdu(0.0f),
			dpdv(0.0f), n(0.0f), dndu(0.0f), dndv(0.0f),
			m(0.0f), eta(1.0f), object(NULL) { }

		/// Map a tangent space displacement into world space
		inline Vector map(Float u, Float v) const {
			Vector2 T = Tp * Vector2(u, v);
			return T.x * dpdu + T.y * dpdv;
		}

		std::string toString() const;
	};

	/// Take the specified step and project back onto the manifold
	bool project(const Vector &d);

	/**
	 * \brief Compute the tangent vectors with the specified
	 * components when projected onto the movable endpoint vertex
	 */
	bool computeTangents();

	void check(SimpleVertex *v);
protected:
	const Scene *m_scene;
	Float m_time;
	int m_iterations, m_maxIterations;

	std::vector<SimpleVertex> m_vertices, m_proposal;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_MANIFOLD_H_ */
