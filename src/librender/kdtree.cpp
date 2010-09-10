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

#include <mitsuba/render/kdtree.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/triaccel.h>

MTS_NAMESPACE_BEGIN

KDTree::KDTree() : m_built(false), m_triangles(NULL)  {
	/* SAH constants used by Wald and Havran */
	m_emptyBonus = 0.8f;
	m_traversalCost = 15;
	m_intersectionCost = 20;
	m_clip = true;
	m_maxBadRefines = 3;
	m_triangleCount = 0;
	m_primitiveCount = 0;
	m_nonTriangleCount = 0;
#if defined(MTS_USE_TRIACCEL4)
	m_packedTriangles = NULL;
	m_stopPrims = 8;
#else
	m_stopPrims = 4;
#endif
}

void KDTree::addShape(const Shape *shape) {
	Assert(!m_built);
	if (shape->isCompound())
		Log(EError, "Cannot add compound shapes to a KD-tree - expand them first!");
	if (shape->getClass()->derivesFrom(TriMesh::m_theClass))
		m_triangleCount += (int) static_cast<const TriMesh *>(shape)->getTriangleCount();
	else
		m_nonTriangleCount += 1;
	m_shapes.push_back(shape);
}

KDTree::~KDTree() {
	if (m_triangles != NULL)
		freeAligned(m_triangles);
#if defined(MTS_USE_TRIACCEL4)
	if (m_packedTriangles != NULL)
		freeAligned(m_packedTriangles);
#endif
}

MTS_IMPLEMENT_CLASS(KDTree, false, Object)
MTS_NAMESPACE_END
