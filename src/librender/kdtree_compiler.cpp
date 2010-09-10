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
#include <mitsuba/render/triaccel.h>
#include <mitsuba/core/timer.h>

MTS_NAMESPACE_BEGIN

void KDTree::build() {
	AssertEx(!m_built, "The kd-tree has already been built!");
	AssertEx(sizeof(KDNode) == 8, "The kd-tree node size should be 8 bytes!");

#if defined(MTS_USE_TRIACCEL)
	if (sizeof(TriAccel) % 16 != 0) 
		Log(EWarn, "The TriAccel data structure won't align properly.");
#endif

#if defined(MTS_USE_TRIACCEL4)
	if (sizeof(TriAccel4) % 16 != 0) 
		Log(EWarn, "The TriAccel4 data structure won't align properly.");
	m_packedTriangleCount = 0;
#endif

	m_primitiveCount = m_triangleCount + m_nonTriangleCount;
	if (m_primitiveCount == 0) {
		Log(EWarn, "The kd-tree does not contain any geometry!");
		m_nodes.clear();
		m_nodes.push_back(KDNode());
		m_nodes.push_back(KDNode());
		m_nodes[1].setLeaf(0,0);
		m_triangles = static_cast<KDTriangle *>(allocAligned(1));

#if defined(MTS_USE_TRIACCEL4)
		m_packedTriangles = static_cast<TriAccel4 *>(allocAligned(1));
#endif
		return;
	}

	/* Allocate space for triangle data and fill it with 
	   preliminary indexing information */
	m_triangles = static_cast<KDTriangle *>(allocAligned(
		sizeof(KDTriangle) * (m_primitiveCount)));
	unsigned int triIndex = 0;
	for (unsigned int i = 0; i < m_shapes.size(); ++i) {
		const Shape *shape = m_shapes[i];
		if (shape->getClass()->derivesFrom(TriMesh::m_theClass)) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			for (unsigned int j=0; j<mesh->getTriangleCount(); ++j) {
				KDTriangle &ta = m_triangles[triIndex++];
				ta.index = j;
				ta.shapeIndex = i;
				ta.k = 0;
			}
		} else {
			KDTriangle &ta = m_triangles[triIndex++];
			ta.index = KNoTriangleFlag;
			ta.shapeIndex = i;
			ta.k = 0;
		}
	}

	/* Initialize the statistics conters */
	ref<Timer> timer = new Timer();
	m_innerNodes = m_leafNodes = m_leafPrims = 0;
	m_totalLeafDepth = m_maxLeafPrims = m_failures = 0;
	m_badRefines = m_actualMaxDepth = m_clippedAway = 0;
	m_minLeafPrims = INT_MAX;
	m_bucketCount = 2*m_stopPrims;
	m_triBuckets = new int[m_bucketCount];
	for (int i=0; i<m_bucketCount; i++)
		m_triBuckets[i] = 0;

	/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
	m_maxDepth = std::min((int) (8 + 1.3f * log2i(m_primitiveCount)),
		MTS_KD_MAXDEPTH);

	/* Allocate space for AABB edge events. The construction uses a
	   normal and a "right" edge event list during the recursion. The
	   reason for this is that it is usually not possible to cleanly 
	   split a list of primitives. To accommodate the growth caused by
	   straddeling primitives, a per-depth right recursion list is 
	   required. The following fragment tries to conservatively
	   reserve memory in advance. */
	long long eventStorage = 0;
	m_rightEvents = new EdgeEventVec3[m_maxDepth];
	for (int axis=0; axis<3; axis++) {
		Float estimatedEvents = 2.0f * (Float) m_primitiveCount;
 		for (int i=0; i<m_maxDepth; i++) {
			estimatedEvents *= 0.6f;
			int nEvents = (int) estimatedEvents;
			if (nEvents > 0) {
				m_rightEvents[i][axis].reserve(nEvents);
				eventStorage += nEvents;
			}
		}
		m_events[axis].reserve(2 * m_primitiveCount);
		eventStorage += 2 * m_primitiveCount;
	}

	Log(EDebug, "Allocated %.2f MB of memory for kd-tree construction", 
		eventStorage * sizeof(EdgeEvent) / (Float) (1024 * 1024));

	Log(EDebug, "kd-tree configuration:");
	Log(EDebug, "   Primitive clipping     : %s", m_clip ? "on" : "off");
	Log(EDebug, "   Traversal cost         : %.2f", m_traversalCost);
	Log(EDebug, "   Intersection cost      : %.2f", m_intersectionCost);
	Log(EDebug, "   Max. tree depth        : %i", m_maxDepth);
	Log(EDebug, "   Min. leaf prims        : %i", m_stopPrims);
	Log(EDebug, "   Max. bad refines       : %i", m_maxBadRefines);
	Log(EDebug, "   Empty space bonus      : %.2f", m_emptyBonus);
#if defined(MTS_USE_TRIACCEL4)
	Log(EDebug, "   Intersection algorithm : TriAccel4 (SSE)");
#elif defined(MTS_USE_TRIACCEL)
	Log(EDebug, "   Intersection algorithm : TriAccel (scalar)");
#elif defined(MTS_USE_MT)
	Log(EDebug, "   Intersection algorithm : Moeller-Trumbore");
#else
	#error No intersection algorithm has been selected!
#endif
#if defined(MTS_HAS_COHERENT_RT)
	Log(EDebug, "   Coherent ray tracing   : enabled");
#endif

	Log(EDebug, "Constructing a SAH kd-tree (%i primitives) ..", m_primitiveCount);

	/* Generate edge events and calculate the root bounds */
	for (unsigned int i=0; i<m_primitiveCount; ++i) {
		const KDTriangle &ta = m_triangles[i];
		AABB aabb;

		if (ta.index == KNoTriangleFlag) {
			aabb = m_shapes[ta.shapeIndex]->getAABB();
		} else {
			const TriMesh *mesh = static_cast<const TriMesh *>(m_shapes[ta.shapeIndex]);
			const Triangle &tri = mesh->getTriangles()[ta.index];
			aabb = tri.getAABB(mesh->getVertexBuffer());
		}

		m_rootBounds.expandBy(aabb);

		for (int axis=0; axis<3; axis++) {
			if (aabb.max[axis] == aabb.min[axis]) {
				/* The triangle is planar wrt. the current axis */
				m_events[axis].push_back(EdgeEvent(EdgeEvent::EEdgePlanar,
					aabb.min[axis], i));
			} else {
				/* The triangle covers a non-zero interval. Insert
					start+end events */
				m_events[axis].push_back(EdgeEvent(EdgeEvent::EEdgeStart,
					aabb.min[axis], i));
				m_events[axis].push_back(EdgeEvent(EdgeEvent::EEdgeEnd,
					aabb.max[axis], i));
			}
		}
	}
			
	/* Generate a bounding sphere */
	m_bsphere.center = m_rootBounds.getCenter();
	for (unsigned int i=0; i<m_primitiveCount; ++i) {
		const KDTriangle &ta = m_triangles[i];
		const Shape *shape = m_shapes[ta.shapeIndex];
		if (ta.index == KNoTriangleFlag) {
			AABB aabb = shape->getAABB();
			for (int i=0; i<8; ++i) 
				m_bsphere.expandBy(aabb.getCorner(i));
		} else {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			const Vertex *vb = mesh->getVertexBuffer();
			const Triangle &tri = mesh->getTriangles()[ta.index];
			for (int j=0; j<3; j++)
				m_bsphere.expandBy(vb[tri.idx[j]].v);
		}
	}

	/* Sort events for each axis */
	for (int axis=0; axis<3; axis++) {
		std::sort(m_events[axis].begin(), m_events[axis].end(), 
			EdgeEventSorter());
	}

	/* Allocate the first node and a (dummy) alignment node */
	m_nodes.clear();
	m_nodes.push_back(KDNode());
	m_nodes.push_back(KDNode());

	/* Slightly enlarge */
	m_rootBounds.max += (m_rootBounds.max-m_rootBounds.min) * 0.01f + Vector(Epsilon, Epsilon, Epsilon);
	m_rootBounds.min -= (m_rootBounds.max-m_rootBounds.min) * 0.01f + Vector(Epsilon, Epsilon, Epsilon);

	/* Recursively build the tree */
	buildTree(1, 0, 0, m_primitiveCount, m_rootBounds, m_events);

#ifdef MTS_USE_TRIACCEL
	Log(EDebug, "Pre-computing triangle intersection information ..");
	/* Pre-compute constants for triangle intersection */
		for (unsigned int i=0; i<m_primitiveCount; ++i) {
		KDTriangle &ta = m_triangles[i];
		if (ta.index == KNoTriangleFlag)
			continue;
		const TriMesh *mesh = static_cast<const TriMesh *>(m_shapes[ta.shapeIndex]);
		const Vertex *vertices = mesh->getVertexBuffer();
		const Triangle &tri = mesh->getTriangles()[ta.index];
		Point A = vertices[tri.idx[0]].v,
			  B = vertices[tri.idx[1]].v,
			  C = vertices[tri.idx[2]].v;
		m_failures += ta.load(A, B, C);
	}
#endif
	if (m_failures > 0)
		Log(EWarn, "Encountered %i malformed primitives!", m_failures);

	Log(EDebug, "Building time : %i ms", timer->getMilliseconds());

	/* Calculate the used amount of temporary memory */
	eventStorage = m_eventsL.capacity() + m_eventsR.capacity()
		+ m_sEventsL.capacity() + m_sEventsR.capacity();
	for (int axis=0; axis<3; axis++) {
		eventStorage += m_events[axis].capacity();
		for (int i=0; i<m_maxDepth; i++)
			eventStorage += m_rightEvents[i][axis].capacity();
	}

	Log(EDebug, "Freeing %.2f MB of temporary memory", 
		eventStorage * sizeof(EdgeEvent) / (Float) (1024 * 1024));

	Log(EDebug, "Detailed kd-tree statistics:");
	Log(EDebug, "   Triangle primitives    : %i", m_triangleCount);
	Log(EDebug, "   Other primitives       : %i", m_nonTriangleCount);
	Log(EDebug, "   Inner nodes            : %i", m_innerNodes);
	Log(EDebug, "   Leaf nodes             : %i", m_leafNodes);
	Log(EDebug, "   Total leaf primitives  : %i", m_leafPrims);
	if (m_clip)
		Log(EDebug, "   Pruned primitives      : %i", m_clippedAway);
	Log(EDebug, "   Min leaf primitives    : %i", m_minLeafPrims);
	Log(EDebug, "   Max leaf primitives    : %i", m_maxLeafPrims);
	Log(EDebug, "   Avg. leaf primitives   : %.2f", 
		(Float) m_leafPrims / (Float) m_leafNodes);
	Log(EDebug, "   Total bad refines      : %i", m_badRefines);
	Log(EDebug, "   Max. leaf depth        : %i", m_actualMaxDepth);
	Log(EDebug, "   Avg. leaf depth        : %.2f", 
		(Float) m_totalLeafDepth / (Float) m_leafNodes);
	size_t nBytes = sizeof(KDNode) * m_nodes.size();

#if defined(MTS_USE_TRIACCEL4)
	nBytes += sizeof(TriAccel4) * m_packedTriangleCount;
#endif

#if defined(MTS_USE_TRIACCEL) || defined(MTS_USE_MT)
	nBytes += sizeof(int) * m_indices.size()
			+ sizeof(KDTriangle) * m_primitiveCount;
#endif

	Log(EDebug, "   Total memory usage     : %.2f MB", nBytes / (Float) (1024*1024));
	std::ostringstream oss;
	oss << "   Leaf node histogram    : ";
	for (int i=0; i<m_bucketCount; i++) {
		oss << i << "(" << m_triBuckets[i] << ") ";
		if ((i+1)%4==0 && i+1<m_bucketCount) {
			Log(EDebug, "%s", oss.str().c_str());
			oss.str("");
			oss << "                            ";
		}
	}
	delete[] m_triBuckets;
	Log(EDebug, "%s", oss.str().c_str());

	/* Free up memory */
	for (int axis=0; axis<3; axis++) 
		EdgeEventVec().swap(m_events[axis]);
	EdgeEventVec().swap(m_eventsL);
	EdgeEventVec().swap(m_eventsR);
	EdgeEventVec().swap(m_sEventsL);
	EdgeEventVec().swap(m_sEventsR);
	delete[] m_rightEvents;
#if defined(MTS_USE_TRIACCEL4)
	#if !defined(MTS_USE_TRIACCEL) && !defined(MTS_USE_MT)
	/* Indirection not needed anymore, since everything is stored explicitly */
	freeAligned(m_triangles);
	m_triangles = NULL;
	#endif

	/* Re-allocate using aligned memory and free unused memory */
	m_packedTriangles = static_cast<TriAccel4 *>(allocAligned(
		sizeof(TriAccel4) * m_packedTriangleCount));
	for (unsigned int i=0; i<m_packedTriangleCount; i++)
		m_packedTriangles[i] = m_tempPackedTriangles[i];
	std::vector<TriAccel4>().swap(m_tempPackedTriangles);
#endif
	Vector diff = m_rootBounds.max - m_rootBounds.min;
	m_rootBounds.min -= diff * 0.001f;
	m_rootBounds.max += diff * 0.001f;
	m_built = true;
}

void KDTree::buildTree(int nodeIndex, int depth, int badRefines,
		int numPrims, const AABB &aabb, EdgeEventVec3 &allEvents) {
	Float invSA = 1.0f / aabb.getSurfaceArea();
	Float nodeCost = m_intersectionCost * numPrims;
	Score bestSplit;
	
	if (depth >= m_maxDepth || numPrims <= m_stopPrims) {
		createLeaf(nodeIndex, depth, numPrims, allEvents);
		return;
	}

	/* ==================================================================== */
    /*                        Split candidate search                        */
    /* ==================================================================== */

	/* First, find the optimal splitting plane according to the
	   surface area heuristic. To do this in O(n), the search is
	   implemented as a sweep over the edge events */
	for (int axis=0; axis<3; axis++) {
		/* Initially, the split plane is placed left of the scene
		   and thus all geometry is on its right side */
		int numLeft = 0, numPlanar = 0, numRight = numPrims;
		const EdgeEventVec &events = allEvents[axis];
		EdgeEventVec::const_iterator it = events.begin();

		/* Iterate over all events on the current axis */
		while (it != events.end()) {
			/* Record the current position and count all 
			   other events, which are also here */
			Float t = (*it).t;
			int numStart = 0, numEnd = 0;

			/* Count "end" events */
			while (it != events.end() && (*it).t == t 
				&& (*it).type == EdgeEvent::EEdgeEnd) {
				++numEnd; ++it;
			}

			/* Count "planar" events */
			while (it != events.end() && (*it).t == t 
				&& (*it).type == EdgeEvent::EEdgePlanar) {
				++numPlanar; ++it;
			}

			/* Count "start" events */
			while (it != events.end() && (*it).t == t 
				&& (*it).type == EdgeEvent::EEdgeStart) {
				++numStart; ++it;
			}

			/* The split plane is now moved onto 't'. Therefore,
			   all planar primitives and ending primitives are
			   removed from the right side */
			numRight -= numPlanar; numRight -= numEnd;

			/* Calculate a score using the surface area heuristic */
			if (t >= aabb.min[axis] && t <= aabb.max[axis]) {
				Score score = SAH(axis, invSA, aabb, t, numLeft, 	
					numRight, numPlanar);
				if (score < bestSplit) {
					bestSplit = score;
				}
			} else {
				/* Clipped AABBs cannot be outside of the voxel except
				   if this is a non-triangle primitive */
				AssertEx(m_clip == false || 
					m_triangles[(*it).index].index == KNoTriangleFlag, "Internal error");
			}

			/* The split plane is moved past 't'. All prims,
			   which were planar on 't', are moved to the left
			   side. Also, starting prims are now also left of
			   the split plane. */
			numLeft += numStart; numLeft += numPlanar;
			numPlanar = 0;
		}
		/* Sanity checks. Everything should now be on the 
		   left side of the split plane */
		Assert(numRight == 0 && numLeft == numPrims);
	}

	/* "Bad refines" heuristic from PBRT */
	if (bestSplit.score > nodeCost) {
		if ((bestSplit.score > 4 * nodeCost && numPrims < 16)
				|| badRefines >= m_maxBadRefines
				|| bestSplit.score == std::numeric_limits<Float>::infinity()) {
			createLeaf(nodeIndex, depth, numPrims, allEvents);
			return;
		}
		++badRefines; ++m_badRefines;
	}

	AABB leftAABB = aabb, rightAABB = aabb;
	leftAABB.max[bestSplit.axis] = bestSplit.t;
	rightAABB.min[bestSplit.axis] = bestSplit.t;

	/* ==================================================================== */
    /*                      Primitive Classification                        */
    /* ==================================================================== */

	/* Initially mark all prims as being located on both sides */
	EdgeEventVec &events = allEvents[bestSplit.axis];
	for (EdgeEventVec::const_iterator it = events.begin(); 
		it != events.end(); ++it)
		m_triangles[(*it).index].k = EBothSides;

	int primsLeft = 0, primsRight = 0, primsBoth = numPrims;
	/* Sweep over all edge events and classify the primitives wrt. the split */
	for (EdgeEventVec::const_iterator it = events.begin(); 
		it != events.end(); ++it) {
		const EdgeEvent &event = *it;

		if (event.type == EdgeEvent::EEdgeEnd && event.t <= bestSplit.t) {
			/* The primitive's interval ends before or on the split plane
			   -> Completely move to the left side */
			Assert(m_triangles[event.index].k == EBothSides);
			m_triangles[event.index].k = ELeftSide;
			primsBoth--;
			primsLeft++;
		} else if (event.type == EdgeEvent::EEdgeStart 
			&& event.t >= bestSplit.t) {
			/* The primitive's interval starts after or on the split plane
			   -> Completely move to the right side */
			Assert(m_triangles[event.index].k == EBothSides);
			m_triangles[event.index].k = ERightSide;
			primsBoth--;
			primsRight++;
		} else if (event.type == EdgeEvent::EEdgePlanar) {
			/* If the planar primitive is not on the split plane, the
			   classification is easy. Otherwise, place it on the side with
			   the better SAH score */
			Assert(m_triangles[event.index].k == EBothSides);
			if (event.t < bestSplit.t || (event.t == bestSplit.t && 
				bestSplit.planarSide == Score::EPlanarLeft)) {
				m_triangles[event.index].k = ELeftSide;
				primsBoth--;
				primsLeft++;
			} else if (event.t > bestSplit.t || (event.t == bestSplit.t && 
				bestSplit.planarSide == Score::EPlanarRight)) {
				m_triangles[event.index].k = ERightSide;
				primsBoth--;
				primsRight++;
			} else {
				AssertEx(false, "Internal error!");
			}
		}
	}

	/* Some sanity checks */
	Assert(primsLeft + primsRight + primsBoth == numPrims);
	Assert(primsLeft + primsBoth == bestSplit.nLeft);
	Assert(primsRight + primsBoth == bestSplit.nRight);

	/* ==================================================================== */
    /*                              Splicing                                */
    /* ==================================================================== */

	int clippedL = 0, clippedR = 0;

	for (int axis=0; axis<3; axis++) {
		EdgeEventVec &events = allEvents[axis];
		
		m_eventsL.clear();
		m_eventsR.clear();
		m_sEventsL.clear();
		m_sEventsR.clear();

		for (EdgeEventVec::const_iterator it = events.begin(); 
			it != events.end(); ++it) {
			const EdgeEvent &event = *it;

			if (m_triangles[event.index].k == ELeftSide) {
				/* Left-only primitive. Move to the left list and advance */
				m_eventsL.push_back(event);
			} else if (m_triangles[event.index].k == ERightSide) {
				/* Right-only primitive. Move to the right list and advance */
				m_eventsR.push_back(event);
			} else if (m_triangles[event.index].k == EBothSides) {
				/* The primitive overlaps the split plane. Its edge events
				   must be added to both lists. */
				if (!m_clip) {
					m_eventsL.push_back(event);
					m_eventsR.push_back(event);
					continue;
				}

				/* This primitive has now been processed */
				m_triangles[event.index].k = EBothSidesProcessed;

				const KDTriangle &kdTri = m_triangles[event.index];
				const Shape *shape = m_shapes[kdTri.shapeIndex];
				AABB primAABBLeft, primAABBRight;

				/* If primitive clipping is enabled, calculate special AABBs for each
				   of the child voxels. This can make a big difference in terms 
				   of overall kd-tree quality */
				if (kdTri.index != KNoTriangleFlag) {
					const TriMesh *mesh = static_cast<const TriMesh *>(shape);
					const Triangle &tri = mesh->getTriangles()[kdTri.index];

					primAABBLeft = tri.getClippedAABB(
						mesh->getVertexBuffer(), leftAABB);
					primAABBRight = tri.getClippedAABB(
						mesh->getVertexBuffer(), rightAABB);
				} else {
					if (shape->isClippable()) {
						primAABBLeft = shape->getClippedAABB(leftAABB);
						primAABBRight = shape->getClippedAABB(rightAABB);
					} else {
						primAABBLeft = primAABBRight = shape->getAABB();
						primAABBLeft.clip(leftAABB);
						primAABBRight.clip(rightAABB);
					}
				}

				/* Clipping the primitive can lead to an empty AABB, in which 
				   case it won't be added to the child node */
				if (primAABBLeft.isValid() && primAABBLeft.getSurfaceArea() > 0) {
					if (primAABBLeft.max[axis] == primAABBLeft.min[axis]) {
						/* The overlapping primitive is planar wrt. the current
						   axis. Thus, this can't be the split axis */
						// Some tiny rounding errors can cause this to fail
						// Assert(axis != bestSplit.axis);
						m_sEventsL.push_back(EdgeEvent(EdgeEvent::EEdgePlanar, 
							primAABBLeft.min[axis], event.index));
					} else {
						Assert(primAABBLeft.max[axis] > primAABBLeft.min[axis]);
						m_sEventsL.push_back(EdgeEvent(EdgeEvent::EEdgeStart, 
							primAABBLeft.min[axis], event.index));
						m_sEventsL.push_back(EdgeEvent(EdgeEvent::EEdgeEnd, 
							primAABBLeft.max[axis], event.index));
					}
				} else {
					/* The clipped primitive does not overlap the left node */
					clippedL++;
				}

				if (primAABBRight.isValid() && primAABBRight.getSurfaceArea() > 0) {
					if (primAABBRight.max[axis] == primAABBRight.min[axis]) {
						/* The overlapping primitive is planar wrt. the current
						   axis. Thus, this can't be the split axis */
						// Assert(axis != bestSplit.axis); Some tiny rounding errors can cause this
						m_sEventsR.push_back(EdgeEvent(EdgeEvent::EEdgePlanar, 
							primAABBRight.min[axis], event.index));
					} else {
						Assert(primAABBRight.max[axis] > primAABBRight.min[axis]);
						m_sEventsR.push_back(EdgeEvent(EdgeEvent::EEdgeStart, 
							primAABBRight.min[axis], event.index));
						m_sEventsR.push_back(EdgeEvent(EdgeEvent::EEdgeEnd, 
							primAABBRight.max[axis], event.index));
					}
				} else {
					/* The clipped primitive does not overlap the right node */
					clippedR++;
				}
			} else if (m_triangles[event.index].k != EBothSidesProcessed) {
				AssertEx(false, "Internal error");
			}
		}

		/* Remove the 'processed' flag */
		for (EdgeEventVec::const_iterator it = events.begin(); 
			it != events.end(); ++it) {
			if (m_triangles[(*it).index].k == EBothSidesProcessed) 
				m_triangles[(*it).index].k = EBothSides;
		}
	
		/* ================================================================ */
		/*                            Re-sorting                            */
		/* ================================================================ */

		/* The left and right lists now contain sorted event entries
		   from the non-overlapping primitives. The overlapping 
		   primitive events must still be sorted and merged into 
		   those list. */
		EdgeEventVec &rightEvents = m_rightEvents[depth][axis];
	
		/* Sort the events from overlapping prims */
		std::sort(m_sEventsL.begin(), m_sEventsL.end(), EdgeEventSorter());
		std::sort(m_sEventsR.begin(), m_sEventsR.end(), EdgeEventSorter());

		/* Set the target list sizes */
		allEvents[axis].resize(m_eventsL.size() + m_sEventsL.size());
		rightEvents.resize(m_eventsR.size() + m_sEventsR.size());

		/* Merge the left list */
		std::merge(m_eventsL.begin(), m_eventsL.end(),
			m_sEventsL.begin(), m_sEventsL.end(), allEvents[axis].begin(),
			EdgeEventSorter());

		/* Merge the right list */
		std::merge(m_eventsR.begin(), m_eventsR.end(),
			m_sEventsR.begin(), m_sEventsR.end(), rightEvents.begin(),
			EdgeEventSorter());
	}

	/* ==================================================================== */
    /*                              Recursion                               */
    /* ==================================================================== */

	if (m_clip) {
		clippedL /= 3; clippedR /= 3;
		m_clippedAway += clippedL + clippedR;
		primsLeft -= clippedL;
		primsRight -= clippedR;
	}

	m_innerNodes++;

	int leftNodeIndex = (int) m_nodes.size();
	int rightNodeIndex = (int) m_nodes.size()+1;

	/* Allocate room for the right node and its left child */
	m_nodes.push_back(KDNode()); m_nodes.push_back(KDNode());

	/* Update the current node and recurse */
	m_nodes[nodeIndex].setInner(bestSplit.axis, leftNodeIndex-nodeIndex, bestSplit.t);
	buildTree(leftNodeIndex, depth+1, badRefines, primsLeft + primsBoth, 
		leftAABB, allEvents);
	buildTree(rightNodeIndex, depth+1, badRefines, primsRight + primsBoth,
		rightAABB, m_rightEvents[depth]);
}

void KDTree::createLeaf(int nodeIndex, int depth, int numPrims, EdgeEventVec3 &allEvents) {
	int numAdded = 0;

	/* Keep some tree statistics */
	if (depth > m_actualMaxDepth)
		m_actualMaxDepth = depth;
	m_leafNodes++;
	m_leafPrims += numPrims;
	m_totalLeafDepth += depth;
	m_maxLeafPrims = std::max(m_maxLeafPrims, numPrims);
	m_minLeafPrims = std::min(m_minLeafPrims, numPrims);
	if (numPrims < m_bucketCount)
		m_triBuckets[numPrims]++;

#if defined(MTS_USE_TRIACCEL4)
	/* Directly pack and store the triangles */
	m_nodes[nodeIndex].setLeaf(
		(unsigned int) m_tempPackedTriangles.size(), 
		(int) std::ceil(numPrims/4.0f));
	Point A[4], B[4], C[4];
	uint32_t packIndex = 0, shapeIndex[4], index[4];
	bool nonTriFlag = false;
	for (unsigned int i=0; i<allEvents[0].size(); i++) {
		const EdgeEvent &event = allEvents[0][i];
		if (event.type == EdgeEvent::EEdgeStart 
				|| event.type == EdgeEvent::EEdgePlanar) {
			if (packIndex == 0) {
				/* Starting a new TriAccel4 */
				memset(A, 0, sizeof(Point)*4);
				memset(B, 0, sizeof(Point)*4);
				memset(C, 0, sizeof(Point)*4);
				memset(index, 0, sizeof(int)*4);
				memset(shapeIndex, 0, sizeof(int)*4);
				nonTriFlag = false;
			}
			const KDTriangle &t = m_triangles[event.index];
			if (t.index != KNoTriangleFlag) {
				const TriMesh *mesh = static_cast<const TriMesh *>(m_shapes[t.shapeIndex]);
				const Vertex *vertices = mesh->getVertexBuffer();
				const Triangle &tri = mesh->getTriangles()[t.index];
				A[packIndex] = vertices[tri.idx[0]].v;
				B[packIndex] = vertices[tri.idx[1]].v;
				C[packIndex] = vertices[tri.idx[2]].v;
			} else {
				nonTriFlag = true;
			}
			shapeIndex[packIndex] = t.shapeIndex;
			index[packIndex] = t.index;
			packIndex++; numAdded++;
			if (packIndex == 4) {
				/* Full, add to the list */
				packIndex = 0;
				TriAccel4 ta;
				m_failures += ta.load(A, B, C, shapeIndex, index);
				ta.nonTriFlag = nonTriFlag;
				#if defined(MTS_USE_TRIACCEL)
					/* Both TriAccel and TriAccel4 are enabled (for 
					   fast coherent+non-coherent RT). */
					ta.indirectionIndex = (unsigned int) m_indices.size();
					ta.indirectionCount = numPrims;
				#endif
				m_tempPackedTriangles.push_back(ta);
				m_packedTriangleCount++;
			}
		}
	}
	if (packIndex > 0) {
		/* Partially filled */
		TriAccel4 ta;
		m_failures += ta.load(A, B, C, shapeIndex, index);
		#if defined(MTS_USE_TRIACCEL)
			ta.indirectionIndex = (unsigned int) m_indices.size();
			ta.indirectionCount = numPrims;
		#endif
		ta.nonTriFlag = nonTriFlag;
		m_tempPackedTriangles.push_back(ta);
		m_packedTriangleCount++;
	}
	#if defined(MTS_USE_TRIACCEL)
		for (EdgeEventVec::const_iterator it = allEvents[0].begin(); 
			it != allEvents[0].end(); ++it) {
			const EdgeEvent &event = *it;
			if (event.type == EdgeEvent::EEdgeStart 
					|| event.type == EdgeEvent::EEdgePlanar) {
				m_indices.push_back(event.index);
			}
		}
	#endif
#else
	/* Allocate redirection indices */
	m_nodes[nodeIndex].setLeaf(m_indices.size(), numPrims);
	for (EdgeEventVec::const_iterator it = allEvents[0].begin(); 
		it != allEvents[0].end(); ++it) {
		const EdgeEvent &event = *it;
		if (event.type == EdgeEvent::EEdgeStart 
				|| event.type == EdgeEvent::EEdgePlanar) {
			m_indices.push_back(event.index);
			numAdded++;
		}
	}
#endif
	Assert(numAdded == numPrims);
}

MTS_NAMESPACE_END
