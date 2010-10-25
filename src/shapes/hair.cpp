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

#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/gkdtree.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

#define MTS_HAIR_USE_FANCY_CLIPPING 1

MTS_NAMESPACE_BEGIN

/**
 * \brief Space-efficient acceleration structure for cylindrical hair
 * segments with miter joints.
 */
class HairKDTree : public GenericKDTree<HairKDTree> {
	friend class GenericKDTree<HairKDTree>;
public:
	HairKDTree(std::vector<Point> &vertices, 
			std::vector<bool> &vertexStartsFiber, Float radius)
			: m_radius(radius) {
		/* Take the supplied vertex & start fiber arrays (without copying) */
		m_vertices.swap(vertices);
		m_vertexStartsFiber.swap(vertexStartsFiber);

		/* Compute the index of the first vertex in each segment. */
		m_segIndex.reserve(m_vertices.size());
		for (size_t i=0; i<m_vertices.size()-1; i++)
			if (!m_vertexStartsFiber[i+1])
				m_segIndex.push_back(i);
		m_segmentCount = m_segIndex.size();

		Log(EDebug, "Building a kd-tree for " SIZE_T_FMT " hair vertices, "
			SIZE_T_FMT " segments,", m_vertices.size(), m_segmentCount);

		/* Ray-cylinder intersections are expensive. Use only the
		   SAH cost as the tree subdivision stopping criterion, 
		   not the number of primitives */
		setStopPrims(0);
		setTraversalCost(10);
		setIntersectionCost(30);
		buildInternal();

		/* Optimization: replace all primitive indices by the
		   associated vertex indices (this avoids an extra 
		   indirection during traversal later on) */
		for (size_type i=0; i<m_indexCount; ++i)
			m_indices[i] = m_segIndex[m_indices[i]];

		/* Free the segIndex array, it is not needed anymore */
		std::vector<index_type>().swap(m_segIndex);
	}

	/// Return the AABB of the hair kd-tree
	inline const AABB &getAABB() const {
		return m_aabb;
	}

	/// Return the list of vertices underlying the hair kd-tree
	inline const std::vector<Point> &getVertices() const {
		return m_vertices;
	}

	/**
	 * Return a boolean list specifying whether a vertex 
	 * marks the beginning of a new fiber
	 */
	inline const std::vector<bool> &getStartFiber() const {
		return m_vertexStartsFiber;
	}

	/// Return the radius of the hairs stored in the kd-tree
	inline Float getRadius() const {
		return m_radius;
	}

	/// Return the total number of segments
	inline size_t getSegmentCount() const {
		return m_segmentCount;
	}

	/// Intersect a ray with all segments stored in the kd-tree
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt, 
			Float &t, void *temp) const {
		Float tempT = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<false>(ray, mint, maxt, tempT, temp)) {
					t = tempT;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * \brief Intersect a ray with all segments stored in the kd-tree
	 * (Visiblity query version)
	 */
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt) const {
		Float tempT = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL)) 
					return true;
			}
		}
		return false;
	}

#if defined(MTS_HAIR_USE_FANCY_CLIPPING)
	/**
	 * Compute the ellipse created by the intersection of an infinite
	 * cylinder and a plane. Returns false in the degenerate case.
	 * Based on:
	 * www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
	 */
	bool intersectCylPlane(Point planePt, Normal planeNrml,
			Point cylPt, Vector cylD, Float radius, Point &center,
			Vector *axes, Float *lengths) const {
		if (absDot(planeNrml, cylD) < Epsilon)
			return false;

		Assert(std::abs(planeNrml.length()-1) <Epsilon);
		Vector B, A = cylD - dot(cylD, planeNrml)*planeNrml;

		Float length = A.length();
		if (length > Epsilon && planeNrml != cylD) {
			A /= length;
			B = cross(planeNrml, A);
		} else {
			coordinateSystem(planeNrml, A, B);
		}

		Vector delta = planePt - cylPt,
			   deltaProj = delta - cylD*dot(delta, cylD);

		Float aDotD = dot(A, cylD);
		Float bDotD = dot(B, cylD);
		Float c0 = 1-aDotD*aDotD;
		Float c1 = 1-bDotD*bDotD;
		Float c2 = 2*dot(A, deltaProj);
		Float c3 = 2*dot(B, deltaProj);
		Float c4 = dot(delta, deltaProj) - radius*radius;

		Float lambda = (c2*c2/(4*c0) + c3*c3/(4*c1) - c4)/(c0*c1);

		Float alpha0 = -c2/(2*c0),
			  beta0 = -c3/(2*c1);

		lengths[0] = std::sqrt(c1*lambda),
		lengths[1] = std::sqrt(c0*lambda);

		center = planePt + alpha0 * A + beta0 * B;
		axes[0] = A;
		axes[1] = B;
		return true;
	}

	AABB intersectCylFace(int axis,
			const Point &min, const Point &max,
			const Point &cylPt, const Vector &cylD) const {
		int axis1 = (axis + 1) % 3;
		int axis2 = (axis + 2) % 3;

		Normal planeNrml(0.0f);
		planeNrml[axis] = 1;

		Point ellipseCenter;
		Vector ellipseAxes[2];
		Float ellipseLengths[2];

		AABB aabb;
		if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius, 
			ellipseCenter, ellipseAxes, ellipseLengths)) {
			/* Degenerate case -- return an invalid AABB. This is
			   not a problem, since one of the other faces will provide
			   enough information to arrive at a correct clipped AABB */
			return aabb;
		}

		/* Intersect the ellipse against the sides of the AABB face */
		for (int i=0; i<4; ++i) {
			Point p1, p2;
			p1[axis] = p2[axis] = min[axis];
			p1[axis1] = ((i+1) & 2) ? min[axis1] : max[axis1];
			p1[axis2] = ((i+0) & 2) ? min[axis2] : max[axis2];
			p2[axis1] = ((i+2) & 2) ? min[axis1] : max[axis1];
			p2[axis2] = ((i+1) & 2) ? min[axis2] : max[axis2];

			Point2 p1l(
				dot(p1 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
				dot(p1 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);
			Point2 p2l(
				dot(p2 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
				dot(p2 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);

			Vector2 rel = p2l-p1l;
			Float A = dot(rel, rel);
			Float B = 2*dot(Vector2(p1l), rel);
			Float C = dot(Vector2(p1l), Vector2(p1l))-1;

			Float x0, x1;
			if (solveQuadratic(A, B, C, x0, x1)) {
				if (x0 >= 0 && x0 <= 1)
					aabb.expandBy(p1+(p2-p1)*x0);
				if (x1 >= 0 && x1 <= 1)
					aabb.expandBy(p1+(p2-p1)*x1);
			}
		}

		ellipseAxes[0] *= ellipseLengths[0];
		ellipseAxes[1] *= ellipseLengths[1];
		AABB faceBounds(min, max);

		/* Find the componentwise maxima of the ellipse */
		for (int i=0; i<2; ++i) {
			int j = (i==0) ? axis1 : axis2;
			Float alpha = ellipseAxes[0][j];
			Float beta = ellipseAxes[1][j];
			Float ratio = beta/alpha, tmp = std::sqrt(1+ratio*ratio);
			Float cosTheta = 1/tmp, sinTheta = ratio/tmp;
			Point p1 = ellipseCenter + cosTheta*ellipseAxes[0] + sinTheta*ellipseAxes[1];
			Point p2 = ellipseCenter - cosTheta*ellipseAxes[0] - sinTheta*ellipseAxes[1];

			if (faceBounds.contains(p1)) 
				aabb.expandBy(p1);
			if (faceBounds.contains(p2)) 
				aabb.expandBy(p2);
		}

		return aabb;
	}

	AABB getAABB(index_type index) const {
		index_type iv = m_segIndex.at(index);
		Point center;
		Vector axes[2];
		Float lengths[2];

		bool success = intersectCylPlane(firstVertex(iv), firstMiterNormal(iv), 
			firstVertex(iv), tangent(iv), m_radius, center, axes, lengths);
		Assert(success);

		AABB result;
		axes[0] *= lengths[0]; axes[1] *= lengths[1];
		for (int i=0; i<3; ++i) {
			Float range = std::sqrt(axes[0][i]*axes[0][i] + axes[1][i]*axes[1][i]);
			result.min[i] = std::min(result.min[i], center[i]-range);
			result.max[i] = std::max(result.max[i], center[i]+range);
		}

		success = intersectCylPlane(secondVertex(iv), secondMiterNormal(iv), 
			secondVertex(iv), tangent(iv), m_radius, center, axes, lengths);
		Assert(success);

		axes[0] *= lengths[0]; axes[1] *= lengths[1];
		for (int i=0; i<3; ++i) {
			Float range = std::sqrt(axes[0][i]*axes[0][i] + axes[1][i]*axes[1][i]);
			result.min[i] = std::min(result.min[i], center[i]-range);
			result.max[i] = std::max(result.max[i], center[i]+range);
		}
		return result;
	}

	AABB getClippedAABB(index_type index, const AABB &box) const {
		/* Compute a base bounding box */
		AABB base(getAABB(index));
		base.clip(box);

		index_type iv = m_segIndex.at(index);

		Point cylPt = firstVertex(iv);
		Vector cylD = tangent(iv);

		/* Now forget about the cylinder ends and 
		   intersect an infinite cylinder with each AABB face */
		AABB clippedAABB;
		clippedAABB.expandBy(intersectCylFace(0, 
				Point(base.min.x, base.min.y, base.min.z),
				Point(base.min.x, base.max.y, base.max.z),
				cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(0,
				Point(base.max.x, base.min.y, base.min.z),
				Point(base.max.x, base.max.y, base.max.z),
				cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(1, 
				Point(base.min.x, base.min.y, base.min.z),
				Point(base.max.x, base.min.y, base.max.z),
				cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(1,
				Point(base.min.x, base.max.y, base.min.z),
				Point(base.max.x, base.max.y, base.max.z),
				cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(2, 
				Point(base.min.x, base.min.y, base.min.z),
				Point(base.max.x, base.max.y, base.min.z),
				cylPt, cylD));

		clippedAABB.expandBy(intersectCylFace(2,
				Point(base.min.x, base.min.y, base.max.z),
				Point(base.max.x, base.max.y, base.max.z),
				cylPt, cylD));

		clippedAABB.clip(base);
		return clippedAABB;
	}
#else
	/// Compute the AABB of a segment (only used during tree construction)
	AABB getAABB(int index) const {
		index_type iv = m_segIndex.at(index);

		// cosine of steepest miter angle
		const Float cos0 = dot(firstMiterNormal(iv), tangent(iv));
		const Float cos1 = dot(secondMiterNormal(iv), tangent(iv));
		const Float maxInvCos = 1.0 / std::min(cos0, cos1);
		const Vector expandVec(m_radius * maxInvCos);

		const Point a = firstVertex(iv);
		const Point b = secondVertex(iv);

		AABB aabb;
		aabb.expandBy(a - expandVec);
		aabb.expandBy(a + expandVec);
		aabb.expandBy(b - expandVec);
		aabb.expandBy(b + expandVec);
		return aabb;
	}

	/// Compute the clipped AABB of a segment (only used during tree construction)
	AABB getClippedAABB(int index, const AABB &box) const {
		AABB aabb(getAABB(index));
		aabb.clip(box);
		return aabb;
	}
#endif

	/// Return the total number of segments
	inline int getPrimitiveCount() const {
		return m_segIndex.size();
	}

	inline EIntersectionResult intersect(const Ray &ray, index_type iv, 
		Float mint, Float maxt, Float &t, void *tmp) const {
		/* First compute the intersection with the infinite cylinder */
		Float nearT, farT;
		Vector axis = tangent(iv);

		// Projection of ray onto subspace normal to axis
		Vector relOrigin = ray.o - firstVertex(iv);
		Vector projOrigin = relOrigin - dot(axis, relOrigin) * axis;
		Vector projDirection = ray.d - dot(axis, ray.d) * axis;

		// Quadratic to intersect circle in projection
		const Float A = projDirection.lengthSquared();
		const Float B = 2 * dot(projOrigin, projDirection);
		const Float C = projOrigin.lengthSquared() - m_radius*m_radius;

		if (!solveQuadratic(A, B, C, nearT, farT))
			return ENever;

		if (nearT > maxt || farT < mint)
			return ENo;

		/* Next check the intersection points against the miter planes */
		Point pointNear = ray(nearT);
		Point pointFar = ray(farT);
		if (dot(pointNear - firstVertex(iv), firstMiterNormal(iv)) >= 0 &&
			dot(pointNear - secondVertex(iv), secondMiterNormal(iv)) <= 0 &&
			nearT >= mint) {
			t = nearT;
		} else if (dot(pointFar - firstVertex(iv), firstMiterNormal(iv)) >= 0 &&
				dot(pointFar - secondVertex(iv), secondMiterNormal(iv)) <= 0) {
			if (farT > maxt)
				return ENo;
			t = farT;
		} else {
			return ENo;
		}

		index_type *storage = static_cast<index_type *>(tmp);
		if (storage)
			*storage = iv;

		return EYes;
	}
	
	inline EIntersectionResult intersect(const Ray &ray, index_type iv, 
		Float mint, Float maxt) const {
		Float tempT;
		return intersect(ray, iv, mint, maxt, tempT, NULL);
	}

	/* Some utility functions */
	inline Point firstVertex(index_type iv) const { return m_vertices[iv]; }
	inline Point secondVertex(index_type iv) const { return m_vertices[iv+1]; }
	inline Point prevVertex(index_type iv) const { return m_vertices[iv-1]; }
	inline Point nextVertex(index_type iv) const { return m_vertices[iv+2]; }

	inline bool prevSegmentExists(index_type iv) const { return !m_vertexStartsFiber[iv]; }
	inline bool nextSegmentExists(index_type iv) const { return !m_vertexStartsFiber[iv+2]; }

	inline Vector tangent(index_type iv) const { return normalize(secondVertex(iv) - firstVertex(iv)); }
	inline Vector prevTangent(index_type iv) const { return normalize(firstVertex(iv) - prevVertex(iv)); }
	inline Vector nextTangent(index_type iv) const { return normalize(nextVertex(iv) - secondVertex(iv)); }

	inline Vector firstMiterNormal(index_type iv) const {
		if (prevSegmentExists(iv))
			return normalize(prevTangent(iv) + tangent(iv));
		else
			return tangent(iv);
	}

	inline Vector secondMiterNormal(index_type iv) const {
		if (nextSegmentExists(iv))
			return normalize(tangent(iv) + nextTangent(iv));
		else
			return tangent(iv);
	}

	MTS_DECLARE_CLASS()
protected:
	std::vector<Point> m_vertices;
	std::vector<bool> m_vertexStartsFiber;
	std::vector<index_type> m_segIndex;
	size_t m_segmentCount;
	Float m_radius;
};

class Hair : public Shape {
public:
	Hair(const Properties &props) : Shape(props) {
		fs::path path = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		Float radius = (Float) props.getFloat("radius", 0.05f);

		/* Object-space -> World-space transformation */
		Transform objectToWorld = props.getTransform("toWorld", Transform());

		Log(EInfo, "Loading hair geometry from \"%s\" ..", path.leaf().c_str());

		fs::ifstream is(path);
		if (is.fail())
			Log(EError, "Could not open \"%s\"!", path.file_string().c_str());

		std::string line;
		bool newFiber = true;
		Point p, lastP(0.0f);
		std::vector<Point> vertices;
		std::vector<bool> vertexStartsFiber;
		size_t nDegenerate = 0;

		while (is.good()) {
			std::getline(is, line);
			if (line.length() > 0 && line[0] == '#') {
				newFiber = true;
				continue;
			}
			std::istringstream iss(line);
			iss >> p.x >> p.y >> p.z;
			if (!iss.fail()) {
				if (newFiber || p != lastP) {
					vertices.push_back(objectToWorld(p));
					vertexStartsFiber.push_back(newFiber);
					lastP = p;
				} else {
					nDegenerate++;
				}
				newFiber = false;
			} else {
				newFiber = true;
			}
		}

		if (nDegenerate > 0)
			Log(EInfo, "Encountered " SIZE_T_FMT 
				" degenerate segments!", nDegenerate);

		vertexStartsFiber.push_back(true);

		m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
	}

	Hair(Stream *stream, InstanceManager *manager) 
		: Shape(stream, manager) {
		Float radius = stream->readFloat();
		size_t vertexCount = (size_t) stream->readUInt();

		std::vector<Point> vertices(vertexCount);
		std::vector<bool> vertexStartsFiber(vertexCount+1);
		stream->readFloatArray((Float *) &vertices[0], vertexCount * 3);

		for (size_t i=0; i<vertexCount; ++i) 
			vertexStartsFiber[i] = stream->readBool();
		vertexStartsFiber[vertexCount] = true;

		m_kdtree = new HairKDTree(vertices, vertexStartsFiber, radius);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		const std::vector<Point> &vertices = m_kdtree->getVertices();
		const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();

		stream->writeFloat(m_kdtree->getRadius());
		stream->writeUInt((uint32_t) vertices.size());
		stream->writeFloatArray((Float *) &vertices[0], vertices.size() * 3);
		for (size_t i=0; i<vertices.size(); ++i)
			stream->writeBool(vertexStartsFiber[i]);
	}

	bool rayIntersect(const Ray &ray, Float mint, 
			Float maxt, Float &t, void *temp) const {
		return m_kdtree->rayIntersect(ray, mint, maxt, t, temp);
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		return m_kdtree->rayIntersect(ray, mint, maxt);
	}

	void fillIntersectionRecord(const Ray &ray, Float t, 
		const void *temp, Intersection &its) const {
		its.p = ray(t);

		/* No UV coordinates for now */
		its.uv = Point2(0,0);
		its.dpdu = Vector(0,0,0);
		its.dpdv = Vector(0,0,0);

		const HairKDTree::index_type *storage = 
			static_cast<const HairKDTree::index_type *>(temp);
		HairKDTree::index_type iv = *storage;

		its.geoFrame.s = m_kdtree->tangent(iv);
		const Vector relHitPoint = its.p - m_kdtree->firstVertex(iv);
		const Vector axis = m_kdtree->tangent(iv);
		its.geoFrame.n = Normal(relHitPoint - dot(axis, relHitPoint) * axis);
		its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);
		its.shFrame = its.geoFrame;
		its.wi = its.toLocal(-ray.d);
		its.hasUVPartials = false;
		its.shape = this;
	}

	ref<TriMesh> createTriMesh() {
		/// Choice of discretization
		const size_t phiSteps = 10;
		const Float dPhi   = (2*M_PI) / phiSteps;

		size_t nSegments = m_kdtree->getSegmentCount();
		ref<TriMesh> mesh = new TriMesh("Hair mesh approximation",
			phiSteps*2*nSegments, phiSteps*2*nSegments, true, false, false);

		Point *vertices = mesh->getVertexPositions();
		Normal *normals = mesh->getVertexNormals();
		Triangle *triangles = mesh->getTriangles();
		size_t triangleIdx = 0, vertexIdx = 0;
		
		const std::vector<Point> &hairVertices = m_kdtree->getVertices();
		const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();
		const Float radius = m_kdtree->getRadius();
		Float *cosPhi = new Float[phiSteps];
		Float *sinPhi = new Float[phiSteps];
		for (size_t i=0; i<phiSteps; ++i) {
			sinPhi[i] = std::sin(i*dPhi);
			cosPhi[i] = std::cos(i*dPhi);
		}

		size_t hairIdx = 0;
		for (size_t iv=0; iv<hairVertices.size()-1; iv++) {
			if (!vertexStartsFiber[iv+1]) {
				for (size_t phi=0; phi<phiSteps; ++phi) {
					Vector tangent = m_kdtree->tangent(iv);
					Vector dir = Normal(Frame(tangent).toWorld(
							Vector(cosPhi[phi], sinPhi[phi], 0)));
					Normal miterNormal1 = m_kdtree->firstMiterNormal(iv);
					Normal miterNormal2 = m_kdtree->secondMiterNormal(iv);
					Float t1 = dot(miterNormal1, radius*dir) / dot(miterNormal1, tangent);
					Float t2 = dot(miterNormal2, radius*dir) / dot(miterNormal2, tangent);

					normals[vertexIdx] = Normal(dir);
					vertices[vertexIdx++] = m_kdtree->firstVertex(iv) + radius*dir - tangent*t1;
					normals[vertexIdx] = Normal(dir);
					vertices[vertexIdx++] = m_kdtree->secondVertex(iv) + radius*dir - tangent*t2;

					int idx0 = 2*(phi + hairIdx*phiSteps), idx1 = idx0+1;
					int idx2 = (2*phi+2) % (2*phiSteps) + 2*hairIdx*phiSteps, idx3 = idx2+1;
					triangles[triangleIdx].idx[0] = idx0;
					triangles[triangleIdx].idx[1] = idx2;
					triangles[triangleIdx].idx[2] = idx1;
					triangleIdx++;
					triangles[triangleIdx].idx[0] = idx1;
					triangles[triangleIdx].idx[1] = idx2;
					triangles[triangleIdx].idx[2] = idx3;
					triangleIdx++;
				}
				hairIdx++;
			}
		}
		Assert(triangleIdx == phiSteps*2*nSegments);
		Assert(vertexIdx == phiSteps*2*nSegments);

		delete[] cosPhi;
		delete[] sinPhi;

		mesh->setBSDF(m_bsdf);
		mesh->setLuminaire(m_luminaire);
		mesh->configure();

		return mesh.get();
	}

	const AbstractKDTree *getKDTree() const {
		return m_kdtree.get();
	}

	AABB getAABB() const {
		return m_kdtree->getAABB();
	}

	Float getSurfaceArea() const {
		Log(EError, "Hair::getSurfaceArea(): Not implemented.");
		return -1;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Hair[" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<HairKDTree> m_kdtree;
};

MTS_IMPLEMENT_CLASS(HairKDTree, false, GenericKDTree)
MTS_IMPLEMENT_CLASS_S(Hair, false, Shape)
MTS_EXPORT_PLUGIN(Hair, "Hair intersection primitive");
MTS_NAMESPACE_END
