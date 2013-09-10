/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/timer.h>

#define MTS_QTREE_MAXDEPTH 50

MTS_NAMESPACE_BEGIN

namespace {
	/// Find the smallest t >= 0 such that a*t + b is a multiple of c
	inline Float nextMultiple(Float a, Float b, Float c) {
		Float tmp     = b/c,
			  rounded = (a > 0 ? std::ceil(tmp) : std::floor(tmp)) * c,
			  diff    = rounded - b;

		if (diff == 0)
			diff = signum(a) * c;

		return diff / a;
	}

	/// Temporary storage for patch-ray intersections
	struct PatchIntersectionRecord {
		Point p;
		int x, y;
	};

	/// Stack entry for recursive quadtree traversal
	struct StackEntry {
		int level, x, y;
	};
};

class Heightfield : public Shape {
public:
	Heightfield(const Properties &props) : Shape(props), m_data(NULL), m_normals(NULL) {
		m_sizeHint = Vector2i(
			props.getInteger("width", -1),
			props.getInteger("height", -1)
		);

		m_objectToWorld = props.getTransform("toWorld", Transform());
		m_shadingNormals = props.getBoolean("shadingNormals", true);
	}

	Heightfield(Stream *stream, InstanceManager *manager)
		: Shape(stream, manager), m_data(NULL), m_normals(NULL) {
	}

	~Heightfield() {
		if (m_data) {
			freeAligned(m_data);
			for (int i=0; i<m_levelCount; ++i)
				freeAligned(m_minmax[i]);
			delete[] m_minmax;
			delete[] m_levelSize;
			delete[] m_numChildren;
			delete[] m_blockSize;
		}
		if (m_normals)
			freeAligned(m_normals);
    }

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
	}

	AABB getAABB() const {
		AABB result;
		for (int i=0; i<8; ++i)
			result.expandBy(m_objectToWorld(m_dataAABB.getCorner(i)));

		return result;
	}

	Float getSurfaceArea() const {
		return m_surfaceArea; /// XXX transformed surface area? ...
	}

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return (size_t) m_levelSize[0].x * (size_t) m_levelSize[0].y;
	}

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *tmp) const {
		StackEntry stack[MTS_QTREE_MAXDEPTH];

		/* Transform ray into object space */
		Ray ray;
		m_objectToWorld.inverse()(_ray, ray);

		/* Ray length to cross a single cell along the X or Y axis */
		Float tDeltaXSingle = std::abs(ray.dRcp.x),
		      tDeltaYSingle = std::abs(ray.dRcp.y);

		/* Cell coordinate increments for steps along the ray */
		int iDeltaX = signumToInt(ray.d.x),
			iDeltaY = signumToInt(ray.d.y);

		if (iDeltaX == 0 && iDeltaY == 0) {
			/* TODO: special case for perpendicular rays */
			return false;
		}

		int stackIdx = 0;
		stack[stackIdx].level = m_levelCount-1;
		stack[stackIdx].x = 0;
		stack[stackIdx].y = 0;

		while (stackIdx >= 0) {
			StackEntry entry         = stack[stackIdx--];
			const Interval &interval = m_minmax[entry.level][entry.x + entry.y * m_levelSize[entry.level].x];
			const Vector2 &blockSize = m_blockSize[entry.level];
			Ray localRay(Point(ray.o.x - entry.x*blockSize.x, ray.o.y - entry.y*blockSize.y, ray.o.z), ray.d, 0);

			/* Intersect against the current min-max quadtree node */
			AABB aabb(
				Point3(0, 0, interval.min),
				Point3(blockSize.x, blockSize.y, interval.max)
			);

			Float nearT = mint, farT = maxt;
			Point enterPt, exitPt;

			if (!aabb.rayIntersect(localRay, nearT, farT, enterPt, exitPt))
				continue;

			Float tMax = farT - nearT;

			if (entry.level > 0) {
				/* Inner node -- push child nodes in 2D DDA order */
				const Vector2i &numChildren = m_numChildren[entry.level];
				const Vector2 &subBlockSize = m_blockSize[--entry.level];
				entry.x *= numChildren.x; entry.y *= numChildren.y;

				int x = (exitPt.x >= subBlockSize.x) ? numChildren.x-1 : 0;
				int y = (exitPt.y >= subBlockSize.y) ? numChildren.y-1 : 0;

				Float tDeltaX = tDeltaXSingle * subBlockSize.x,
				      tDeltaY = tDeltaYSingle * subBlockSize.y,
					  tNextX  = nextMultiple(-ray.d.x, exitPt.x, subBlockSize.x),
					  tNextY  = nextMultiple(-ray.d.y, exitPt.y, subBlockSize.y),
					  t       = 0;

				while ((uint32_t) x < (uint32_t) numChildren.x &&
				       (uint32_t) y < (uint32_t) numChildren.y && t <= tMax) {
					stack[++stackIdx].level = entry.level;
					stack[stackIdx].x = entry.x + x;
					stack[stackIdx].y = entry.y + y;

					if (tNextX < tNextY) {
						t = tNextX;
						tNextX += tDeltaX;
						x -= iDeltaX;
					} else {
						t = tNextY;
						tNextY += tDeltaY;
						y -= iDeltaY;
					}
				}
			} else {
				Float
					f00 = m_data[entry.y * m_dataSize.x + entry.x],
					f01 = m_data[(entry.y + 1) * m_dataSize.x + entry.x],
					f10 = m_data[entry.y * m_dataSize.x + entry.x + 1],
					f11 = m_data[(entry.y + 1) * m_dataSize.x + entry.x + 1];

				Float A = ray.d.x * ray.d.y * (f00 - f01 - f10 + f11);
				Float B = ray.d.y * (f01 - f00 + enterPt.x * (f00 - f01 - f10 + f11))
				        + ray.d.x * (f10 - f00 + enterPt.y * (f00 - f01 - f10 + f11))
						- ray.d.z;
				Float C = (enterPt.x - 1) * (enterPt.y - 1) * f00
				        + enterPt.y * f01 + enterPt.x * (f10 - enterPt.y * (f01 + f10 - f11))
				        - enterPt.z;

				Float t0, t1;
				if (!solveQuadratic(A, B, C, t0, t1))
					continue;

				if (t0 >= -Epsilon && t0 <= tMax + Epsilon)
					t = t0;
				else if (t1 >= -Epsilon && t1 <= tMax + Epsilon)
					t = t1;
				else
					continue;

				if (tmp) {
					PatchIntersectionRecord &temp = *((PatchIntersectionRecord *) tmp);
					Point pLocal = enterPt + ray.d * t;
					temp.x = entry.x;
					temp.y = entry.y;
					temp.p = pLocal;
					t += nearT;
				}

				return true;
			}
		}

		return false;
	}

	void fillIntersectionRecord(const Ray &ray,
			const void *tmp, Intersection &its) const {
		PatchIntersectionRecord &temp = *((PatchIntersectionRecord *) tmp);

		int x = temp.x, y = temp.y, width = m_dataSize.x;
		Float
			f00 = m_data[y     * width + x],
			f01 = m_data[(y+1) * width + x],
			f10 = m_data[y     * width + x + 1],
			f11 = m_data[(y+1) * width + x + 1];

		its.p  = Point(temp.p.x + temp.x, temp.p.y + temp.y, temp.p.z);

		its.uv = Point2(its.p.x / m_levelSize[0].x, its.p.y / m_levelSize[0].y);
		its.dpdu = Vector(1, 0, (1.0f - temp.p.y) * (f10 - f00) + temp.p.y * (f11 - f01)) * m_levelSize[0].x;
		its.dpdv = Vector(0, 1, (1.0f - temp.p.x) * (f01 - f00) + temp.p.x * (f11 - f10)) * m_levelSize[0].y;

		its.geoFrame.s = normalize(its.dpdu);
		its.geoFrame.t = normalize(its.dpdv - dot(its.dpdv, its.geoFrame.s) * its.geoFrame.s);
		its.geoFrame.n = cross(its.geoFrame.s, its.geoFrame.t);

		if (m_shadingNormals) {
			Normal
				n00 = m_normals[y     * width + x],
				n01 = m_normals[(y+1) * width + x],
				n10 = m_normals[y     * width + x + 1],
				n11 = m_normals[(y+1) * width + x + 1];

			its.shFrame.n = normalize(
				(1 - temp.p.x) * ((1-temp.p.y) * n00 + temp.p.y * n01)
				   + temp.p.x  * ((1-temp.p.y) * n10 + temp.p.y * n11));

			its.shFrame.s = normalize(its.geoFrame.s - dot(its.geoFrame.s, its.shFrame.n) * its.shFrame.n);
			its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
		} else {
			its.shFrame = its.geoFrame;
		}
		its.shape = this;

 		its.wi = its.toLocal(-ray.d);
 		its.hasUVPartials = false;
		its.instance = NULL;
		its.time = ray.time;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Float t;
		return rayIntersect(ray, mint, maxt, t, NULL);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(Texture::m_theClass)) {
			if (m_data != NULL)
				Log(EError, "Attempted to attach multiple textures to a height field shape!");

			ref<Bitmap> bitmap = static_cast<Texture *>(child)->getBitmap(m_sizeHint);

			m_dataSize = bitmap->getSize();
			if (m_dataSize.x < 2) m_dataSize.x = 2;
			if (m_dataSize.y < 2) m_dataSize.y = 2;
			if (!isPowerOfTwo(m_dataSize.x - 1)) m_dataSize.x = (int) roundToPowerOfTwo((uint32_t) m_dataSize.x - 1) + 1;
			if (!isPowerOfTwo(m_dataSize.y - 1)) m_dataSize.y = (int) roundToPowerOfTwo((uint32_t) m_dataSize.y - 1) + 1;

			if (bitmap->getSize() != m_dataSize) {
				Log(EInfo, "Resampling heightfield texture from %ix%i to %ix%i..",
					bitmap->getWidth(), bitmap->getHeight(), m_dataSize.x, m_dataSize.y);

				bitmap = bitmap->resample(NULL, ReconstructionFilter::EClamp,
					ReconstructionFilter::EClamp, m_dataSize,
					-std::numeric_limits<Float>::infinity(),
					std::numeric_limits<Float>::infinity());
			}

			size_t size = (size_t) m_dataSize.x * (size_t) m_dataSize.y * sizeof(Float),
			       storageSize = size;
			m_data = (Float *) allocAligned(size);

			if (m_shadingNormals) {
				size *= 3;
				m_normals = (Normal *) allocAligned(size);
				memset(m_normals, 0, size);
				storageSize += size;
			}

			bitmap->convert(m_data, Bitmap::ELuminance, Bitmap::EFloat);

			Log(EInfo, "Building acceleration data structure for %ix%i height field ..", m_dataSize.x, m_dataSize.y);

			ref<Timer> timer = new Timer();
			m_levelCount = (int) std::max(log2i((uint32_t) m_dataSize.x-1), log2i((uint32_t) m_dataSize.y-1)) + 1;

			m_levelSize = new Vector2i[m_levelCount];
			m_numChildren = new Vector2i[m_levelCount];
			m_blockSize = new Vector2[m_levelCount];
			m_minmax = new Interval*[m_levelCount];

			m_levelSize[0]  = Vector2i(m_dataSize.x - 1, m_dataSize.y - 1);
			m_blockSize[0] = Vector2(1, 1);
			m_surfaceArea = 0;
			size = (size_t) m_levelSize[0].x * (size_t) m_levelSize[0].y * sizeof(Interval);
			m_minmax[0] = (Interval *) allocAligned(size);
			storageSize += size;

			/* Build the lowest MIP layer directly from the heightfield data */
			Interval *bounds = m_minmax[0];
			for (int y=0; y<m_levelSize[0].y; ++y) {
				for (int x=0; x<m_levelSize[0].x; ++x) {
					Float f00 = m_data[y * m_dataSize.x + x];
					Float f10 = m_data[y * m_dataSize.x + x + 1];
					Float f01 = m_data[(y + 1) * m_dataSize.x + x];
					Float f11 = m_data[(y + 1) * m_dataSize.x + x + 1];
					Float fmin = std::min(std::min(f00, f01), std::min(f10, f11));
					Float fmax = std::max(std::max(f00, f01), std::max(f10, f11));
					*bounds++ = Interval(fmin, fmax);

					/* Estimate the total surface area (this is approximate) */
					Float diff0 = f01-f10, diff1 = f00-f11;
					m_surfaceArea += std::sqrt(1.0f + .5f * (diff0*diff0 + diff1*diff1));
				}
			}

			/* Propagate height bounds upwards to the other layers */
			for (int level=1; level<m_levelCount; ++level) {
				Vector2i &cur  = m_levelSize[level],
				         &prev = m_levelSize[level-1];

				/* Calculate size of this layer */
				cur.x = prev.x > 1 ? (prev.x / 2) : 1;
				cur.y = prev.y > 1 ? (prev.y / 2) : 1;

				m_numChildren[level].x = prev.x > 1 ? 2 : 1;
				m_numChildren[level].y = prev.y > 1 ? 2 : 1;
				m_blockSize[level] = Vector2(
					m_levelSize[0].x / cur.x,
					m_levelSize[0].y / cur.y
				);

				/* Allocate memory for interval data */
				Interval *prevBounds = m_minmax[level-1], *curBounds;
				size_t size = (size_t) cur.x * (size_t) cur.y * sizeof(Interval);
				m_minmax[level] = curBounds = (Interval *) allocAligned(size);
				storageSize += size;

				/* Build by querying the previous layer */
				for (int y=0; y<cur.y; ++y) {
					int y0 = std::min(2*y,   prev.y-1),
						y1 = std::min(2*y+1, prev.y-1);
					for (int x=0; x<cur.x; ++x) {
						int x0 = std::min(2*x,   prev.x-1),
							x1 = std::min(2*x+1, prev.x-1);
						const Interval &f00 = prevBounds[y0 * prev.x + x0], &f01 = prevBounds[y0 * prev.x + x1];
						const Interval &f10 = prevBounds[y1 * prev.x + x0], &f11 = prevBounds[y1 * prev.x + x1];
						Interval combined(f00);
						combined.expandBy(f01);
						combined.expandBy(f10);
						combined.expandBy(f11);
						*curBounds++ = combined;
					}
				}
			}

			if (m_shadingNormals) {
				Log(EInfo, "Precomputing shading normals ..");

				for (int y=0; y<m_levelSize[0].y; ++y) {
					for (int x=0; x<m_levelSize[0].x; ++x) {
						Float f00 = m_data[y * m_dataSize.x + x];
						Float f10 = m_data[y * m_dataSize.x + x + 1];
						Float f01 = m_data[(y + 1) * m_dataSize.x + x];
						Float f11 = m_data[(y + 1) * m_dataSize.x + x + 1];

						m_normals[y       * m_dataSize.x + x]     += normalize(Normal(f00 - f10, f00 - f01, 1));
						m_normals[y       * m_dataSize.x + x + 1] += normalize(Normal(f00 - f10, f10 - f11, 1));
						m_normals[(y + 1) * m_dataSize.x + x]     += normalize(Normal(f01 - f11, f00 - f01, 1));
						m_normals[(y + 1) * m_dataSize.x + x + 1] += normalize(Normal(f01 - f11, f10 - f11, 1));
					}
				}

				for (int y=0; y<m_dataSize.y; ++y) {
					for (int x=0; x<m_dataSize.x; ++x) {
						Normal &normal = m_normals[x + y * m_dataSize.x];
						normal /= normal.length();
					}
				}
			}

			Log(EInfo, "Done (took %i ms, uses %s of memory)", timer->getMilliseconds(),
					memString(storageSize).c_str());

			m_dataAABB = AABB(
				Point3(0, 0, m_minmax[m_levelCount-1][0].min),
				Point3(m_levelSize[0].x, m_levelSize[0].y, m_minmax[m_levelCount-1][0].max)
			);
		} else {
			Shape::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Heightfield[" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_objectToWorld;
	Vector2i m_sizeHint;
	AABB m_dataAABB;
	bool m_shadingNormals;

	/* Height field data */
	Float *m_data;
	Normal *m_normals;
	Vector2i m_dataSize;
	Float m_surfaceArea;

	/* Min-max quadtree data */
	int m_levelCount;
	Vector2i *m_levelSize;
	Vector2i *m_numChildren;
	Vector2 *m_blockSize;
	Interval **m_minmax;
};

MTS_IMPLEMENT_CLASS_S(Heightfield, false, Shape)
MTS_EXPORT_PLUGIN(Heightfield, "Height field intersection primitive");
MTS_NAMESPACE_END

