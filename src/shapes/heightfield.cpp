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


#define MTS_QTREE_MAXDEPTH 24 /* Up to 16M x 16M */

		/* TODO: special case for perpendicular rays */
MTS_NAMESPACE_BEGIN

namespace {
	/// Find the smallest t >= 0 such that a*t + b is a multiple of c
	inline Float nextMultiple(Float a, Float b, Float c) {
		if (a == 0)
			return std::numeric_limits<Float>::infinity();
		else if (a > 0)
			return (std::ceil(b/c)*c - b) / a;
		else
			return (std::floor(b/c)*c - b) / a;
	}
};

class Heightfield : public Shape {
public:
	Heightfield(const Properties &props) : Shape(props), m_data(NULL) {
		m_sizeHint = Vector2i(
			props.getInteger("width", -1),
			props.getInteger("height", -1)
		);

		m_objectToWorld = props.getTransform("toWorld", Transform());
	}

	Heightfield(Stream *stream, InstanceManager *manager)
		: Shape(stream, manager), m_data(NULL) {
	}

	~Heightfield() {
		if (m_data) {
			freeAligned(m_data);
			for (int i=0; i<m_levelCount; ++i)
				freeAligned(m_minmax[i]);
			delete[] m_minmax;
			delete[] m_levelSize;
			delete[] m_blockSizeF;
		}
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
		return 0; /// XXX
	}

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return (size_t) m_levelSize[0].x * (size_t) m_levelSize[0].y;
	}
	struct StackEntry {
		uint32_t level;
		uint32_t x, y;

		inline std::string toString() const {
			std::ostringstream oss;
			oss << "StackEntry[level=" << level << ", x=" << x << ", y=" << y << "]" << endl;
			return oss.str();
		}
	};

	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *tmp) const {
		StackEntry stack[MTS_QTREE_MAXDEPTH];

		/* Transform ray into object space */
		Ray ray;
		m_objectToWorld.inverse()(_ray, ray);

		/* Clamp ray to the bounding box */
	    Float nearT, farT;
	    if (!m_dataAABB.rayIntersect(ray, nearT, farT))
			return false;
		if (nearT > mint) mint = nearT;
		if (farT  < maxt) maxt = farT;
		if (maxt <= mint) return false;

		/* Rescale the ray so that it has unit length in 2D */
//		Float length2D = hypot2(ray.d.x * (maxt-mint), ray.d.y * (maxt-mint));
//		ray.d /= length2D; ray.dRcp *= length2D;

		/* Ray length to cross a single cell along the X or Y axis */
		Float tDeltaXSingle  = std::abs(ray.dRcp.x),
		      tDeltaYSingle  = std::abs(ray.dRcp.y);

		/* Cell coordinate increments for steps along the ray */
		int iDeltaX = (int) signum(ray.d.x),
			iDeltaY = (int) signum(ray.d.y);

		int stackIdx = 0;
		stack[stackIdx].level = m_levelCount-1;
		stack[stackIdx].x = 0;
		stack[stackIdx].y = 0;

		while (stackIdx >= 0) {
			StackEntry entry         = stack[stackIdx];
			const Vector2i &size     = m_levelSize[entry.level];
			const Interval &interval = m_minmax[entry.level][entry.x + entry.y * m_levelSize[entry.level].x];
			const Vector2 &blockSize = m_blockSizeF[entry.level];

			/* Intersect against the current min-max quadtree node */
			AABB aabb(
				Point3(entry.x       * blockSize.x, entry.y       * blockSize.y, interval.min),
				Point3((entry.x + 1) * blockSize.x, (entry.y + 1) * blockSize.y, interval.max)
			);

			Float nearT, farT;
			bool match = aabb.rayIntersect(ray, nearT, farT);

			if (!match || farT < mint || nearT > maxt) {
				/* Miss -> pop the stack */
				--stackIdx;
				continue;
			}

			if (entry.level > 0) {
				/* Inner node -- push child nodes in 2D DDA order */
				const Vector2 &subBlockSize = m_blockSizeF[entry.level+1];
				Point2 p(ray.o.x + ray.d.x * nearT,
				         ray.o.y + ray.d.y * nearT);

				int nextLevel = entry.level + 1;
				int x = p.x - aabb.min.x >= subBlockSize.x ? 1 : 0;
				int y = p.y - aabb.min.y >= subBlockSize.y ? 1 : 0;

				floorToInt(p.x / subBlockSize.x); /// Uh oh, this can't work
				int y = floorToInt(p.y / subBlockSize.y);

				Float tDeltaX = tDeltaXSingle * subBlockSize.x,
				      tDeltaY = tDeltaYSingle * subBlockSize.y,
					  tNextX  = nextMultiple(ray.d.x, p.x, subBlockSize.x),
					  tNextY  = nextMultiple(ray.d.y, p.y, subBlockSize.y);

				Float t = 0;
				while (t < maxt) {
					if (tNextX < tNextY) {
						tNextX += tDeltaX;
						x += iDeltaX;
					} else {
						tNextY += tDeltaY;
						y += iDeltaY;
					}

					

					stack[stackIdx+1].level = entry.level + 1;
					stack[stackIdx+1].x = x;
					stack[stackIdx+1].y = y;
					++stackIdx;
				}
			} else {
				cout << "X(" << entry.x + 1 << ", " << entry.y+1 << ") = 1;" << endl;
				--stackIdx;
			}
		}

		return false;
	}

#if 0
	struct PatchRecord {
		Point2 uv;
		Point p;
		int x, y;
	};

	void fillIntersectionRecord(const Ray &ray,
			const void *tmp, Intersection &its) const {
		PatchRecord &temp = *((PatchRecord *) tmp);

		int x = temp.x, y = temp.y, width = m_levelSize[0].x;
		Float
			f00 = m_data[x   + y     * width],
			f01 = m_data[x   + (y+1) * width],
			f10 = m_data[x+1 + y     * width],
			f11 = m_data[x+1 + (y+1) * width];

		its.uv = Point2((x+temp.uv.x) * m_invSize.x, (y+temp.uv.y) * m_invSize.y);
		its.p  = temp.p;

		its.dpdu = Vector(1, 0, (1.0f - its.uv.y) * (f10 - f00) + its.uv.y * (f11 - f01));
		its.dpdv = Vector(0, 1, (1.0f - its.uv.x) * (f01 - f00) + its.uv.x * (f11 - f10));
		its.geoFrame.n = cross(its.dpdu, its.dpdv);
		its.geoFrame.s = normalize(its.dpdu);
		its.geoFrame.t = normalize(its.dpdv);

		its.shFrame = its.geoFrame;
		its.shape = this;

 		its.wi = its.toLocal(-ray.d);
 		its.hasUVPartials = false;
		its.instance = NULL;
		its.time = ray.time;
	}
#endif

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Float t;
		return rayIntersect(ray, mint, maxt, t, NULL);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(Texture::m_theClass)) {
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
			bitmap->convert(m_data, Bitmap::ELuminance, Bitmap::EFloat);

			Log(EInfo, "Building acceleration data structure for %ix%i height field ..", m_dataSize.x, m_dataSize.y);

			ref<Timer> timer = new Timer();
			m_levelCount = (int) std::max(log2i((uint32_t) m_dataSize.x-1), log2i((uint32_t) m_dataSize.y-1)) + 1;

			m_levelSize = new Vector2i[m_levelCount];
			m_blockSizeF = new Vector2[m_levelCount];
			m_minmax = new Interval*[m_levelCount];

			m_levelSize[0]  = Vector2i(m_dataSize.x - 1, m_dataSize.y - 1);
			m_blockSizeF[0] = Vector2(1, 1);
			size = (size_t) m_levelSize[0].x * (size_t) m_levelSize[0].y * sizeof(Interval);
			m_minmax[0] = (Interval *) allocAligned(size);
			storageSize += size;

			/* Build the lowest MIP layer directly from the heightfield data */
			Interval *bounds = m_minmax[0];
			for (int y=0; y<m_levelSize[0].y; ++y) {
				for (int x=0; x<m_levelSize[0].x; ++x) {
					Float v00 = m_data[y * m_levelSize[0].x + x];
					Float v01 = m_data[y * m_levelSize[0].x + x + 1];
					Float v10 = m_data[(y + 1) * m_levelSize[0].x + x];
					Float v11 = m_data[(y + 1) * m_levelSize[0].x + x + 1];
					Float vmin = std::min(std::min(v00, v01), std::min(v10, v11));
					Float vmax = std::max(std::max(v00, v01), std::max(v10, v11));
					*bounds++ = Interval(vmin, vmax);
				}
			}

			/* Propagate height bounds upwards to the other layers */
			for (int level=1; level<m_levelCount; ++level) {
				Vector2i &cur  = m_levelSize[level],
				         &prev = m_levelSize[level-1];

				/* Calculate size of this layer */
				cur.x = prev.x > 1 ? (prev.x / 2) : 1;
				cur.y = prev.y > 1 ? (prev.y / 2) : 1;
				m_blockSizeF[level] = Vector2(
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
						const Interval &v00 = prevBounds[y0 * prev.x + x0], &v01 = prevBounds[y0 * prev.x + x1];
						const Interval &v10 = prevBounds[y1 * prev.x + x0], &v11 = prevBounds[y1 * prev.x + x1];
						Interval combined(v00);
						combined.expandBy(v01);
						combined.expandBy(v10);
						combined.expandBy(v11);
						*curBounds++ = combined;
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

	/* Height field data */
	Float *m_data;
	Vector2i m_dataSize;

	/* Min-max quadtree data */
	int m_levelCount;
	Vector2i *m_levelSize;
	Vector2 *m_blockSizeF;
	Interval **m_minmax;
};

MTS_IMPLEMENT_CLASS_S(Heightfield, false, Shape)
MTS_EXPORT_PLUGIN(Heightfield, "Height field intersection primitive");
MTS_NAMESPACE_END

