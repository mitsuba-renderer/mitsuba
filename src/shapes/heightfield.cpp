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

MTS_NAMESPACE_BEGIN

/*!\plugin{heightfield}{Height field}
 * \order{12}
 *
 * Developed by Milo\^{s} Ha\^{s}an
 */
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
		}
    }

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
	}

	AABB getAABB() const {
		AABB aabb(
			Point3(0, 0, m_minmax[m_levelCount-1][0].min),
			Point3(m_levelSize[0].x, m_levelSize[0].y, m_minmax[m_levelCount-1][0].max)
		);

		AABB result;
		for (int i=0; i<8; ++i)
			result.expandBy(m_objectToWorld(aabb.getCorner(i)));

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
		int level;
		int x, y;
	};

#if 0
	bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *tmp) const {
		StackEntry stack[MTS_KD_MAXDEPTH];
		Ray ray;
		m_objectToWorld(_ray, ray);
		int idx = 0;

		stack[idx].level = m_levelCount-1;
		stack[idx].x = 0;
		stack[idx].y = 0;

		while (idx >= 0) {
			const StackEntry &entry    = stack[idx];
			const Vector2i &size       = m_levelSize[entry.level];
			const Vector2i &blockSize  = m_levelBlockSize[entry.level];
			const Interval &interval   = m_minmax[entry.level][entry.x + entry.y * m_levelSize[entry.level].x];

			AABB aabb(
				Point3(entry.x * blockSize.x, entry.y * blockSize.y, interval.min),
				Point3((entry.x + 1) * blockSize.x, (entry.y + 1) * blockSize.y, interval.max)
			);

			Float nearT, farT;
			bool match = aabb.rayIntersect(ray, nearT, farT);

			if (!match || farT < mint || nearT > maxt) {
				--idx;
				continue;
			}

			if (level == 0) {
				t = mint;
				return true;
			} else {
			}
		}

		return false;
	}
#endif

	// =============================================================
	//  Ray <-> Bilinear patch intersection computation
	//  Based on the paper "Ray Bilinear Patch Intersections"
    //  by Ramsey et al., JGT Volume 9, 3, 2004
	// =============================================================

	struct PatchRecord {
		Point2 uv;
		Point p;
		int x, y;
	};

	inline Float computeU(Float v, Float A1, Float A2, Float B1, Float B2,
			Float C1, Float C2, Float D1, Float D2) const {
		Float a = v*A2+B2;
		Float b = v*(A2-A1)+B2-B1;

		if (std::abs(b) >= std::abs(a))
			return (v*(C1-C2)+D1-D2)/b;
		else
			return -(v*C2+D2)/a;
	}


	bool rayIntersectPatch(const Ray &ray, int x, int y, Float mint, Float maxt, Float &t, void *tmp) const {
		int width = m_levelSize[0].x;

		const Float
			f00 = m_data[x   + y     * width],
			f01 = m_data[x   + (y+1) * width],
			f10 = m_data[x+1 + y     * width],
			f11 = m_data[x+1 + (y+1) * width];

		/* Variables for substitution */
		const Float
			a  = f01 + f10 - f11 - f00,
			b  = f00 - f10,
			c  = f00 - f01,
		    dx = x - ray.o.x,
		    dy = y - ray.o.y,
		    dz = f00 - ray.o.z,
			A1 = a*ray.d.x,
			A2 = a*ray.d.y,
			B1 = ray.d.z + b*ray.d.x,
			B2 = b*ray.d.y,
			C1 = c*ray.d.x,
			C2 = ray.d.z + c*ray.d.y,
			D1 = dx*ray.d.z - dz*ray.d.x,
			D2 = dy*ray.d.z - dz*ray.d.y;

		/* Quadratic equation coefficients */
		Float A = A2*C1 - A1*C2;
		Float B = A2*D1 - A1*D2 + B2*C1 - B1*C2;
		Float C = B2*D1 - B1*D2;

		const Float inf = std::numeric_limits<Float>::infinity();
		Float v0, v1;
		if (!solveQuadratic(A, B, C, v0, v1))
			return false;

		Float u0=inf, u1=inf, t0=inf, t1=inf;
		Point p0, p1;

		/* Validate & find coordinates of first potential intersection */
		if (v0 >= 0 && v0 <= 1) {
			u0 = computeU(v0, A1, A2, B1, B2, C1, C2, D1, D2);
			if (u0 >= 0 && u0 <= 1) {
				p0 = Point(x+u0, y+v0, (1.0f - u0) * ((1.0f - v0) * f00 + v0 * f01)
					+ u0 * ((1.0f - v0) * f10 + v0 * f11));
				t0 = dot(p0-ray.o, ray.d);
				if (t0 < mint || t0 > maxt)
					t0 = inf;
			}
		}

		/* Validate & find coordinates of second potential intersection */
		if (v1 >= 0 && v1 <= 1 && v1 != v0) {
			u1 = computeU(v1, A1, A2, B1, B2, C1, C2, D1, D2);
			if (u1 >= 0 && u1 <= 1) {
				p1 = Point(x+u1, y+v1, (1.0f - u1) * ((1.0f - v1) * f00 + v1 * f01)
					+ u1 * ((1.0f - v1) * f10 + v1 * f11));
				t1 = dot(p1-ray.o, ray.d);
				if (t1 < mint || t1 > maxt)
					t1 = inf;
			}
		}

		if (t0 == inf && t1 == inf)
			return false;

		if (tmp) {
			PatchRecord &temp = *((PatchRecord *) tmp);
			if (t0 <= t1) {
				t = t0;
				temp.uv = Point2(u0, v0);
				temp.p  = p0;
				temp.x  = x;
				temp.y  = y;
			} else {
				t = t1;
				temp.uv = Point2(u1, v1);
				temp.p  = p1;
				temp.x  = x;
				temp.y  = y;
			}
		}

		return true;
	}

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

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &t, void *tmp) const {
		for (int y=0; y<m_levelSize[0].y-1; ++y) {
			for (int x=0; x<m_levelSize[0].x-1; ++x) {
				if (rayIntersectPatch(ray, x, y, mint, maxt, t, tmp))
					return true;
			}
		}
		return false;
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		Float t;
		return rayIntersect(ray, mint, maxt, t, NULL);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(Texture::m_theClass)) {
			ref<Bitmap> bitmap = static_cast<Texture *>(child)->getBitmap(m_sizeHint);

			Vector2i size = bitmap->getSize();
			if (size.x < 2) size.x = 2;
			if (size.y < 2) size.y = 2;
			if (!isPowerOfTwo(size.x)) size.x = (int) roundToPowerOfTwo((uint32_t) size.x);
			if (!isPowerOfTwo(size.y)) size.y = (int) roundToPowerOfTwo((uint32_t) size.y);
			m_sizeF = Vector2(size);
			m_invSizeF = Vector2(1/size.x, 1/size.y);

			if (bitmap->getSize() != size) {
				Log(EWarn, "Heightfield texture size should be at least 2x2 and "
					"a power of two. Resampling from %ix%i to %ix%i..",
					bitmap->getWidth(), bitmap->getHeight(), size.x, size.y);

				bitmap = bitmap->resample(NULL, ReconstructionFilter::EClamp,
					ReconstructionFilter::EClamp, size,
					-std::numeric_limits<Float>::infinity(),
					std::numeric_limits<Float>::infinity());
			}

			m_data = (Float *) allocAligned((size_t) size.x * (size_t) size.y * sizeof(Float));
			bitmap->convert(m_data, Bitmap::ELuminance, Bitmap::EFloat);

			Log(EInfo, "Building acceleration data structure for %ix%i height field ..", size.x, size.y);

			ref<Timer> timer = new Timer();
			m_levelCount = (int) std::max(log2i((uint32_t) size.x), log2i((uint32_t) size.y)) + 1;

			m_levelSize = new Vector2i[m_levelCount];
			m_minmax = new Interval*[m_levelCount];
			m_levelSize[0] = size;
			m_minmax[0] = NULL;
			size_t totalStorage = (size_t) size.x * (size_t) size.y * sizeof(float);
			for (int i=0; m_levelSize[i].x > 1 || m_levelSize[i].y > 1; ++i) {
				m_levelSize[i+1].x = m_levelSize[i].x > 1 ? (m_levelSize[i].x / 2) : 1;
				m_levelSize[i+1].y = m_levelSize[i].y > 1 ? (m_levelSize[i].y / 2) : 1;
				size_t size = (size_t) m_levelSize[i+1].x * (size_t) m_levelSize[i+1].y * sizeof(Interval);
				m_minmax[i+1] = (Interval *) allocAligned(size);
				totalStorage += size;
			}

			/* Build the lowest MIP layer directly from the heightfield data */
			Interval *bounds = m_minmax[1];
			for (int y=0; y<m_levelSize[1].y; ++y) {
				int y0 = std::min(2*y,   m_levelSize[0].y-1),
					y1 = std::min(2*y+1, m_levelSize[0].y-1);
				for (int x=0; x<m_levelSize[1].x; ++x) {
					int x0 = std::min(2*x,   m_levelSize[0].x-1),
						x1 = std::min(2*x+1, m_levelSize[0].x-1);
					Float v00 = m_data[y0 * m_levelSize[0].x + x0];
					Float v01 = m_data[y0 * m_levelSize[0].x + x1];
					Float v10 = m_data[y1 * m_levelSize[0].x + x0];
					Float v11 = m_data[y1 * m_levelSize[0].x + x1];
					Float vmin = std::min(std::min(v00, v01), std::min(v10, v11));
					Float vmax = std::max(std::max(v00, v01), std::max(v10, v11));
					*bounds++ = Interval(vmin, vmax);
				}
			}

			/* Propagate height bounds upwards to the other layers */
			for (int level=2; level<m_levelCount; ++level) {
				const Vector2i &cur  = m_levelSize[level],
					           &prev = m_levelSize[level-1];
				Interval *curBounds  = m_minmax[level];
				Interval *prevBounds = m_minmax[level-1];
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
			Log(EInfo, "Done (took %i ms, %s)", timer->getMilliseconds(),
					memString(totalStorage).c_str());
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
	Vector2i m_sizeHint, *m_levelSize, *m_levelBlockSize;
	Vector2 m_sizeF, m_invSizeF;
	int m_levelCount;
	Float *m_data;
	Interval **m_minmax;
	Transform m_objectToWorld;
};

MTS_IMPLEMENT_CLASS_S(Heightfield, false, Shape)
MTS_EXPORT_PLUGIN(Heightfield, "Height field intersection primitive");
MTS_NAMESPACE_END

