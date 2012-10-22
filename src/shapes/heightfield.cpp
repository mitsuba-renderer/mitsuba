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
#include <mitsuba/core/properties.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{heightfield}{Height field}
 * \order{12}
 *
 * Developed by Milo\^{s} Ha\^{s}an
 */
class Heightfield : public Shape {
private:
    inline Float getIJ(int i, int j) const {
        i = std::max(0, std::min(m_rows-1, i));
        j = std::max(0, std::min(m_cols-1, j));
        return m_data[i*m_cols + j];
    }

    inline Normal getNormalIJ(int i, int j) const {
        i = std::max(0, std::min(m_rows-1, i));
        j = std::max(0, std::min(m_cols-1, j));
        return m_normals[i*m_cols + j];
    }

    inline Float getXY(int x, int y) const {
        return getIJ(m_rows - y - 1, x);
    }

    inline Normal getNormalXY(int x, int y) const {
        return getNormalIJ(m_rows - y - 1, x);
    }

    inline Float getBilinear(Float x, Float y) const {
        // cout << "getBilinear: " << x << ' ' << y << ' ';
        int xi = (int) std::floor(x);
        int yi = (int) std::floor(y);
        Float u = x - xi, v = y - yi;
        Float z00 = getXY(xi, yi);
        Float z01 = getXY(xi, yi+1);
        Float z10 = getXY(xi+1, yi);
        Float z11 = getXY(xi+1, yi+1);
        Float result = (1-u)*(1-v)*z00 + (1-u)*v*z01 +  u*(1-v)*z10 +  u*v*z11;
        // cout << result << '\n';
        return result;
    }

    inline Vector2 getBilinearGradient(Float x, Float y) const {
        int xi = (int) std::floor(x);
        int yi = (int) std::floor(y);
        Float u = x - xi, v = y - yi;
        Float z00 = getXY(xi, yi);
        Float z01 = getXY(xi, yi+1);
        Float z10 = getXY(xi+1, yi);
        Float z11 = getXY(xi+1, yi+1);
        Float du = (1-v) * (z10 - z00) + v * (z11 - z01);
        Float dv = (1-u) * (z01 - z00) + u * (z11 - z10);
        return Vector2(du, dv);
    }

    inline Normal getInterpolatedNormal(Float x, Float y) const {
        int xi = (int) std::floor(x);
        int yi = (int) std::floor(y);
        Float u = x - xi, v = y - yi;
        Normal z00 = getNormalXY(xi, yi);
        Normal z01 = getNormalXY(xi, yi+1);
        Normal z10 = getNormalXY(xi+1, yi);
        Normal z11 = getNormalXY(xi+1, yi+1);
        Normal result = (1-u)*(1-v)*z00 + (1-u)*v*z01 +  u*(1-v)*z10 +  u*v*z11;
        return normalize(result);
    }

public:
	Heightfield(const Properties &props) : Shape(props) {
        m_pixelSize = props.getFloat("pixelSize", 0.01f);
        Float multiplier = props.getFloat("multiplier", 1);
        m_flipNormal = props.getBoolean("flipNormal", false);

        // read bitmap; use only red channel
        fs::path hfPath = props.getString("hf");
		ref<FileStream> stream = new FileStream(hfPath);
        ref<Bitmap> hf = new Bitmap(Bitmap::EOpenEXR, stream);
        m_rows = hf->getHeight();
        m_cols = hf->getWidth();
        Vector4* hfData = (Vector4*) hf->getFloatData();

        m_data = new float[m_rows * m_cols];
        for (int i = 0; i < m_rows * m_cols; i++)
            m_data[i] = hfData[i].x * multiplier;

        // compute bounding box
        Float xMax = (m_cols-1) * m_pixelSize / 2;
        Float yMax = (m_rows-1) * m_pixelSize / 2;
        Float zMin = std::numeric_limits<Float>::infinity(), zMax = -zMin;
        for (int i = 0; i < m_rows * m_cols; i++) {
            Float v = m_data[i];
            zMin = std::min(zMin, v);
            zMax = std::max(zMax, v);
        }
        zMin -= Epsilon; zMax += Epsilon;
        m_aabb = AABB(Point(-xMax, -yMax, zMin), Point(xMax, yMax, zMax));
        // cout << m_aabb.toString() << '\n';

        // compute normals
        m_normals = 0;
        if (props.getBoolean("interpolateNormals", false)) {
            m_normals = new Normal[m_rows * m_cols];
            int index = 0;

            for (int i = 0; i < m_rows; i++)
                for (int j = 0; j < m_cols; j++) {
                    Float dx = getIJ(i, j+1) - getIJ(i, j-1);
                    Float dy = getIJ(i-1, j) - getIJ(i+1, j); // reversal between i and y
                    Normal n(-dx / m_pixelSize, -dy / m_pixelSize, 2);
                    m_normals[index++] = normalize(n);
                }
        }
	}

	Heightfield(Stream *stream, InstanceManager *manager)
			: Shape(stream, manager) {
		Assert(false);
	}

	~Heightfield() {
	    delete[] m_data;
	    if (m_normals) delete[] m_normals;
    }

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		Assert(false);
	}

	AABB getAABB() const {
		return m_aabb;
	}

	Float getSurfaceArea() const {
		return 0;
	}

	/// find smallest t >= 0 such that a*t + b is integer
	inline static Float nextInt(Float a, Float b) {
		if (a == 0) return std::numeric_limits<Float>::infinity();
		else if (a > 0) return (std::ceil(b) - b) / a;
		else return (std::floor(b) - b) / a;
	}

    /// clip ray against bounding box
	bool clipRay(const Ray& ray, Float& mint, Float& maxt) const {
	    Float tmp1, tmp2;
	    if (!m_aabb.rayIntersect(ray, tmp1, tmp2)) return false;
        mint = std::max(mint, tmp1);
        maxt = std::min(maxt, tmp2);
        if (maxt <= mint) return false;
        return true;
	}

	/// transform from world to grid space
	Point world2grid(Point p) const {
        p.x /= m_pixelSize;
        p.y /= m_pixelSize;
        p.x += (m_cols - 1) / 2.0f;
        p.y += (m_rows - 1) / 2.0f;
        return p;
	}

    /// Subtract ray, fit quadratic function and solve.
    /// x0, xMid and x1 assumed to be 0, 0.5 and 1.
    /// Return -1 if no solution in [0,1].
    inline Float fitAndSolve(Float f0, Float fMid, Float f1, Float z0, Float z1) const {
        // subtract the "z" line of the ray
        f0 -= z0; f1 -= z1; fMid -= (z0 + z1) / 2;

        // fit quadratic function to values
        Float a = 2*f0 - 4*fMid + 2*f1;
        Float b = -3*f0 + 4*fMid - f1;
        Float c = f0;

        // solve and return the smaller root within [0,1] (if any)
        Float t0, t1;
        if (!solveQuadratic(a, b, c, t0, t1)) return -1;
        if (t0 >= 0 && t0 <= 1) return t0;
        if (t1 >= 0 && t1 <= 1) return t1;
        return -1;
    }

	/// find nearest intersection of ray with bilinear interpolant
	Float intersectPatch(Float x0, Float y0, Float z0,
	    Float x1, Float y1, Float z1) const {
        Float f0 = getBilinear(x0, y0);
        Float fMid = getBilinear((x0 + x1) / 2, (y0 + y1) / 2);
        Float f1 = getBilinear(x1, y1);
        return fitAndSolve(f0, fMid, f1, z0, z1);
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt, Float &tResult, void *tmp) const {
        if (!clipRay(ray, mint, maxt)) return false;

	    // ray segment endpoints in grid space
        Point start = world2grid(ray(mint));
        Point end = world2grid(ray(maxt));

        // the 2D vector from start to end of ray
        Vector d = end - start;
        if (d.lengthSquared() == 0) return false;

        // handle trivial case of perpendicular ray
        if (d.x == 0 && d.y == 0) {
            Assert(d.z != 0);
            Float z = getBilinear(start.x, start.y);

            if ((z - start.z) * (z - end.z) > 0) return false;

            Float t = (z - start.z) / (end.z - start.z);
            tResult = (1-t) * mint + t * maxt;
            return true;
        }

        // the 2D ray parameter and its deltas (can be infinite)
        Float t = 0;
        Float tMax = std::sqrt(d.x*d.x + d.y*d.y);
        Float tDeltaX = std::abs(tMax / d.x);
        Float tDeltaY = std::abs(tMax / d.y);
        d /= tMax;
        Float tMaxX = nextInt(d.x, start.x);
        Float tMaxY = nextInt(d.y, start.y);
        Float tNext;

        while (t < tMax) {
            if (tMaxX < tMaxY) {
                // advance in x
                tNext = std::min(tMaxX, tMax);
                tMaxX += tDeltaX;
            } else {
                // advance in y
                tNext = std::min(tMaxY, tMax);
                tMaxY += tDeltaY;
            }

            // get intersection point (shifted into [0,1])
            Point p = start + t * d;
            Point q = start + tNext * d;
            Float tPatch = intersectPatch(p.x, p.y, p.z, q.x, q.y, q.z);

            if (tPatch >= 0) {
                // found intersection, shift back from [0,1]
                tPatch = t + tPatch * (tNext - t);
                tResult = mint + tPatch * (maxt - mint) / tMax;
                return true;
            }

            t = tNext;
        }

		return false;
	}

	bool rayIntersect(const Ray &r, Float mint, Float maxt) const {
	    Float t;
		return rayIntersect(r, mint, maxt, t, 0);
	}

	void fillIntersectionRecord(const Ray &ray,
			const void *temp, Intersection &its) const {
        Float t = its.t;
        memset(&its, 0, sizeof(Intersection));

        // position
        its.t = t;
		its.p = ray(t);

		// normals
		Point pGrid = world2grid(its.p);
        Vector2 grad = getBilinearGradient(pGrid.x, pGrid.y);
		Normal gn = normalize(Vector(-grad.x, -grad.y, m_pixelSize));
		Normal sn = m_normals ? getInterpolatedNormal(pGrid.x, pGrid.y) : gn;

		// primitive uv-s
        its.uv = Point2(pGrid.x, pGrid.y);
        its.dpdu = Vector(m_pixelSize, 0, grad.x);
        its.dpdv = Vector(0, m_pixelSize, grad.y);

		if (m_flipNormal) { gn *= -1; sn *= -1; }
    	its.geoFrame = Frame(gn);
    	its.shFrame = Frame(sn);

 		its.wi = its.toLocal(-ray.d);
		its.shape = this;
	}

	void getNormalDerivative(const Intersection &its,
        Vector &dndu, Vector &dndv, bool shadingFrame) const {
    	if (!m_normals) {
    		dndu = Vector(0.0f);
    		dndv = Vector(0.0f);
            return;
        }

        Float x = its.uv.x, y = its.uv.y;
        int xi = (int) std::floor(x);
        int yi = (int) std::floor(y);
        Float u = x - xi, v = y - yi;
        Normal n00 = getNormalXY(xi, yi);
        Normal n01 = getNormalXY(xi, yi+1);
        Normal n10 = getNormalXY(xi+1, yi);
        Normal n11 = getNormalXY(xi+1, yi+1);

    	Normal N((1-u)*(1-v)*n00 + (1-u)*v*n01 +  u*(1-v)*n10 +  u*v*n11);
    	Float il = 1.0f / N.length(); N *= il;

        dndu = (1-v) * (n10 - n00) + v * (n11 - n01);
    	dndu *= il; dndu -= N * dot(N, dndu);
        dndv = (1-u) * (n01 - n00) + u * (n11 - n10);
    	dndv *= il; dndv -= N * dot(N, dndv);
    }

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return 1;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Heightfield[" << endl
			<< "  pixelSize = " << m_pixelSize << ", " << endl
			<< "  rows = " << m_rows << ", " << endl
			<< "  cols = " << m_cols << ", " << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
    AABB m_aabb;
    Float m_pixelSize;
    int m_rows, m_cols;
    float* m_data;
    Normal* m_normals;
    bool m_flipNormal;
};

MTS_IMPLEMENT_CLASS_S(Heightfield, false, Shape)
MTS_EXPORT_PLUGIN(Heightfield, "Height field intersection primitive");
MTS_NAMESPACE_END

