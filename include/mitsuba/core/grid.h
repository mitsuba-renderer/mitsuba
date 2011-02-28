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

#if !defined(__GRID_H)
#define __GRID_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Uniform 3D grid for storing and manipulating arbitrary quantities
 *
 * \ingroup libcore
 */
template <typename ValueType> class Grid {
public:
	/// Construct a new grid with the given resolution (initialized to zero)
	Grid(const Vector3i &res, const AABB &aabb) : m_res(res), m_aabb(aabb) {
		m_slab = res.x*res.y;
		m_numCells = m_slab*res.z;
		m_cells = new ValueType[m_numCells];
		for (int i=0; i<3; ++i)
			m_cellWidth[i] = m_aabb.getExtents()[i] / (Float) res[i];
		clear();
	}

	/// Unserialize a grid from a binary data stream
	Grid(Stream *stream) {
		m_res = Vector3i(stream);
		m_aabb = AABB(stream);
		m_slab = m_res.x*m_res.y;
		m_numCells = m_slab*m_res.z;
		m_cells = new ValueType[m_numCells];
		stream->read(m_cells, sizeof(ValueType)*m_numCells);
		for (int i=0; i<3; ++i)
			m_cellWidth[i] = m_aabb.getExtents()[i] / (Float) m_res[i];
	}

	/// Serialize a grid to a binary data stream
	inline void serialize(Stream *stream) const {
		m_res.serialize(stream);
		m_aabb.serialize(stream);
		for (int i=0; i<3; ++i)
			stream->writeSingle((float) m_aabb.min[i]);
		for (int i=0; i<3; ++i)
			stream->writeSingle((float) m_aabb.max[i]);
		stream->write(m_cells, sizeof(ValueType)*m_numCells);
	}

	/// Reset everything to zero
	inline void clear() {
		memset(m_cells, 0, sizeof(ValueType) * m_numCells);
	}

	/// Add the values from another grid of identical shape and type
	inline void operator+=(const Grid<ValueType> &grid) {
		SAssert(grid.m_numCells == m_numCells);
		for (size_t i=0; i<m_numCells; ++i)
			m_cells[i] += grid.m_cells[i];
	}

	/// Multiply all entries by a constant
	inline void operator*=(Float value) {
		for (size_t i=0; i<m_numCells; ++i)
			m_cells[i] *= value;
	}

	/// Return the grid AABB
	inline const AABB &getAABB() const { return m_aabb; }

	/// Return the grid resolution
	inline const Vector3i &getResolution() const { return m_res; }

	/// Return the cell size
	inline const Vector &getCellWidth() const { return m_cellWidth; }

	/// Perform a lookup
	inline ValueType &operator()(int x, int y, int z) {
		return m_cells[x + y*m_res.x + z*m_slab];
	}

	/// Perform a lookup (const version)
	inline const ValueType &operator()(int x, int y, int z) const {
		return m_cells[x + y*m_res.x + z*m_slab];
	}
	/// Return a pointer to the underlying array
	inline ValueType *getData() const { return m_cells; }

	/// Return a string representation
	std::string toString() const {
		std::ostringstream oss;
		
		oss << m_aabb.min.x << " " << m_aabb.max.x << " " << m_res.x << endl;
		oss << m_aabb.min.y << " " << m_aabb.max.y << " " << m_res.y << endl;
		oss << m_aabb.min.z << " " << m_aabb.max.z << " " << m_res.z << endl;

		for (int z=0; z<m_res.z; ++z) 
			for (int y=0; y<m_res.y; ++y) 
				for (int x=0; x<m_res.x; ++x) 
					oss << m_cells[x + y*m_res.x + z*m_slab] << " ";
		return oss.str();
	}

	/// Modify the voxel containing a certain point
	void setValue(const Point &p, ValueType value) {
		/* Intersect with the voxel grid */
		if (!m_aabb.contains(p)) {
			std::ostringstream oss;
			oss << "The grid " << m_aabb.toString() << " does not contain the "
				"position " << p.toString() << "!" << endl;
			throw std::runtime_error(oss.str());
		}
		Vector pos = p - m_aabb.min;
		Vector3i dpos(
			std::max(0, std::min((int) (pos.x / m_cellWidth.x), m_res.x-1)),
			std::max(0, std::min((int) (pos.y / m_cellWidth.y), m_res.y-1)),
			std::max(0, std::min((int) (pos.z / m_cellWidth.z), m_res.z-1))
		);
		m_cells[dpos.x + dpos.y * m_res.x + dpos.z * m_slab] = value;
	}

	/// Apply the given functor to a grid cell
	template <typename Functor> void apply(const Point &p, const Functor &functor) {
		/* Intersect with the voxel grid */
		if (!m_aabb.contains(p)) {
			std::ostringstream oss;
			oss << "The grid " << m_aabb.toString() << " does not contain the "
				"position " << p.toString() << "!" << endl;
			throw std::runtime_error(oss.str());
		}
		Vector pos = p - m_aabb.min;
		Vector3i dpos(
			std::max(0, std::min((int) (pos.x / m_cellWidth.x), m_res.x-1)),
			std::max(0, std::min((int) (pos.y / m_cellWidth.y), m_res.y-1)),
			std::max(0, std::min((int) (pos.z / m_cellWidth.z), m_res.z-1))
		);
		int index = dpos.x + dpos.y * m_res.x + dpos.z * m_slab;
		m_cells[index] = functor(m_cells[index]);
	}

	/**
	 * Apply the given functor to a grid cell - don't throw an error if 
	 * the point is not part of the grid. Returns 'true' upon success
	 */
	template <typename Functor> bool applyIfContained(const Point &p, const Functor &functor) {
		if (!m_aabb.contains(p))
			return false;

		Vector pos = p - m_aabb.min;
		Vector3i dpos(
			std::max(0, std::min((int) (pos.x / m_cellWidth.x), m_res.x-1)),
			std::max(0, std::min((int) (pos.y / m_cellWidth.y), m_res.y-1)),
			std::max(0, std::min((int) (pos.z / m_cellWidth.z), m_res.z-1))
		);
		int index = dpos.x + dpos.y * m_res.x + dpos.z * m_slab;
		m_cells[index] = functor(m_cells[index]);
		return true;
	}

	/**
	 * \brief Rasterize a ray to the grid and apply the functor to 
	 * every traversed cell
	 */
	template <typename Functor> void rasterize(const Ray &ray, Functor &functor) {
		Float mint, maxt, t;

		/* Intersect with the voxel grid */
		if (!m_aabb.rayIntersect(ray, mint, maxt))
			return;

		/* Find the covered range in the ray space */
		mint = std::max(mint, ray.mint); maxt = std::max(mint, std::min(maxt, ray.maxt));
		if (mint == maxt)
			return;

		/* Compute the discrete coordinates of the first intersected voxel */
		Vector pos = ray(mint) - m_aabb.min;
		Vector3i dpos(
			std::max(0, std::min((int) (pos.x / m_cellWidth.x), m_res.x-1)),
			std::max(0, std::min((int) (pos.y / m_cellWidth.y), m_res.y-1)),
			std::max(0, std::min((int) (pos.z / m_cellWidth.z), m_res.z-1))
		);
		t = mint;

		/* Precompute useful traversal information */
		Vector corner1 = Vector(dpos.x * m_cellWidth.x, 
			dpos.y * m_cellWidth.y, dpos.z * m_cellWidth.z);
		Vector corner2 = Vector((dpos.x+1) * m_cellWidth.x, 
			(dpos.y+1) * m_cellWidth.y, (dpos.z+1) * m_cellWidth.z);
		Vector delta, next;
		Vector3i step, bounds;
		for (int i=0; i<3; ++i) {
			if (std::abs(ray.d[i]) < Epsilon) {
				delta[i] = 0; step[i] = 0; next[i] =
					std::numeric_limits<Float>::infinity();
			} else if (ray.d[i] > 0) {
				delta[i] = m_cellWidth[i]/ray.d[i];
				next[i] = mint + (corner2[i] - pos[i]) / ray.d[i];
				step[i] = 1;
				bounds[i] = m_res[i];
			} else {
				delta[i] = -m_cellWidth[i]/ray.d[i];
				next[i] = mint + (corner1[i] - pos[i]) / ray.d[i];
				step[i] = -1;
				bounds[i] = -1;
			}
		}

		/* Walk through the voxel grid (3D DDA) */
		while (true) {
			Float nextT = std::min(std::min(std::min(next.x, next.y), next.z), maxt);

			const int arrayIdx = dpos.x + dpos.y * m_res.x + dpos.z * m_slab;
			m_cells[arrayIdx] = functor(m_cells[arrayIdx], nextT - t);
			t = nextT;

			if (next.x <= next.y && next.x <= next.z) {
				if (next.x > maxt)
					break;
				dpos.x += step.x;
				if (dpos.x == bounds.x)
					break;
				next.x += delta.x;
			} else if (next.y <= next.x && next.y <= next.z) {
				if (next.y > maxt)
					break;
				dpos.y += step.y;
				if (dpos.y == bounds.y)
					break;
				next.y += delta.y;
			} else {
				if (next.z > maxt)
					break;
				dpos.z += step.z;
				if (dpos.z == bounds.z)
					break;
				next.z += delta.z;
			}
		}
	}

	/// Return a string representation (MATLAB format)
	std::string toStringMATLAB() const {
		std::ostringstream oss;

		oss << "v=reshape([";
		for (int z=0; z<m_res.z; ++z) // intentional order (matlab array conventions..)
			for (int x=0; x<m_res.x; ++x)
				for (int y=0; y<m_res.y; ++y)
					oss << m_cells[x + y*m_res.x + z*m_slab] << " ";

		oss << "], " << m_res.x << "," << m_res.y << "," << m_res.z << ");" << endl;
		oss << "[X Y Z]=meshgrid( ..." << endl;
		oss << "\tlinspace(" << m_aabb.min.x << ", " << m_aabb.max.x << ", " << m_res.x << "), ..." << endl;
		oss << "\tlinspace(" << m_aabb.min.y << ", " << m_aabb.max.y << ", " << m_res.y << "), ..." << endl;
		oss << "\tlinspace(" << m_aabb.min.z << ", " << m_aabb.max.z << ", " << m_res.z << ")  ..." << endl;
		oss << ");" << endl;
		return oss.str();
	}

	/// Do a lookup using trilinear interpolation
	ValueType lookup(const Point &_p) const {
		Point p  = Point(_p - m_aabb.min);
		
		Float x = p.x / m_cellWidth.x - .5f;
		Float y = p.y / m_cellWidth.y - .5f;
		Float z = p.z / m_cellWidth.z - .5f;

		const int i = (int) x, j = (int) y, k = (int) z, 
			  pos=i+j*m_res.x+k*m_slab;

		if (i < 0 || j < 0 || k < 0 || i >= m_res.x-1 || 
			j >= m_res.y || k >= m_res.z)
			return ValueType();

		const Float alpha = x-i,
					beta = y-j,
					gamma = z-k;

		const ValueType
					A1 = m_cells[pos],
					B1 = (i+1<m_res.x) ? m_cells[pos+1] : ValueType(0),
					C1 = (j+1<m_res.y) ? m_cells[pos+m_res.x] : ValueType(0),
					D1 = (i+1<m_res.x && j+1<m_res.y) ? m_cells[pos+m_res.x+1] : ValueType(0);

		ValueType A2, B2, C2, D2;
		if (k + 1 < m_res.z) {
			A2 = m_cells[pos+m_slab];
			B2 = (i+1<m_res.x) ? m_cells[pos+1+m_slab] : ValueType(0);
			C2 = (j+1<m_res.y) ? m_cells[pos+m_res.x+m_slab] : ValueType(0);
			D2 = (i+1<m_res.x && j+1<m_res.y) ? m_cells[pos+m_res.x+m_slab+1] : ValueType(0);
		} else {
			A2 = B2 = C2 = D2 = ValueType(0);
		}

		return  (A1 * ((1-alpha) * (1-beta))
			  +  B1 * (   alpha  * (1-beta))
			  +  C1 * ((1-alpha) *	beta)
			  +  D1 *	 alpha   *	beta) * (1-gamma)
			  + (A2 * ((1-alpha) * (1-beta))
			  +  B2 * (   alpha  * (1-beta))
			  +  C2 * ((1-alpha) *	beta)
			  +  D2 *	 alpha   *	beta) * gamma;
	}


	/// Release all memory
	~Grid() {
		delete[] m_cells;
	}
private:
	Vector3i m_res;
	AABB m_aabb;
	size_t m_numCells;
	size_t m_slab;
	Vector m_cellWidth;
	ValueType *m_cells;
};

MTS_NAMESPACE_END

#endif /* __GRID_H */
