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

#include <mitsuba/render/mipmap3d.h>
#include <mitsuba/core/timer.h>

MTS_NAMESPACE_BEGIN

SparseMipmap3D::SparseMipmap3D(const AABB &aabb, size_t size, const float *data, 
		Float maxError, Float offset)  : m_aabb(aabb), m_size(size) {
	Assert(isPowerOfTwo(m_size));
	m_levels = 1 + log2i(m_size);
	m_aabbSum = Vector(m_aabb.min + m_aabb.max);

	float **pyramid = new float*[m_levels];

	// Stored efficiently due to a partial template specialization
	std::vector<bool> *bitPyramid = new std::vector<bool>[m_levels];

	Log(EDebug, "Pass 1: Building a dense %ix%ix%i mipmap", size, size, size);
	ref<Timer> timer = new Timer();

	size_t nEntries = size*size*size;
	pyramid[0] = new float[nEntries];
	bitPyramid[0].resize(nEntries);

	pyramid[0] = const_cast<float *>(data);
	for (size_t i=0; i<nEntries; ++i)
		bitPyramid[0][i] = true;

	for (size_t level=1; level<m_levels; ++level) {
		size_t res = size >> level, res2 = res << 1,
			   slab = res*res, slab2 = slab << 2;
		float *children = pyramid[level-1], *currentLevel;
		pyramid[level] = currentLevel = new float[res*res*res];
		std::vector<bool> &used = bitPyramid[level];
		used.resize(res*res*res);

		for (size_t z2=0, z=0; z<res; z2 = ++z << 1) {
			for (size_t y2=0, y=0; y<res; y2 = ++y << 1) {
				for (size_t x2=0, x=0; x<res; x2 = ++x << 1) {
					used[z*slab + y*res + x] = true;
					currentLevel[z*slab + y*res + x] = 0.125f * (
						children[ z2*slab2    +  y2    * res2 +  x2   ] +
						children[ z2*slab2    +  y2    * res2 + (x2+1)] +
						children[ z2*slab2    + (y2+1) * res2 +  x2   ] +
						children[ z2*slab2    + (y2+1) * res2 + (x2+1)] +
						children[(z2+1)*slab2 +  y2    * res2 +  x2   ] +
						children[(z2+1)*slab2 +  y2    * res2 + (x2+1)] +
						children[(z2+1)*slab2 + (y2+1) * res2 +  x2   ] +
						children[(z2+1)*slab2 + (y2+1) * res2 + (x2+1)]);
				}
			}
		}
	}
	Log(EDebug, "Pass 1: Done. (%i ms)", timer->getMilliseconds());
	timer->reset();

	Log(EDebug, "Pass 2: Compressing ...");
	size_t innerCount = 0, leafCount = 0, discardedCount = 0;
	for (size_t level=0; level<m_levels-1; ++level) {
		size_t res = size >> level, halfRes = res >> 1, res2 = res << 1, 
				slab = res*res, halfSlab = slab >> 2, slab2 = slab << 2;
		float *parentLevel = pyramid[level+1];
		for (size_t z2=0, z=0; z<res; z2 = ++z << 1) {
			for (size_t y2=0, y=0; y<res; y2 = ++y << 1) {
				for (size_t x2=0, x=0; x<res; x2 = ++x << 1) {
					if (level > 0) {
						const std::vector<bool> &childUsed = bitPyramid[level-1];
						/* To be able to discard this node, all of its descendants 
						   must be unused */
						if (childUsed[ z2*slab2    +  y2    * res2 +  x2   ] ||
							childUsed[ z2*slab2    +  y2    * res2 + (x2+1)] ||
							childUsed[ z2*slab2    + (y2+1) * res2 +  x2   ] ||
							childUsed[ z2*slab2    + (y2+1) * res2 + (x2+1)] ||
							childUsed[(z2+1)*slab2 +  y2    * res2 +  x2   ] ||
							childUsed[(z2+1)*slab2 +  y2    * res2 + (x2+1)] ||
							childUsed[(z2+1)*slab2 + (y2+1) * res2 +  x2   ] ||
							childUsed[(z2+1)*slab2 + (y2+1) * res2 + (x2+1)]) {
							innerCount++;	
							continue;
						}
					}

					/* Test whether discarding this node would produce an
					   excessive relative error */
					float parentValue = parentLevel[(z>>1)*halfSlab+(y>>1)*halfRes+(x>>1)];
					int support = 1<<level,
						xs = x * support, xe = xs + support,
						ys = y * support, ye = ys + support,
						zs = z * support, ze = zs + support,
						fullSlab = size*size;

					Float error = 0;
					for (int zp=zs; zp<ze; ++zp) {
						for (int yp=ys; yp<ye; ++yp) {
							for (int xp=xs; xp<xe; ++xp) {
								float referenceValue = data[zp*fullSlab + yp*size + xp];
								error = std::max(error,
									std::abs((parentValue-referenceValue)/(referenceValue+offset)));
								if (error >= maxError)
									goto stop; /* Argh :) */
							}
						}
					}
					stop:

					if (error >= maxError) {
						leafCount++;
						continue;
					}

					bitPyramid[level][z*slab+y*res + x] = false;
					++discardedCount;
				}
			}
		}
	}
	Log(EDebug, "Pass 2: Done. (%i ms)", timer->getMilliseconds());
	timer->reset();

	Log(EDebug, "Pass 3: Building an octree (%i inner nodes, %i leaf nodes, %i discarded, %i KiB/%i KiB) ...", 
		innerCount, leafCount, discardedCount, (leafCount+innerCount*sizeof(Node)) / 1024,
			(nEntries * sizeof(float)) / 1024);

	m_nodes.reserve(innerCount+1);

	int nodeID = build(m_levels-1, Point3i(0, 0, 0), pyramid, bitPyramid);
	Assert(nodeID == 0);
	
	Log(EDebug, "Pass 3: Done. (%i ms)", timer->getMilliseconds());

	for (size_t i=1; i<m_levels; ++i) 
		delete[] pyramid[i];
	delete[] pyramid;
	delete[] bitPyramid;
}

SparseMipmap3D::SparseMipmap3D(Stream *stream, InstanceManager *manager) {
	m_aabb = AABB(stream);
	m_size = (size_t) stream->readUInt();
	
	size_t nodeCount = stream->readSize();
	m_nodes.resize(nodeCount);
	for (size_t i=0; i<nodeCount; ++i) {
		stream->readIntArray(m_nodes[i].child, 8);
		m_nodes[i].value = stream->readFloat();
	}

	m_levels = 1 + log2i(m_size);
	m_aabbSum = Vector(m_aabb.min + m_aabb.max);
}

void SparseMipmap3D::serialize(Stream *stream, InstanceManager *manager) const {
	m_aabb.serialize(stream);
	stream->writeUInt(m_size);
	stream->writeSize(m_nodes.size());

	for (size_t i=0; i<m_nodes.size(); ++i) {
		stream->writeIntArray(m_nodes[i].child, 8);
		stream->writeFloat(m_nodes[i].value);
	}
}

uint32_t SparseMipmap3D::build(int level, const Point3i &p, float **pyramid, 
	std::vector<bool> *bitPyramid) {
	bool leaf = (level == 0);

	int x2 = p.x << 1, y2 = p.y << 1, z2 = p.z << 1, 
		res = m_size >> level, res2 = res << 1,
		slab = res*res, slab2 = slab << 2;

	Assert(p.x < res);
	Assert(p.y < res);
	Assert(p.z < res);

	if (level > 0) {
		const std::vector<bool> &childUsed = bitPyramid[level-1];
		if (!childUsed[ z2*slab2    +  y2    * res2 +  x2   ] &&
			!childUsed[ z2*slab2    +  y2    * res2 + (x2+1)] &&
			!childUsed[ z2*slab2    + (y2+1) * res2 +  x2   ] &&
			!childUsed[ z2*slab2    + (y2+1) * res2 + (x2+1)] &&
			!childUsed[(z2+1)*slab2 +  y2    * res2 +  x2   ] &&
			!childUsed[(z2+1)*slab2 +  y2    * res2 + (x2+1)] &&
			!childUsed[(z2+1)*slab2 + (y2+1) * res2 +  x2   ] &&
			!childUsed[(z2+1)*slab2 + (y2+1) * res2 + (x2+1)])
			leaf = true;
	}

	float value = pyramid[level][p.z * slab + p.y * res + p.x] ;

	if (leaf) {
		union {
			float floatVal;
			uint32_t intVal;
		} tmp;
		tmp.floatVal = value;

		/* Sacrifice 1 bit of accuracy to avoid creating an extra node */
		return (tmp.intVal >> 1) | 0x80000000;
	}

	if (m_nodes.size()+1 > m_nodes.capacity())
		Log(EError, "Overflow!");
	m_nodes.push_back(Node(value));
	int entry = m_nodes.size()-1;

	Node &node = m_nodes[entry];

	for (int i=0; i<8; ++i) 
		node.child[i] = build(level-1, 
			Point3i(x2 + ((i&4)>>2), y2 + ((i&2)>>1), z2 + (i&1)),
			pyramid, bitPyramid);

	return entry;
}

/*
  Octree traversal algorithm presented in
  "An Efficient Parametric Algorithm for Octree Traversal"
  by J. Revelles, C. Urena, M. Lastra
*/

/* Tables 1 and 2 */
inline int first_node(Float tx0, Float ty0, Float tz0, Float txm, Float tym, Float tzm) {
	uint8_t result = 0;

	if (tx0 > ty0 && tx0 > tz0) {
		if (tym < tx0) result |= 2;
		if (tzm < tx0) result |= 1;
	} else if (ty0 > tz0) {
		if (txm < ty0) result |= 4;
		if (tzm < ty0) result |= 1;
	} else {
		if (txm < tz0) result |= 4;
		if (tym < tz0) result |= 2;
	}

	return result;
}

/* Table 3 */
inline int new_node(Float t1, int a, Float t2, int b, Float t3, int c) {
	return ((t1 < t2 && t1 < t3) ? a : (t2 < t3 ? b : c));
}

Float SparseMipmap3D::lineIntegral(const Ray &r) const {
	Float length = r.d.length();

	Ray ray(r(r.mint), r.d/length, 0, (r.maxt-r.mint)*length, 0.0f);

	uint8_t a = 0;
	if (ray.d.x < 0) {
		ray.o.x = m_aabbSum.x-ray.o.x;
		ray.d.x = -ray.d.x;
		ray.dRcp.x = -ray.dRcp.x;
		a |= 4;
	}
	if (ray.d.y < 0) {
		ray.o.y = m_aabbSum.y-ray.o.y;
		ray.d.y = -ray.d.y;
		ray.dRcp.y = -ray.dRcp.y;
		a |= 2;
	}
	if (ray.d.z < 0) {
		ray.o.z = m_aabbSum.z-ray.o.z;
		ray.d.z = -ray.d.z;
		ray.dRcp.z = -ray.dRcp.z;
		a |= 1;
	}

	Float tx0 = (m_aabb.min.x-ray.o.x)*ray.dRcp.x,
		  ty0 = (m_aabb.min.y-ray.o.y)*ray.dRcp.y,
		  tz0 = (m_aabb.min.z-ray.o.z)*ray.dRcp.z,
		  tx1 = (m_aabb.max.x-ray.o.x)*ray.dRcp.x,
		  ty1 = (m_aabb.max.y-ray.o.y)*ray.dRcp.y,
		  tz1 = (m_aabb.max.z-ray.o.z)*ray.dRcp.z,
		  mint = std::max(std::max(tx0, ty0), tz0),
		  maxt = std::min(ray.maxt, std::min(std::min(tx1, ty1), tz1));

	if (mint >= maxt)
		return 0.0f;

	const QueryContext ctx(a, maxt);
	return lineIntegral(0, tx0, ty0, tz0, tx1, ty1, tz1, ctx);
}

bool SparseMipmap3D::invertLineIntegral(const Ray &r, Float desiredDensity,
		Float &accumDensity, Float &samplePos, Float &sampleDensity) const {
	Float length = r.d.length();

	Ray ray(r(r.mint), r.d/length, 0, (r.maxt-r.mint)*length, 0.0f);

	uint8_t a = 0;
	if (ray.d.x < 0) {
		ray.o.x = m_aabbSum.x-ray.o.x;
		ray.d.x = -ray.d.x;
		ray.dRcp.x = -ray.dRcp.x;
		a |= 4;
	}
	if (ray.d.y < 0) {
		ray.o.y = m_aabbSum.y-ray.o.y;
		ray.d.y = -ray.d.y;
		ray.dRcp.y = -ray.dRcp.y;
		a |= 2;
	}
	if (ray.d.z < 0) {
		ray.o.z = m_aabbSum.z-ray.o.z;
		ray.d.z = -ray.d.z;
		ray.dRcp.z = -ray.dRcp.z;
		a |= 1;
	}

	Float tx0 = (m_aabb.min.x-ray.o.x)*ray.dRcp.x,
		  ty0 = (m_aabb.min.y-ray.o.y)*ray.dRcp.y,
		  tz0 = (m_aabb.min.z-ray.o.z)*ray.dRcp.z,
		  tx1 = (m_aabb.max.x-ray.o.x)*ray.dRcp.x,
		  ty1 = (m_aabb.max.y-ray.o.y)*ray.dRcp.y,
		  tz1 = (m_aabb.max.z-ray.o.z)*ray.dRcp.z,
		  mint = std::max(std::max(tx0, ty0), tz0),
		  maxt = std::min(ray.maxt, std::min(std::min(tx1, ty1), tz1));

	if (mint >= maxt) {
		accumDensity = 0.0f;
		return false;
	}

	QueryContext ctx(a, maxt);
	ctx.samplePos = std::numeric_limits<Float>::infinity();
	ctx.sampleDensity = -std::numeric_limits<Float>::infinity();
	ctx.remaining = desiredDensity;

	invertLineIntegral(0, tx0, ty0, tz0, tx1, ty1, tz1, ctx);

	if (ctx.remaining == 0) {
		accumDensity = desiredDensity;
		sampleDensity = ctx.sampleDensity;
		samplePos = r.mint*length + ctx.samplePos;
		return true;
	} else {
		accumDensity = desiredDensity - ctx.remaining;
		return false;
	}
}

Float SparseMipmap3D::lineIntegral(int32_t idx,
	Float tx0, Float ty0, Float tz0,
	Float tx1, Float ty1, Float tz1, const QueryContext &ctx) const {

	/* Stop if we're outside of the ray bounds */
	if (tx1 < 0 || ty1 < 0 || tz1 < 0 || tx0 > ctx.maxt || ty0 > ctx.maxt || tz0 > ctx.maxt)
		return 0.0f;

	if (idx < 0) {
		/* Arrived at a leaf node */
		Float overlap = std::min(std::min(std::min(tx1, ty1), tz1), ctx.maxt)
			- std::max((Float) 0, std::max(std::max(tx0, ty0), tz0));
		union {
			float floatValue;
			int intValue;
		} tmp;
		tmp.intValue = idx << 1;
		return tmp.floatValue * overlap;
	}

	const Node &node = m_nodes[idx];
	const Float txm = .5f * (tx0+tx1),
				tym = .5f * (ty0+ty1),
				tzm = .5f * (tz0+tz1);
	unsigned int curChild = first_node(tx0, ty0, tz0, txm, tym, tzm);
	Float result = 0;

	do {
		switch (curChild) {
			case 0:
				result += lineIntegral(node.child[ctx.a], tx0, ty0, tz0, txm, tym, tzm, ctx);
				curChild = new_node(txm, 4, tym, 2, tzm, 1);
				break;

			case 1:
				result += lineIntegral(node.child[1^ctx.a], tx0, ty0, tzm, txm, tym, tz1, ctx);
				curChild = new_node(txm, 5, tym, 3, tz1, 8);
				break;

			case 2:
				result += lineIntegral(node.child[2^ctx.a], tx0, tym, tz0, txm, ty1, tzm, ctx);
				curChild = new_node(txm, 6, ty1, 8, tzm, 3);
				break;

			case 3:
				result += lineIntegral(node.child[3^ctx.a], tx0, tym, tzm, txm, ty1, tz1, ctx);
				curChild = new_node(txm, 7, ty1, 8, tz1, 8);
				break;

			case 4:
				result += lineIntegral(node.child[4^ctx.a], txm, ty0, tz0, tx1, tym, tzm, ctx);
				curChild = new_node(tx1, 8, tym, 6, tzm, 5);
				break;

			case 5:
				result += lineIntegral(node.child[5^ctx.a], txm, ty0, tzm, tx1, tym, tz1, ctx);
				curChild = new_node(tx1, 8, tym, 7, tz1, 8);
				break;

			case 6:
				result += lineIntegral(node.child[6^ctx.a], txm, tym, tz0, tx1, ty1, tzm, ctx);
				curChild = new_node(tx1, 8, ty1, 8, tzm, 7);
				break;

			case 7:
				result += lineIntegral(node.child[7^ctx.a], txm, tym, tzm, tx1, ty1, tz1, ctx);
				curChild = 8;
				break;
		}
	} while (curChild < 8);

	return result;
}

void SparseMipmap3D::invertLineIntegral(int32_t idx,
	Float tx0, Float ty0, Float tz0,
	Float tx1, Float ty1, Float tz1, QueryContext &ctx) const {

	/* Stop if we're outside of the ray bounds */
	if (tx1 < 0 || ty1 < 0 || tz1 < 0 || tx0 > ctx.maxt || ty0 > ctx.maxt || tz0 > ctx.maxt)
		return;

	if (idx < 0) {
		/* Arrived at a leaf node */
		Float overlap = std::min(std::min(std::min(tx1, ty1), tz1), ctx.maxt)
			- std::max((Float) 0, std::max(std::max(tx0, ty0), tz0));
		union {
			float floatValue;
			int intValue;
		} tmp;

		tmp.intValue = idx << 1;
		const Float dist = ctx.remaining/tmp.floatValue;

		if (EXPECT_NOT_TAKEN(dist <= overlap)) {
			ctx.samplePos = std::max(std::max(tx0, ty0), tz0) + dist;
			ctx.sampleDensity = tmp.floatValue;
			ctx.remaining = 0;
		} else {
			ctx.remaining -= tmp.floatValue * overlap;
		}
		return;
	}

	const Node &node = m_nodes[idx];
	const Float txm = .5f * (tx0+tx1),
				tym = .5f * (ty0+ty1),
				tzm = .5f * (tz0+tz1);
	unsigned int curChild = first_node(tx0, ty0, tz0, txm, tym, tzm);

	do {
		switch (curChild) {
			case 0:
				invertLineIntegral(node.child[ctx.a], tx0, ty0, tz0, txm, tym, tzm, ctx);
				curChild = new_node(txm, 4, tym, 2, tzm, 1);
				break;

			case 1:
				invertLineIntegral(node.child[1^ctx.a], tx0, ty0, tzm, txm, tym, tz1, ctx);
				curChild = new_node(txm, 5, tym, 3, tz1, 8);
				break;

			case 2:
				invertLineIntegral(node.child[2^ctx.a], tx0, tym, tz0, txm, ty1, tzm, ctx);
				curChild = new_node(txm, 6, ty1, 8, tzm, 3);
				break;

			case 3:
				invertLineIntegral(node.child[3^ctx.a], tx0, tym, tzm, txm, ty1, tz1, ctx);
				curChild = new_node(txm, 7, ty1, 8, tz1, 8);
				break;

			case 4:
				invertLineIntegral(node.child[4^ctx.a], txm, ty0, tz0, tx1, tym, tzm, ctx);
				curChild = new_node(tx1, 8, tym, 6, tzm, 5);
				break;

			case 5:
				invertLineIntegral(node.child[5^ctx.a], txm, ty0, tzm, tx1, tym, tz1, ctx);
				curChild = new_node(tx1, 8, tym, 7, tz1, 8);
				break;

			case 6:
				invertLineIntegral(node.child[6^ctx.a], txm, tym, tz0, tx1, ty1, tzm, ctx);
				curChild = new_node(tx1, 8, ty1, 8, tzm, 7);
				break;

			case 7:
				invertLineIntegral(node.child[7^ctx.a], txm, tym, tzm, tx1, ty1, tz1, ctx);
				curChild = 8;
				break;
		}

		if (EXPECT_NOT_TAKEN(ctx.remaining == 0))
			break;
	} while (curChild < 8);
}

std::string SparseMipmap3D::toString() const {
	std::ostringstream oss;
	oss << "SparseMipmap3D[" << endl
		<< "  res = " << m_size << "," << endl
		<< "  levels = " << m_levels << "," << endl
		<< "  aabb = " << m_aabb.toString() << "," << endl
		<< "  mem = " << m_nodes.size() * sizeof(Node) / 1024 << " KiB" << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(SparseMipmap3D, false, SerializableObject)
MTS_NAMESPACE_END
