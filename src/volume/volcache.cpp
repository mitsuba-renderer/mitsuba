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

#include <mitsuba/render/volume.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/lrucache.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

static StatsCounter statsHitRate("Volume cache", "Cache hit rate", EPercentage);
static StatsCounter statsCreate("Volume cache", "Block creations");
static StatsCounter statsDestruct("Volume cache", "Block destructions");
static StatsCounter statsEmpty("Volume cache", "Empty blocks", EPercentage);

/* Lexicographic ordering for Vector3i */
struct Vector3iKeyOrder : public std::binary_function<Vector3i, Vector3i, bool> {
	inline bool operator()(const Vector3i &v1, const Vector3i &v2) const {
		if (v1.x < v2.x) return true;
		else if (v1.x > v2.x) return false;
		if (v1.y < v2.y) return true;
		else if (v1.y > v2.y) return false;
		if (v1.z < v2.z) return true;
		else if (v1.z > v2.z) return false;
		return false;
	}
};

/*!\plugin{volcache}{Caching volume data source}
 * \parameters{
 *     \parameter{blockSize}{\Integer}{
 *         Size of the individual cache blocks
 *         \default{8, i.e. $8\times8\times 8$}
 *     }
 *     \parameter{voxelWidth}{\Float}{
 *         Width of a voxel (in a cache block) expressed in
 *         world-space units. \default{set to the ray marching
 *         step size of the nested medium}
 *     }
 *     \parameter{memoryLimit}{\Integer}{
 *         Maximum allowed memory usage in MiB. \default{1024, i.e. 1 GiB}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *         Optional linear transformation that should be applied
 *         to the volume data
 *     }
 *     \parameter{\Unnamed}{\Volume}{
 *         A nested volume data source
 *     }
 * }
 *
 * This plugin can be added between the renderer and another
 * data source, for which it caches all data lookups using a
 * LRU scheme. This is useful when the nested volume data source
 * is expensive to evaluate.
 *
 * The cache works by performing on-demand rasterization of subregions
 * of the nested volume into blocks ($8\times 8 \times 8$ by default).
 * These are kept in memory until a user-specifiable threshold is exeeded,
 * after which point a \emph{least recently used} (LRU) policy removes
 * records that haven't been accessed in a long time.
 */
class CachingDataSource : public VolumeDataSource {
public:
	typedef LRUCache<Vector3i, Vector3iKeyOrder, float *> BlockCache;

	CachingDataSource(const Properties &props)
		: VolumeDataSource(props) {
		/// Size of an individual block (must be a power of 2)
		m_blockSize = props.getInteger("blockSize", 8);

		if (!math::isPowerOfTwo(m_blockSize))
			Log(EError, "Block size must be a power of two!");

		/* Width of an individual voxel. Will use the step size of the
		   nested medium by default */
		m_voxelWidth = props.getFloat("voxelWidth", -1);

		/* Permissible memory usage in MiB. Default: 1GiB */
		m_memoryLimit = (size_t) props.getLong("memoryLimit", 1024) * 1024 * 1024;

		m_stepSizeMultiplier = (Float) props.getFloat("stepSizeMultiplier", 1.0f);

		m_volumeToWorld = props.getTransform("toWorld", Transform());
	}

	CachingDataSource(Stream *stream, InstanceManager *manager)
	: VolumeDataSource(stream, manager) {
		m_nested = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		configure();
	}

	virtual ~CachingDataSource() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);
		manager->serialize(stream, m_nested.get());
	}

	void configure() {
		if (m_nested == NULL)
			Log(EError, "A nested volume data source is needed!");
		m_aabb = m_nested->getAABB();
		if (!m_aabb.isValid())
			Log(EError, "Nested axis-aligned bounding box was invalid!");

		if (m_voxelWidth == -1)
			m_voxelWidth = m_nested->getStepSize();

		size_t memoryLimitPerCore = m_memoryLimit
			/ std::max((size_t) 1, Scheduler::getInstance()->getLocalWorkerCount());

		Vector totalCells  = m_aabb.getExtents() / m_voxelWidth;
		for (int i=0; i<3; ++i)
			m_cellCount[i] = (int) std::ceil(totalCells[i]);

		if (m_nested->supportsFloatLookups())
			m_channels = 1;
		else if (m_nested->supportsVectorLookups())
			m_channels = 1;
		else if (m_nested->supportsSpectrumLookups())
			m_channels = SPECTRUM_SAMPLES;
		else
			Log(EError, "Nested volume offers no access methods!");

		m_blockRes = m_blockSize+1;
		int blockMemoryUsage = (int) std::pow((Float) m_blockRes, 3) * m_channels * sizeof(float);
		m_blocksPerCore = memoryLimitPerCore / blockMemoryUsage;

		m_worldToVolume = m_volumeToWorld.inverse();
		m_worldToGrid = Transform::scale(Vector(1/m_voxelWidth))
			 * Transform::translate(-Vector(m_aabb.min)) * m_worldToVolume;
		m_voxelMask = m_blockSize-1;
		m_blockMask = ~(m_blockSize-1);
		m_blockShift = math::log2i((uint32_t) m_blockSize);

		Log(EInfo, "Volume cache configuration");
		Log(EInfo, "   Block size in voxels      = %i", m_blockSize);
		Log(EInfo, "   Voxel width               = %f", m_voxelWidth);
		Log(EInfo, "   Memory usage of one block = %s", memString(blockMemoryUsage).c_str());
		Log(EInfo, "   Memory limit              = %s", memString(m_memoryLimit).c_str());
		Log(EInfo, "   Memory limit per core     = %s", memString(memoryLimitPerCore).c_str());
		Log(EInfo, "   Max. blocks per core      = %i", m_blocksPerCore);
		Log(EInfo, "   Effective resolution      = %s", totalCells.toString().c_str());
		Log(EInfo, "   Effective storage         = %s", memString((size_t)
			(totalCells[0]*totalCells[1]*totalCells[2]*sizeof(float)*m_channels)).c_str());
	}

	Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		int x = (int) p.x, y = (int) p.y, z = (int) p.z;

		if (EXPECT_NOT_TAKEN(
			x < 0 || x >= m_cellCount.x ||
			y < 0 || y >= m_cellCount.y ||
			z < 0 || z >= m_cellCount.z))
			return 0.0f;

		BlockCache *cache = m_cache.get();
		if (EXPECT_NOT_TAKEN(cache == NULL)) {
			cache = new BlockCache(m_blocksPerCore,
				boost::bind(&CachingDataSource::renderBlock, this, _1),
				boost::bind(&CachingDataSource::destroyBlock, this, _1));
			m_cache.set(cache);
		}

#if defined(VOLCACHE_DEBUG)
		if (cache->isFull()) {
			/* For debugging: when the cache is full, dump locations
			   of all cache records into an OBJ file and exit */
			std::vector<Vector3i> keys;
			cache->get_keys(std::back_inserter(keys));

			std::ofstream os("keys.obj");
			os << "o Keys" << endl;
			for (size_t i=0; i<keys.size(); i++) {
				Vector3i key = keys[i];
				key = key * m_blockSize + Vector3i(m_blockSize/2);

				Point p(key.x * m_voxelWidth + m_aabb.min.x,
					key.y * m_voxelWidth + m_aabb.min.y,
					key.z * m_voxelWidth + m_aabb.min.z);

				os << "v " << p.x << " " << p.y << " " << p.z << endl;
			}

			/// Need to generate some fake geometry so that blender will import the points
			for (size_t i=3; i<=keys.size(); i++)
				os << "f " << i << " " << i-1 << " " << i-2 << endl;
			os.close();
			_exit(-1);
		}
#endif

		bool hit = false;
		float *blockData = cache->get(Vector3i(
			(x & m_blockMask) >> m_blockShift,
			(y & m_blockMask) >> m_blockShift,
			(z & m_blockMask) >> m_blockShift), hit);

		statsHitRate.incrementBase();
		if (hit)
			++statsHitRate;

		if (blockData == NULL)
			return 0.0f;

		const int x1 = x & m_voxelMask, y1 = y & m_voxelMask, z1 = z & m_voxelMask,
				x2 = x1 + 1, y2 = y1 + 1, z2 = z1 + 1;

		const float fx = (float) p.x - x, fy = (float) p.y - y, fz = (float) p.z - z,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		const float
			&d000 = blockData[(z1*m_blockRes + y1)*m_blockRes + x1],
			&d001 = blockData[(z1*m_blockRes + y1)*m_blockRes + x2],
			&d010 = blockData[(z1*m_blockRes + y2)*m_blockRes + x1],
			&d011 = blockData[(z1*m_blockRes + y2)*m_blockRes + x2],
			&d100 = blockData[(z2*m_blockRes + y1)*m_blockRes + x1],
			&d101 = blockData[(z2*m_blockRes + y1)*m_blockRes + x2],
			&d110 = blockData[(z2*m_blockRes + y2)*m_blockRes + x1],
			&d111 = blockData[(z2*m_blockRes + y2)*m_blockRes + x2];

		float result = ((d000*_fx + d001*fx)*_fy +
				 (d010*_fx + d011*fx)*fy)*_fz +
				((d100*_fx + d101*fx)*_fy +
				 (d110*_fx + d111*fx)*fy)*fz;

		return result;
	}

	Spectrum lookupSpectrum(const Point &_p) const {
		return Spectrum(0.0f);
	}

	Vector lookupVector(const Point &_p) const {
		return Vector(0.0f);
	}

	bool supportsFloatLookups() const {
		return m_nested->supportsFloatLookups();
	}

	bool supportsSpectrumLookups() const {
		return m_nested->supportsSpectrumLookups();
	}

	bool supportsVectorLookups() const {
		return m_nested->supportsVectorLookups();
	}

	Float getStepSize() const {
		return m_voxelWidth * m_stepSizeMultiplier;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(VolumeDataSource::m_theClass)) {
			Assert(m_nested == NULL);
			m_nested = static_cast<VolumeDataSource *>(child);
		} else {
			VolumeDataSource::addChild(name, child);
		}
	}

	float *renderBlock(const Vector3i &blockIdx) const {
		float *result = new float[m_blockRes*m_blockRes*m_blockRes];
		Point offset = m_aabb.min + Vector(
			blockIdx.x * m_blockSize * m_voxelWidth,
			blockIdx.y * m_blockSize * m_voxelWidth,
			blockIdx.z * m_blockSize * m_voxelWidth);

		int idx = 0;
		bool nonempty = false;
		for (int z = 0; z<m_blockRes; ++z) {
			for (int y = 0; y<m_blockRes; ++y) {
				for (int x = 0; x<m_blockRes; ++x) {
					Point p = offset + Vector((Float) x, (Float) y, (Float) z) * m_voxelWidth;
					float value = (float) m_nested->lookupFloat(p);
					result[idx++] = value;
					nonempty |= (value != 0);
				}
			}
		}

		++statsCreate;
		statsEmpty.incrementBase();

		if (nonempty) {
			return result;
		} else {
			++statsEmpty;
			delete[] result;
			return NULL;
		}
	}

	void destroyBlock(float *ptr) const {
		++statsDestruct;
		delete[] ptr;
	}

	Float getMaximumFloatValue() const {
		return m_nested->getMaximumFloatValue();
	}

	MTS_DECLARE_CLASS()
protected:
	ref<VolumeDataSource> m_nested;
	Transform m_volumeToWorld;
	Transform m_worldToVolume;
	Transform m_worldToGrid;
	Float m_voxelWidth;
	Float m_stepSizeMultiplier;
	size_t m_memoryLimit;
	size_t m_blocksPerCore;
	int m_channels;
	int m_blockSize, m_blockRes;
	int m_blockMask, m_voxelMask, m_blockShift;
	Vector3i m_cellCount;
	mutable ThreadLocal<BlockCache> m_cache;
};

MTS_IMPLEMENT_CLASS_S(CachingDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(CachingDataSource, "Caching data source");
MTS_NAMESPACE_END
