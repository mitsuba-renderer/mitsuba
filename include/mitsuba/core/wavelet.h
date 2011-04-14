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

#if !defined(__WAVELET_H)
#define __WAVELET_H

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/serialization.h>

//#define USE_GOOGLE_SPARSE_HASHMAP 1
//#define USE_GOOGLE_DENSE_HASHMAP 1
//#define USE_STL_HASHMAP 1

#if defined(USE_GOOGLE_SPARSE_HASHMAP)
#include <google/sparse_hash_map>
#include <backward/hash_fun.h>
#elif defined(USE_GOOGLE_DENSE_HASHMAP)
#include <google/dense_hash_map>
#include <backward/hash_fun.h>
#elif defined(USE_STL_HASHMAP)
#include <tr1/unordered_map>
#endif

MTS_NAMESPACE_BEGIN

/**
 * \brief Performs non-standard 2D Haar wavelet transformations.
 *
 * Based on "Wavelets for computer graphics: A primer, part 1" by
 * Eric J. Stollnitz, Tony D. DeRose, and David H. Salesin
 * (IEEE Computer Graphics and Applications, May 1995)
 * \ingroup libcore
 */
class MTS_EXPORT_CORE Wavelet2D : public Object {
public:
	/**
	 * \brief Create a wavelet representation from a given bitmap.
	 *
	 * Only one color channel is supported for this encoding, so
	 * the desired channel must be selected using the `colorChannel'
	 * parameter.
	 */
	Wavelet2D(const Bitmap *pBitmap, int colorChannel = 0);
	
	/// Create the wavelet from a sparse representation
	Wavelet2D(const SparseWavelet2D *sw);

	/**
	 * \brief Turn the wavelet representation back into an image.
	 *
	 * Optionally, scale+offset factors can be supplied to
	 * map the bitmap to a desired brightness
	 */
	void decode(Bitmap *pBitmap, float pOffset = 0, float pScale = 1);

	/// Discard a given fraction of wavelet coefficients (in [0,1])
	void discard(Float fraction);

	/**
	 * \brief Discard components such that the relative L^2-error is below the 
	 * given bound. Returns the achieved compression ratio.
	 */
	Float compress(Float maxError);

	/// Turn the wavelet into a sparse representation
	SparseWavelet2D *toSparseWavelet() const;

	/// Return the wavelet's size
	inline size_t getSize() const { return m_size; }

	MTS_DECLARE_CLASS()
protected:
	/// Horizontal 1D decomposition step
	void forwardStepX(size_t y, size_t size);

	/// Vertical 1D decomposition step
	void forwardStepY(size_t x, size_t size);

	/// Horizontal 1D decoding step
	void backwardStepX(size_t y, size_t size);

	/// Vertical 1D decoding step
	void backwardStepY(size_t x, size_t size);

	/// Perform the non-standard wavelet decomposition
	void nonstandardDecomposition();

	/// Virtual destructor
	virtual ~Wavelet2D();
protected:
	/// \cond
	struct coeff_t {
		inline coeff_t() { } 
		inline coeff_t(size_t index, float value)
		 : index(index), value(std::abs(value)) {
		}
		inline bool operator<(const coeff_t &coeff) const {
			return value < coeff.value;
		}
		size_t index;
		float value;
	};
	/// \endcond
protected:
	float *m_data;
	float *m_temp;
	size_t m_size;
};

/**
 * \brief Implements the non-standard 3D wavelet transform using 
 * Haar basis functions.
 */
class MTS_EXPORT_CORE Wavelet3D : public Object {
public:
	/**
	 * Create a wavelet representation of the given volume data.
	 * 'resolution' specifies the side-length of the cube, which
	 * must be a power of two.
	 */
	Wavelet3D(const float *data, size_t resolution);

	/**
	 * \brief Turn the wavelet representation back into a dense format
	 */
	void decode(float *target);

	/// Discard a given fraction of wavelet coefficients (in [0,1])
	void discard(Float fraction);

	/**
	 * \brief Discard components such that the relative L^2-error is below the 
	 * given bound. Returns the achieved compression ratio.
	 */
	Float compress(Float maxRelError);

	/// Turn the wavelet into a sparse octree representation
	SparseWaveletOctree *toOctree() const;

	/// Return the side length of the encoded cube
	inline size_t getSize() const { return m_size; }

	MTS_DECLARE_CLASS()
protected:
	/// Forward transform steps
	void forwardStepX(size_t y, size_t z, size_t size);
	void forwardStepY(size_t x, size_t z, size_t size);
	void forwardStepZ(size_t x, size_t y, size_t size);

	/// Backward transform steps
	void backwardStepX(size_t y, size_t z, size_t size);
	void backwardStepY(size_t x, size_t z, size_t size);
	void backwardStepZ(size_t x, size_t y, size_t size);

	/// Perform the non-standard tensor product wavelet transform
	void nonstandardDecomposition();

	/// Virtual destructor
	virtual ~Wavelet3D();
protected:
	/// \cond
	struct coeff_t {
		inline coeff_t() { } 
		inline coeff_t(size_t index, float value)
		 : index(index), value(std::abs(value)) {
		}

		inline bool operator<(const coeff_t &coeff) const {
			return value < coeff.value;
		}

		size_t index;
		float value;
	};
	/// \endcond
protected:
	float *m_data;
	float *m_temp;
	size_t m_size, m_slab;
};

/**
 * \brief Sparse 2D wavelet representation using the Haar basis
 * \ingroup libcore
 */
class MTS_EXPORT_CORE SparseWavelet2D : public SerializableObject {
public:
	/// 2D wavelet coefficient index for \ref SparseWavelet2D
	struct Key {
		uint16_t empty;
		uint8_t level; // level in the hierarchy
		uint8_t type;  // wavelet type (0..2)
		uint16_t i, j; // horizontal and vertical offset

		/// Create a new wavelet key
		inline static Key create(uint8_t level, uint8_t type, uint16_t i, uint16_t j) {
			Key key;
			key.empty = 0;
			key.level = level;
			key.type = type;
			key.i = i;
			key.j = j;
			return key;
		}

		/// Unpack a wavelet key from a 64bit representation
		inline static Key unpack(uint64_t packed) {
			union {
				Key a;
				uint64_t b;
			};

			b = packed;
			return a;
		}

		/// Turn a wavelet key into a 64bit representation
		inline uint64_t pack() const {
			union {
				Key a;
				uint64_t b;
			};

			a = *this;
			return b;
		}

		/// Return the sign of a quadrant within this wavelet
		inline float quadrantSign(int x, int y) const {
			const float signs[3][4] = {
				{1, -1, 1, -1}, /* Horizontal differencing */
				{1, 1, -1, -1}, /* Vertical differencing */
				{1, -1, -1, 1}  /* Diagonal differencing */
			};

			return signs[type][x+y*2];
		}

		/// Return a string representation
		inline std::string toString() const {
			std::ostringstream oss;
			oss << "Key[level=" << (int) level << ", type=" << (int) type 
				<< ", i=" << (int) i << ", j=" << (int) j << "]";
			return oss.str();
		}
	};

#if defined(USE_GOOGLE_SPARSE_HASHMAP)
	typedef google::sparse_hash_map<uint64_t, float> CoefficientMap;
#elif defined(USE_GOOGLE_DENSE_HASHMAP)
	typedef google::dense_hash_map<uint64_t, float> CoefficientMap;
#else
	typedef std::map<uint64_t, float, std::less<uint64_t> > CoefficientMap;
#endif
	typedef CoefficientMap::const_iterator CoefficientIterator;

public:
	/// Construct a sparse wavelet representation with all-zero coefficients
	SparseWavelet2D(size_t size);

	/// Copy constructor
	SparseWavelet2D(const SparseWavelet2D *sw);

	/// Unserialize from a binary data stream
	SparseWavelet2D(Stream *stream, InstanceManager *manager);

	/// Set the value of the scaling function
	inline void setScalingFunction(float value) { m_scalingFunction = value; }

	/// Return the value of the scaling function
	inline float getScalingFunction() const { return m_scalingFunction; }

	/// Return the side length of the encoded data
	inline size_t getSize() const { return m_size; }

	/// Write one of the wavelet coefficients
	inline void put(const Key &key, float value) {
		m_data[key.pack()] = value;
	}

	/// Read one of the wavelet coefficients
	inline float get(const Key &key) const {
		CoefficientIterator it = m_data.find(key.pack());
		if (it == m_data.end())
			return 0.0f;
		return (*it).second;
	}

	/// Evaluate the sparse representation at the given pixel position
	Float getPixel(const Point2i &pt) const;

	/**
	 * \brief Compute a line integral in 2D wavelet space.
	 *
	 * Coordinates are expected as pixel coordinates in [0,0]-[size,size],
	 * but are allowed to be fractional
	 */
	Float lineIntegral(Point2 start, Point2 end) const;

	/// Set the whole function to zero
	void clear();

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~SparseWavelet2D() { }
protected:
	CoefficientMap m_data;
	float m_scalingFunction;
	size_t m_size;
	int m_maxLevel;
};

/**
 * \brief Sparse 3D wavelet representation using the Haar basis and an 
 * octree structure
 * \ingroup libcore
 */
class MTS_EXPORT_CORE SparseWaveletOctree : public Object {
public:
	/// Construct a sparse wavelet representation with all-zero coefficients
	SparseWaveletOctree(size_t size, float scalingFunction);

	/// Return the side length of the encoded data
	inline size_t getSize() const { return m_size; }

	/// Set one of the wavelet coefficients (all 7 types at once)
	void put(int level, int i, int j, int k, float coeff[7]);

	/**
	 * Compute a line integral in Octree wavelet space. Coordinates are
	 * expected as pixel coordinates in [0,0,0]-[size,size,size],
	 * but are allowed to be fractional
	 */
	Float lineIntegral(Point start, Point end) const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~SparseWaveletOctree() { }
private:
	struct Node {
		inline Node(float value) : value(value) {
			for (int i=0; i<8; ++i)
				child[i] = 0;
		}

		int32_t child[8];
		float value;
	};

	Float lineIntegral(int32_t idx,
		Float tx0, Float ty0, Float tz0,
		Float tx1, Float ty1, Float tz1, uint8_t a) const;
private:
	Node *m_root;
	std::vector<Node> m_nodes;
	size_t m_size;
	int m_maxLevel;
};

MTS_NAMESPACE_END

#endif /* __WAVELET_H */
