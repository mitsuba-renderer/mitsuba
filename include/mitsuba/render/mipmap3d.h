#if !defined(__MIPMAP3D_H)
#define __MIPMAP3D_H

#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/**
 * Sparse mipmap data structure based on an adaptive octree representation
 */
class MTS_EXPORT_RENDER SparseMipmap3D : public SerializableObject {
public:
	/**
	 * Construct a new mipmap from the given volume data
	 *
	 * @param aabb
	 *    Axis-aligned bounding box specifying the size of the volume region
	 * @param size
	 *    Size-length of the volume data in pixels (must be a power of two)
	 * @param maxError
	 *    Maximum permitted "relative error"
	 * @param offset
	 *    Fudge factor for the denominator, meant to avoid the relative
	 *    error singularity at zero
	 */
	SparseMipmap3D(const AABB &aabb, size_t size, const float *data, 
		Float maxError, Float offset);

	/// Unserialize from a binary data stream
	SparseMipmap3D(Stream *stream, InstanceManager *manager);

	/**
	 * Compute a line integral through the octree. Coordinates are
	 * expected to be in the same coordinate system as the provided AABB.
	 */
	Float lineIntegral(const Ray &ray) const;

	/**
	 * Invert a line integral through the octree
	 *
	 * @param ray
	 *    Specifies the ray, along which the integral should be inverted.
	 *    Coordinates are expected to be in the same coordinate system as the
	 *    provided AABB.
	 * @param desiredDensity
	 *    Try to integrate, until this much density has been accumulated
	 * @param accumDensity
	 *    When the inversion fails and the function returns false, this 
	 *    parameter records the actual amount of accumulated density.
	 * @param samplePos
	 *    Upon sucess, this parameter returns the sampled position.
	 * @param sampleDensity
	 *    Upon success, this parameter stores the density value at the 
	 *    right integration domain boundary.
	 */
	bool invertLineIntegral(const Ray &ray, Float desiredDensity,
			Float &accumDensity, Float &samplePos, Float &sampleDensity) const;

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	struct QueryContext {
		uint8_t a;  // ray directionality flag
		Float maxt;	
		Float remaining;
		Float samplePos;
		Float sampleDensity;

		inline QueryContext(uint8_t a, Float maxt)
		  : a(a), maxt(maxt) { }
	};

	struct Node {
		inline Node() { }
		inline Node(float value) : value(value) { }
		int32_t child[8];
		float value;
	};

	uint32_t build(int level, const Point3i &p, float **pyramid, 
		std::vector<bool> *bitPyramid);

	Float lineIntegral(int32_t idx,
		Float tx0, Float ty0, Float tz0,
		Float tx1, Float ty1, Float tz1, 
		const QueryContext &ctx) const;
	
	void invertLineIntegral(int32_t idx,
		Float tx0, Float ty0, Float tz0,
		Float tx1, Float ty1, Float tz1, 
		QueryContext &ctx) const;

	/// Virtual destructor
	virtual ~SparseMipmap3D() { }
private:
	AABB m_aabb;
	std::vector<Node> m_nodes;
	size_t m_size, m_levels;
	Vector m_aabbSum;
};

MTS_NAMESPACE_END

#endif /* __MIPMAP3D_H */
