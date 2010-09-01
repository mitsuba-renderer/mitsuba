#ifndef MITERSEG_H_
#define MITERSEG_H_

#include <mitsuba/render/shape.h>
#include <mitsuba/core/ref.h>

#include "hair.h"

MTS_NAMESPACE_BEGIN

class MiterHairSegment : public Shape {
private:
	ref<Hair> m_hair;  // the hair array to which this segment belongs
	int m_iv;  // the index of this hair segment within the parent hair array

public:

	MiterHairSegment(ref<Hair> parent, int index);

	MiterHairSegment(const Properties &props) : Shape(props) {
		Log(EError, "Miter Hair Segments cannot be created directly; they are created by Hair instances.");
	}

	MiterHairSegment(Stream *stream, InstanceManager *manager)
		: Shape(stream, manager) {
		AssertEx(false, "Hair Segments do not support serialization.");
	}

	bool isClippable() const {
		return true;
	}

	void configure();

	AABB getWorldAABB(Float start, Float end) const;

	AABB getClippedAABB(const AABB &aabb) const;

	bool rayIntersect(const Ray &_ray, Float start, Float end, Float &t) const;

	bool rayIntersect(const Ray &ray, Intersection &its) const;

#if defined(MTS_SSE)
	__m128 rayIntersectPacket(const RayPacket4 &, const __m128, __m128, __m128, Intersection4 &) const;
#endif

	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const;

	void serialize(Stream *stream, InstanceManager *manager) const {
		AssertEx(false, "Hair Segments do not support serialization.");
	}

	std::string toString() const;

	inline Point firstVertex() const { return m_hair->vertex(m_iv); }
	inline Point secondVertex() const { return m_hair->vertex(m_iv+1); }
	inline Vector tangent() const { return normalize(secondVertex() - firstVertex()); }

	inline bool prevSegmentExists() const { return !m_hair->vertexStartsFiber(m_iv); }
	inline Point prevVertex() const { return m_hair->vertex(m_iv-1); }
	inline Vector prevTangent() const { return normalize(firstVertex() - prevVertex()); }

	inline bool nextSegmentExists() const { return !m_hair->vertexStartsFiber(m_iv+2); }
	inline Point nextVertex() const { return m_hair->vertex(m_iv+2); }
	inline Vector nextTangent() const { return normalize(nextVertex() - secondVertex()); }

	inline Vector firstMiterNormal() const {
		if (prevSegmentExists())
			return normalize(prevTangent() + tangent());
		else
			return tangent();
	}

inline Vector secondMiterNormal() const {
		if (nextSegmentExists())
			return normalize(tangent() + nextTangent());
		else
			return tangent();
	}

	MTS_DECLARE_CLASS()
};

MTS_NAMESPACE_END


#endif /* MITERSEG_H_ */
