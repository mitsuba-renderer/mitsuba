#if !defined(__BSPHERE_H)
#define __BSPHERE_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Bounding sphere class
 */
struct BSphere {
	Point center;
	Float radius;

	/// Construct an empty bounding sphere
	inline BSphere() {
		radius = 0.0f;
	}
	
	/// Unserialize an AABB from a stream
	inline BSphere(Stream *stream) {
		center = Point(stream);
		radius = stream->readFloat();
	}

	/// Create a bounding sphere from a given center point and a radius
	inline BSphere(const Point &pCenter, Float pRadius)
		: center(pCenter), radius(pRadius) {
	}

	/// Copy-constructor
	inline BSphere(const BSphere &boundingSphere) 
		: center(boundingSphere.center), radius(boundingSphere.radius) {
	}

	/// Return whether this bounding sphere is empty
	inline bool isEmpty() const {
		return radius <= 0.0f;
	}

	/** \brief Expands the bounding sphere to contain another point.
	 * Does not move the center point
	 */
	inline void expandBy(const Point p) {
		Vector dir = p - center;
		radius = std::max(radius, (p-center).length());
	}

	/// Comparison operator
	inline bool operator==(const BSphere &boundingSphere) const {
		return center == boundingSphere.center && radius == boundingSphere.radius;
	}

	/// Comparison operator
	inline bool operator!=(const BSphere &boundingSphere) const {
		return !operator==(boundingSphere);
	}

	/// Calculate the intersection points with the given ray
	inline bool rayIntersect(const Ray &ray, Float &nearHit, Float &farHit) const {
		Vector originToCenter = center - ray.o;
		Float distToRayClosest = dot(originToCenter, ray.d);
		Float tmp1 = originToCenter.lengthSquared() - radius*radius;

		if (tmp1 <= 0.0f) {
			/* Inside the sphere */
			nearHit = farHit = 
				std::sqrt(distToRayClosest * distToRayClosest - tmp1)
					+ distToRayClosest;
			return true;
		}

		/* Points in different direction */
		if (distToRayClosest < 0.0f)
			return false;

		Float sqrOriginToCenterLength = originToCenter.lengthSquared();
		Float sqrHalfChordDist = radius * radius - sqrOriginToCenterLength
			+ distToRayClosest * distToRayClosest;

		if (sqrHalfChordDist < 0) // Miss
			return false;

		// Hit
		Float hitDistance = std::sqrt(sqrHalfChordDist);
		nearHit = distToRayClosest - hitDistance;
		farHit = distToRayClosest + hitDistance;

		if (nearHit == 0)
			nearHit = farHit;

		return true;
	}

	/// Serialize this AABB to a stream
	inline void serialize(Stream *stream) const {
		center.serialize(stream);
		stream->writeFloat(radius);
	}

	/// Returns a string representation of the bounding sphere
	inline std::string toString() const {
		std::ostringstream oss;
		oss << "BSphere[center = " << center.toString()
			<< ", radius = " << radius << "]";
		return oss.str();
	}
};

MTS_NAMESPACE_END

#endif /* __BSPHERE_H */
