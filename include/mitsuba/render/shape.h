#if !defined(__SHAPE_H)
#define __SHAPE_H

#include <mitsuba/core/bsphere.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/serialization.h>
#include <mitsuba/render/records.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/subsurface.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract base class of all shapes
 */
class MTS_EXPORT_RENDER Shape : public ConfigurableObject {
public:
	/// Is this a compound shape consisting of several sub-objects?
	virtual bool isCompound() const;

	/// Return a sub-element of a compound shape
	virtual Shape *getElement(int i);

	/// Does this shape support primitive clipping?
	virtual bool isClippable() const;

	/**
	 * Returns the minimal axis-aligned bounding box of this shape 
	 * after it has clipped to the extends of another AABB.
	 * This is extremely useful to construct better kd-trees.
	 */
	virtual AABB getClippedAABB(const AABB &aabb) const;

	/// Fast ray intersection test (only calculates the intersection distance)
	virtual bool rayIntersect(const Ray &ray, Float start, Float end, Float &t) const;

	/// More detailed ray intersection test, which stores local geometry information in 'its'
	virtual bool rayIntersect(const Ray &ray, Intersection &its) const;

#if defined(MTS_SSE)
	/// Perform 4 simultaneous intersection tests using SSE
	virtual __m128 rayIntersectPacket(const RayPacket4 &packet, const
        __m128 mint, __m128 maxt, __m128 inactive, Intersection4 &its) const;
#endif

	/// Return a bounding sphere containing the shape
	inline const BSphere &getBSphere() const { return m_bsphere; }
	
	/// Return a bounding box containing the shape
	inline const AABB &getAABB() const { return m_aabb; }

	/**
	 * Sample a point on the shape - should be uniform
	 * wrt. surface area. Returns the associated probability
	 * density
	 */
	virtual Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const;

	/**
	 * Return the probability density of sampling the 
	 * given point using sampleArea()
	 */
	virtual Float pdfArea(const ShapeSamplingRecord &sRec) const;

	/**
	 * Sample a point on the shape - should be uniform
	 * wrt. solid angle. Returns the associated probability
	 * density
	 */
	virtual Float sampleSolidAngle(ShapeSamplingRecord &sRec, const Point &from, const Point2 &sample) const;

	/**
	 * Return the probability density of sampling the 
	 * given point using sampleSolidAngle()
	 */
	virtual Float pdfSolidAngle(const ShapeSamplingRecord &sRec, const Point &from) const;

	/// Return the shape's surface area
	inline Float getSurfaceArea() const { return m_surfaceArea; }

	/// Return the shape's BSDF
	inline const BSDF *getBSDF() const { return m_bsdf.get(); }
	/// Return the shape's BSDF
	inline BSDF *getBSDF() { return m_bsdf.get(); }

	/// Return the name of this shape
	inline const std::string &getName() const { return m_name; }

	/// Does this shape have a sub-surface integrator?
	inline bool hasSubsurface() const { return m_subsurface.get() != NULL; }
	/// Return the associated sub-surface integrator
	inline Subsurface *getSubsurface() { return m_subsurface; }
	/// Return the associated sub-surface integrator 
	inline const Subsurface *getSubsurface() const { return m_subsurface.get(); }

	/// Is this shape also an area luminaire?
	inline bool isLuminaire() const { return m_luminaire.get() != NULL; }
	/// Return the associated luminaire (if any)
	inline Luminaire *getLuminaire() { return m_luminaire; }
	/// Return the associated luminaire (if any)
	inline const Luminaire *getLuminaire() const { return m_luminaire.get(); }

	/// Called once after parsing
	virtual void configure();
	/// Serialize this shape to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Add a child (e.g. a luminaire/sub surface integrator) to this shape
	void addChild(const std::string &name, ConfigurableObject *child);

	MTS_DECLARE_CLASS()
protected:
	/// Create a new shape
	Shape(const Properties &props);

	/// Unserialize a shape
	Shape(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Shape();
protected:
	Transform m_worldToObject, m_objectToWorld;
	std::string m_name;
	AABB m_aabb;
	ref<BSDF> m_bsdf;
	ref<Subsurface> m_subsurface;
	ref<Luminaire> m_luminaire;
	Float m_surfaceArea, m_invSurfaceArea;
	BSphere m_bsphere;
};

MTS_NAMESPACE_END

#endif /* __SHAPE_H */
