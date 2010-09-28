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

#if !defined(__SHAPE_H)
#define __SHAPE_H

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/** \brief Data record, which holds sampling-related information
 *  for a shape.
 */
struct MTS_EXPORT_RENDER ShapeSamplingRecord {
public:
	/// Create a sampling record (does no initialization!)
	inline ShapeSamplingRecord() { }

	/// Create a sampling record with the specified position and normal
	inline ShapeSamplingRecord(const Point &p, const Normal &n)
		: p(p), n(n) { }

	/// Return a string representation
	std::string toString() const;
public:
	/// Sampled surface position
	Point p;

	/// Sampled surface normal
	Normal n;
};

/** \brief Container for all information related to
 * a surface intersection
 */
struct MTS_EXPORT_RENDER Intersection {
public:
	inline Intersection() : t(std::numeric_limits<Float>::infinity()), shape(NULL) {
	}

	/* Convert a vector expressed inside the shading frame into world
	   coordinates */
	inline Vector toWorld(const Vector &v) const {
		return shFrame.toWorld(v);
	}

	/* Convert a vector expressed inside world coordinates frame into 
	   shading frame coordinates */
	inline Vector toLocal(const Vector &v) const {
		return shFrame.toLocal(v);
	}

	/// Is the current intersection valid?
	inline bool isValid() const {
		return t != std::numeric_limits<Float>::infinity();
	}

	/// Is the intersected shape also a luminaire?
	inline bool isLuminaire() const;

	/// Does the intersected shape have a subsurface integrator?
	inline bool hasSubsurface() const;

	/**
	 * Returns the BSDF of the intersected shape. The
	 * parameter ray must match the one used to create
	 * the intersection record. Computes texture coordinate
	 * partials if this is required by the BSDF.
	 * Should only be called if there is a valid
	 * intersection!
	 */
	inline const BSDF *getBSDF(const RayDifferential &ray);

	/**
	 * Returns radiance emitted into direction d.
	 * Should only be called if the intersected
	 * shape is indeed a luminaire!
	 */
	inline Spectrum Le(const Vector &d) const;

	/**
	 * Returns radiance from a subsurface integrator
	 * emitted into direction d.
	 * Should only be called if the intersected
	 * shape does indeed have a subsurface integrator!
	 */
	inline Spectrum LoSub(const Scene *scene, const Vector &d) const;

	/// Computes texture coordinate partials
	void computePartials(const RayDifferential &ray);

	/// Return a string representation
	std::string toString() const;
public:
	/// Incident direction in the local frame
	Vector wi;

	/// Distance traveled along the ray
	Float t;

	/* Intersection point in 3D coordinates */
	Point p;

	/// Geometry frame
	Frame geoFrame;

	/// Shading frame
	Frame shFrame;

	/// UV surface coordinates
	Point2 uv;

	/// Position partials wrt. to changes in texture-space
	Vector dpdu, dpdv;

	/// Texture coordinate mapping partials wrt. changes in screen-space
	Float dudx, dudy, dvdx, dvdy;

#if defined(MTS_HAS_VERTEX_COLORS)
	/// Interpolated vertex color (if enabled)
	Spectrum color;
#endif

	/// Affected shape
	const Shape *shape;

	/// Have texture coordinate partials been computed
	bool hasUVPartials;
};

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
	 * wrt. solid angle as seen from \a x. Returns the associated probability
	 * density
	 */
	virtual Float sampleSolidAngle(ShapeSamplingRecord &sRec, const Point &x, const Point2 &sample) const;

	/**
	 * Return the probability density of sampling the given point using sampleSolidAngle()
	 */
	virtual Float pdfSolidAngle(const ShapeSamplingRecord &sRec, const Point &x) const;

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


