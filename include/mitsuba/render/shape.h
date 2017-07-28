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

#pragma once
#if !defined(__MITSUBA_RENDER_SHAPE_H_)
#define __MITSUBA_RENDER_SHAPE_H_

#include <mitsuba/render/common.h>
#include <mitsuba/core/cobject.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/** \brief Container for all information related to
 * a surface intersection
 * \ingroup librender
 * \ingroup libpython
 */
struct MTS_EXPORT_RENDER Intersection {
public:
    inline Intersection() :
        shape(NULL), t(std::numeric_limits<Float>::infinity()) { }

    /// Convert a local shading-space vector into world space
    inline Vector toWorld(const Vector &v) const {
        return shFrame.toWorld(v);
    }

    /// Convert a world-space vector into local shading coordinates
    inline Vector toLocal(const Vector &v) const {
        return shFrame.toLocal(v);
    }

    /// Is the current intersection valid?
    inline bool isValid() const {
        return t != std::numeric_limits<Float>::infinity();
    }

    /// Is the intersected shape also a emitter?
    inline bool isEmitter() const;

    /// Is the intersected shape also a sensor?
    inline bool isSensor() const;

    /// Does the intersected shape have a subsurface integrator?
    inline bool hasSubsurface() const;

    /// Does the surface mark a transition between two media?
    inline bool isMediumTransition() const;

    /**
     * \brief Determine the target medium
     *
     * When \c isMediumTransition() = \c true, determine the medium that
     * contains the ray (\c this->p, \c d)
     */
    inline const Medium *getTargetMedium(const Vector &d) const;

    /**
     * \brief Determine the target medium based on the cosine
     * of the angle between the geometric normal and a direction
     *
     * Returns the exterior medium when \c cosTheta > 0 and
     * the interior medium when \c cosTheta <= 0.
     */
    inline const Medium *getTargetMedium(Float cosTheta) const;

    /**
     * \brief Returns the BSDF of the intersected shape.
     *
     * The parameter ray must match the one used to create the intersection
     * record. This function computes texture coordinate partials if this is
     * required by the BSDF (e.g. for texture filtering).
     *
     * \remark This function should only be called if there is a valid
     * intersection!
     */
    inline const BSDF *getBSDF(const RayDifferential &ray);

    /// Returns the BSDF of the intersected shape
    inline const BSDF *getBSDF() const;

    /**
     * \brief Returns radiance emitted into direction d.
     *
     * \remark This function should only be called if the
     * intersected shape is actually an emitter.
     */
    inline Spectrum Le(const Vector &d) const;

    /**
     * \brief Returns radiance from a subsurface integrator
     * emitted into direction d.
     *
     * \remark Should only be called if the intersected
     * shape is actually a subsurface integrator.
     */
    inline Spectrum LoSub(const Scene *scene, Sampler *sampler,
            const Vector &d, int depth=0) const;

    /// Computes texture coordinate partials
    void computePartials(const RayDifferential &ray);

    /// Move the intersection forward or backward through time
    inline void adjustTime(Float time);

    /// Calls the suitable implementation of \ref Shape::getNormalDerivative()
    inline void getNormalDerivative(Vector &dndu, Vector &dndv,
        bool shadingFrame = true) const;

    /// Return a string representation
    std::string toString() const;
public:
    /// Pointer to the associated shape
    const Shape *shape;

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

    /// Position partials wrt. the UV parameterization
    Vector dpdu, dpdv;

    /// UV partials wrt. changes in screen-space
    Float dudx, dudy, dvdx, dvdy;

    /// Time value associated with the intersection
    Float time;

    /// Interpolated vertex color
    Spectrum color;

    /// Incident direction in the local shading frame
    Vector wi;

    /// Have texture coordinate partials been computed
    bool hasUVPartials : 1;

    /// Primitive index, e.g. the triangle ID (if applicable)
    uint32_t primIndex : 31;

    /// Stores a pointer to the parent instance, if applicable
    const Shape *instance;
};

/** \brief Abstract base class of all shapes
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER Shape : public ConfigurableObject {
public:
    // =============================================================
    //! @{ \name Query functions to be implemented in subclasses
    // =============================================================

    /// Return the name of this shape (e.g. the filename)
    virtual std::string getName() const;

    /// Is this a compound shape consisting of several sub-objects?
    virtual bool isCompound() const;

    /**
     * \brief Return a sub-element of a compound shape.
     *
     * When expanding shapes, the scene will repeatedly call this
     * function with increasing indices. Returning \a NULL indicates
     * that no more elements are available.
     */
    virtual Shape *getElement(int i);

    /**
     * \brief Return the shape's surface area
     *
     * Assumes that the object is not undergoing some kind of
     * time-dependent scaling.
     *
     * The default implementation throws an exception
     */
    virtual Float getSurfaceArea() const;

    /// Return a bounding box containing the shape
    virtual AABB getAABB() const = 0;

    /**
     * \brief Returns the minimal axis-aligned bounding box
     * of this shape when clipped to another bounding box.
     *
     * This is extremely important to construct decent kd-trees.
     * The default implementation just takes the bounding box
     * returned by \ref getAABB() and clips it to \a box.
     */
    virtual AABB getClippedAABB(const AABB &box) const;

    /**
     * \brief Create a triangle mesh approximation of this shape
     *
     * This function is used by the realtime preview and
     * certain integrators, which rely on hardware rasterization.
     *
     * The default implementation simply returns \a NULL.
     */
    virtual ref<TriMesh> createTriMesh();

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Ray tracing routines
    // =============================================================

    /**
     * \brief Fast ray intersection test
     *
     * Check whether the shape is intersected by the given ray. Some
     * temporary space (\ref MTS_KD_INTERSECTION_TEMP-4 bytes) is,
     * supplied which can be used to cache information about the
     * intersection. The function \ref fillIntersectionRecord()
     * can later use this information to fill in a detailed
     * intersection record.
     *
     * \remark In Python, this function also calls \c fillIntersectionRecord
     * and has the signature
     * <tt>intersection = shape.rayIntersect(ray, mint, maxt)</tt>
     */
    virtual bool rayIntersect(const Ray &ray, Float mint,
            Float maxt, Float &t, void *temp) const;

    /**
     * \brief Fast ray intersection test for visibility queries
     *
     * Check whether the shape is intersected by the given ray.
     * No details about the intersection are returned, hence the
     * function is only useful for visibility queries. For most
     * shapes, this will simply call forward the call to \ref
     * rayIntersect. When the shape actually contains a nested
     * kd-tree, some optimizations are possible.
     *
     * \remark This function is not exposed in Python
     */
    virtual bool rayIntersect(const Ray &ray, Float mint, Float maxt) const;

    /**
     * \brief Given that an intersection has been found, create a
     * detailed intersection record
     *
     * \remark This function is not directly exposed in Python.
     * It is implicitly called as part of \c rayIntersect.
     */
    virtual void fillIntersectionRecord(const Ray &ray,
            const void *temp, Intersection &its) const;

    /**
     * \brief Return the derivative of the normal vector with
     * respect to the UV parameterization
     *
     * This can be used to compute Gaussian and principal curvatures,
     * amongst other things.
     *
     * \param its
     *     Intersection record associated with the query
     * \param dndu
     *     Parameter used to store the partial derivative of the
     *     normal vector with respect to \c u
     * \param dndv
     *     Parameter used to store the partial derivative of the
     *     normal vector with respect to \c v
     * \param shadingFrame
     *     Specifies whether to compute the derivative of the
     *     geometric normal \a or the shading normal of the surface
     *
     * \remark In Python, the signature of this function is
     * <tt>dndu, dndv = shape.getNormalDerivative(its, shadingFrame)</tt>
     */
    virtual void getNormalDerivative(const Intersection &its,
        Vector &dndu, Vector &dndv, bool shadingFrame = true) const;

    /**
     * \brief Compute the Gaussian and mean curvature at the given
     * surface intersection.
     *
     * \param its
     *     Intersection record associated with the query
     * \param H
     *     Parameter used to store the mean curvature
     * \param K
     *     Parameter used to store the Gaussian curvature
     * \param shadingFrame
     *     Specifies whether to compute the curvature based on the
     *     geometric normal \a or the shading normal of the surface
     *
     * \remark In Python, the signature of this function is
     * <tt>H, K = shape.getCurvature(its, shadingFrame)</tt>
     */
    void getCurvature(const Intersection &its, Float &H, Float &K,
        bool shadingFrame = true) const;

    /**
     * Adjust an intersection record to a different time value
     */
    virtual void adjustTime(Intersection &its, Float time) const;

    /**
     * \brief Return the internal kd-tree of this shape (if any)
     *
     * This function is used by the kd-tree visualization in
     * the interactive walkthrough. The default implementation
     * simply returns NULL.
     *
     * \remark This function is not exposed in Python
     */
    virtual const KDTreeBase<AABB> *getKDTree() const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Sampling routines
    // =============================================================

    /**
     * \brief Sample a point on the surface of this shape instance
     * (with respect to the area measure)
     *
     * The returned sample density will be uniform over the surface.
     *
     * \param pRec
     *     A position record, which will be used to return the sampled
     *     position, as well as auxilary information about the sample.
     *
     * \param sample
     *     A uniformly distributed 2D vector
     */
    virtual void samplePosition(PositionSamplingRecord &pRec,
            const Point2 &sample) const;

    /**
     * \brief Query the probability density of \ref samplePosition() for
     * a particular point on the surface.
     *
     * This method will generally return the inverse of the surface area.
     *
     * \param pRec
     *     A position record, which will be used to return the sampled
     *     position, as well as auxilary information about the sample.
     */

    virtual Float pdfPosition(const PositionSamplingRecord &pRec) const;

    /**
     * \brief Sample a point on the surface of this shape instance
     * (with respect to the solid angle measure)
     *
     * The sample density should ideally be uniform in direction as seen from
     * the reference point \c dRec.p.
     *
     * This general approach for sampling positions is named "direct" sampling
     * throughout Mitsuba motivated by direct illumination rendering techniques,
     * which represent the most important application.
     *
     * When no implementation of this function is supplied, the \ref Shape
     * class will revert to the default approach, which piggybacks on
     * \ref sampleArea(). This usually results in a a suboptimal sample
     * placement, which can manifest itself in the form of high variance
     *
     * \param dRec
     *    A direct sampling record that specifies the reference point and a
     *    time value. After the function terminates, it will be populated
     *    with the position sample and related information
     *
     * \param sample
     *     A uniformly distributed 2D vector
     */
    virtual void sampleDirect(DirectSamplingRecord &dRec,
            const Point2 &sample) const;

    /**
     * \brief Query the probability density of \ref sampleDirect() for
     * a particular point on the surface.
     *
     * \param dRec
     *    A direct sampling record, which specifies the query
     *    location. Note that this record need not be completely
     *    filled out. The important fields are \c p, \c n, \c ref,
     *    \c dist, \c d, \c measure, and \c uv.
     *
     * \param p
     *     An arbitrary point used to define the solid angle measure
     */
    virtual Float pdfDirect(const DirectSamplingRecord &dRec) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Miscellaneous
    // =============================================================

    /// Does the surface of this shape mark a medium transition?
    inline bool isMediumTransition() const { return m_interiorMedium.get() || m_exteriorMedium.get(); }
    /// Return the medium that lies on the interior of this shape (\c NULL == vacuum)
    inline Medium *getInteriorMedium() { return m_interiorMedium; }
    /// Return the medium that lies on the interior of this shape (\c NULL == vacuum, const version)
    inline const Medium *getInteriorMedium() const { return m_interiorMedium.get(); }
    /// Return the medium that lies on the exterior of this shape (\c NULL == vacuum)
    inline Medium *getExteriorMedium() { return m_exteriorMedium; }
    /// Return the medium that lies on the exterior of this shape (\c NULL == vacuum, const version)
    inline const Medium *getExteriorMedium() const { return m_exteriorMedium.get(); }

    /// Does this shape have a sub-surface integrator?
    inline bool hasSubsurface() const { return m_subsurface.get() != NULL; }
    /// Return the associated sub-surface integrator
    inline Subsurface *getSubsurface() { return m_subsurface; }
    /// Return the associated sub-surface integrator
    inline const Subsurface *getSubsurface() const { return m_subsurface.get(); }

    /// Is this shape also an area emitter?
    inline bool isEmitter() const { return m_emitter.get() != NULL; }
    /// Return the associated emitter (if any)
    inline Emitter *getEmitter() { return m_emitter; }
    /// Return the associated emitter (if any)
    inline const Emitter *getEmitter() const { return m_emitter.get(); }
    /// Set the emitter of this shape
    inline void setEmitter(Emitter *emitter) { m_emitter = emitter; }

    /// Is this shape also an area sensor?
    inline bool isSensor() const { return m_sensor.get() != NULL; }
    /// Return the associated sensor (if any)
    inline Sensor *getSensor() { return m_sensor; }
    /// Return the associated sensor (if any)
    inline const Sensor *getSensor() const { return m_sensor.get(); }

    /// Does the shape have a BSDF?
    inline bool hasBSDF() const { return m_bsdf.get() != NULL; }
    /// Return the shape's BSDF
    inline const BSDF *getBSDF() const { return m_bsdf.get(); }
    /// Return the shape's BSDF
    inline BSDF *getBSDF() { return m_bsdf.get(); }
    /// Set the BSDF of this shape
    inline void setBSDF(BSDF *bsdf) { m_bsdf = bsdf; }

    /**
     * \brief Return the number of primitives (triangles, hairs, ..)
     * contributed to the scene by this shape
     *
     * Does not include instanced geometry
     */
    virtual size_t getPrimitiveCount() const = 0;

    /**
     * \brief Return the number of primitives (triangles, hairs, ..)
     * contributed to the scene by this shape
     *
     * Includes instanced geometry
     */
    virtual size_t getEffectivePrimitiveCount() const = 0;

    /// Copy attachments (BSDF, Emitter, ..) from another shape
    void copyAttachments(Shape *shape);

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name ConfigurableObject interface
    // =============================================================

    /// Called once after constructing the object
    virtual void configure();

    /// Serialize this shape to a stream
    virtual void serialize(Stream *stream, InstanceManager *manager) const;

    /// Add a child (e.g. a emitter/sub surface integrator) to this shape
    void addChild(const std::string &name, ConfigurableObject *child);

    /// Add an unnamed child
    inline void addChild(ConfigurableObject *child) { addChild("", child); }

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Create a new shape
    Shape(const Properties &props);

    /// Unserialize a shape
    Shape(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~Shape();
protected:
    std::string m_name;
    ref<BSDF> m_bsdf;
    ref<Subsurface> m_subsurface;
    ref<Emitter> m_emitter;
    ref<Sensor> m_sensor;
    ref<Medium> m_interiorMedium;
    ref<Medium> m_exteriorMedium;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_SHAPE_H_ */


