#if !defined(__RECORDS_H)
#define __RECORDS_H

#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

class Shape;
class Luminaire;
class Medium;
class Scene;
class Sampler;
class BSDF;

/**
 * Specifies the transported quantity when sampling a BSDF
 */
enum ETransportQuantity {
	ERadiance = 1,
	EImportance = 2
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

	/* Return a string representation */
	std::string toString() const;
public:
	/* Incident direction in the local frame */
	Vector wi;

	/* Distance traveled along the ray */
	Float t;

	/* Intersection point in 3D coordinates */
	Point p;

	/* Geometry frame */
	Frame geoFrame;

	/* Shading frame */
	Frame shFrame;

	/* UV surface coordinates */
	Point2 uv;

	/* Position partials wrt. to changes in texture-space */
	Vector dpdu, dpdv;

	/* Texture coordinate mapping partials wrt. changes in screen-space */
	Float dudx, dudy, dvdx, dvdy;

	/* Affected shape */
	const Shape *shape;

	/* Have texture coordinate partials been computed */
	bool hasUVPartials;
};

struct MTS_EXPORT_RENDER ShapeSamplingRecord {
public:
	inline ShapeSamplingRecord() { }

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

/**
 * Data structure used to record information associated with
 * sampled shadow rays
 */
struct MTS_EXPORT_RENDER LuminaireSamplingRecord {
public:
	/// Create an invalid record
	inline LuminaireSamplingRecord() : luminaire(NULL) { }
	
	/**
	 * When a ray strikes a luminaire that is part of the scene,
	 * the associated intersection record can be converted into
	 * a luminaire sampling record in order to query the luminaire
	 * for emitted radiance. (defined in shape.h)
	 */
	inline LuminaireSamplingRecord(const Intersection &its, const Vector &direction);

	/// Return a string representation
	std::string toString() const;
public:
	/// Associated luminaire
	const Luminaire *luminaire;

	/// Data record of the associated shape sample
	ShapeSamplingRecord sRec;

	/// Direction vector pointing away from the light source
	Vector d;

	/// Probability density of the sampled point on the luminaire
	Float pdf;

	/**
	 * Emitted radiance at 'p' into direction 'd' divided by the associated
	 * probability. Already contains the geometric term and optionally 
	 * attenuation when generated via Scene::sampleLuminaireAttenuated.
	 */
	Spectrum Le;
};

struct MTS_EXPORT_RENDER EmissionRecord {
public:
	enum ESamplingType {
		ENormal,
		EPreview
	};

	/// Construct a luminaire sampling record that can be used to query a luminaire
	inline EmissionRecord(const Luminaire *luminaire, 
			const ShapeSamplingRecord &sRec, const Vector &d) 
		: luminaire(luminaire), type(ENormal), sRec(sRec), d(d) { }

	inline EmissionRecord() : luminaire(NULL), type(ENormal) { }

	/// Return a string representation
	std::string toString() const;
public:
	/// Associated luminaire
	const Luminaire *luminaire;

	ESamplingType type;

	/// Data record of the associated shape sample
	ShapeSamplingRecord sRec;

	/// Direction vector pointing away from the light source
	Vector d;

	/**
	 * Radiant emittance at the sampled point. When this 
	 * record was populated using Scene::sampleEmission(), 'P'
	 * has already been multiplied by the directional 
	 * scattering distribution and divided by the associated 
	 * sampling densities.
	 */
	Spectrum P;

	/// Area probability density
	Float pdfArea;

	/// Directional probability density (wrt. projected solid angles)
	Float pdfDir;
};

/**
 * Data record associated with the sampling procedure responsible for
 * choosing a point on the in-scattering line integral (while solving 
 * the radiative transfer equation using Monte Carlo methods).
 */
struct MediumSamplingRecord {
public:
	inline MediumSamplingRecord() : medium(NULL) { }

	/// Return a string representation
	std::string toString() const;
public:
	/* Traveled distance */
	Float t;

	/* Interaction point */
	Point p;

	/* Local particle orientation */
	Vector orientation;

	/* Reference to the associated medium */
	const Medium *medium;

	/* Specifies the attenuation along the segment [mint, t].
	   When sampling a distance fails, this contains the 
	   attenuation along the whole ray.
	*/
	Spectrum attenuation;

	/* The medium's absorption coefficient at that point */
	Spectrum sigmaA;

	/* The medium's scattering coefficient at that point */
	Spectrum sigmaS;

	/**
	 * Can contain two things:
	 * If a medium interaction occurred, this records the probability 
	 * of having sampled the point p. Otherwise, it contains the
	 * probability of moving through the medium without an interaction.
	 */
	Float pdf;

	/// Max. albedo over all spectral samples
	Float albedo;

	/// Multiple importance sampling weight
	Float miWeight;
};

struct MTS_EXPORT_RENDER RadianceQueryRecord {
public:
	/**
	 * List of suported query types. These can be combined by a binary OR.
	 */
	enum ERadianceQuery {
		/* Emitted radiance from a luminaire intersected by the ray */
		EEmittedRadiance     = 0x0001,

		/* Emitted radiance from a subsurface integrator */
		ESubsurfaceRadiance  = 0x0002,

		/* Direct (surface) radiance */
		EDirectRadiance      = 0x0004,

		/* Indirect (surface) radiance, where the last bounce did not go
		   through a Dirac delta BSDF */
		EIndirectRadiance    = 0x0008,

		/* Indirect (surface) radiance, where the last bounce went
		   through a Dirac delta BSDF */
		ECausticRadiance     = 0x0010,

		/* In-scattered radiance due to volumetric scattering (direct) */
		EInscatteredDirectRadiance = 0x0020,

		/* In-scattered radiance due to volumetric scattering (indirect) */
		EInscatteredIndirectRadiance = 0x0040,

		/* Distance to the next surface intersection */
		EDistance            = 0x0080,

		/* Opacity value: 1 when a surface was hit, 0 when the ray leads
		   into empty space. When there is a participating medium,
		   this can also take on fractional values. */
		EOpacity             = 0x0100,

		/* A ray intersection may need to be performed. This can be set to 
		   zero if the caller has already provided the intersection */
		EIntersection        = 0x0200,

		/* Radiance from volumes */
		EVolumeRadiance      = EInscatteredDirectRadiance | EInscatteredIndirectRadiance,

		/* Radiance query without emitted radiance, ray intersection required */
		ERadianceNoEmission  = ESubsurfaceRadiance | EDirectRadiance | EIndirectRadiance
			| ECausticRadiance | EInscatteredDirectRadiance | EInscatteredIndirectRadiance | EIntersection,

		/* Default radiance query, ray intersection required */
		ERadiance = ERadianceNoEmission | EEmittedRadiance,

		/* Radiance + opacity */
		ECameraRay = ERadiance | EOpacity
	};

	/// Construct an invalid radiance query record
	inline RadianceQueryRecord() 
	 : type(0), scene(NULL), sampler(NULL),
	   depth(0), alpha(0), dist(-1), wavelength(-1), extra(0) {
	}

	/// Construct a radiance query record for the given scene and sampler
	inline RadianceQueryRecord(const Scene *scene, Sampler *sampler) 
	 : type(0), scene(scene), sampler(sampler), 
	   depth(0), alpha(0), dist(-1), wavelength(-1), extra(0) {
	}
	
	/// Copy constructor
	inline RadianceQueryRecord(const RadianceQueryRecord &rRec) 
	 : type(rRec.type), scene(rRec.scene), sampler(rRec.sampler), 
	   depth(rRec.depth), alpha(rRec.alpha), dist(rRec.dist),
	   wavelength(rRec.wavelength), extra(rRec.extra) {
	}

	/// Begin a new query of the given type
	inline void newQuery(int _type) {
		type = _type;
		depth = 1;
		wavelength = -1;
		extra = 0;
	}

	/// Initialize the query record for a recursive query
	inline void recursiveQuery(const RadianceQueryRecord &parent, int _type) {
		type = _type;
		scene = parent.scene;
		sampler = parent.sampler;
		depth = parent.depth+1;
		wavelength = parent.wavelength;
		extra = 0;
	}

	/**
	 * Search for a ray intersection (in records.inl). This
	 * does several things at once - if the intersection has 
	 * already been provided, the function returns.
	 * Otherwise, it
	 * - performs the ray intersection
	 * - computes the attenuation due to participating media
	 *   and stores it in <tt>attenuation</tt>.
	 * - sets the alpha value (if <tt>EAlpha</tt> is set in <tt>type</tt>)
	 * - sets the distance value (if <tt>EDistance</tt> is set in <tt>type</tt>)
	 * - clears the <tt>EIntersection</tt> flag in <tt>type</tt>
	 * Returns true if there is a valid intersection.
	 */
	inline bool rayIntersect(const RayDifferential &ray);

	/// Retrieve a 2D sample
	inline Point2 nextSample2D();

	/// Retrieve a 1D sample
	inline Float nextSample1D();

	/// Return a string representation
	std::string toString() const;
public:
	// An asterisk (*) marks entries, which may be overwritten
	// by the callee.

	/// Query type (*)
	int type;

	/// Pointer to the associated scene
	const Scene *scene;

	/// Sample generator
	Sampler *sampler;

	/// Current depth value (# of light bounces) (*)
	int depth;

	/// Surface interaction data structure (*)
	Intersection its;

	/// Attenuation along the current ray (*)
	Spectrum attenuation;

	/// Opacity value of the associated pixel (*)
	Float alpha;

	/**
	 * Ray distance to the first surface interaction 
	 * (if requested by the query type EDistance) (*)
	 */
	Float dist;

	/**
	 * In some cases, the integrator may be forced to restrict
	 * radiance computations to one wavelength (e.g. when intersecting
	 * a dielectric material with dispersion or while path
	 * tracing through a highly scattering medium with a non-constant
	 * scattering coefficient). This attribute is used to store the
	 * chosen wavelength. (*)
	 */
	int wavelength;

	/**
	 * Internal flag, which can be used to pass additional information 
	 * amonst recursive calls inside an integrator. The use
	 * is dependent on the particular integrator implementation. (*)
	 */
	int extra;
};

/**
 * Data structure, which contains all information required to
 * sample or query a BSDF. 
 */
struct BSDFQueryRecord {
public:
	/**
	 * Given a surface interaction and an incident direction 
	 * construct a query record which can be used to sample 
	 * an outgoing direction.
	 * For convenience, this function uses the local incident direction 
	 * vector contained in the supplied intersection record.
	 */
	explicit inline BSDFQueryRecord(RadianceQueryRecord &rRec, const Intersection &its, Point2 sample) 
	  : rRec(&rRec), its(its), wi(its.wi), sample(sample), quantity(ERadiance),
	    typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
	}

	/**
	 * Given a surface interaction and an incident direction, 
	 * construct a query record which can be used to sample
	 * an outgoing direction.
	 * For convenience, this function uses the local incident direction 
	 * vector contained in the supplied intersection record.
	 */
	explicit inline BSDFQueryRecord(const Intersection &its, Point2 sample) 
	  : rRec(NULL), its(its), wi(its.wi), sample(sample), quantity(ERadiance),
	    typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
	}

	/**
	 * Given a surface interaction an an incident/exitant direction 
	 * pair (wi, wo), create a BSDF query record to evaluate f(wi, wo).
	 * For convenience, this function uses the local incident direction 
	 * vector contained in the supplied intersection record.
	 */
	explicit inline BSDFQueryRecord(RadianceQueryRecord &rRec, 
		const Intersection &its, const Vector &wo) 
	  : rRec(&rRec), its(its), wi(its.wi), wo(wo), sample(sample), quantity(ERadiance),
	    typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
	}

	/**
	 * Given a surface interaction an an incident/exitant direction 
	 * pair (wi, wo), create a BSDF query record to evaluate f(wi, wo). 
	 * For convenience, this function uses the local incident direction 
	 * vector contained in the supplied intersection record.
	 */
	explicit inline BSDFQueryRecord(const Intersection &its, const Vector &wo) 
	  : rRec(NULL), its(its), wi(its.wi), wo(wo), sample(sample), quantity(ERadiance),
	    typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
	}

	/**
	 * Given a surface interaction an an incident/exitant direction 
	 * pair (wi, wo), create a BSDF query record to evaluate f(wi, wo).
	 */
	explicit inline BSDFQueryRecord(const Intersection &its, const Vector &wi, const Vector &wo) 
	  : rRec(NULL), its(its), wi(wi), wo(wo), sample(sample), quantity(ERadiance),
	    typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
	}

	inline BSDFQueryRecord(const BSDFQueryRecord &r) :
		rRec(r.rRec), its(r.its), wi(r.wi), wo(r.wo), sample(r.sample), quantity(r.quantity),
		typeMask(r.typeMask), sampledType(r.sampledType), component(r.component),
		sampledComponent(r.sampledComponent) {
	}

	/// Return a string representation
	std::string toString() const;
public:
	/* Pointer to the associated radiance query record (or NULL) */
	RadianceQueryRecord *rRec;

	/* Surface interaction */
	const Intersection &its;

	/* Incident direction */
	Vector wi;

	/* Outgoing/sampled direction */
	Vector wo;

	/* Random sample used to generate the new direction */
	Point2 sample;

	/* Transported quantity (radiance or importance) -- required for 
	   non-reciprocal BSDFs such as transmission through a dielectric
	   material */
	ETransportQuantity quantity;

	/* Bit mask containing the component types, which may be sampled.
	   After sampling has been performed, the component type is stored
	   inside 'sampledType'. */
	unsigned int typeMask, sampledType;

	/* To sample a specific BSDF component, this entry must be non-negative.
	   After sampling has been performed, the component index is stored
	   inside 'sampledComponent' */
	int component, sampledComponent;
};

extern void operator<<(const ETransportQuantity &quantity, std::ostream &os);

MTS_NAMESPACE_END

#endif /* __RECORDS_H */
