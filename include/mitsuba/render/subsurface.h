#if !defined(__SUBSURFACE_H)
#define __SUBSURFACE_H

#include <mitsuba/render/records.h>
#include <mitsuba/core/netobject.h>

MTS_NAMESPACE_BEGIN

class Scene;
class RenderQueue;
class RenderJob;

/**
 * Abstract subsurface integrator -- can be attached to an arbitrary 
 * shape to compute exitant radiance due to internal scattering.
 */
class MTS_EXPORT_RENDER Subsurface : public NetworkedObject {
public:
	/**
	 * Possibly perform a pre-process task. The last three parameters are
	 * resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to all local and remote workers.
	 */
	virtual void preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) = 0;

	/// Return the list of shapes associated with this subsurface integrator
	inline const std::vector<Shape *> getShapes() const { return m_shapes; }

	/// Get the exitant radiance for a point on the surface
	virtual Spectrum Lo(const Scene *scene, const Intersection &its, const Vector &d) const = 0;

	/// Serialize this subsurface integrator to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Set the parent node of the subsurface integrator
	void setParent(ConfigurableObject *parent);

	MTS_DECLARE_CLASS()
protected:
	/// Create a new subsurface scattering class
	Subsurface(const Properties &props);

	/// Unserialize a subsurface integrator from a binary data stream
	Subsurface(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Subsurface();
protected:
	Spectrum m_sigmaS;
	Spectrum m_sigmaA;
	Spectrum m_sigmaT;
	Float m_eta, m_sizeMultiplier;
	std::vector<Shape *> m_shapes;
};

MTS_NAMESPACE_END

#endif /* __SUBSURFACE_H */
