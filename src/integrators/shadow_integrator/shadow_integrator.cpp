#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

class shadow_integrator : public SamplingIntegrator {
public:
	MTS_DECLARE_CLASS()

public:
	shadow_integrator(const Properties& props) :SamplingIntegrator(props) {
		Spectrum default_color;
		default_color.fromLinearRGB(0.2f, 0.5f, 0.2f);
		m_color = props.getSpectrum("color", default_color);
	}

	// unserialize, in order
	shadow_integrator(Stream* stream, InstanceManager* manager) :
		SamplingIntegrator(stream, manager) {
		m_color = Spectrum(stream);
		m_max_dist = stream->readFloat();
	}

	// serialize
	void serialize(Stream* stream, InstanceManager* manager) const {
		SamplingIntegrator::serialize(stream, manager);
		m_color.serialize(stream);
		stream->writeFloat(m_max_dist);
	}

	/// Query for an unbiased estimate of the radiance along r
	Spectrum Li(const RayDifferential& r, RadianceQueryRecord& rRec) const {
		typedef TVector3<Float> vec3;
		if (rRec.rayIntersect(r)) {
			Float distance = rRec.its.t;
			vec3 normal_value = vec3(rRec.its.shFrame.n.x, rRec.its.shFrame.n.y, rRec.its.shFrame.n.z);
			Log(EInfo, (std::string("before: ") + normal_value.toString()).c_str());
			normal_value = 0.5f * normalize(normal_value) + vec3(0.5f);
			Log(EInfo, (std::string("after: ") + normal_value.toString()).c_str());

			Spectrum ret;
			ret.fromLinearRGB(normal_value.x, normal_value.y, normal_value.z);

			return ret;
		}

		return Spectrum(0.0f);
	}


	bool preprocess(const Scene* scene, RenderQueue* queue, const RenderJob* job, int sceneResID, int cameraResID, int samplerResID) {
		SamplingIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

		const AABB& scene_aabb = scene->getAABB();
		Point camera_position = scene->getSensor()->getWorldTransform()->eval(0).transformAffine(Point(0.0f));
		m_max_dist = -std::numeric_limits<Float>::infinity();

		for (int i = 0; i < 8; ++i) {
			m_max_dist = std::max(m_max_dist, (camera_position - scene_aabb.getCorner(i)).length());
		}

		return true;
	}

private:
	Spectrum m_color;
	Float m_max_dist;
};

MTS_IMPLEMENT_CLASS_S(shadow_integrator, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(shadow_integrator, "A contrived integrator");
MTS_NAMESPACE_END