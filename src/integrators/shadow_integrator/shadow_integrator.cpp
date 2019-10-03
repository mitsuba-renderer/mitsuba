#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

class shadow_integrator: public SamplingIntegrator {
public:
    MTS_DECLARE_CLASS()

public:
	shadow_integrator(const Properties &props):SamplingIntegrator(props) {
        Spectrum default_color;
        default_color.fromLinearRGB(0.2f, 0.5f, 0.2f);
        m_color = props.getSpectrum("color", default_color);
    }

    // unserialize
	shadow_integrator(Stream *stream, InstanceManager *manager):
        SamplingIntegrator(stream, manager) {
        m_color = Spectrum(stream);
    }

    // serialize
    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        m_color.serialize(stream);
    }

    /// Query for an unbiased estimate of the radiance along r
    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        return m_color;
    }

private:
    Spectrum m_color;
};

MTS_IMPLEMENT_CLASS_S(shadow_integrator, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(shadow_integrator, "A contrived integrator");
MTS_NAMESPACE_END