/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

class MarschnerShader : public Subsurface {
public:
	MarschnerShader(const Properties &props)
		: Subsurface(props) {
	}

	MarschnerShader(Stream *stream, InstanceManager *manager)
	 : Subsurface(stream, manager) {
		configure();
	}

	virtual ~MarschnerShader() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID) {
		if (!scene->getIntegrator()->getClass()->derivesFrom(MTS_CLASS(SampleIntegrator)))
			Log(EError, "Must be used with a SampleIntegrator!");
		return true;
	}

	void cancel() {
	}

	Spectrum Lo(const Scene *scene, const Intersection &its, const Vector &d, int depth) const {
//		Vector wiLocal = its.wi;
	//	Vector wiWorld = d;

		/// Compute scattered direction
		Vector wo = -d;
		const Shape *shape = its.shape;

		RadianceQueryRecord rRec(scene, const_cast<Sampler *>(scene->getSampler()));
		rRec.newQuery(RadianceQueryRecord::ERadiance, NULL);
		rRec.depth = depth + 1;
		Spectrum recursiveRadiance = static_cast<const SampleIntegrator *>(scene->getIntegrator())->Li(RayDifferential(its.p, wo, its.time), rRec);

		Spectrum modelResult(0.5f);

		return modelResult * recursiveRadiance;
	}

	void configure() {
	}


	MTS_DECLARE_CLASS()
private:
};

MTS_IMPLEMENT_CLASS_S(MarschnerShader, false, Subsurface)
MTS_EXPORT_PLUGIN(MarschnerShader, "Marschner model");
MTS_NAMESPACE_END
