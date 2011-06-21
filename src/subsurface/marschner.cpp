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
#include "../shapes/hair.h"

MTS_NAMESPACE_BEGIN

class MarschnerModel : public Subsurface {
public:
	MarschnerModel(const Properties &props)
		: Subsurface(props) {
	}

	MarschnerModel(Stream *stream, InstanceManager *manager)
	 : Subsurface(stream, manager) {
		configure();
	}

	virtual ~MarschnerModel() {
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

	void configure() {
		/* Precompute certain things if necessary */
	}

	/// Set the parent object
	void setParent(ConfigurableObject *parent) {
		/// Can't use derivesFrom() here for subtle linker/shared
		/// library reasons on windows
		if (parent->getClass()->getName() != "HairShape")
			Log(EError, "Can only be attached to a HairShape!");
	}

	Spectrum Lo(const Scene *scene, Sampler *sampler, 
			const Intersection &its, const Vector &d, int depth) const {
		Spectrum result(0.0f);

		const HairShape *shape = static_cast<const HairShape *>(its.shape);

		Vector wiLocal = its.wi, wiWorld = d;

		LuminaireSamplingRecord lRec;
		if (scene->sampleLuminaire(its.p, its.time, lRec, sampler->next2D())) {
			/* Do something with lRec */
		}

		/* Recursively gather radiance, but don't include emission */
		RadianceQueryRecord rRec(scene, sampler);
		rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission, NULL);
		rRec.depth = depth + 1;
		
		/// Compute scattered direction
		Vector wo = -d;

		Spectrum recursiveRadiance = static_cast<const SampleIntegrator *>(
			scene->getIntegrator())->Li(RayDifferential(its.p, wo, its.time), rRec);

		return result;
	}

	MTS_DECLARE_CLASS()
private:
};

MTS_IMPLEMENT_CLASS_S(MarschnerModel, false, Subsurface)
MTS_EXPORT_PLUGIN(MarschnerModel, "Marschner hair scattering model");
MTS_NAMESPACE_END
