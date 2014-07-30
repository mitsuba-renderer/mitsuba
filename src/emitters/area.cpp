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

#include <mitsuba/render/emitter.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{area}{Area light}
 * \icon{emitter_area}
 * \order{2}
 * \parameters{
 *     \parameter{radiance}{\Spectrum}{
 *         Specifies the emitted radiance in units of
 *         power per unit area per unit steradian.
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This plugin implements an area light, i.e. a light source that emits
 * diffuse illumination from the exterior of an arbitrary shape.
 * Since the emission profile of an area light is completely diffuse, it
 * has the same apparent brightness regardless of the observer's viewing
 * direction. Furthermore, since it occupies a nonzero amount of space, an
 * area light generally causes scene objects to cast soft shadows.
 *
 * When modeling scenes involving area lights, it is preferable
 * to use spheres as the emitter shapes, since they provide a
 * particularly good direct illumination sampling strategy (see
 * the \pluginref{sphere} plugin for an example).
 *
 * To create an area light source, simply instantiate the desired
 * emitter shape and specify an \code{area} instance as its child:
 *
 * \vspace{4mm}
 * \begin{xml}
 * <!-- Create a spherical light source at the origin -->
 * <shape type="sphere">
 *     <emitter type="area">
 *         <spectrum name="radiance" value="1"/>
 *     </emitter>
 * </shape>
 * \end{xml}
 */

class AreaLight : public Emitter {
public:
	AreaLight(const Properties &props) : Emitter(props) {
		m_type |= EOnSurface;

		if (props.hasProperty("toWorld"))
			Log(EError, "Found a 'toWorld' transformation -- this is not "
				"allowed -- the area light inherits this transformation from "
				"its parent shape");

		m_radiance = props.getSpectrum("radiance", Spectrum::getD65());
		m_power = Spectrum(0.0f); /// Don't know the power yet
	}

	AreaLight(Stream *stream, InstanceManager *manager)
		: Emitter(stream, manager) {
		m_radiance = Spectrum(stream);
		m_power = Spectrum(stream);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Emitter::serialize(stream, manager);
		m_radiance.serialize(stream);
		m_power.serialize(stream);
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec,
			const Point2 &sample, const Point2 *extra) const {
		m_shape->samplePosition(pRec, sample);
		return m_power;
	}

	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
		return m_radiance * M_PI;
	}

	Spectrum eval(const Intersection &its, const Vector &d) const {
		if (dot(its.shFrame.n, d) <= 0)
			return Spectrum(0.0f);
		else
			return m_radiance;
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_shape->pdfPosition(pRec);
	}

	Spectrum sampleDirection(DirectionSamplingRecord &dRec,
			PositionSamplingRecord &pRec,
			const Point2 &sample, const Point2 *extra) const {
		Vector local = warp::squareToCosineHemisphere(sample);
		dRec.d = Frame(pRec.n).toWorld(local);
		dRec.pdf = warp::squareToCosineHemispherePdf(local);
		dRec.measure = ESolidAngle;
		return Spectrum(1.0f);
	}

	Spectrum evalDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {
		Float dp = dot(dRec.d, pRec.n);

		if (dRec.measure != ESolidAngle || dp < 0)
			dp = 0.0f;

		return Spectrum(INV_PI * dp);
	}

	Float pdfDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {
		Float dp = dot(dRec.d, pRec.n);

		if (dRec.measure != ESolidAngle || dp < 0)
			dp = 0.0f;

		return INV_PI * dp;
	}

	Spectrum sampleRay(Ray &ray,
			const Point2 &spatialSample,
			const Point2 &directionalSample,
			Float time) const {
		PositionSamplingRecord pRec(time);
		m_shape->samplePosition(pRec, spatialSample);
		Vector local = warp::squareToCosineHemisphere(directionalSample);
		ray.setTime(time);
		ray.setOrigin(pRec.p);
		ray.setDirection(Frame(pRec.n).toWorld(local));
		return m_power;
	}

	Spectrum sampleDirect(DirectSamplingRecord &dRec,
			const Point2 &sample) const {
		m_shape->sampleDirect(dRec, sample);

		/* Check that the emitter and reference position are oriented correctly
		   with respect to each other. Note that the >= 0 check
		   for 'refN' is intentional -- those sampling requests that specify
		   a reference point within a medium or on a transmissive surface
		   will set dRec.refN = 0, hence they should always be accepted. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0 && dRec.pdf != 0) {
			return m_radiance / dRec.pdf;
		} else {
			dRec.pdf = 0.0f;
			return Spectrum(0.0f);
		}
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		/* Check that the emitter and receiver are oriented correctly
		   with respect to each other. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0) {
			return m_shape->pdfDirect(dRec);
		} else {
			return 0.0f;
		}
	}

	void setParent(ConfigurableObject *parent) {
		Emitter::setParent(parent);

		if (parent->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Shape *shape = static_cast<Shape *>(parent);
			if (m_shape == shape || shape->isCompound())
				return;

			if (m_shape != NULL)
				Log(EError, "An area light cannot be parent of multiple shapes");

			m_shape = shape;
			m_shape->configure();
			m_power = m_radiance * M_PI * m_shape->getSurfaceArea();
		} else {
			Log(EError, "An area light must be child of a shape instance");
		}
	}

	AABB getAABB() const {
		return m_shape->getAABB();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AreaLight[" << endl
			<< "  radiance = " << m_radiance.toString() << "," << endl
			<< "  samplingWeight = " << m_samplingWeight << "," << endl
			<< "  surfaceArea = ";
		if (m_shape)
			oss << m_shape->getSurfaceArea();
		else
			oss << "<no shape attached!>";
		oss << "," << endl
		    << "  medium = " << indent(m_medium.toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_radiance, m_power;
};

// ================ Hardware shader implementation ================

class AreaLightShader : public Shader {
public:
	AreaLightShader(Renderer *renderer, const Spectrum &radiance)
		: Shader(renderer, EEmitterShader), m_radiance(radiance) {
	}

	void resolve(const GPUProgram *program, const std::string &evalName,
			std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_radiance", false));
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_radiance;" << endl
			<< endl
			<< "vec3 " << evalName << "_area(vec2 uv) {" << endl
			<< "    return " << evalName << "_radiance * pi;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "    if (cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return vec3(inv_pi);" << endl
			<< "}" << endl;
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_radiance);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_radiance;
};

Shader *AreaLight::createShader(Renderer *renderer) const {
	return new AreaLightShader(renderer, m_radiance);
}

MTS_IMPLEMENT_CLASS(AreaLightShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(AreaLight, false, Emitter)
MTS_EXPORT_PLUGIN(AreaLight, "Area light");
MTS_NAMESPACE_END
