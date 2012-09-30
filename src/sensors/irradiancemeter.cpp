/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/render/sensor.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{irradiancemeter}{Irradiance meter}
 * \order{6}
 * \parameters{
 *     \parameter{shutterOpen, shutterClose}{\Float}{
 *         Specifies the time interval of the measurement---this
 *         is only relevant when the scene is in motion.
 *         \default{0}
 *     }
 * }
 *
 * This sensor plugin implements a simple irradiance meter, which 
 * measures the average incident power per unit area over a 
 * provided surface. 
 * Such a sensor is useful for conducting virtual experiments and 
 * testing the renderer for correctness.
 *
 * To create an irradiance meter, instantiate the desired measurement 
 * shape and specify the sensor as its child. Note that when the
 * sensor's film resolution is larger than $1\times 1$, each pixel
 * will record the average irradiance over a rectangular part of the
 * shape's UV parameterization.
 *
 * \vspace{4mm}
 * \begin{xml}
 * <scene version=$\MtsVer$>
 *     <!-- Measure the average irradiance arriving on
 *          a unit radius sphere located at the origin -->
 *     <shape type="sphere">
 *         <sensor type="irradiancemeter">
 *             <!-- Write the output to a MATLAB M-file. The output file will
 *                  contain a 1x1 matrix storing an estimate of the average 
 *                  irradiance over the surface of the sphere. -->
 *             <film type="mfilm"/>
 *
 *             <!-- Use 1024 samples for the measurement -->
 *             <sampler type="independent">
 *                 <integer name="sampleCount" value="1024"/>
 *             </sampler>
 *         </sensor>
 *     </shape>
 *
 *     <!-- ... other scene declarations ... -->
 * </scene>
 * \end{xml}
 */

class IrradianceMeter : public Sensor {
public:
	IrradianceMeter(const Properties &props) : Sensor(props) {
		m_type |= ENeedsApertureSample | EOnSurface;

		if (props.hasProperty("toWorld"))
			Log(EError, "Found a 'toWorld' transformation -- this is not "
				"allowed -- the irradiance meter inherits this transformation from "
				"its parent shape");
	}

	IrradianceMeter(Stream *stream, InstanceManager *manager) 
		: Sensor(stream, manager) {
		m_shape = static_cast<Shape *>(manager->getInstance(stream));
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sensor::serialize(stream, manager);
		manager->serialize(stream, m_shape);
	}

	void configure() {
		Sensor::configure();
		m_invResolution = Vector2(
			1.0f / m_film->getCropSize().x,
			1.0f / m_film->getCropSize().y
		);
	}

	Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
			const Point2 &otherSample, Float timeSample) const {
		ray.time = sampleTime(timeSample);
		ray.mint = Epsilon;
		ray.maxt = std::numeric_limits<Float>::infinity();

		PositionSamplingRecord pRec(ray.time);
		m_shape->samplePosition(pRec, Point2(
			pixelSample.x * m_invResolution.x,
			pixelSample.y * m_invResolution.y));

		ray.setOrigin(pRec.p);
		ray.setDirection(Frame(pRec.n).toWorld(
			Warp::squareToCosineHemisphere(otherSample)));

		return Spectrum(1.0f);
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec,
			const Point2 &sample, const Point2 *extra) const {
		Point2 samplePos(sample);

		if (extra) {
			/* The caller wants to condition on a specific pixel position */
			samplePos.x = (extra->x + sample.x) * m_invResolution.x;
			samplePos.y = (extra->y + sample.y) * m_invResolution.y;
		}

		pRec.uv = Point2(samplePos.x * m_resolution.x,
			samplePos.y * m_resolution.y);

		m_shape->samplePosition(pRec, samplePos);
		return Spectrum(1.0f);
	}

	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
		return Spectrum(m_shape->pdfPosition(pRec));
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_shape->pdfPosition(pRec);
	}

	Spectrum sampleDirection(DirectionSamplingRecord &dRec, 
			PositionSamplingRecord &pRec,
			const Point2 &sample, const Point2 *extra) const {
		Vector local = Warp::squareToCosineHemisphere(sample);
		dRec.d = Frame(pRec.n).toWorld(local);
		dRec.pdf = Warp::squareToCosineHemispherePdf(local);
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

	Spectrum sampleDirect(DirectSamplingRecord &dRec, 
			const Point2 &sample) const {
		m_shape->sampleDirect(dRec, sample);

		/* Check that the sensor and reference position are oriented correctly 
		   with respect to each other. Note that the >= 0 check
		   for 'refN' is intentional -- those sampling requests that specify
		   a reference point within a medium or on a transmissive surface
		   will set dRec.refN = 0, hence they should always be accepted. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0 && dRec.pdf != 0.0f) {
			dRec.uv = Point2(0.5f);
			return Spectrum(1.0f / dRec.pdf);
		} else {
			dRec.pdf = 0.0f;
			return Spectrum(0.0f);
		}
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		/* Check that the sensor and reference position are oriented correctly 
		   with respect to each other. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0) {
			return m_shape->pdfDirect(dRec);
		} else {
			return 0.0f;
		}
	}

	void setParent(ConfigurableObject *parent) {
		Sensor::setParent(parent);

		if (parent->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Shape *shape = static_cast<Shape *>(parent);
			if (m_shape == shape || shape->isCompound())
				return;

			if (m_shape != NULL)
				Log(EError, "An irradiance sensor cannot be parent of multiple shapes");

			m_shape = shape;
		} else {
			Log(EError, "An irradiance sensor must be child of a shape instance");
		}
	}

	AABB getAABB() const {
		return m_shape->getAABB();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "IrradianceMeter[" << endl
			<< "  sampler = " << indent(m_sampler->toString()) << "," << endl
			<< "  film = " << indent(m_film->toString()) << "," << endl
			<< "  medium = " << indent(m_medium.toString()) << "," << endl
			<< "  shutterOpen = " << m_shutterOpen << "," << endl
			<< "  shutterOpenTime = " << m_shutterOpenTime << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(IrradianceMeter, false, Sensor)
MTS_EXPORT_PLUGIN(IrradianceMeter, "Irradiance meter");
MTS_NAMESPACE_END
