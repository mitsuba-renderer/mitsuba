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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/core/autodiff.h>
#include <mitsuba/core/statistics.h>
#include <boost/algorithm/string.hpp>

DECLARE_DIFFSCALAR_BASE();

MTS_NAMESPACE_BEGIN

static StatsCounter statsConverged(
		" manifold", "Converged manifold walks", EPercentage);

/*!\plugin{motion}{Motion and specular motion vector integrator}
 * \parameters{
 *     \parameter{time}{\Float}{
 *       Denotes the time stamp of the target frame of the motion vectors. 
 *       The current frame is specified via the sensor's \code{shutterOpen}
 *       and \code{shutterClose} parameters, which should both be set to 
 *       the same value. \default{0}
 *     }
 *     \parameter{time}{\String}{
 *        Path configuration of the desired motion vectors:
 *        \begin{enumerate}[(i)]
 *           \item \textbf{d}: Primary (non-specular) hit points\vspace{-1mm}
 *           \item \textbf{rd}: A non-specular surface seen through a specular reflection\vspace{-1mm}
 *           \item \textbf{ttd}: A non-specular surface seen through a pair of specular refractions\vspace{-1mm}
 *           \item \textbf{trtd}: A non-specular surface seen through a sequence of refraction, reflection, and refraction events.\vspace{-1mm}
 *        \end{enumerate}
 *        etc.
 *     }
 *     \parameter{derivativesOnly}{\Boolean}{
 *       By default, the Manifold Exploration technique is used to accurately
 *       solve the underlying specular flow problem. When this parameter is set to \code{true},
 *       the nonlinear solver is deactivated, and only first-order extrapolations are provided.
 *       \default{\code{false}}
 *     }
 *     \parameter{glossyThreshold}{\Float}{
 *        Threshold on the roughness parameter of reflectance models
 *        to be classified as specular.
 *        \default{\texttt{0}, i.e.~only perfectly specular materials are classified as specular}
 *     }
 *     \parameter{maxTimeSteps}{\Integer}{
 *        Maximum number of temporal sub-steps \default{5}
 *     }
 *     \parameter{maxSpaceSteps}{\Integer}{
 *        Maximum number of spatial sub-steps \default{10}
 *     }
 * }
 * This integrator extracts motion vectors for animated input scenes, alternatively
 * at primary hit points or at hit points observed through sequences of reflective
 * and refractive objects. The first two color components (R and G) of the resulting
 * rendering specify the screen-space motion in 2D pixel coordinates, and the last
 * component (B) denotes the change in distance of the observed 3D point to the
 * camera position. Sometimes a specular path cannot be tracked from one frame to the other,
 * e.g. because it does not exist, or because the solver did not converge. In this case,
 * the pixel color is set to infinity. The images on the following page show
 * motion vectors obtained for a sphere that is moving from the left to the right.
 *
 * \renderings{
 *     \rendering{Input scene at time $t=0$}{integrator_motion_sphere_1}
 *     \rendering{Input scene at time $t=1$}{integrator_motion_sphere_2}
 * }
 * \renderings{
 *     \medrendering{\code{config="d"}}{integrator_motion_path_d}
 *     \medrendering{\code{config="rd"}}{integrator_motion_path_rd}
 *     \medrendering{\code{config="ttd"}}{integrator_motion_path_ttd}
 * }
 * \renderings{
 *     \rendering{\code{config="trtd"}}{integrator_motion_path_trtd}
 *     \rendering{\code{config="trrtd"}}{integrator_motion_path_trrtd}
 * }
 * \clearpage
 *
 * \begin{xml}[caption={Exemplary scene configuration for computing specular motion vectors}]
 * <scene>
 *     <integrator type="motion">
 *         <string name="config" value="ttd"/>
 *         <float name="time" value="1"/>
 *     </integrator>
 *     
 *     <shape type="serialized">
 *         <string name="filename" value="..."/>
 *         <animation name="toWorld">
 *             <transform time="0">
 *                 <!-- Transformation at time zero -->
 *             </transform>
 *             <transform time="1">
 *                 <!-- Transformation at time one -->
 *             </transform>
 *         </animation>
 *          <bsdf type="dielectric"/>
 *     </shape>
 *  
 *     <sensor type="perspective">
 *         <float name="shutterOpen" value="0"/>
 *         <float name="shutterClose" value="0"/>
 *  
 *         <sampler type="ldsampler">
 *             <integer name="sampleCount" value="1"/>
 *             <boolean name="pixelCenters" value="true"/>
 *         </sampler>
 *          
 *         <film type="hdrfilm" id="film">
 *             <string name="pixelFormat" value="rgb"/>
 *             <boolean name="banner" value="false"/>
 *             <rfilter type="box"/>
 *         </film>
 *     </sensor>
 * </scene>
 * \end{xml}
 */

class MotionIntegrator : public SamplingIntegrator {
public:
	typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> EMatrix;
	typedef Eigen::Matrix<Float, Eigen::Dynamic, 1> EVector;
	typedef Eigen::Matrix<Float, 1, 7> Gradient;
	typedef DScalar1<Float, Gradient> DScalar;
	typedef DScalar::DVector3 DVector;

	MotionIntegrator(const Properties &props) : SamplingIntegrator(props) {
		m_time = props.getFloat("time");
		m_config = boost::to_lower_copy(props.getString("config", "d"));
		if (m_config.length() == 0)
			Log(EError, "Path configuration string must have at least one entry!");
		if (m_config[m_config.length()-1] != 'd')
			Log(EError, "Configuration string must end with 'd'!");

		m_derivativesOnly = props.getBoolean("derivativesOnly", false);
		m_maxTimeSteps = props.getInteger("maxTimeSteps", 5);
		m_maxSpaceSteps = props.getInteger("maxSpaceSteps", 10);
		m_glossyThreshold = props.getFloat("glossyThreshold", 0);
		m_subSteps = props.getInteger("subSteps", 1);
	}

	MotionIntegrator(Stream *stream, InstanceManager *manager)
	 : SamplingIntegrator(stream, manager) {
		 m_time = stream->readFloat();
		 m_config = stream->readString();
		 m_derivativesOnly = stream->readBool();
		 m_maxTimeSteps = stream->readInt();
		 m_maxSpaceSteps = stream->readInt();
		 m_glossyThreshold = stream->readFloat();
		 m_subSteps = stream->readInt();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		stream->writeFloat(m_time);
		stream->writeString(m_config);
		stream->writeBool(m_derivativesOnly);
		stream->writeInt(m_maxTimeSteps);
		stream->writeInt(m_maxSpaceSteps);
		stream->writeFloat(m_glossyThreshold);
		stream->writeInt(m_subSteps);
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		const Point2 apertureSample(0.5f);
		Point p0, p1;

		if (m_config.length() == 1) {
			Intersection &its = rRec.its;

			// compute motion of environment
			if (!rRec.rayIntersect(r)) {
				BSphere sphere = rRec.scene->getBSphere();
				Float nearT, farT;
				sphere.radius = 1e6;
				sphere.rayIntersect(r, nearT, farT);
				p0 = r(nearT);
			} else {
				p0 = its.p;
			}

			p1 = p0;

			// reproject point according to the motion of the intersected shape
			if (its.isValid() && its.instance != NULL)  {
				Intersection its2(its);
				its2.adjustTime(m_time);
				p1 = its2.p;
			}
		} else {
			std::vector<Intersection> source, target, temp, temp2;
			Ray ray(r);

			/* Trace an initial light path with the given configuration*/
			if (!tracePath(rRec, ray, source))
				return Spectrum(0.0f);

			m_scene = rRec.scene;
			p0 = source[1].p;

			int timeIteration = 0;

			Float stepSizeReduction = 1.0f;
			while (true) {
				if (++timeIteration > m_maxTimeSteps)
					return Spectrum(std::numeric_limits<Float>::infinity());

				Float maxStepSize = (m_time-r.time) / m_subSteps;
				Float timeStepSize = std::min((Float) 1.0f, maxStepSize / (m_time-source[0].time));
				timeStepSize *= stepSizeReduction;

				/* Compute updated intersection records for time 'm_time' */
				adjustTime(rRec, apertureSample, source, target, timeStepSize);

				if (timeIteration == 1) {
					bool moved = false;
					for (size_t i=0; i<source.size(); ++i)
						if ((source[i].p-target[i].p).length() > 1e-4f)
							moved = true;

					if (!moved)
						return Spectrum(0.0f);

					if (m_derivativesOnly) {
						p1 = extrapolateTimePoint(source, target);
						break;
					} else {
						statsConverged.incrementBase();
					}
				}

				if (!timeStep(rRec, source, target, temp, temp2)) {
					stepSizeReduction *= 0.5f;
				} else {
					p1 = source[1].p;
					stepSizeReduction = std::min((Float) 1.0f, stepSizeReduction * 2);

					if (std::abs(source[0].time-m_time) < 1e-5f) {
						++statsConverged;
						break;
					}
				}
			}
		}

		const Sensor *sensor = rRec.scene->getSensor();
		DirectSamplingRecord dRec0(p0, r.time), dRec1(p1, m_time);
		sensor->sampleDirect(dRec0, apertureSample);
		sensor->sampleDirect(dRec1, apertureSample);

		/* Step 4: Compute depth difference */
		Float dDelta = dRec1.dist - dRec0.dist;
		Spectrum result(0.0f);
		result.fromLinearRGB(dRec1.uv.x-dRec0.uv.x,
				dRec1.uv.y-dRec0.uv.y,
				std::isfinite(dDelta) ? dDelta : (Float) 0);
		return result;
	}

	bool timeStep(RadianceQueryRecord &rRec, std::vector<Intersection> &source, const std::vector<Intersection> &target, std::vector<Intersection> &temp, std::vector<Intersection> &temp2) const {
		Ray ray = extrapolateTimeRay(source, target);

		if (!tracePath(rRec, ray, temp))
			return false;

		Float error = computeError(temp, target);

		Float spaceStepSize = 1.0f;
		int spaceIteration = 0;
		while (error > 1e-5f) {
			++spaceIteration;

			if (spaceIteration > m_maxSpaceSteps) {
				return false;
			}

			Ray candidateRay = extrapolateSpaceRay(temp, target, spaceStepSize);

			Float candidateError = 0;
			if (!tracePath(rRec, candidateRay, temp2))
				candidateError = std::numeric_limits<Float>::infinity();
			else
				candidateError = computeError(temp2, target);

			if (candidateError < error) {
				temp = temp2;
				error = candidateError;
				spaceStepSize = std::min((Float) 1.0f, spaceStepSize * 2);
			} else {
				spaceStepSize *= 0.5f;
			}
		}

		source = temp;
		return true;
	}

	bool tracePath(RadianceQueryRecord &rRec, Ray ray, std::vector<Intersection> &intersections) const {
		int depth = 0;

		Intersection its;
		memset(&its, 0, sizeof(Intersection));
		its.p = ray.o;
		its.time = ray.time;

		intersections.clear();
		intersections.push_back(its);

		while (depth < (int) m_config.size()) {
			rRec.scene->rayIntersect(ray, its);

			char interactionType = m_config[depth++];
			if (interactionType == 'd') {
				if (!its.isValid()) {
					BSphere sphere = rRec.scene->getBSphere();
					Float nearT, farT;
					sphere.radius *= 1000;
					bool success = sphere.rayIntersect(ray, nearT, farT);
					Assert(success && nearT < 0 && farT > 0);
					its.p = ray(farT);
					its.geoFrame = its.shFrame = Frame(-ray.d);
					its.dpdu = its.geoFrame.s;
					its.dpdv = its.geoFrame.t;
					its.time = ray.time;
					its.uv = Point2(0.0f);
					its.shape = its.instance = NULL;
				}

				if (its.isValid() && !(its.shape->getBSDF()->getType() & BSDF::EDiffuseReflection) && !its.isEmitter())
					return false;

				intersections.push_back(its);

				break;
			}

			if (!its.isValid())
				return false;

			intersections.push_back(its);

			/* Sample BSDF * cos(theta) */
			const BSDF *bsdf = its.shape->getBSDF();

			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);

			if (interactionType == 'r') {
				if (bsdf->getType() & BSDF::EDeltaReflection) {
					bRec.typeMask = BSDF::EDeltaReflection;
				} else if (bsdf->getType() & BSDF::EGlossyReflection) {
					for (int i=0; i<bsdf->getComponentCount(); ++i)
						if ((bsdf->getType(i) & BSDF::EGlossyReflection) && bsdf->getRoughness(its, i) < m_glossyThreshold)
							bRec.typeMask = BSDF::EGlossyReflection;
				}
			} else if (interactionType == 't') {
				if (bsdf->getType() & BSDF::EDeltaTransmission) {
					bRec.typeMask = BSDF::EDeltaTransmission;
				} else if (bsdf->getType() & BSDF::EGlossyTransmission) {
					for (int i=0; i<bsdf->getComponentCount(); ++i)
						if ((bsdf->getType(i) & BSDF::EGlossyTransmission) && bsdf->getRoughness(its, i) < m_glossyThreshold)
							bRec.typeMask = BSDF::EGlossyTransmission;
				}
			}

			if (bRec.typeMask == BSDF::EAll)
				return false;

			Spectrum bsdfWeight = bsdf->sample(bRec, Point2(0.0f));
			if (bsdfWeight.isZero())
				return false;

			const Vector wo = its.toWorld(bRec.wo);

			ray = Ray(its.p, wo, ray.time);
		}

		return true;
	}

	void adjustTime(const RadianceQueryRecord &rRec, const Point2 &apertureSample, const std::vector<Intersection> &source, std::vector<Intersection> &target, Float timeStepSize) const {
		target = source;

		Float targetTime = (1-timeStepSize) * source[0].time + timeStepSize * m_time;

		const Sensor *sensor = rRec.scene->getSensor();
		DirectSamplingRecord dRec0(Point(0.0f), source[0].time),
							 dRec1(Point(0.0f), targetTime);
		sensor->sampleDirect(dRec0, apertureSample);
		sensor->sampleDirect(dRec1, apertureSample);
		Assert((source[0].p-dRec0.p).lengthSquared() < 1e-6);

		target[0].p = dRec1.p;
		target[0].time = targetTime;

		for (size_t i=1; i<target.size(); ++i)
			target[i].adjustTime(targetTime);
	}

	DVector getVertexPosition(const std::vector<Intersection> &source, const std::vector<Intersection> &target, int i, int rel) const {
		DScalar u(rel*2, 0), v(rel*2+1, 0), time(6, 0);

		return DScalar::vector(source[i].p)
		     + DScalar::vector(source[i].dpdu) * u
		     + DScalar::vector(source[i].dpdv) * v
		     + DScalar::vector(target[i].p-source[i].p) * time;
	}

	void getVertexFrame(const std::vector<Intersection> &source, const std::vector<Intersection> &target, int i, DVector &s, DVector &t, DVector &n) const {
		DScalar u(2, 0), v(3, 0), time(6, 0);
		const BSDF *bsdf = source[i].shape->getBSDF();

		Frame frame, du, dv;
		frame = bsdf->getFrame(source[i]);
		bsdf->getFrameDerivative(source[i], du, dv);

		DVector n0 = DScalar::vector(frame.n)
			+ DScalar::vector(du.n) * u
			+ DScalar::vector(dv.n) * v;
		DVector s0 = DScalar::vector(frame.s)
			+ DScalar::vector(du.s) * u
			+ DScalar::vector(dv.s) * v;
		DVector t0 = DScalar::vector(frame.t)
			+ DScalar::vector(du.t) * u
			+ DScalar::vector(dv.t) * v;

		frame = bsdf->getFrame(target[i]);
		bsdf->getFrameDerivative(target[i], du, dv);

		DVector n1 = DScalar::vector(frame.n)
			+ DScalar::vector(du.n) * u
			+ DScalar::vector(dv.n) * v;
		DVector s1 = DScalar::vector(frame.s)
			+ DScalar::vector(du.s) * u
			+ DScalar::vector(dv.s) * v;
		DVector t1 = DScalar::vector(frame.t)
			+ DScalar::vector(du.t) * u
			+ DScalar::vector(dv.t) * v;

		n = (1-time)*n0 + time*n1;
		s = (1-time)*s0 + time*s1;
		t = (1-time)*t0 + time*t1;
	}

	void assembleMatrix(const std::vector<Intersection> &source, const std::vector<Intersection> &target, EMatrix &M) const {
		DScalar::setVariableCount(7);
		M.resize(2*(source.size()-2), 2*source.size()+1);
		M.setZero();

		for (int i=1; i < (int) source.size()-1; ++i) {
			DVector p_pred = getVertexPosition(source, target, i-1, 0);
			DVector p_cur  = getVertexPosition(source, target, i, 1);
			DVector p_succ = getVertexPosition(source, target, i+1, 2);

			DVector s, t, n;
			getVertexFrame(source, target, i, s, t, n);

			DVector wi = p_pred - p_cur, wo = p_succ - p_cur;
			wi *= inverse(sqrt(wi.dot(wi)));
			wo *= inverse(sqrt(wo.dot(wo)));

			Float eta = source[i].shape->getBSDF()->getEta();

			if (m_config[i-1] == 'r')
				eta = 1;
			else if (wi.dot(n).getValue() < 0)
				eta = 1/eta;

			DVector H = wi + wo * DScalar(eta);
			H *= inverse(sqrt(H.dot(H)));

			Gradient H_s = H.dot(s).getGradient(),
			         H_t = H.dot(t).getGradient();

			M.block<1, 6>((i-1)*2+0, (i-1)*2) = H_s.segment<6>(0);
			M.block<1, 6>((i-1)*2+1, (i-1)*2) = H_t.segment<6>(0);
			M((i-1)*2+0, M.cols()-1) = H_s(6);
			M((i-1)*2+1, M.cols()-1) = H_t(6);
		}
	}

	Point extrapolateTimePoint(const std::vector<Intersection> &source, const std::vector<Intersection> &target) const {
		EMatrix M;
		assembleMatrix(source, target, M);
		EVector b = -M.block(0, 2, (source.size()-2)*2, (source.size()-2)*2).lu().solve(M.col(M.cols()-1));
		return target[1].p + target[1].dpdu * b[0] + target[1].dpdv * b[1];
	}

	Ray extrapolateTimeRay(const std::vector<Intersection> &source, const std::vector<Intersection> &target) const {
		EMatrix M;
		assembleMatrix(source, target, M);

		EVector b = -M.block(0, 2, (source.size()-2)*2, (source.size()-2)*2).lu().solve(M.col(M.cols()-1));
		Point rayOrigin = target[0].p;
		Point rayTarget = target[1].p + target[1].dpdu * b[0] + target[1].dpdv * b[1];
		return Ray(rayOrigin, normalize(rayTarget-rayOrigin), target[1].time);
	}

	Float computeError(const std::vector<Intersection> &source, const std::vector<Intersection> &target) const {
		int last = source.size()-1;
		Float scale = std::max(Epsilon,std::max(std::abs(target[last].p.x), std::max(std::abs(target[last].p.y), std::abs(target[last].p.z))));
		return (target[last].p-source[last].p).length() / scale;
	}

	Ray extrapolateSpaceRay(const std::vector<Intersection> &source, const std::vector<Intersection> &target, Float stepSize) const {
		EMatrix M;
		assembleMatrix(source, target, M);

		Eigen::PartialPivLU<EMatrix> lu = M.block(0, 2, (source.size()-2)*2, (source.size()-2)*2).lu();

		int last = source.size()-1;
		Vector rel = target[last].p-source[last].p,
		       dpdu = source[last].dpdu,
			   dpdv = source[last].dpdv;

		Float b1 = dot(rel, dpdu),
			  b2 = dot(rel, dpdv),
			  a11 = dot(dpdu, dpdu), a12 = dot(dpdu, dpdv),
			  a22 = dot(dpdv, dpdv),
			  det = a11 * a22 - a12 * a12;

		Float invDet = 1.0f / det,
		      du = ( a22 * b1 - a12 * b2) * invDet,
		      dv = (-a12 * b1 + a11 * b2) * invDet;

		EVector b = -(du*lu.solve(M.col(M.cols()-3)) + dv*lu.solve(M.col(M.cols()-2)));

		Point rayTarget = source[1].p + stepSize * (source[1].dpdu * b[0] + source[1].dpdv * b[1]);

		return Ray(source[0].p, normalize(rayTarget-source[0].p), source[1].time);
	}

	std::string toString() const {
		return "MotionIntegrator[]";
	}

	MTS_DECLARE_CLASS()
private:
	Float m_time;
	std::string m_config;
	bool m_derivativesOnly;
	int m_maxSpaceSteps;
	int m_maxTimeSteps;
	int m_subSteps;
	Float m_glossyThreshold;
	mutable const Scene *m_scene;
};

MTS_IMPLEMENT_CLASS_S(MotionIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MotionIntegrator, " motion vector integrator");
MTS_NAMESPACE_END
