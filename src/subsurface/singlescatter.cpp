/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

	This plugin computes single scattering for objects with refractive triangle
	mesh boundaries and was developed by Nicolas Holzschuch

    It is the implementation of the paper
	"Accurate computation of single scattering in participating media with
	refractive boundaries" by Nicolas Holzschuch
	Computer Graphics Forum, Wiley-Blackwell, 2014, pp.1-12. <10.1111/cgf.12517>

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
#include <mitsuba/core/sse.h>
#include <mitsuba/core/ssemath.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN


//////////////////////////////////////////////////////////////////////////////
/// \brief Evaluate the Henyey-Greenstein phase function.
///
Spectrum hg(Float cosTheta, const Spectrum &g) {
	Spectrum temp = Spectrum(1) + g * g + 2 * g * cosTheta;
	return INV_FOURPI * (Spectrum(1) - g * g) / (temp * temp.sqrt());
}

static ref<Mutex> mutex = new Mutex;

/*!\plugin{singlescatter}{Single scattering in participating media}
 * \parameters{
 *     \parameter{material}{\String}{
 *         Name of a material preset, see
 *         \tblref{medium-coefficients}. \default{\texttt{skin1}}
 *     }
 *     \parameter{sigmaA, sigmaS}{\Spectrum}{
 *         Absorption and scattering
 *         coefficients of the medium in inverse scene units.
 *         These parameters are mutually exclusive with \code{sigmaT} and \code{albedo}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{sigmaT, albedo}{\Spectrum}{
 *         Extinction coefficient in inverse scene units
 *         and a (unitless) single-scattering albedo.
 *         These parameters are mutually exclusive with \code{sigmaA} and \code{sigmaS}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{scale}{\Float}{
 *         Optional scale factor that will be applied to the \code{sigma*} parameters.
 *         It is provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}}
 *     
 * }
 *
 * \renderings{
 *    \rendering{The bumpy sphere test scene rendered with amber material}{bumpy_sphere.jpg}
 *    \rendering{The Stanford Bunny rendered with a translucent material}{bunny_single.jpg}
 * }

 * This plugin implements the single scattering model in participating media
 * \cite{Holzschuch2015}.  
 *
 * There are two different ways of configuring the medium properties.
 * One possibility is to load a material preset
 * using the \code{material} parameter---see \tblref{medium-coefficients}
 * for details. Alternatively, when specifying parameters by hand, they
 * can either be provided using the scattering and absorption coefficients,
 * or by declaring the extinction coefficient and single scattering albedo
 * (whichever is more convenient). Mixing these parameter initialization
 * methods is not allowed.
 *
 * All scattering parameters (named \code{sigma*}) should
 * be provided in inverse scene units. For instance, when a world-space
 * distance of 1 unit corresponds to a meter, the scattering coefficents must
 * be in units of inverse meters. For convenience, the \code{scale}
 * parameter can be used to correct this. For instance, when the scene is
 * in meters and the coefficients are in inverse millimeters, set
 * \code{scale=1000}.
 *
 * Note that a subsurface integrator can be associated with an \code{id}
 * and shared by several shapes using the reference mechanism introduced in
 * \secref{format}. This can be useful when an object is made up of many
 * separate sub-shapes.
 *
 *
 * \remarks{
 *    \item This plugin implements only the single scattering 
 *    component.
 *    Single scattering comes in two flavors: fast, the default straight line 
 *    approximation (with multiple samples along the line), and slow, that is 
 *    the accurate version described in \cite{Holzschuch2015}.
 *
 *   \item It is quite important that the \code{sigma*} parameters have the right units.
 *   For instance: if the \code{sigmaT} parameter is accidentally set to a value that
 *   is too small by a factor of 1000, the plugin will attempt to create
 *   one million times as many irradiance samples, which will likely cause
 *   the rendering process to crash with an ``out of memory'' failure.
 * }
 */

class SingleScatter : public Subsurface {
public:
	SingleScatter(const Properties &props) : Subsurface(props) {
		/* Single scattering strategy: use fast single scatter? (Jensen) */
		m_fastSingleScatter = props.getBoolean("fastSingleScatter", true);

		/* Single scattering: number of samples along the inside ray ? */
		m_fastSingleScatterSamples = props.getInteger("fssSamples", 2);

		/* Single scattering: use shadow rays? */
		/* This flag only makes sense when debugging, i.e. it should always be
		 * true */
		m_singleScatterShadowRays =
			props.getBoolean("singleScatterShadowRays", true);

		/* Single scattering: compute transmittance? */
		/* This flag only makes sense when debugging, i.e. it should always be
		 * true */
		m_singleScatterTransmittance =
			props.getBoolean("singleScatterTransmittance", true);

		/* Single scattering: number of total internal reflexion? */
		m_singleScatterDepth = props.getInteger("singleScatterDepth", 4);

		/* Get the material parameters: */
		lookupMaterial(props, m_sigmaS, m_sigmaA, m_g);

		/* Test for old stuff, unsupported anymore */
		if (props.hasProperty("specularReflectance"))
			Log(EError, "specularReflectance is not supported. Use a BSDF "
						"plugin instead.");
		if (props.hasProperty("intIOR") || props.hasProperty("extIOR"))
			Log(EError, "intIOR or extIOR are not supported. Use a BSDF plugin "
						"instead.");
	}

	//------------------------------------------------------------------------
	virtual void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF)))
			if (!m_BSDF.get())
				m_BSDF = dynamic_cast<BSDF *>(child);
			else
				Log(EError, "SingleScatter nodes should have a unique BSDF child.");
		else
			Subsurface::addChild(name, child);
	}

	SingleScatter(Stream *stream, InstanceManager *manager)
		: Subsurface(stream, manager) {
		m_BSDF = static_cast<BSDF *>(manager->getInstance(stream));
		m_sigmaS = Spectrum(stream);
		m_sigmaA = Spectrum(stream);
		m_g = Spectrum(stream);
		m_eta = stream->readFloat();
		// Additions for single scatter
		m_fastSingleScatter = stream->readBool();
		m_fastSingleScatterSamples = stream->readInt();
		m_singleScatterShadowRays = stream->readBool();
		m_singleScatterTransmittance = stream->readBool();
		m_singleScatterDepth = stream->readInt();
		configure();
	}

	virtual ~SingleScatter() {}

	void bindUsedResources(ParallelProcess *proc) const {}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
		manager->serialize(stream, m_BSDF.get());
		m_sigmaS.serialize(stream);
		m_sigmaA.serialize(stream);
		m_g.serialize(stream);
		stream->writeFloat(m_eta);
		// Additions for single scatter
		stream->writeBool(m_fastSingleScatter);
		stream->writeInt(m_fastSingleScatterSamples);
		stream->writeBool(m_singleScatterShadowRays);
		stream->writeBool(m_singleScatterTransmittance);
		stream->writeInt(m_singleScatterDepth);
	}

	//---------------- Begin set of functions for single scattering --------------------
	Spectrum attenuation(const Spectrum &muT, Float negDistance) const {
		Spectrum result(1.0);
		for (int c = 0; c < SPECTRUM_SAMPLES; ++c)
			if (m_sigmaT[c])
				result[c] = math::fastexp(muT[c] * negDistance);
		return result;
	}

	//------------------------------------------------------------------------
	void boundingConeNormals(const Triangle &tri, Vector &axis, Float &angle,
							 const Vector *normals) const {
		// bounding sphere for the triangle:
		const Vector tN[3] = { normals[tri.idx[0]], normals[tri.idx[1]],
							   normals[tri.idx[2]] };
		const Vector a(tN[1] - tN[0]);
		const Vector b(tN[2] - tN[0]);
		const Float a2 = dot(a, a);
		const Float b2 = dot(b, b);
		const Float da = std::sqrt(a2);
		Float db = std::sqrt(b2);
		const Vector axb = cross(a, b);
		const Float axb2 = dot(axb, axb);
		const Float daxb = std::sqrt(axb2);
		if (axb2 != 0.0f) {
			angle = math::safe_asin(da * db * (a - b).length() / (2 * daxb));
			axis = normalize(tN[0] + cross(a2 * b - b2 * a, axb) / (2 * axb2));
		} else {
			angle = 0;
			axis = tN[0];
		}
	}

	//------------------------------------------------------------------------
	bool aabbSegmentTest(const AABB &aabb, const Point &L, const Point &V1,
						 const Point &V2) const {
		// Is there a ray from anywhere on segment [V1, V2] through any triangle
		// inside the AABB connecting to L?
		// Basic splindle test, but with cones.
		// Bounding sphere of the aabb:
		const BSphere aabbSphere = aabb.getBSphere();
		// Bounding cone for omega_L:
		Vector omegaL = L - aabbSphere.center;
		const Float domegaL = omegaL.length();
		omegaL /= domegaL;
		// Cone must be tangent to bounding sphere.
		if (domegaL < aabbSphere.radius)
			return true; // thetaL = M_PI. All tests will send true.

		const Float thetaL = math::safe_asin(aabbSphere.radius / domegaL);
		Vector omegaV1 = (V1 - aabbSphere.center);
		Vector omegaV2 = (V2 - aabbSphere.center);
		Vector V1V2 = omegaV2 - omegaV1;
		const Float dV1V2 = V1V2.length();
		V1V2 /= dV1V2;
		// Shortest distance between C and [V1 V2]
		Float dV1 = omegaV1.length();
		omegaV1 /= dV1;
		Float dV2 = omegaV2.length();
		omegaV2 /= dV2;
		Float minDist = 0, cosTheta = 0;
		// Numerical stability issues. We can have dV1 >> dV2
		// I pick the smallest
		if (dV1 < dV2) {
			minDist = dV1;
			cosTheta = dot(-omegaV1, V1V2);
			if (dV1 * cosTheta > dV1V2)
				minDist = dV2;
			else if (cosTheta > 0)
				minDist = dV1 * std::sqrt(1 - cosTheta * cosTheta);
		} else {
			minDist = dV2;
			cosTheta = dot(-omegaV2, -V1V2);
			if (dV2 * cosTheta > dV1V2)
				minDist = dV1;
			else if (cosTheta > 0)
				minDist = dV2 * std::sqrt(1 - cosTheta * cosTheta);
		}
		if (minDist < aabbSphere.radius)
			return true; // thetaH = M_PI. All tests will send true.

		const Float thetaV = math::safe_asin(aabbSphere.radius / minDist);
		// omegaV is included in a cone whose direction is varying from omegaV1
		// to omegaV2, angle thetaV
		// omegaL is in the cone whose direction is omegaL, angle thetaL.
		// Spindle test:
		// If omegaV is inside the cone of axis -omegaL, spindleAngle,
		// then there can be a spindle compatible intersection
		const Float spindleAngle =
			Float(0.5f * M_PI - math::safe_asin(m_invEta) + thetaL + thetaV);
		if (spindleAngle < M_PI) {
			const Float cosSpindle = std::cos(spindleAngle);
			const Float value0 = dot(omegaV2, -omegaL);
			const Float value1 = dot(omegaV1, -omegaL);
			// We're looking for the smallest angle, thus the largest cosine
			Float cosMax = (value0 > value1 ? value0 : value1);
			if (cosMax < cosSpindle) {
				// same sign as the derivative of the cosine
				const Float deriv0 =
					dot(V1V2, cross(omegaV2, cross(omegaL, omegaV2)));
				const Float deriv1 =
					dot(V1V2, cross(omegaV1, cross(omegaL, omegaV1)));
				if ((deriv0 * deriv1 < 0) && (deriv0 > 0)) {
					// cosine has a maximum in [0,1].
					const Vector n = normalize(cross(omegaV2, omegaV1));
					const Vector u = normalize(-omegaL - dot(n, -omegaL) * n);
					cosMax = dot(u, -omegaL) > cosMax ? dot(u, -omegaL) : cosMax;
				}
				if (cosMax < cosSpindle)
					return false;
			}
		}
		return true;
	}

	//------------------------------------------------------------------------
	Float minDistanceV(Float dV1, const Vector &omegaV1, Float dV2,
					   const Vector &omegaV2, Float dV1V2,
					   const Vector &V1V2) const {
		// Numerical stability issues. We can have dV1 >> dV2
		// I pick the smallest
		Float minDist = 0;
		if (dV1 < dV2) {
			minDist = dV1;
			Float cosTheta = dot(-omegaV1, V1V2);
			if (dV1 * cosTheta > dV1V2)
				minDist = dV2;
			else if (cosTheta > 0)
				minDist = dV1 * std::sqrt(1 - cosTheta * cosTheta);
		} else {
			minDist = dV2;
			Float cosTheta = dot(-omegaV2, -V1V2);
			if (dV2 * cosTheta > dV1V2)
				minDist = dV1;
			else if (cosTheta > 0)
				minDist = dV2 * std::sqrt(1 - cosTheta * cosTheta);
		}
		return minDist;
	}

	//------------------------------------------------------------------------
	bool triangleSegmentTest(const Triangle &tri, const Point &L,
							 const Point &V2, const Point &V1,
							 const Point *positions, const Vector *normals,
							 Float &alphaMin, Float &alphaMax) const {
		// Is there a ray from anywhere on segment [V1, V2] through triangle tri
		// connecting to L?
		// Bounding sphere of the triangle:
		alphaMin = 0;
		alphaMax = 1;
		const BSphere triSphere = tri.getBSphere(positions);

		// Bounding cone for omega_L:
		Vector omegaL = L - triSphere.center;
		const Float domegaL = omegaL.length();
		omegaL /= domegaL;

		// Cone must be tangent to bounding sphere.
		if (domegaL < triSphere.radius)
			return true; // thetaL = M_PI. All tests will send true.

		const Float thetaL = math::safe_asin(triSphere.radius / domegaL);
		Vector omegaV1 = (V1 - triSphere.center);
		const Float dV1 = omegaV1.length();
		omegaV1 /= dV1;
		Vector omegaV2 = (V2 - triSphere.center);
		const Float dV2 = omegaV2.length();
		omegaV2 /= dV2;
		Vector V1V2 = V2 - V1;
		Float dV1V2 = V1V2.length();
		V1V2 /= dV1V2;

		// Shortest distance between C and [V1 V2]
		Float minDist = minDistanceV(dV1, omegaV1, dV2, omegaV2, dV1V2, V1V2);
		if (minDist < triSphere.radius)
			return true; // thetaH = M_PI. All tests will send true.

		Float thetaV = math::safe_asin(triSphere.radius / minDist);
		// omegaV is included in a cone whose direction is varying from omegaV1
		// to omegaV2, angle thetaV
		// omegaL is in the cone whose direction is omegaL, angle thetaL.
		// Spindle test:
		// If omegaV is inside the cone of axis -omegaL, spindleAngle,
		// then there is a spindle intersection
		const Float spindleAngle =
			Float(0.5f * M_PI - math::safe_asin(m_invEta) + thetaL + thetaV);
		if (spindleAngle < M_PI) {
			const Float cosSpindle = std::cos(spindleAngle);
			const Float value0 = dot(omegaV2, -omegaL);
			const Float value1 = dot(omegaV1, -omegaL);

			// We're looking for the smallest angle, thus the largest cosine
			Float cosMax = (value0 > value1 ? value0 : value1);
			if (cosMax < cosSpindle) {
				// same sign as the derivative of the cosine
				const Float deriv0 =
					dot(V1V2, cross(omegaV2, cross(omegaL, omegaV2)));
				const Float deriv1 =
					dot(V1V2, cross(omegaV1, cross(omegaL, omegaV1)));
				if ((deriv0 * deriv1 < 0) && (deriv0 > 0)) {
					// cosine has a maximum in [0,1].
					const Vector n = normalize(cross(omegaV2, omegaV1));
					const Vector u = normalize(-omegaL - dot(n, -omegaL) * n);
					cosMax = dot(u, -omegaL) > cosMax ? dot(u, -omegaL) : cosMax;
				}
				if (cosMax < cosSpindle)
					return false;
			}
		}

		// H vector cone.
		// omegaH = omegaV + (1/eta) omegaL
		const Vector omegaL_eta = omegaL * m_invEta;
		Vector omegaH1 = omegaV1 + omegaL_eta;
		Vector omegaH2 = omegaV2 + omegaL_eta;
		// 2 spheres, radius =
		const Float rH = 2 * std::sin(thetaV / 2) + (2 * m_invEta) * std::sin(thetaL / 2);
		// Shortest distance between C and [H1 H2]
		// Again
		// V1V2 hasn't changed.
		Float value1 = omegaH1.length();
		omegaH1 /= value1;
		Float value0 = omegaH2.length();
		omegaH2 /= value0;
		Float minDistH = (value0 < value1) ? value0 : value1;
		if (minDistH < rH)
			return true;

		const Float derivH0 = dot(V1V2, cross(omegaV2, cross(omegaL, omegaV2)));
		const Float derivH1 = dot(V1V2, cross(omegaV1, cross(omegaL, omegaV1)));
		Vector n = normalize(cross(omegaV2, omegaV1));
		if ((derivH0 * derivH1 < 0) && (derivH0 > 0)) {
			// omegaH.length() has a minimum in [0,1]:
			const Vector u = normalize(-omegaL + dot(n, omegaL) * n);
			const Vector h = u + omegaL_eta;
			if (h.length() < minDistH)
				minDistH = h.length();
		}
		// vector H is included in the sweeping cone, angle thetaH, vector
		// omegaV1-omegaV2 + 1/eta omegaL
		if (minDistH < rH)
			return true;
		const Float thetaH = math::safe_asin(rH / minDistH);

		// Now we get the cone bounding the normals:
		Vector omegaN;
		Float thetaN;
		boundingConeNormals(tri, omegaN, thetaN, normals);

		// Now, is there an intersection between (-omegaN, thetaN) and the
		// sweeping cone for omegaH? Equivalent to knowing whether the axis for
		// the H cone is inside the cone (-omegaN, thetaN + thetaH).
		if (thetaH + thetaN > 0.5 * M_PI)
			return true;

		const Float cosCone = std::cos(thetaH + thetaN);
		const Float sinCone = std::sin(thetaH + thetaN);
		value0 = dot(omegaH2, -omegaN);
		value1 = dot(omegaH1, -omegaN);
		const Vector perp0 = cross(omegaL_eta, cross(n, omegaV2)) + n;
		const Vector perp1 = cross(omegaL_eta, cross(n, omegaV1)) + n;
		Float deriv0 = dot(-omegaN, cross(omegaL_eta + omegaV2, perp0));
		Float deriv1 = dot(-omegaN, cross(omegaL_eta + omegaV1, perp1));

		// We are only interested by one branch of the hyperbola (Pa . (-omegaN)
		// > 0)
		// We must cut the other branch.
		Float K = dot(omegaL_eta, -omegaN);
		K = K * K;
		const Float a0 = dV1V2 * dot(-V1V2, -omegaN);
		const Float a1 = dV2 * dot(omegaV2, -omegaN);
		const Float a = a0 * a0 - K * dV1V2 * dV1V2;
		const Float b = a0 * a1 - K * dV2 * dV1V2 * dot(omegaV2, -V1V2);
		const Float c = a1 * a1 - K * dV2 * dV2;
		Float delta = b * b - a * c;
		Float al0 = 0, al1 = 1;
		if (delta > 0) {
			delta = math::safe_sqrt(delta);
			if (a > 0) {
				al0 = (-b - delta) / a;
				al1 = (-b + delta) / a;
			} else {
				al1 = (-b - delta) / a;
				al0 = (-b + delta) / a;
			}
			bool btw0 =
				(al0 > Epsilon) && (al0 < 1 - Epsilon) &&
				(fabs(dot(normalize(omegaL_eta +
									normalize(al0 * dV1 * omegaV1 +
											  (1 - al0) * dV2 * omegaV2)), -omegaN)) < Epsilon);
			bool btw1 =
				(al1 > Epsilon) && (al1 < 1 - Epsilon) &&
				(fabs(dot(normalize(omegaL_eta +
									normalize(al1 * dV1 * omegaV1 +
											  (1 - al1) * dV2 * omegaV2)), -omegaN)) < Epsilon);
			if (btw0 && btw1) {
				// Difficult case. Would require splitting the interval in 2.
				// First, see if we can get rid of one interval:
				if ((value0 < cosCone) && (deriv0 > 0)) {
					btw0 = false;
					value0 = 0;
					deriv0 = -1;
				}
				if ((value1 < cosCone) && (deriv1 < 0)) {
					btw1 = false;
					value1 = 0;
					deriv1 = +1;
				}
			}
			if (btw0) {
				if (btw1) {
					// 2 roots = 1 or 2 intervals.
					if (value0 < 0) {
						alphaMin = al0;
						value0 = 0;
						deriv0 = -1;
						alphaMax = al1;
						value1 = 0;
						deriv1 = +1;
					} else {
						// Really two intervals of interest.
						alphaMin = 0;
						alphaMax = 1;
						return true;
					}
				} else {
					// 1 root, btw0
					if (value0 > 0) {
						alphaMax = al0;
						value1 = 0;
						deriv1 = +1;
					} else {
						alphaMin = al0;
						value0 = 0;
						deriv0 = -1;
					}
				}
			} else if (btw1) {
				// 1 root, btw1
				if (value0 > 0) {
					alphaMax = al1;
					value1 = 0;
					deriv1 = +1;
				} else {
					alphaMin = al1;
					value0 = 0;
					deriv0 = -1;
				}
			} else {
				// 0 roots between 0 and 1.
				if ((value0 < 0) && (value1 < 0))
					return false;
			}
		} else {
			// only one branch. Negative or positive ?
			if ((value0 < 0) && (value1 < 0))
				return false;
		}

		// Is the origin inside this branch of the hyperbola?
		const Float lambda = dot(omegaL_eta, n) / dot(omegaN, n);
		const Vector vn = omegaN * lambda - omegaL_eta;
		bool originInside = (lambda <= 0) && (dot(vn, vn) <= 1.0);
		// Done splitting
		if (value0 > cosCone) {
			return true;
		}
		if (value1 > cosCone) {
			return true;
		}
		// at least one end below threshold
		if ((deriv0 * deriv1 > 0) || (deriv0 > 0))
			return false;
		do {
			// searching for the minimum
			const Float alpha = 0.5f * (alphaMax + alphaMin);
			const Vector omegaVa =
				normalize(alpha * dV1 * omegaV1 + (1 - alpha) * dV2 * omegaV2);
			Vector Ha = omegaL_eta + omegaVa;
			Float dHa = Ha.length();
			Vector Han = Ha / dHa;
			if (dot(-omegaN, Han) > cosCone) {
				return true;
			}
			const Vector perpa = cross(Ha, cross(n, omegaVa));
			const Float deriv = dot(-omegaN, cross(Ha, perpa));
			if (!originInside) {
				const Vector projPerpa =
					normalize(perpa - dot(perpa, omegaN) * omegaN);
				const Float minProjectedDistance =
					(cosCone / dot(Ha, omegaN)) * dot(Ha, projPerpa);
				if (minProjectedDistance > sinCone)
					return false;
			}
			if (deriv * deriv0 < 0) {
				deriv1 = deriv;
				alphaMax = alpha;
			} else {
				deriv0 = deriv;
				alphaMin = alpha;
			}
		} while (alphaMax - alphaMin > 0.01);
		return false;
	}

	//------------------------------------------------------------------------
	Vector fAndJacobian(const Vector &params, const Matrix3x3 &JP,
						const Matrix3x3 &JN, const Vector &PV0,
						const Vector &PL0, const Vector &N0,
						Matrix3x3 &J) const {
		// Given parameters, compute the function (h+n) and its Jacobian (J)
		// Input parameters: params: x = distance along the ray, (y,z) = homog.
		// coords on the triangle
		// JP = (dInternal, P2 - P0, P2 -P1) Jacobian of PV (also of PL)
		// JN = (0, N0 - N2, N1 - N2) Jacobian of N
		Matrix3x3 JL = JP;
		JL.m[0][0] = 0;
		JL.m[1][0] = 0;
		JL.m[2][0] = 0;
		Vector omegaV = PV0 + JP * params;
		Vector omegaL = PL0 + JL * params;
		Vector omegaN = N0 + JN * params;
		Float dL = omegaL.length();
		omegaL /= dL;
		Float dN = omegaN.length();
		omegaN /= dN;
		Float dV = omegaV.length();
		if (dV < 1e-7f) {
			const Float bp = 2 * dot(omegaN, omegaL) * m_invEta;
			const Float delta = bp * bp + 1 - (m_invEta * m_invEta);
			Float x = -bp + math::safe_sqrt(delta);
			omegaV = -omegaL * m_invEta - x * omegaN;
			dV = 1.0f;
		} else
			omegaV /= dV;
		Vector lineVec = JN.preMult(omegaN);
		// Jacobian of unit vector N
		const Matrix3x3 Jon =
			(JN - Matrix3x3(omegaN * lineVec.x, omegaN * lineVec.y,
							omegaN * lineVec.z)) /
			dN;
		// Jacobian of unit vector omegaL
		lineVec = JL.preMult(omegaL);
		const Matrix3x3 Jol =
			(JL - Matrix3x3(omegaL * lineVec.x, omegaL * lineVec.y,
							omegaL * lineVec.z)) / dL;
		// Jacobian of unit vector omegaV
		lineVec = JP.preMult(omegaV);
		const Matrix3x3 Jov =
			(JP - Matrix3x3(omegaV * lineVec.x, omegaV * lineVec.y,
							omegaV * lineVec.z)) / dV;
		Vector H = m_eta * omegaV + omegaL;
		Float dH = H.length();
		H /= dH;
		const Matrix3x3 JH = m_eta * Jov + Jol;
		lineVec = JH.preMult(H);
		J = (JH - Matrix3x3(H * lineVec.x, H * lineVec.y, H * lineVec.z)) / dH +
			Jon;
		return H + omegaN;
	}

	//------------------------------------------------------------------------
	Spectrum contributionFromThatPoint(const Vector &paramP, const Vector &dPdu,
									   const Vector &dPdv, const Vector &dNsdu,
									   const Vector &dNsdv, const Point &L,
									   const Point &V0, const Vector &dInternal,
									   const Point &P0, const Vector tN[3],
									   const Vector &Ng, Float a11, Float a12,
									   Float a22, const Spectrum &inputSpectrum,
									   const Scene *scene,
									   Float time = 0.0f) const {
		const Point Pc = P0 + paramP[1] * dPdu + paramP[2] * dPdv;
		const Point V = V0 + paramP[0] * dInternal;
		Vector omegaV = V - Pc;
		Float domegaV = omegaV.length();
		omegaV /= domegaV;
		Vector omegaL = L - Pc;
		Float domegaL = omegaL.length();
		omegaL /= domegaL;
		Spectrum result = inputSpectrum;

		if (m_singleScatterShadowRays) {
			// Shadow test 1: is the outgoing point visible from the light
			// source?
			const Ray shadow1 = Ray(Pc, omegaL, ShadowEpsilon,
									domegaL * (1 - ShadowEpsilon), time);
			if (scene->rayIntersect(shadow1)) {
				return Spectrum(0.0f);
			}
		}

		Vector Ns = (tN[0] + paramP[1] * dNsdu + paramP[2] * dNsdv);
		Float idNs = 1.0f / Ns.length();
		Ns *= idNs;
		const Float cosThetaL = dot(omegaL, Ns);
		const Float cosThetaV = dot(omegaV, Ns);

		/* Fresnel transmittance at the new position */
		const Float F = fresnelDielectricExt(cosThetaL, m_eta);

		/* Evaluate the Henyey-Greenstein model */
		const Float cosThetaInternal = dot(omegaV, dInternal);
		Spectrum phase = hg(cosThetaInternal,
							m_g); // reproduces results with +cosThetaInternal.
		result *= (1 - F) * phase;
		result *= m_sigmaS * attenuation(m_sigmaT, -(paramP[0] + domegaV));
		// computing D with ray differentials as in [Walter 2009]
		// u_p and u_s are omega'_V
		const Float mu = cosThetaL + m_eta * cosThetaV;
		// u_p : perpendicular vector
		const Vector u_p = normalize(cross(omegaV, Ns));
		const Vector dPdu_p =
			domegaV * (u_p - (dot(u_p, Ng) / dot(omegaV, Ng)) * omegaV);

		// Normal derivatives are OK
		const Float dudu_p =
			(a22 * dot(dPdu_p, dPdu) - a12 * dot(dPdu_p, dPdv));
		const Float dvdu_p =
			(-a12 * dot(dPdu_p, dPdu) + a11 * dot(dPdu_p, dPdv));
		const Float dwdu_p = -dudu_p - dvdu_p;
		const Vector dNdu_p = dwdu_p * tN[0] + dudu_p * tN[1] + dvdu_p * tN[2];
		const Vector dNndu_p = dNdu_p * idNs - dot(Ns, dNdu_p * idNs) * Ns;
		const Float dmudu_p =
			m_eta * (mu / cosThetaL) * (dot(-u_p, Ns) + dot(omegaV, dNndu_p));
		const Vector domegaLdu_p = m_eta * u_p + dmudu_p * Ns + mu * dNndu_p;
		const Vector L_p =
			dPdu_p - dot(dPdu_p, omegaL) * omegaL + domegaL * domegaLdu_p;

		// u_s : parallel vector
		const Vector u_s = normalize(cross(u_p, omegaV));
		const Vector dPdu_s =
			domegaV * (u_s - (dot(u_s, Ng) / dot(omegaV, Ng)) * omegaV);
		// Normal derivatives
		const Float dudu_s =
			(a22 * dot(dPdu_s, dPdu) - a12 * dot(dPdu_s, dPdv));
		const Float dvdu_s =
			(-a12 * dot(dPdu_s, dPdu) + a11 * dot(dPdu_s, dPdv));
		const Float dwdu_s = -dudu_s - dvdu_s;
		const Vector dNdu_s = dwdu_s * tN[0] + dudu_s * tN[1] + dvdu_s * tN[2];
		const Vector dNndu_s = dNdu_s * idNs - dot(Ns, dNdu_s * idNs) * Ns;
		const Float dmudu_s =
			m_eta * (mu / cosThetaL) * (dot(-u_s, Ns) + dot(omegaV, dNndu_s));
		const Vector domegaLdu_s = m_eta * u_s + dmudu_s * Ns + mu * dNndu_s;
		const Vector L_s =
			dPdu_s - dot(dPdu_s, omegaL) * omegaL + domegaL * domegaLdu_s;
		const Float D = cross(L_p, L_s).length();

		// For debug/explanation only: result without the ray-differentials
		// D = (domegaV + m_eta * domegaL) * (std::abs(cosThetaL/cosThetaV)*domegaV
		// + std::abs(cosThetaV/cosThetaL)*m_eta*domegaL);
		return result / D;
	}

	//------------------------------------------------------------------------
	bool findZero(Vector &params, const Matrix3x3 &JP, const Matrix3x3 &JN,
				  const Vector &PV0, const Vector &PL0, const Vector &N0,
				  Matrix3x3 &J, int axis, const Vector2 &limits) const {
		// axis = 0: params.y & params.z moving, params.x fixed
		// axis = 1: params.x & params.z moving, params.y fixed (most likely = 0)
		// axis = 2: params.x & params.y moving, params.z fixed (most likely = 0)
		// axis = 3: params.x moving, 1 - params.y - params.z = 0;
		// This function works with Pseudo-Inverse
		// if ((params[0] < limits.x) || (params[0] > limits.y)) Log(EInfo, "#
		// Outside of limits : %f < %f < %f", limits.x, params[0], limits.y);
		// TODO: lower back to 10 + 1e-3. See if visible degrade in quality /
		// improve in time
		// 40 + 1e-6 = 32 mn
		// 10 + 1e-3 = 27 mn
		// stepMax = 0.8 -->
		int maxNumTests = 10; //  low value = fast computations; 20 = 4.6 mn, 90
							  //  = 6.7 mn. Try 10?
		Float limit_df = 1e-3; // precision for searching interval limits
		if (axis == 0)
			limit_df = 1e-6; // We need more accuracy for actual points

		const Float stepMaxMax =
			1.5f; // Given that the triangle is [0,1]^2, 1.5 is *huge*
		Float stepMax =
			0.8f; // Variable. Anything in the range 0.5 - 0.9 is faster.
		int numTests = 0;
		Vector oldParams = params;
		Float old_df = 2.0;
		bool foundP = false;
		do {
			const Vector f = fAndJacobian(params, JP, JN, PV0, PL0, N0, J);
			const Float df = f.length();
			foundP = df < limit_df;
			if (df > old_df) {
				// We did not improve. Go back
				params = oldParams;
				stepMax /= 2.0;
			} else if (!foundP) {
				// Jacobian is the 3x2 matrix (dfdx, dfdv)
				// We compute its pseudo-inverse (J^tJ)^{-1} J^t
				Float a = 0, b = 0, c = 0;
				Float Jm0 = 0, Jm1 = 0, Jm2 = 0;
				switch (axis) {
					case 0:
						a = J.m[0][1] * J.m[0][1] + J.m[1][1] * J.m[1][1] +
							J.m[2][1] * J.m[2][1];
						b = J.m[0][1] * J.m[0][2] + J.m[1][1] * J.m[1][2] +
							J.m[2][1] * J.m[2][2];
						c = J.m[0][2] * J.m[0][2] + J.m[1][2] * J.m[1][2] +
							J.m[2][2] * J.m[2][2];
						break;
					case 1:
						a = J.m[0][0] * J.m[0][0] + J.m[1][0] * J.m[1][0] +
							J.m[2][0] * J.m[2][0];
						b = J.m[0][0] * J.m[0][2] + J.m[1][0] * J.m[1][2] +
							J.m[2][0] * J.m[2][2];
						c = J.m[0][2] * J.m[0][2] + J.m[1][2] * J.m[1][2] +
							J.m[2][2] * J.m[2][2];
						break;
					case 2:
						a = J.m[0][0] * J.m[0][0] + J.m[1][0] * J.m[1][0] +
							J.m[2][0] * J.m[2][0];
						b = J.m[0][0] * J.m[0][1] + J.m[1][0] * J.m[1][1] +
							J.m[2][0] * J.m[2][1];
						c = J.m[0][1] * J.m[0][1] + J.m[1][1] * J.m[1][1] +
							J.m[2][1] * J.m[2][1];
						break;
					case 3:
						Jm0 = J.m[0][1] - J.m[0][2];
						Jm1 = J.m[1][1] - J.m[1][2];
						Jm2 = J.m[2][1] - J.m[2][2];
						a = J.m[0][0] * J.m[0][0] + J.m[1][0] * J.m[1][0] +
							J.m[2][0] * J.m[2][0];
						b = J.m[0][0] * Jm0 + J.m[1][0] * Jm1 + J.m[2][0] * Jm2;
						c = Jm0 * Jm0 + Jm1 * Jm1 + Jm2 * Jm2;
						break;
					default:
						return false;
				}
				const Matrix2x2 JtJ(a, b, b, c);
				Vector2 Jtf;
				switch (axis) {
					case 0:
						Jtf.x = J.m[0][1] * f[0] + J.m[1][1] * f[1] +
								J.m[2][1] * f[2];
						Jtf.y = J.m[0][2] * f[0] + J.m[1][2] * f[1] +
								J.m[2][2] * f[2];
						break;
					case 1:
						Jtf.x = J.m[0][0] * f[0] + J.m[1][0] * f[1] +
								J.m[2][0] * f[2];
						Jtf.y = J.m[0][2] * f[0] + J.m[1][2] * f[1] +
								J.m[2][2] * f[2];
						break;
					case 2:
						Jtf.x = J.m[0][0] * f[0] + J.m[1][0] * f[1] +
								J.m[2][0] * f[2];
						Jtf.y = J.m[0][1] * f[0] + J.m[1][1] * f[1] +
								J.m[2][1] * f[2];
						break;
					case 3:
						Jtf.x = J.m[0][0] * f[0] + J.m[1][0] * f[1] +
								J.m[2][0] * f[2];
						Jtf.y = Jm0 * f[0] + Jm1 * f[1] + Jm2 * f[2];
						break;
					default:
						return false;
				}
				Matrix2x2 JtJinv(0.0);
				JtJ.invert(JtJinv);
				Vector2 step = JtJinv * Jtf;
				const Float dStep = step.length();
				if (axis > 0) {
					// Sometimes the function behaves erratically near the limit
					// (because it's 2 branches of a hyperbola)
					// Must make sure we don't cross the limit, even if it slows
					// us down.
					Float dist = -1;
					if (params[0] - step[0] > limits.y) {
						dist = 0.75f * (limits.y - params[0]); // /(-step[0]);
						if (dist < stepMax)
							stepMax = dist;
					}
					if (params[0] - step[0] < limits.x) {
						dist = 0.75f * (params[0] - limits.x); // /step[0];
						if (dist < stepMax)
							stepMax = dist;
					}
					if (stepMax > stepMaxMax)
						stepMax = stepMaxMax;
				}
				if (dStep > stepMax)
					step *= (stepMax / dStep);
				oldParams = params;
				old_df = df;
				switch (axis) {
					case 0:
						params[1] -= step[0];
						params[2] -= step[1];
						break;
					case 1:
						params[0] -= step[0];
						params[2] -= step[1];
						break;
					case 2:
						params[0] -= step[0];
						params[1] -= step[1];
						break;
					case 3:
						params[0] -= step[0];
						params[1] -= step[1];
						params[2] += step[1];
						break;
					default:
						return false;
				}
				// Try to augment step...
				stepMax *= 1.2;
				if (stepMax > stepMaxMax)
					stepMax = stepMaxMax;
			}
			numTests++;
		} while (!foundP && (stepMax > 1e-6) && (numTests < maxNumTests) &&
				 (params.y < 2) && (params.y > -1) && (params.z < 2) &&
				 (params.z > -1));
		return foundP;
	}

	//------------------------------------------------------------------------
	void findPointOnTriangleBoundary(int axis, Vector &params,
									 Vector startingPoint[3], bool toTest[3],
									 bool &foundP, bool &isMax, bool &isInside,
									 Vector &pmin, Vector &pmax,
									 const Matrix3x3 &JP, const Matrix3x3 &JN,
									 const Vector &PV0, const Vector &PL0,
									 const Vector &N0, Matrix3x3 &J,
									 Vector2 &limits) const {
		params = startingPoint[axis - 1];
		foundP = findZero(params, JP, JN, PV0, PL0, N0, J, axis, limits);
		if (!foundP) {
			isInside = false;
			return;
		}
		Vector gradient;
		{
			Vector w1(J.m[0][0], J.m[0][1], J.m[0][2]);
			Vector w2(J.m[1][0], J.m[1][1], J.m[1][2]);
			Vector w3(J.m[2][0], J.m[2][1], J.m[2][2]);
			Vector w13 = cross(w1, w3);
			const Float dw13 = w13.length();
			Vector w23 = cross(w2, w3);
			const Float dw23 = w23.length();
			Vector w12 = cross(w1, w2);
			const Float dw12 = w12.length();
			if ((dw13 > dw23) && (dw13 > dw12))
				gradient = w13 / dw13;
			else if ((dw12 > dw23) && (dw12 > dw13))
				gradient = w12 / dw12;
			else
				gradient = w23 / dw23;
		}

		// Occasionnally, isMax gives the wrong answer because g.x == 0
		// We will catch this later.
		switch (axis) {
			case 1:
				isInside = (params[2] >= 0.0) && (params[2] <= 1.0);
				isMax = (-gradient.y * gradient.x > 0);
				break;
			case 2:
				isInside = (params[1] >= 0.0) && (params[1] <= 1.0);
				isMax = (-gradient.z * gradient.x > 0);
				break;
			case 3:
				isInside = (params[1] >= 0.0) && (params[1] <= 1.0);
				isMax = ((gradient.y + gradient.z) * gradient.x > 0);
				break;
			default:
				return;
		}
		if (isMax) {
			if (isInside)
				pmax = params;
			if (params.x < limits.y)
				limits.y = params.x;
		} else {
			if (isInside)
				pmin = params;
			if (params.x > limits.x)
				limits.x = params.x;
		}
		// Now use the gradient to find new starting points
		if (axis != 1) {
			const Float lam = -params.y / gradient.y;
			const Vector sp = params + lam * gradient;
			if ((sp[2] >= -0.5f) && (sp[2] <= 1.5f)) {
				toTest[0] = true;
				startingPoint[0] = sp;
				if (startingPoint[0][0] > limits.y)
					startingPoint[0][0] =
						limits.x + 0.8f * (limits.y - limits.x);
				if (startingPoint[0][0] < limits.x)
					startingPoint[0][0] =
						limits.x + 0.2f * (limits.y - limits.x);
			} else {
				// We didn't find a starting point, but we update the depth
				// parameter on the existing one
				if ((startingPoint[0][0] > limits.y) ||
					(startingPoint[0][0] < limits.x))
					startingPoint[0][0] = 0.5f * (limits.x + limits.y);
			}
		}
		if (axis != 2) {
			const Float lam = -params.z / gradient.z;
			const Vector sp = params + lam * gradient;
			if ((sp[1] >= -0.5f) && (sp[1] <= 1.5f)) {
				toTest[1] = true;
				startingPoint[1] = sp;
				if (startingPoint[1][0] > limits.y)
					startingPoint[1][0] =
						limits.x + 0.8f * (limits.y - limits.x);
				if (startingPoint[1][0] < limits.x)
					startingPoint[1][0] =
						limits.x + 0.2f * (limits.y - limits.x);
			} else {
				// We didn't find a starting point, but we update the depth
				// parameter on the existing one
				if ((startingPoint[1][0] > limits.y) ||
					(startingPoint[1][0] < limits.x))
					startingPoint[1][0] = 0.5f * (limits.x + limits.y);
			}
		}
		if (axis != 3) {
			const Float lam =
				(1 - params.y - params.z) / (gradient.y + gradient.z);
			const Vector sp = params + lam * gradient;
			if ((sp[2] >= -0.5f) && (sp[2] <= 1.5f)) {
				toTest[2] = true;
				startingPoint[2] = sp;
				if (startingPoint[2][0] > limits.y)
					startingPoint[2][0] =
						limits.x + 0.8f * (limits.y - limits.x);
				if (startingPoint[2][0] < limits.x)
					startingPoint[2][0] =
						limits.x + 0.2f * (limits.y - limits.x);
			} else {
				// We didn't find a starting point, but we update the depth
				// parameter on the existing one
				if ((startingPoint[2][0] > limits.y) ||
					(startingPoint[2][0] < limits.x))
					startingPoint[2][0] = 0.5f * (limits.x + limits.y);
			}
		}
	}

	//------------------------------------------------------------------------
	Spectrum testThisTriangle(const Triangle &tri, const Point &L,
							  const Point &V0, const Vector &dInternal,
							  Float xmin, Float xmax, const Point *positions,
							  const Vector *normals,
							  const Spectrum &inputSpectrum, const Scene *scene,
							  Float time = 0.) const {

		const Point tP[3] = { positions[tri.idx[0]], positions[tri.idx[1]],
							  positions[tri.idx[2]] };
		const Vector dPdu = tP[1] - tP[0];
		const Vector dPdv = tP[2] - tP[0];
		const Matrix3x3 JP(dInternal, -dPdu, -dPdv);
		Vector Ng = cross(dPdu, dPdv);
		const Float lNg = Ng.length();
		if (lNg < 1e-7f)
			return Spectrum(0.0f);
		Ng /= lNg;

		// Sidedness agreement, vs. geometric normal.
		const Vector PL0 = L - tP[0];
		const Vector PV0 = V0 - tP[0];
		// Latest try
		if (dot(PL0, Ng) < 0)
			return Spectrum(0.0f);

		// Triangle has passed all the obvious tests.
		// End points inside the triangle?
		const Float dI = dot(dInternal, Ng);
		const Vector tN[3] = { normals[tri.idx[0]], normals[tri.idx[1]],
							   normals[tri.idx[2]] };
		const Vector dNsdv = (tN[2] - tN[0]);
		const Vector dNsdu = (tN[1] - tN[0]);
		const Matrix3x3 JN(Vector(0, 0, 0), dNsdu, dNsdv);
		int numInside = 0;
		bool extremitiesInside = false;
		Matrix3x3 JPinv(0.0);
		JP.invert(JPinv);
		const Vector paramP = -(JPinv * PV0); // paramP = coord of intersection
											  // ray/triangle. x = coord along
											  // the ray, yz = barycentric.
		if (dI >= 0) {
			if (paramP.x < xmin)
				return Spectrum(0.0f);
			if (paramP.x < xmax)
				xmax = paramP.x;
		} else {
			if (paramP.x > xmax)
				return Spectrum(0.0f);
			if (paramP.x > xmin)
				xmin = paramP.x;
		}
		Vector startingPoint[3]; // Ideal starting points, according to gradient
		bool toTest[3] = { false, false,
						   false }; // should we test this segment?
		if ((paramP.y < 5) && (paramP.y > -4) && (paramP.z < 5) &&
			(paramP.z > -4)) {
			// if we're not too far from the triangle, let's use the gradient at
			// the point of entry
			// for our starting points on each interval.
			const Vector omegaL = normalize(L - V0);
			const Vector omegaN = normalize(tN[0] + JN * paramP);
			const Float bp = 2 * dot(omegaN, omegaL) * m_invEta;
			const Float delta = bp * bp + 1 - (m_invEta * m_invEta);
			const Float Hnorm = -bp + std::sqrt(delta);
			const Vector omegaV = -omegaL * m_invEta - Hnorm * omegaN;
			const Vector entryGradient = normalize(JPinv * omegaV);
			// Ray touches triangle plane at paramP, direction entryGradient in
			// homog coord.
			Float lambda = -paramP.y / entryGradient.y;
			toTest[0] = (lambda > 0);
			if (toTest[0]) {
				startingPoint[0] = paramP + lambda * entryGradient;
				if ((startingPoint[0][2] < -0.5) || (startingPoint[0][2] > 1.5))
					toTest[0] = false;
				else {
					if (startingPoint[0][0] > xmax)
						startingPoint[0][0] = xmin + 0.8f * (xmax - xmin);
					if (startingPoint[0][0] < xmin)
						startingPoint[0][0] = xmin + 0.2f * (xmax - xmin);
				}
			}
			lambda = -paramP.z / entryGradient.z;
			toTest[1] = (lambda > 0);
			if (toTest[1]) {
				startingPoint[1] = paramP + lambda * entryGradient;
				if ((startingPoint[1][1] < -0.5) || (startingPoint[1][1] > 1.5))
					toTest[1] = false;
				else {
					if (startingPoint[1][0] > xmax)
						startingPoint[1][0] = xmin + 0.8f * (xmax - xmin);
					if (startingPoint[1][0] < xmin)
						startingPoint[1][0] = xmin + 0.2f * (xmax - xmin);
				}
			}
			lambda =
				(1 - paramP.y - paramP.z) / (entryGradient.y + entryGradient.z);
			toTest[2] = (lambda > 0);
			if (toTest[2]) {
				startingPoint[2] = paramP + lambda * entryGradient;
				if ((startingPoint[2][1] < -0.5) || (startingPoint[2][1] > 1.5))
					toTest[2] = false;
				else {
					if (startingPoint[2][0] > xmax)
						startingPoint[2][0] = xmin + 0.8f * (xmax - xmin);
					if (startingPoint[2][0] < xmin)
						startingPoint[2][0] = xmin + 0.2f * (xmax - xmin);
				}
			}
		}
		// otherwise, we'll just use the middle of the segment as a wild guess
		if (!toTest[0]) {
			toTest[0] = true;
			startingPoint[0] = Vector(0.5f * (xmin + xmax), 0, 0.5f);
		}
		if (!toTest[1]) {
			toTest[1] = true;
			startingPoint[1] = Vector(0.5f * (xmin + xmax), 0.5f, 0);
		}
		if (!toTest[2]) {
			toTest[2] = true;
			startingPoint[2] = Vector(0.5f * (xmin + xmax), 0.5f, 0.5f);
		}
		Vector2 limits(xmin, xmax);
		Vector pmin(xmin, 0.3, 0.3);
		Vector pmax(xmax, 0.3, 0.3);
		if ((paramP.y >= 0) && (paramP.y <= 1) && (paramP.z >= 0) &&
			(paramP.z <= 1) && (1 - paramP.y - paramP.z >= 0)) {
			numInside++;
			extremitiesInside = true; // if (dI > 0) ?
			if (dI >= 0)
				pmax = paramP;
			else
				pmin = paramP;
		}
		Float a11 = dot(dPdu, dPdu);
		Float a12 = dot(dPdu, dPdv);
		Float a22 = dot(dPdv, dPdv);
		const Float det = a11 * a22 - a12 * a12;
		if (det == 0)
			return Spectrum(0.0f);

		{
			const Float invDet = 1.0f / det;
			a11 *= invDet;
			a12 *= invDet;
			a22 *= invDet;
		}
		bool foundP[3] = { false, false, false };
		bool isMax[3] = { false, false, false };
		bool isInside[3] = { false, false, false };
		Vector params[3];
		// 1. 1st pass, with starting points from enter point
		Matrix3x3 J;
		for (int axis = 0; axis < 3; axis++) {
			if (toTest[axis]) {
				findPointOnTriangleBoundary(
					axis + 1, params[axis], startingPoint, toTest, foundP[axis],
					isMax[axis], isInside[axis], pmin, pmax, JP, JN, PV0, PL0,
					tN[0], J, limits);
				toTest[axis] = false;
				if (limits.x > limits.y)
					return Spectrum(0.0f);
				if (isInside[axis])
					numInside++;
			}
		}
		// 2. 2nd pass, with starting points found from previous pass
		if (numInside < 2) {
			for (int axis = 0; axis < 3; axis++) {
				if (toTest[axis] && !isInside[axis]) {
					findPointOnTriangleBoundary(
						axis + 1, params[axis], startingPoint, toTest,
						foundP[axis], isMax[axis], isInside[axis], pmin, pmax,
						JP, JN, PV0, PL0, tN[0], J, limits);
					toTest[axis] = false;
					if (limits.x > limits.y)
						return Spectrum(0.0f);
					if (isInside[axis])
						numInside++;
				}
			}
		}
		if (numInside == 0)
			return Spectrum(0.0f);

		// If required, a 3rd pass could go there
		if ((numInside < 2) && !extremitiesInside) {
			return Spectrum(0.0f); // should not happen.
		}

		if ((numInside == 2) && ((fabs(limits.x - pmin.x) > 1e-5) ||
								 (fabs(limits.y - pmax.x) > 1e-5))) {
			// Occasionally, a point is identified as a max, even though it's
			// actually a min.
			// We can't detect it in the computation because gradient.x <
			// epsilon. So we fix it here.
			int u1 = 0, u2 = 0;
			if (!isInside[2]) {
				u1 = 0;
				u2 = 1;
			} else if (!isInside[1]) {
				u1 = 0;
				u2 = 2;
			} else if (!isInside[0]) {
				u1 = 1;
				u2 = 2;
			}
			if (params[u1].x < params[u2].x) {
				pmin = params[u1];
				pmax = params[u2];
			} else {
				pmin = params[u2];
				pmax = params[u1];
			}
			limits.x = pmin.x;
			limits.y = pmax.x;
		}
		// Sampling.
		int numSamples =
			1; // for small triangle models, 1 sample per triangle is enough
		Spectrum sum(0.0f);
		if (extremitiesInside && (numInside < 2)) {
			// The ray enters (or leaves) through the triangle, we didn't see
			// the exit, we march along the ray
			// Expensive, but better than nothing
			Float startX = 0, endX = 0, deltaX = 0;
			if (limits.x == 0.0f) {
				deltaX = 0.1f * m_radius;
				startX = limits.x;
				endX = limits.y;
			} else {
				deltaX = -0.1f * m_radius;
				startX = limits.y;
				endX = limits.x;
			}
			Vector testP = paramP; // entry point
			testP.x = startX + 0.5f * deltaX;
			for (Float x = startX + 0.5f * deltaX;
				 (x - endX) * (x - startX) < 0; x += deltaX) {
				if (findZero(testP, JP, JN, PV0, PL0, tN[0], J, 0, limits)) {
					if ((testP.y >= 0) && (testP.y <= 1) && (testP.z >= 0) &&
						(testP.z <= 1) && (1 - testP.y - testP.z >= 0)) {
						Spectrum cfp = contributionFromThatPoint(
							testP, dPdu, dPdv, dNsdu, dNsdv, L, V0, dInternal,
							tP[0], tN, Ng, a11, a12, a22, inputSpectrum, scene,
							time);
						sum += fabs(deltaX) * cfp;
					} else
						break;
				} else
					break;
				testP.x = x + deltaX; // advance along ray by x
			}
		} else {
			Vector deltaP = (pmax - pmin) / Float(numSamples);
			for (int i = 0; i < numSamples; i++) {
				Vector paramP = pmin + deltaP * (i + 0.5f);
				if (findZero(paramP, JP, JN, PV0, PL0, tN[0], J, 0, limits)) {
					if ((paramP.y >= 0) && (paramP.y <= 1) && (paramP.z >= 0) &&
						(paramP.z <= 1) && (1 - paramP.y - paramP.z >= 0)) {
						Spectrum cfp = contributionFromThatPoint(
							paramP, dPdu, dPdv, dNsdu, dNsdv, L, V0, dInternal,
							tP[0], tN, Ng, a11, a12, a22, inputSpectrum, scene,
							time);
						sum += deltaP.x * cfp;
					}
				}
			}
		}
		return sum;
	}

	//------------------------------------------------------------------------
	Spectrum LoSingle(const Scene *scene, Sampler *sampler,
					  const Intersection &its, const Vector &dInternal,
					  int depth, Float z0) const {

		Spectrum result(0.0f);
		if (depth >= m_singleScatterDepth) {
			return result;
		}

		//---- Look at the intersection
		Ray forwardRay(its.p, dInternal, its.time);
		Intersection its2;
		if (EXPECT_NOT_TAKEN(!scene->rayIntersect(forwardRay, its2))) {
			return result; // starting point seems to be outside the object
		}

		// How large is the object?
		const Float thickness = its2.t;

		// 2014-04-22 jDG: Beware of intersecting mesh: we do have a (possible)
		//                 exit transmittance ONLY when hitting the SAME object
		//                 again.
		if (m_singleScatterTransmittance && its2.shape == getShapes()[0]) {
			// 2014-03-10 jDG: Better to use the <<canonic>> way to get
			//                 both the refracted direction AND the attenuation.
			BSDFSamplingRecord bRec(its2, sampler);
			bRec.typeMask = BSDF::EDeltaTransmission;
			Spectrum bsdfAtt = m_BSDF->sample(bRec, sampler->next2D());
			sampler->advance();
			Vector dOutgoing = its2.toWorld(bRec.wo);

			// compute radiance through the object
			if (!bsdfAtt.isZero()) {
				// There is transmitted light ?
				// 2014-04-22 jDG: added indirect surface radiance (so take ALL
				// of them)
				RadianceQueryRecord query(scene, sampler);
				query.newQuery(
					RadianceQueryRecord::ERadiance // All types of radiances
												   // (including us) !
						| RadianceQueryRecord::EIntersection, // Make sure to
															  // compute
															  // intersect
															  // first.
					its.shape->getExteriorMedium());
				query.depth = depth;
				Spectrum refracted = m_integrator->Li(
					RayDifferential(its2.p, dOutgoing, its.time), query);
				refracted *= bsdfAtt;
				refracted *= attenuation(m_sigmaT, -thickness);
				result += refracted;
			}
		}

		// Internal reflection (total or quasi total). Reccursivly call LoSingle
		// again.
		{
			BSDFSamplingRecord bRec(its2, sampler);
			const BSDF *bsdf = 0;
			if (its2.shape == getShapes()[0]) {
				bRec.typeMask =
					BSDF::EReflection; // Evrything except exit rays.
				bsdf = m_BSDF.get();
			} else {
				bRec.typeMask = BSDF::EAll; // Everything.
				bsdf = its2.getBSDF();
			}

			Spectrum bsdfAtt = bsdf->sample(bRec, sampler->next2D());
			sampler->advance();
			Vector dReflect = its2.toWorld(bRec.wo);

			if (!bsdfAtt.isZero()) {
				Spectrum reflected = LoSingle(scene, sampler, its2, dReflect,
											  depth + 1, z0 + thickness);
				reflected *= bsdfAtt;
				// 2014-02-24 jDG: Attenuation on the whole path, at once.
				reflected *= attenuation(m_sigmaT, -thickness);
				result += reflected;
			}
		}

		/* Sample a point on a light source */
		DirectSamplingRecord dRec(its.p, its.time);
		// 2014-04-30 NH: Added the m_eta^2 coefficient, for the light ray enter
		// the material.
		const Spectrum value =
			m_eta * m_eta *
			scene->sampleEmitterDirect(dRec, sampler->next2D(), false);
		sampler->advance();
		if (value.isZero())
			return result;

		const Point L = dRec.p;
		if (m_fastSingleScatter) {
			// Classical SS approximation. Shoot one ray back to the light
			// source.
			// We allow more than one sample along the ray, though.
			const Float sMax = 1 - exp(-(thickness / m_radius));
			const Float dSamples =
				m_fastSingleScatterSamples
					? sMax / Float(m_fastSingleScatterSamples)
					: sMax;

			const Spectrum weight0 =
				(dSamples * m_radius * (dRec.dist * dRec.dist)) * m_sigmaS;

			for (int s = 0; s < m_fastSingleScatterSamples; ++s) {
				Float sample = sampler->next1D() * sMax;
				sampler->advance();
				const Float dist = -math::fastlog(1 - sample) * m_radius;
				const Point V = its.p + dist * dInternal;
				if (EXPECT_NOT_TAKEN(dist > thickness))
					continue;

				/* Sampling weight, because sampling should have been
				 * chanel-dependant */
				Spectrum weight = weight0 * math::fastexp(m_invRadius * dist);

				/* First, connect to the light source */
				Vector VL = L - V;
				Float dVL = VL.length();
				VL /= dVL;
				Ray toTheLight(V, VL, Epsilon, dVL * (1 - ShadowEpsilon), its.time);
				if (!scene->rayIntersect(toTheLight, its2))
					continue;

				const Point PWorld = its2.p;
				/* Make sure that the light source is not occluded from this
				 * position */
				Vector omegaL = L - PWorld;
				Float dL = omegaL.length();
				omegaL /= dL;

				/* shadow ray */
				Ray ray(PWorld, omegaL, Epsilon, dL * (1 - ShadowEpsilon),
						its2.time);
				if (scene->rayIntersect(ray)) {
					continue;
				}

				Vector omegaV = V - PWorld;
				Float dV = omegaV.length();
				omegaV /= dV;

				/* Account for importance sampling wrt. transmittance */
				const Float cosThetaL = dot(omegaL, its2.shFrame.n);
				const Float cosThetaV = dot(omegaV, its2.shFrame.n);
				if (cosThetaL == 0 || cosThetaV == 0)
					continue;

				/* Fresnel transmittance at the new position */
				const Float F = fresnelDielectricExt(cosThetaL, m_eta);

				/* Evaluate the Henyey-Greenstein model */
				Float cosThetaInternal = dot(omegaV, dInternal);
				Spectrum phase;
				phase = hg(cosThetaInternal, m_g) *
						attenuation(m_sigmaT, -(dist + dV));

				const Float D = (dV + m_eta * dL) *
								(std::abs(cosThetaL / cosThetaV) * dV +
								 std::abs(cosThetaV / cosThetaL) * m_eta * dL);
				result += ((1 - F) / D) * phase * value * weight;
			}
		} else {
			// "Slow" single scatter: find all intersections, sample them
			const TriMesh *triMesh = static_cast<const TriMesh *>(its.shape);
			const Point *positions = triMesh->getVertexPositions();
			const Vector *normals = triMesh->getVertexNormals();

			size_t numTriangles = triMesh->getTriangleCount();
			bool *doneThisTriangleBefore = new bool[numTriangles];
			for (size_t i = 0; i < numTriangles; i++) {
				doneThisTriangleBefore[i] = false;
			}
			// scan the kd-tree, cull all nodes not intersecting with segment,
			// keep going.
			struct spindleStackEntry {
				const ShapeKDTree::KDNode *node;
				AABB aabb;
			};
			spindleStackEntry stack[MTS_KD_MAXDEPTH];
			int stackPos = 0;
			const ShapeKDTree::KDNode *currNode = scene->getKDTree()->getRoot();
			AABB currNodeAABB = scene->getKDTree()->getTightAABB();
			while (currNode != NULL) {
				// a) test if bounding sphere of AABB makes it possible for the
				// spindle test
				// for *any* point on the segment.
				if (aabbSegmentTest(currNodeAABB, L, its.p, its2.p)) {
					// b) if yes:
					// b1) is it a leaf?
					if (currNode->isLeaf()) {
						/// Index number format (max 2^32 prims)
						const uint32_t primEnd = currNode->getPrimEnd();
						for (uint32_t entry = currNode->getPrimStart();
							 entry < primEnd; entry++) {
							uint32_t primIdx = scene->getKDTree()->getIndices()[entry];
							uint32_t shapeIdx = scene->getKDTree()->findShape(primIdx);
							if (its.shape ==
								scene->getKDTree()->getShapes()[shapeIdx]) {
								// it belongs to the original shape
								if (EXPECT_TAKEN(scene->getKDTree()->m_triangleFlag[shapeIdx])) {
									if (!doneThisTriangleBefore[primIdx]) {
										doneThisTriangleBefore[primIdx] = true;
										Float alphaMin, alphaMax;
										if (triangleSegmentTest(triMesh->getTriangles()[primIdx],
												L, its.p, its2.p, positions,
												normals, alphaMin, alphaMax)) {
											result += testThisTriangle(triMesh->getTriangles()[primIdx],
												L, its.p, dInternal,
												alphaMin * thickness,
												alphaMax * thickness, positions,
												normals,
												value * (dRec.dist * dRec.dist),
												scene, its.time);
										}
									}
								}
							}
						}
						// Pop from the stack:
						if (stackPos > 0) {
							--stackPos;
							currNode = stack[stackPos].node;
							currNodeAABB = stack[stackPos].aabb;
						} else {
							break;
						}
					} else {
						// b2) if it's not a leaf.
						// add the two children to the stack, and iterate
						Float splitValue = currNode->getSplit();
						int axis = currNode->getAxis();
						stack[stackPos].node = currNode->getLeft();
						stack[stackPos].aabb = currNodeAABB;
						stack[stackPos].aabb.max[axis] = splitValue;
						currNode = currNode->getRight();
						currNodeAABB.min[axis] = splitValue;
						++stackPos;
					}
				} else {
					// Pop from the stack:
					if (stackPos > 0) {
						--stackPos;
						currNode = stack[stackPos].node;
						currNodeAABB = stack[stackPos].aabb;
					} else {
						break;
					}
				}
			}
			delete[] doneThisTriangleBefore;
		}
		return result;
	}
	//---------------- End set of functions for single scattering --------------------

	Spectrum Lo(const Scene *scene, Sampler *sampler, const Intersection &its,
				const Vector &d, int depth) const {
		//---- Initialize intergator and BSDF stuff from the first intersection seen.
		if (!m_integrator) {
			LockGuard lock(mutex);
			if (!m_integrator) {
				SingleScatter *self = const_cast<SingleScatter *>(this);

				self->m_integrator = dynamic_cast<const MonteCarloIntegrator *>(
					scene->getIntegrator());
				if (!m_integrator)
					Log(EError, "Single scatter requires a sampling-based "
								"surface integrator!");
				if (!m_integrator->getClass()->derivesFrom(
						MTS_CLASS(SamplingIntegrator)))
					Log(EError, "Single scatter requires a sampling-based "
								"surface integrator!");
			}
		}

		Spectrum result(0.0f);

		//---- Perform Reflections (if any) ----------------------------------
		{
			BSDFSamplingRecord bRec(its, sampler);
			bRec.typeMask = BSDF::EDeltaReflection;
			Spectrum reflectAttenuation =
				m_BSDF->sample(bRec, sampler->next2D());
			Vector dBounced = its.toWorld(bRec.wo);
			sampler->advance();
			if (!reflectAttenuation.isZero()) {
				RadianceQueryRecord query(scene, sampler);
				query.newQuery(
					RadianceQueryRecord::ERadiance // All radiances sources.
						| RadianceQueryRecord::EIntersection, // Should compute
															  // first
															  // intersect.
					its.shape->getExteriorMedium());
				query.depth = depth + 1;
				RayDifferential ray(its.p, dBounced, its.time);
				result += reflectAttenuation * m_integrator->Li(ray, query);
			}
		}
		//---- Perform refractions (if any) and single scatter ----------------------
		{
			BSDFSamplingRecord bRec(its, sampler);
			bRec.typeMask = BSDF::EDeltaTransmission;
			Spectrum refractAttenuation =
				m_BSDF->sample(bRec, sampler->next2D());
			Vector dInternal = its.toWorld(bRec.wo);
			sampler->advance();

			if (!refractAttenuation.isZero()) {
				result +=
					refractAttenuation *
					LoSingle(scene, sampler, its, dInternal, depth + 1, 0);
			}
		}
		return result;
	}

	void configure() {
		if (!m_BSDF.get()) {
			Log(EError, "Single scatter should have a BSDF child node.");
			m_eta = 1;
		} else
			m_eta = m_BSDF->getEta();
		m_invEta = 1 / m_eta;

		if (m_eta < 1)
			Log(EError, "Unsupported material configuration (intIOR/extIOR < 1)");

		m_sigmaT = m_sigmaA + m_sigmaS;

		/* Find the smallest mean-free path over all wavelengths */
		Spectrum mfp = Spectrum(1.0f) / m_sigmaT;
		m_radius = std::numeric_limits<Float>::max();
		for (int lambda = 0; lambda < SPECTRUM_SAMPLES; lambda++)
			m_radius = std::min(m_radius, mfp[lambda]);
		m_invRadius = 1.0f / m_radius;
	}

	bool preprocess(const Scene *scene, RenderQueue *queue,
					const RenderJob *job, int sceneResID, int cameraResID,
					int _samplerResID) {
		if (!scene->getIntegrator()->getClass()->derivesFrom(
				MTS_CLASS(SamplingIntegrator)))
			Log(EError, "The single scattering pluging requires "
						"a sampling-based surface integrator!");
		return true;
	}

	void wakeup(ConfigurableObject *parent,
				std::map<std::string, SerializableObject *> &params) {}

	void cancel() {}

	MTS_DECLARE_CLASS()
private:
	ref<const MonteCarloIntegrator> m_integrator;

	Float m_radius, m_invRadius;
	Float m_eta, m_invEta;
	Spectrum m_sigmaS, m_sigmaA, m_sigmaT, m_g;
	ref<const BSDF> m_BSDF;

	bool m_fastSingleScatter;
	int m_fastSingleScatterSamples;
	bool m_singleScatterShadowRays;
	bool m_singleScatterTransmittance;
	int m_singleScatterDepth;
};

MTS_IMPLEMENT_CLASS_S(SingleScatter, false, Subsurface)
MTS_EXPORT_PLUGIN(SingleScatter, "Single scattering model");
MTS_NAMESPACE_END
