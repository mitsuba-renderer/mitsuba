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

#include <mitsuba/render/camera.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter cameraRays("General", "Camera ray generations");

/**
 * Simple environment camera model
 * - based on the version in PBRT
 */
class EnvironmentCamera : public Camera {
public:
	EnvironmentCamera(const Properties &props) 
		: Camera(props) { }

	EnvironmentCamera(Stream *stream, InstanceManager *manager) 
	 : Camera(stream, manager) {
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Camera::serialize(stream, manager);
	}

	void configure() {
		Camera::configure();
		Vector2i filmSize = m_film->getSize();
		m_invResolution = Point2(
			1.0f / filmSize.x,
			1.0f / filmSize.y
		);
	}

	/**
     * Cartesian-to-spherical coordinate mapping with
     * v=0 => Y=1
     */
	Vector squareToSphereY(Float u, Float v) const {
		Float cosTheta = std::cos(v * M_PI), 
			  sinTheta = std::sin(v * M_PI),//std::sqrt(1-cosTheta*cosTheta),
              phi = u * 2 * M_PI,
			  cosPhi = std::cos(phi), sinPhi = std::sin(phi);
		return Vector(
			sinTheta * sinPhi, cosTheta, -sinTheta*cosPhi);
	}

	/* Corresponding reverse mapping */
	Point2 sphereToSquareY(const Vector &d) const {
		Float u = std::atan2(d.x,-d.z) * (0.5f * INV_PI),
			  v = std::acos(std::max((Float) -1.0f, 
				  std::min((Float) 1.0f, d.y))) / M_PI;
		if (u < 0)
			u += 1;
		return Point2(u, v);
	}

	void generateRay(const Point2 &dirSample, const Point2 &lensSample, 
		Float timeSample, Ray &ray) const {
		++cameraRays;

		Float u = dirSample.x * m_invResolution.x,
			  v = dirSample.y * m_invResolution.y;

		Vector direction = squareToSphereY(u, v);
		Point2 uvPrime = sphereToSquareY(direction);

		if ((std::abs(uvPrime.x-u) > Epsilon || std::abs(uvPrime.y-v)>Epsilon) && u < 1 && v < 1 && u > 0 && v > 0)
			cout << uvPrime.toString() << " vs " << u << ", " << v << endl;

		/* Construct ray in camera space */
		Ray localRay(Point(0.0f), direction,
			m_shutterOpen + m_shutterOpenTime * timeSample);

		/* Transform into world space */
		m_cameraToWorld(localRay, ray);
	}
	
	MTS_DECLARE_CLASS()
private:
	Point2 m_invResolution;
};

MTS_IMPLEMENT_CLASS_S(EnvironmentCamera, false, Camera)
MTS_EXPORT_PLUGIN(EnvironmentCamera, "Environment camera");
MTS_NAMESPACE_END
