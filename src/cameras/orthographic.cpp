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
 * Simple orthographic camera model
 * - based on the version in PBRT
 */
class OrthographicCamera : public ProjectiveCamera {
public:
	OrthographicCamera(const Properties &props) 
		: ProjectiveCamera(props) {
		/* Should the smaller image plane side be mapped to normalized 
		   device coordinates in [0,1]? */
		m_mapSmallerSide = props.getBoolean("mapSmallerSide", true);
		/* Calculate the perspective transform */
		m_cameraToScreen = Transform::orthographic(m_nearClip, m_farClip);
		m_cameraToScreenGL = Transform::glOrthographic(m_nearClip, m_farClip);
		if (std::abs(1-m_cameraToWorld(Vector(0, 0, 1)).length()) > Epsilon) 
			Log(EError, "The orthographic camera does not allow non-unity scale factors in the 'Z' direction!");
	}

	OrthographicCamera(Stream *stream, InstanceManager *manager) 
	 : ProjectiveCamera(stream, manager) {
		m_worldToScreen = Transform(stream);
		m_rasterToCamera = Transform(stream);
		m_cameraToRaster = Transform(stream);
		m_areaDensity = stream->readFloat();
		m_rasterToScreen = Transform(stream);
		m_screenToRaster = Transform(stream);
		m_mapSmallerSide = stream->readBool();
	}

	void configure() {
		ProjectiveCamera::configure();
		Vector2i size = m_film->getSize();

		bool mapYToNDC01 = (m_aspect >= 1.0f);
		if (!m_mapSmallerSide)
			mapYToNDC01 = !mapYToNDC01;

		/* Maps from the image plane space to raster space. The
		   smaller 2D coordinate axis on the film will have
		   normalized device coordinates in [0, 1]. Also inverts
		   the Y coordinate
		*/
		if (mapYToNDC01) {
			m_screenToRaster = 
				Transform::scale(Vector((Float) size.x, (Float) size.y, 1.0f))
				* Transform::scale(Vector(1/(2*m_aspect), -0.5f, 1.0f))
				* Transform::translate(Vector(m_aspect, -1.0f, 0));
		} else {
			m_screenToRaster = 
				Transform::scale(Vector((Float) size.x, (Float) size.y, 1.0f))
				* Transform::scale(Vector(0.5f, -0.5f * m_aspect, 1.0f))
				* Transform::translate(Vector(1.0f, - 1 / m_aspect, 0));
		}

		m_rasterToScreen = m_screenToRaster.inverse();

		m_worldToScreen = m_cameraToScreen * m_worldToCamera;
		m_rasterToCamera = m_cameraToScreen.inverse() * m_rasterToScreen;
		m_cameraToRaster = m_rasterToCamera.inverse();
		m_areaDensity = std::abs(4*m_cameraToWorld.det3x3());
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		ProjectiveCamera::serialize(stream, manager);

		m_worldToScreen.serialize(stream);
		m_rasterToCamera.serialize(stream);
		m_cameraToRaster.serialize(stream);
		stream->writeFloat(m_areaDensity);
		m_rasterToScreen.serialize(stream);
		m_screenToRaster.serialize(stream);
		stream->writeBool(m_mapSmallerSide);
	}

	bool needsLensSample() const {
		return false;
	}

	void generateRay(const Point2 &dirSample, const Point2 &lensSample, 
		Float timeSample, Ray &ray) const {
		++cameraRays;
		Point rasterCoords(dirSample.x, dirSample.y, 0);
		Point imageCoords;
		m_rasterToCamera(rasterCoords, imageCoords);

		/* Construct ray in camera space */
		Ray localRay(imageCoords, Vector(0, 0, 1),
			m_shutterOpen + m_shutterOpenTime * timeSample);
		localRay.mint = 0;
		localRay.maxt = m_farClip - m_nearClip;

		/* Transform into world space */
		m_cameraToWorld(localRay, ray);
	}
	
	virtual Point getPosition(const Point2 &sample) const {
		Point rasterCoords(sample.x, sample.y, 0);
		Point imageCoords, result;
		m_rasterToCamera(rasterCoords, imageCoords);
		m_cameraToWorld(imageCoords, result);
		return result;
	}

	bool positionToSample(const Point &p, Point2 &sample) const {
		Point localP, imageP;
		m_worldToCamera(p, localP);
		if (localP.z <= 0.0f) {
			return false;
		}
		m_cameraToRaster(localP, imageP);
		sample.x = imageP.x;
		sample.y = imageP.y;
		const Point2i max(m_film->getCropOffset() + m_film->getCropSize());
		return 
			 sample.x >= m_film->getCropOffset().x 
			&& sample.y >= m_film->getCropOffset().y
			&& sample.x < max.x && sample.y < max.y;
	}
	
	Float areaDensity(const Point2 &p) const {
		return m_areaDensity;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "OrthographicCamera[" << std::endl
			<< "  sampler = " << indent(m_sampler->toString()) << "," << std::endl
			<< "  film = " << indent(m_film->toString()) << "," << std::endl
			<< "  aspect = " << m_aspect << "," << std::endl
			<< "  nearClip = " << m_nearClip << "," << std::endl
			<< "  farClip = " << m_farClip << "," << std::endl
			<< "  shutterOpen = " << m_shutterOpen << "," << std::endl
			<< "  shutterClose = " << m_shutterClose << "," << std::endl
			<< "  areaDensity = " << m_areaDensity << "," << std::endl
			<< "  cameraToWorld = " << indent(m_cameraToWorld.toString()) << "," << std::endl
			<< "  cameraToScreen = " << indent(m_cameraToScreen.toString()) << "," << std::endl
			<< "  rasterToCamera = " << indent(m_rasterToCamera.toString()) << std::endl
			<< "]";
		return oss.str();
	}

	Transform getGLProjectionTransform(const Point2 &jitter) const {
		Log(EError, "Not implemented!");
		return m_cameraToScreenGL;
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_worldToScreen;
	Transform m_rasterToCamera, m_cameraToRaster;
	Float m_nearClip, m_farClip;
	Transform m_rasterToScreen, m_screenToRaster;
	bool m_mapSmallerSide;
	Float m_areaDensity;
};

MTS_IMPLEMENT_CLASS_S(OrthographicCamera, false, ProjectiveCamera)
MTS_EXPORT_PLUGIN(OrthographicCamera, "Orthographic camera");
MTS_NAMESPACE_END
