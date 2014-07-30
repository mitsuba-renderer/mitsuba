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

#include <mitsuba/hw/viewer.h>

MTS_NAMESPACE_BEGIN

class CylClip : public Viewer {
public:
	CylClip() : m_red(0.0f), m_blue(0.0f), m_gray(.5f), m_angle(0) {
		m_viewTransform = Transform::lookAt(Point(10*std::sin(m_angle), 0, std::cos(m_angle)*10),
				Point(0, 0, 0), Vector(0, 1, 0));
		m_lineParams = Point2(M_PI/2, 0.28f);
		m_cylPos = Point(0.0f);
		m_device->setFSAA(4);
		m_red[0] = 1.0f;
		m_blue[2] = 1.0f;
		m_showEllipses = false;
		m_showRectangles = false;
		m_showClippedAABB = false;
		m_radius = .2f;
	}

	void windowResized(const DeviceEvent &event) {
		Float aspect = m_device->getAspect();
		m_projTransform = Transform::glPerspective(45, 1e-4, 1e4)
			* Transform::scale(Vector(1.0f, aspect, 1.0f));
	}

	void mouseDragged(const DeviceEvent &event) {
		if (event.getMouseButton() == Device::ELeftButton) {
			m_angle += event.getMouseRelative().x / 100.0f;
			m_viewTransform = Transform::lookAt(Point(10*std::sin(m_angle), 0, std::cos(m_angle)*10),
					Point(0, 0, 0), Vector(0, 1, 0));
		} else if (event.getMouseButton() == Device::ERightButton) {
			m_lineParams += Vector2(
				event.getMouseRelative().x / 500.0f,
				event.getMouseRelative().y / 500.0f
			);
		} else if (event.getMouseButton() == Device::EMiddleButton) {
			m_cylPos += Vector3(
				event.getMouseRelative().x / 300.0f,
				0,
				event.getMouseRelative().y / 300.0f
			);
		}
		redraw();
	}

	/**
	 * Compute the ellipse created by the intersection of an infinite
	 * cylinder and a plane. Returns false in the degenerate case.
	 * Based on:
	 * www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
	 */
	bool intersectCylPlane(Point planePt, Normal planeNrml,
			Point cylPt, Vector cylD, Float radius, Point &center,
			Vector *axes, Float *lengths) const {
		if (absDot(planeNrml, cylD) < Epsilon)
			return false;

		Vector B, A = cylD - dot(cylD, planeNrml)*planeNrml;

		Float length = A.length();
		if (length != 0) {
			A /= length;
			B = cross(planeNrml, A);
		} else {
			coordinateSystem(planeNrml, A, B);
		}

		Vector delta = planePt - cylPt,
			   deltaProj = delta - cylD*dot(delta, cylD);

		Float aDotD = dot(A, cylD);
		Float bDotD = dot(B, cylD);
		Float c0 = 1-aDotD*aDotD;
		Float c1 = 1-bDotD*bDotD;
		Float c2 = 2*dot(A, deltaProj);
		Float c3 = 2*dot(B, deltaProj);
		Float c4 = dot(delta, deltaProj) - radius*radius;

		Float lambda = (c2*c2/(4*c0) + c3*c3/(4*c1) - c4)/(c0*c1);

		Float alpha0 = -c2/(2*c0),
			  beta0 = -c3/(2*c1);

		lengths[0] = std::sqrt(c1*lambda),
		lengths[1] = std::sqrt(c0*lambda);

		center = planePt + alpha0 * A + beta0 * B;
		axes[0] = A;
		axes[1] = B;
		return true;
	}

	AABB intersectCylFace(int axis,
			const Point &min, const Point &max,
			const Point &cylPt, const Vector &cylD) {
		int axis1 = (axis + 1) % 3;
		int axis2 = (axis + 2) % 3;

		Normal planeNrml(0.0f);
		planeNrml[axis] = 1;

		Point ellipseCenter;
		Vector ellipseAxes[2];
		Float ellipseLengths[2];

		AABB aabb;
		if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius,
			ellipseCenter, ellipseAxes, ellipseLengths)) {
			/* Degenerate case -- return an invalid AABB. This is
			   not a problem, since one of the other faces will provide
			   enough information to arrive at a correct clipped AABB */
			return aabb;
		}

		if (m_showEllipses) {
			m_renderer->setColor(m_gray);
			m_renderer->drawEllipse(ellipseCenter,
					ellipseAxes[0]*ellipseLengths[0],
					ellipseAxes[1]*ellipseLengths[1]);
		}

		/* Intersect the ellipse against the sides of the AABB face */
		for (int i=0; i<4; ++i) {
			Point p1, p2;
			p1[axis] = p2[axis] = min[axis];
			p1[axis1] = ((i+1) & 2) ? min[axis1] : max[axis1];
			p1[axis2] = ((i+0) & 2) ? min[axis2] : max[axis2];
			p2[axis1] = ((i+2) & 2) ? min[axis1] : max[axis1];
			p2[axis2] = ((i+1) & 2) ? min[axis2] : max[axis2];

			Point2 p1l(
				dot(p1 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
				dot(p1 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);
			Point2 p2l(
				dot(p2 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
				dot(p2 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);

			Vector2 rel = p2l-p1l;
			Float A = dot(rel, rel);
			Float B = 2*dot(Vector2(p1l), rel);
			Float C = dot(Vector2(p1l), Vector2(p1l))-1;

			Float x0, x1;
			if (solveQuadratic(A, B, C, x0, x1)) {
				m_renderer->setColor(m_red);
				if (x0 >= 0 && x0 <= 1) {
					m_renderer->drawPoint(p1 + (p2-p1) * x0);
					aabb.expandBy(p1+(p2-p1)*x0);
				}
				if (x1 >= 0 && x1 <= 1) {
					m_renderer->drawPoint(p1 + (p2-p1) * x1);
					aabb.expandBy(p1+(p2-p1)*x1);
				}
			}
		}

		ellipseAxes[0] *= ellipseLengths[0];
		ellipseAxes[1] *= ellipseLengths[1];
		m_renderer->setColor(m_blue);

		AABB faceBounds(min, max);
		/* Find the componentwise maxima of the ellipse */
		for (int i=0; i<2; ++i) {
			int j = (i==0) ? axis1 : axis2;
			Float alpha = ellipseAxes[0][j];
			Float beta = ellipseAxes[1][j];
			Float ratio = beta/alpha, tmp = std::sqrt(1+ratio*ratio);
			Float cosTheta = 1/tmp, sinTheta = ratio/tmp;
			Point p1 = ellipseCenter + cosTheta*ellipseAxes[0] + sinTheta*ellipseAxes[1];
			Point p2 = ellipseCenter - cosTheta*ellipseAxes[0] - sinTheta*ellipseAxes[1];

			if (faceBounds.contains(p1)) {
				aabb.expandBy(p1);
				m_renderer->drawPoint(p1);
			}
			if (faceBounds.contains(p2)) {
				aabb.expandBy(p2);
				m_renderer->drawPoint(p2);
			}
		}

		m_renderer->setColor(m_gray);
		if (aabb.isValid() && m_showRectangles)
			m_renderer->drawAABB(aabb);

		return aabb;
	}

	void keyPressed(const DeviceEvent &event) {
		switch (event.getKeyboardKey()) {
			case 'e':
				m_showEllipses = !m_showEllipses;
				break;
			case 'r':
				m_showRectangles = !m_showRectangles;
				break;
			case 'a':
				m_showClippedAABB = !m_showClippedAABB;
				break;
			case '[':
				m_radius *= 1.1;
				break;
			case ']':
				m_radius /= 1.1;
				break;
		}
		redraw();
	}

	bool init(int argc, char **argv) {
		m_renderer->setPointSize(4.0f);
		return true;
	}

	void draw() {
		m_renderer->setDepthTest(true);
		m_renderer->setCamera(m_projTransform.getMatrix(),
			m_viewTransform.inverse().getMatrix());

		AABB aabb(Point(-3, -1, -1), Point(3, 1, 1));
		m_renderer->setColor(Spectrum(0.3f));
		m_renderer->drawAABB(aabb);

		m_renderer->setColor(m_gray);
		Vector cylD(sphericalDirection(m_lineParams.x, m_lineParams.y));

		m_renderer->drawLine(m_cylPos-cylD*1e4, m_cylPos+cylD*1e4);
		AABB clippedAABB;

		clippedAABB.expandBy(intersectCylFace(0,
				Point(aabb.min.x, aabb.min.y, aabb.min.z),
				Point(aabb.min.x, aabb.max.y, aabb.max.z),
				m_cylPos, cylD));

		clippedAABB.expandBy(intersectCylFace(0,
				Point(aabb.max.x, aabb.min.y, aabb.min.z),
				Point(aabb.max.x, aabb.max.y, aabb.max.z),
				m_cylPos, cylD));

		clippedAABB.expandBy(intersectCylFace(1,
				Point(aabb.min.x, aabb.min.y, aabb.min.z),
				Point(aabb.max.x, aabb.min.y, aabb.max.z),
				m_cylPos, cylD));

		clippedAABB.expandBy(intersectCylFace(1,
				Point(aabb.min.x, aabb.max.y, aabb.min.z),
				Point(aabb.max.x, aabb.max.y, aabb.max.z),
				m_cylPos, cylD));

		clippedAABB.expandBy(intersectCylFace(2,
				Point(aabb.min.x, aabb.min.y, aabb.min.z),
				Point(aabb.max.x, aabb.max.y, aabb.min.z),
				m_cylPos, cylD));

		clippedAABB.expandBy(intersectCylFace(2,
				Point(aabb.min.x, aabb.min.y, aabb.max.z),
				Point(aabb.max.x, aabb.max.y, aabb.max.z),
				m_cylPos, cylD));

		m_renderer->setColor(m_gray);

		if (m_showClippedAABB) {
			if (clippedAABB.isValid())
				m_renderer->drawAABB(clippedAABB);
		}

		m_renderer->setDepthTest(false);
		drawHUD(formatString("Cylinder clipping test. LMB-dragging moves the camera, RMB-dragging rotates the cylinder\n"
				"[e] Ellipses: %s\n"
				"[r] Bounding rectangles : %s\n"
				"[a] Clipped AABB : %s",
			m_showEllipses ? "On": "Off",
			m_showRectangles ? "On": "Off",
			m_showClippedAABB ? "On": "Off"
		));
	}

	MTS_DECLARE_UTILITY()
private:
	Transform m_projTransform, m_viewTransform;
	Spectrum m_red, m_blue, m_gray;
	Point2 m_lineParams;
	Point m_cylPos;
	Float m_angle, m_radius;
	bool m_showEllipses;
	bool m_showRectangles;
	bool m_showClippedAABB;
};

MTS_EXPORT_UTILITY(CylClip, "Cylinder clipping test")
MTS_NAMESPACE_END
