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

#include <mitsuba/hw/gizmo.h>

MTS_NAMESPACE_BEGIN

Gizmo::Gizmo() : m_active(false), m_drag(false) { }

Gizmo::~Gizmo() { }

void Gizmo::init(const BSphere &bsphere) {
	m_bsphere = bsphere;
	m_active = true;
	m_drag = false;
}

Transform Gizmo::getTransform() const {
	Vector leg1 = normalize(m_dragStart - m_bsphere.center);
	Vector leg2 = normalize(m_dragEnd - m_bsphere.center);
	Vector rotAxis = normalize(cross(leg1, leg2));
	Float angle = unitAngle(leg1, leg2) * 180/M_PI;
	
	return Transform::translate(Vector(m_bsphere.center)) *
		Transform::rotate(rotAxis, angle) * 
		Transform::translate(-Vector(m_bsphere.center));
}

void Gizmo::startDrag(const Ray &ray) {
	Float nearT, farT;
	if (!m_bsphere.rayIntersect(ray, nearT, farT)
		|| (nearT < 0 && farT < 0)) {
		Log(EWarn, "Gizmo::init(): Internal error (no intersection found)!");
		return;
	} else if (nearT < 0) {
		m_dragStart = ray(farT);
	} else {
		m_dragStart = ray(nearT);
	}
	m_dragEnd = m_dragStart;
	m_drag = true;
}

void Gizmo::dragTo(const Ray &ray, const Camera *camera) {
	Float nearT, farT;
	if (m_bsphere.rayIntersect(ray, nearT, farT)) {
		if (nearT < 0)
			m_dragEnd = ray(farT);
		else
			m_dragEnd = ray(nearT);
	} else {
		Normal n(normalize(m_bsphere.center - camera->getPosition()));
		Float t = (dot(n, Vector(m_bsphere.center)) - dot(n, Vector(ray.o)))/dot(n, ray.d);
		Point closest = ray(t);
		m_dragEnd = m_bsphere.center + normalize(closest - 
				m_bsphere.center) * m_bsphere.radius;
	}
}

void Gizmo::draw(Renderer *renderer, const Camera *camera) {
	if (m_bsphere.isEmpty())
		return;

	/* Compute the tangent circle */
	Vector camToSphere(m_bsphere.center - camera->getPosition());
	const Float length = camToSphere.length(), radius = m_bsphere.radius;
	camToSphere /= length;

	Float tcRadius = std::sqrt(length*length - radius*radius)*radius/length;
	Float tcDist = std::sqrt(radius*radius - tcRadius*tcRadius);

	renderer->setDepthTest(false);
	renderer->drawCircle(m_bsphere.center - camToSphere * tcDist, 
			camToSphere, tcRadius);
	renderer->setDepthTest(true);

	if (m_drag && m_dragStart != m_dragEnd) {
		Spectrum color1, color2;
		color1.fromLinearRGB(0.7f, 0.7f, 1.0f);
		color2.fromLinearRGB(0.3f, 0.3f, 0.3f);
		Transform trafo = getTransform().inverse();
		Point dragEnd = trafo(trafo(m_dragEnd));

		renderer->setColor(color1);
		renderer->drawArc(m_bsphere.center, m_dragStart, dragEnd, true);
		renderer->setColor(color2);
		renderer->drawArc(m_bsphere.center, m_dragStart, dragEnd, false);

		Vector rotAxis = normalize(cross(m_dragStart-m_bsphere.center,
				dragEnd-m_bsphere.center));
		Vector circle2Axis = cross(normalize(m_dragStart-m_bsphere.center), rotAxis);
		renderer->drawCircle(m_bsphere.center, Normal(circle2Axis), radius);
		renderer->setColor(Spectrum(1.0f));
	}
}

MTS_IMPLEMENT_CLASS(Gizmo, false, Object)
MTS_NAMESPACE_END
