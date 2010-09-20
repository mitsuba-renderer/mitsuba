/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>

MTS_NAMESPACE_BEGIN

Shape::Shape(const Properties &props) 
 : ConfigurableObject(props), m_luminaire(NULL) {
	m_objectToWorld = props.getTransform("toWorld", Transform());
	m_worldToObject = m_objectToWorld.inverse();
	m_surfaceArea = 0.0f;
}

Shape::Shape(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_worldToObject = Transform(stream);
	m_aabb = AABB(stream);
	m_bsphere = BSphere(stream);
	m_surfaceArea = stream->readFloat();
	m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));
	m_subsurface = static_cast<Subsurface *>(manager->getInstance(stream));
	m_luminaire = static_cast<Luminaire *>(manager->getInstance(stream));
	m_invSurfaceArea = 1.0f / m_surfaceArea;
	m_objectToWorld = m_worldToObject.inverse();
}

Shape::~Shape() {
}


void Shape::configure() {
	/* Ensure that there is at least some default BSDF */
	if (m_bsdf == NULL) {
		m_bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(BSDF::m_theClass, Properties("lambertian")));
	}
}
	
bool Shape::isCompound() const {
	return false;
}

bool Shape::isClippable() const {
	return false;
}

AABB Shape::getClippedAABB(const AABB &aabb) const {
	AABB result(m_aabb);
	result.clip(aabb);
	return result;
}

Shape *Shape::getElement(int i) {
	return NULL;
}

Float Shape::pdfArea(const ShapeSamplingRecord &sRec) const {
	return m_invSurfaceArea;
}

Float Shape::sampleSolidAngle(ShapeSamplingRecord &sRec, const Point &from, const Point2 &sample) const {
	/* Turns the area sampling routine into one that samples wrt. solid angles */
	Float pdfArea = sampleArea(sRec, sample);
	Vector lumToPoint = from - sRec.p;
	Float distSquared = lumToPoint.lengthSquared(), dp = dot(lumToPoint, sRec.n);
	if (dp > 0) {
		return pdfArea * distSquared * std::sqrt(distSquared) / dp;
	} else {
		return 0.0f;
	}
}

Float Shape::pdfSolidAngle(const ShapeSamplingRecord &sRec, const Point &from) const {
	/* Turns the area sampling routine into one that samples wrt. solid angles */
	Vector lumToPoint = from - sRec.p;
	Float distSquared = lumToPoint.lengthSquared();
	Float invDP = std::max((Float) 0, std::sqrt(distSquared) / dot(lumToPoint, sRec.n));
	return pdfArea(sRec) * distSquared * invDP;
}

void Shape::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();
	if (cClass->derivesFrom(BSDF::m_theClass)) {
		m_bsdf = static_cast<BSDF *>(child);
	} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
		Assert(m_luminaire == NULL);
		m_luminaire = static_cast<Luminaire *>(child);
	} else if (cClass->derivesFrom(Subsurface::m_theClass)) {
		Assert(m_subsurface == NULL);
		m_subsurface = static_cast<Subsurface *>(child);
	} else {
		Log(EError, "Shape: Invalid child node!");
	}
}

/// Ray intersection test
bool Shape::rayIntersect(const Ray &ray, Intersection &its) const {
	Log(EError, "Not implemented!");
	return false;
}

bool Shape::rayIntersect(const Ray &ray, Float start, Float end, Float &t) const {
	Log(EError, "Not implemented!");
	return false;
}

Float Shape::sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
	Log(EError, "Not implemented!");
	return 0;
}

#if defined(MTS_SSE)
__m128 Shape::rayIntersectPacket(const RayPacket4 &packet, const
       __m128 mint, __m128 maxt, __m128 inactive, Intersection4 &its) const {
	Log(EError, "Not implemented!");
	return SSEConstants::zero.ps;
}
#endif

void Shape::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);

	m_worldToObject.serialize(stream);
	m_aabb.serialize(stream);
	m_bsphere.serialize(stream);
	stream->writeFloat(m_surfaceArea);
	manager->serialize(stream, m_bsdf.get());
	manager->serialize(stream, m_subsurface.get());
	manager->serialize(stream, m_luminaire.get());
}

std::string ShapeSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "ShapeSamplingRecord[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  n = " << n.toString() << std::endl
		<< "]";
	return oss.str();
}

std::string Intersection::toString() const {
	if (!isValid())
		return "Intersection[invalid]";
	std::ostringstream oss;
	oss << "Intersection[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  t = " << t << "," << std::endl
		<< "  geoFrame = " << indent(geoFrame.toString()) << "," << std::endl
		<< "  shFrame = " << indent(shFrame.toString()) << "," << std::endl
		<< "  uv = " << uv.toString() << "," << std::endl
		<< "  dpdu = " << dpdu.toString() << "," << std::endl
		<< "  dpdv = " << dpdv.toString() << "," << std::endl
		<< "  shape = " << indent(((Object *)shape)->toString()) << std::endl
		<< "]";
	return oss.str();
}

void Intersection::computePartials(const RayDifferential &ray) {
	/* Compute the texture coordinates partials wrt. 
	   changes in the screen-space position. Based on PBRT */
	if (hasUVPartials)
		return;
	hasUVPartials = true;

	if (!ray.hasDifferentials) {
		dudx = dvdx = dudy = dvdy = 0.0f;
		return;
	}

	/* Offset of the plane passing through the surface */
	const Float d = -dot(geoFrame.n, Vector(p));

	const Float txRecip = dot(geoFrame.n, ray.rx.d),
				tyRecip = dot(geoFrame.n, ray.ry.d);

	if (EXPECT_NOT_TAKEN(txRecip == 0 || tyRecip == 0)) {
		dudx = dvdx = dudy = dvdy = 0.0f;
		return;
	}

	/* Ray distances traveled */
	const Float tx = -(dot(geoFrame.n, Vector(ray.rx.o)) + d) / 
		txRecip;
	const Float ty = -(dot(geoFrame.n, Vector(ray.ry.o)) + d) / 
		tyRecip;

	/* Auxilary intersection point of the adjacent rays */
	Point px = ray.rx(tx), py = ray.ry(ty);

	/* Calculate the U and V partials by solving two out
	   of a set of 3 equations in an overconstrained system */
	Float A[2][2], Bx[2], By[2], x[2];
	int axes[2];

	Float absX = std::abs(geoFrame.n.x),
		  absY = std::abs(geoFrame.n.y),
		  absZ = std::abs(geoFrame.n.z);

	if (absX > absY && absX > absZ) {
		axes[0] = 1; axes[1] = 2;
	} else if (absY > absZ) {
		axes[0] = 0; axes[1] = 2;
	} else {
		axes[0] = 0; axes[1] = 1;
	}

	A[0][0] = dpdu[axes[0]];
	A[0][1] = dpdv[axes[0]];
	A[1][0] = dpdu[axes[1]];
	A[1][1] = dpdv[axes[1]];

	Bx[0] = px[axes[0]] - p[axes[0]];
	Bx[1] = px[axes[1]] - p[axes[1]];
	By[0] = py[axes[0]] - p[axes[0]];
	By[1] = py[axes[1]] - p[axes[1]];

	if (solveLinearSystem2x2(A, Bx, x)) {
		dudx = x[0]; dvdx = x[1];
	} else {
		dudx = 1; dvdx = 0;
	}

	if (solveLinearSystem2x2(A, By, x)) {
		dudy = x[0]; dvdy = x[1];
	} else {
		dudy = 0; dudy = 1;
	}
}

MTS_IMPLEMENT_CLASS(Shape, true, ConfigurableObject)
MTS_NAMESPACE_END
