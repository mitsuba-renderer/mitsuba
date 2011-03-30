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

#include "instance.h"

MTS_NAMESPACE_BEGIN

Instance::Instance(const Properties &props) : Shape(props) {
	m_objectToWorld = props.getTransform("toWorld", Transform());
	m_worldToObject = m_objectToWorld.inverse();
}

Instance::Instance(Stream *stream, InstanceManager *manager) 
	: Shape(stream, manager) {
	m_shapeGroup = static_cast<ShapeGroup *>(manager->getInstance(stream));
	m_objectToWorld = Transform(stream);
	m_worldToObject = m_objectToWorld.inverse();
}

void Instance::serialize(Stream *stream, InstanceManager *manager) const {
	Shape::serialize(stream, manager);
	m_objectToWorld.serialize(stream);
	manager->serialize(stream, m_shapeGroup.get());
}

void Instance::configure() {
	if (!m_shapeGroup)
		Log(EError, "A reference to a 'shapegroup' must be specified!");
}

AABB Instance::getAABB() const {
	const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
	const AABB &aabb = kdtree->getAABB();
	if (!aabb.isValid()) // the geometry group is empty
		return aabb;
	AABB result;
	for (int i=0; i<8; ++i)
		result.expandBy(m_objectToWorld(aabb.getCorner(i)));
	return result;
}

Float Instance::getSurfaceArea() const {
	Log(EError, "Instance::getSurfaceArea(): not supported!");
	return 0.0f;
}

void Instance::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();
	if (cClass->getName() == "ShapeGroup") {
		m_shapeGroup = static_cast<ShapeGroup *>(child);
	} else {
		Shape::addChild(name, child);
	}
}

bool Instance::rayIntersect(const Ray &_ray, Float mint, 
		Float maxt, Float &t, void *temp) const {
	const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
	Ray ray;
	m_worldToObject(_ray, ray);
	return kdtree->rayIntersect(ray, mint, maxt, t, temp);
}

bool Instance::rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
	const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
	Ray ray;
	m_worldToObject(_ray, ray);
	return kdtree->rayIntersect(ray, mint, maxt);
}

void Instance::fillIntersectionRecord(const Ray &ray, 
	const void *temp, Intersection &its) const {
	const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
	kdtree->fillIntersectionRecord<false>(ray, temp, its);
	its.shFrame.n = normalize(m_objectToWorld(its.shFrame.n));
	its.shFrame.s = normalize(m_objectToWorld(its.shFrame.s));
	its.shFrame.t = normalize(m_objectToWorld(its.shFrame.t));
	its.geoFrame = Frame(normalize(m_objectToWorld(its.geoFrame.n)));
	its.wi = its.shFrame.toLocal(-ray.d);
	its.dpdu = m_objectToWorld(its.dpdu);
	its.dpdv = m_objectToWorld(its.dpdv);
}

MTS_IMPLEMENT_CLASS_S(Instance, false, Shape)
MTS_EXPORT_PLUGIN(Instance, "Instanced geometry");
MTS_NAMESPACE_END
