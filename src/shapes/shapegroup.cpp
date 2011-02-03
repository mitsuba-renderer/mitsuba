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

#include "shapegroup.h"

MTS_NAMESPACE_BEGIN

ShapeGroup::ShapeGroup(const Properties &props) : Shape(props) {
	m_kdtree = new ShapeKDTree();
	m_name = props.getID();
}

ShapeGroup::ShapeGroup(Stream *stream, InstanceManager *manager) 
	: Shape(stream, manager) {
	m_kdtree = new ShapeKDTree();
	size_t shapeCount = stream->readUInt();
	for (size_t i=0; i<shapeCount; ++i)
		m_kdtree->addShape(static_cast<Shape *>(manager->getInstance(stream)));
	configure();
}

void ShapeGroup::serialize(Stream *stream, InstanceManager *manager) const {
	Shape::serialize(stream, manager);
	const std::vector<const Shape *> &shapes = m_kdtree->getShapes();
	stream->writeUInt((uint32_t) shapes.size());
	for (size_t i=0; i<shapes.size(); ++i)
		manager->serialize(stream, shapes[i]);
}

void ShapeGroup::configure() {
	if (!m_kdtree->isBuilt())
		m_kdtree->build();
}

AABB ShapeGroup::getAABB() const {
	return AABB();
}

Float ShapeGroup::getSurfaceArea() const {
	return 0.0f;
}

void ShapeGroup::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();
	if (cClass->derivesFrom(ShapeGroup::m_theClass) || cClass->getName() == "ShapeInstance") {
		Log(EError, "Nested instancing is not supported!");
	} else if (cClass->derivesFrom(Shape::m_theClass)) {
		Shape *shape = static_cast<Shape *>(child);
		if (shape->isLuminaire())
			Log(EError, "Instancing of luminaires is not supported");
		if (shape->hasSubsurface())
			Log(EError, "Instancing of subsurface integrators is not supported");
		if (shape->isCompound()) {
			int index = 0;
			do {
				ref<Shape> element = shape->getElement(index++);
				if (element == NULL)
					break;
				addChild("", element);
			} while (true);
		} else {
			m_kdtree->addShape(shape);
		}
	} else {
		Shape::addChild(name, child);
	}
}

bool ShapeGroup::isCompound() const {
	// this shape reduces to nothing (compound, no children)
	return true;
}

std::string ShapeGroup::getName() const {
	return m_name;
}

std::string ShapeGroup::toString() const {
	std::ostringstream oss;
		oss << "ShapeGroup[" << endl
			<< "  name = \"" << m_name << "\", " << endl
			<< "  primCount = " << m_kdtree->getPrimitiveCount() << endl
			<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(ShapeGroup, false, Shape)
MTS_EXPORT_PLUGIN(ShapeGroup, "Grouped geometry for instancing");
MTS_NAMESPACE_END
