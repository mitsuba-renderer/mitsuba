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

#include <mitsuba/core/properties.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

Subsurface::Subsurface(const Properties &props)
 : NetworkedObject(props), m_active(false) { }

Subsurface::Subsurface(Stream *stream, InstanceManager *manager) :
    NetworkedObject(stream, manager) {
    size_t shapeCount = stream->readSize();
    for (size_t i=0; i<shapeCount; ++i) {
        Shape *shape = static_cast<Shape *>(manager->getInstance(stream));
        m_shapes.push_back(shape);
    }
    m_active = false;
}

Subsurface::~Subsurface() { }

void Subsurface::cancel() { }

void Subsurface::setParent(ConfigurableObject *parent) {
    if (parent->getClass()->derivesFrom(MTS_CLASS(Shape))) {
        Shape *shape = static_cast<Shape *>(parent);
        if (shape->isCompound())
            return;
        if (std::find(m_shapes.begin(), m_shapes.end(), shape) == m_shapes.end()) {
            m_shapes.push_back(shape);
        }
    }
}

void Subsurface::serialize(Stream *stream, InstanceManager *manager) const {
    NetworkedObject::serialize(stream, manager);

    stream->writeSize(m_shapes.size());
    for (unsigned int i=0; i<m_shapes.size(); ++i)
        manager->serialize(stream, m_shapes[i]);
}

MTS_IMPLEMENT_CLASS(Subsurface, true, NetworkedObject)
MTS_NAMESPACE_END
