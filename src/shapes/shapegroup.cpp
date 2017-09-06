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

#include "shapegroup.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{shapegroup}{Shape group for geometry instancing}
 * \order{8}
 * \parameters{
 *     \parameter{\Unnamed}{\Shape}{One or more shapes that should be
 *         made available for geometry instancing}
 * }
 *
 * This plugin implements a container for shapes that should be
 * made available for geometry instancing. Any shapes placed in a
 * \pluginref{shapegroup} will not be visible on their own---instead, the
 * renderer will precompute ray intersection acceleration data structures
 * so that they can efficiently be referenced many times using the
 * \pluginref{instance} plugin. This is useful for rendering things like
 * forests, where only a few distinct types of trees have to be kept
 * in memory. An example is given below:
 *
 * \vspace{5mm}
 * \begin{xml}[caption={An example of geometry instancing}, label=lst:instancing]
 * <!-- Declare a named shape group containing two objects -->
 * <shape type="shapegroup" id="myShapeGroup">
 *     <shape type="ply">
 *         <string name="filename" value="data.ply"/>
 *         <bsdf type="roughconductor"/>
 *     </shape>
 *     <shape type="sphere">
 *         <transform name="toWorld">
 *             <scale value="5"/>
 *             <translate y="20"/>
 *         </transform>
 *         <bsdf type="diffuse"/>
 *     </shape>
 * </shape>
 *
 * <!-- Instantiate the shape group without any kind of transformation -->
 * <shape type="instance">
 *     <ref id="myShapeGroup"/>
 * </shape>
 *
 * <!-- Create instance of the shape group, but rotated, scaled, and translated -->
 * <shape type="instance">
 *     <ref id="myShapeGroup"/>
 *     <transform name="toWorld">
 *         <rotate x="1" angle="45"/>
 *         <scale value="1.5"/>
 *         <translate z="10"/>
 *     </transform>
 * </shape>
 * \end{xml}
 */

ShapeGroup::ShapeGroup(const Properties &props) : Shape(props) {
    m_kdtree = new ShapeKDTree();
}

ShapeGroup::ShapeGroup(Stream *stream, InstanceManager *manager)
    : Shape(stream, manager) {
    m_kdtree = new ShapeKDTree();
    size_t shapeCount = stream->readSize();
    for (size_t i=0; i<shapeCount; ++i)
        m_kdtree->addShape(static_cast<Shape *>(manager->getInstance(stream)));
    configure();
}

void ShapeGroup::serialize(Stream *stream, InstanceManager *manager) const {
    Shape::serialize(stream, manager);
    const std::vector<const Shape *> &shapes = m_kdtree->getShapes();
    stream->writeSize(shapes.size());
    for (size_t i=0; i<shapes.size(); ++i)
        manager->serialize(stream, shapes[i]);
}

void ShapeGroup::configure() {
    /* Don't bother showing debug messages if the number
       of triangles is low. This helps loading scenes exported
       from SketchUp which create hundreds of tiny shape groups */
    if (m_kdtree->getPrimitiveCount() < 100*1024)
        m_kdtree->setLogLevel(ETrace);
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
    if (cClass->derivesFrom(MTS_CLASS(ShapeGroup)) || cClass->getName() == "Instance") {
        Log(EError, "Nested instancing is not permitted");
    } else if (cClass->derivesFrom(MTS_CLASS(Shape))) {
        Shape *shape = static_cast<Shape *>(child);
        if (shape->isEmitter())
            Log(EError, "Instancing of emitters is not supported");
        if (shape->isSensor())
            Log(EError, "Instancing of sensors is not supported");
        if (shape->hasSubsurface())
            Log(EError, "Instancing of subsurface scattering models is not supported");
        if (shape->isCompound()) {
            int index = 0;
            do {
                ref<Shape> element = shape->getElement(index++);
                if (element == NULL)
                    break;
                addChild(element);
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

size_t ShapeGroup::getPrimitiveCount() const {
    const std::vector<const Shape *> &shapes = m_kdtree->getShapes();
    size_t result = 0;
    for (size_t i=0; i<shapes.size(); ++i)
        result += shapes[i]->getPrimitiveCount();
    return result;
}

size_t ShapeGroup::getEffectivePrimitiveCount() const {
    return 0;
}

std::string ShapeGroup::toString() const {
    std::ostringstream oss;
        oss << "ShapeGroup[" << endl
            << "  name = \"" << m_name << "\"," << endl
            << "  primCount = " << m_kdtree->getPrimitiveCount() << endl
            << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS_S(ShapeGroup, false, Shape)
MTS_EXPORT_PLUGIN(ShapeGroup, "Grouped geometry for instancing");
MTS_NAMESPACE_END
