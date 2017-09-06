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

#include "instance.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{instance}{Geometry instance}
 * \order{9}
 * \parameters{
 *     \parameter{\Unnamed}{\ShapeGroup}{A reference to a
 *     shape group that should be instantiated}
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional linear instance-to-world transformation.
 *        \default{none (i.e. instance space $=$ world space)}
 *     }
 * }
 * \renderings{
 *    \rendering{Surface viewed from the top}{shape_instance_fractal_top}
 *    \rendering{Surface viewed from the bottom}{shape_instance_fractal_bot}
 *    \caption{
 *       A visualization of a fractal surface by Irving and Segerman.
 *       (a 2D Gospel curve developed up to level 5 along the third
 *       dimension). This scene makes use of instancing to replicate
 *       similar structures to cheaply render a shape that effectively
 *       consists of several hundred millions of triangles.
 *    }
 * }
 *
 * This plugin implements a geometry instance used to efficiently replicate
 * geometry many times. For details on how to create instances, refer to
 * the \pluginref{shapegroup} plugin.
 * \remarks{
 *   \item Note that it is \emph{not} possible to assign a different
 *    material to each instance --- the material assignment specified within
 *    the shape group is the one that matters.
 *   \item Shape groups cannot be used to replicate shapes with
 *   attached emitters, sensors, or subsurface scattering models.
 * }
 */

Instance::Instance(const Properties &props) : Shape(props) {
    m_transform = props.getAnimatedTransform("toWorld", Transform());
}

Instance::Instance(Stream *stream, InstanceManager *manager)
    : Shape(stream, manager) {
    m_shapeGroup = static_cast<ShapeGroup *>(manager->getInstance(stream));
    m_transform = new AnimatedTransform(stream);
}

void Instance::serialize(Stream *stream, InstanceManager *manager) const {
    Shape::serialize(stream, manager);
    manager->serialize(stream, m_shapeGroup.get());
    m_transform->serialize(stream);
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

    std::set<Float> times;
    m_transform->collectKeyframes(times);

    AABB result;
    for (std::set<Float>::iterator it = times.begin(); it != times.end(); ++it) {
        const Transform &trafo = m_transform->eval(*it);

        for (int i=0; i<8; ++i)
            result.expandBy(trafo(aabb.getCorner(i)));
    }

    return result;
}

void Instance::addChild(const std::string &name, ConfigurableObject *child) {
    const Class *cClass = child->getClass();
    if (cClass->getName() == "ShapeGroup") {
        m_shapeGroup = static_cast<ShapeGroup *>(child);
    } else {
        Shape::addChild(name, child);
    }
}

size_t Instance::getPrimitiveCount() const {
    return 0;
}

size_t Instance::getEffectivePrimitiveCount() const {
    return m_shapeGroup->getPrimitiveCount();
}

bool Instance::rayIntersect(const Ray &_ray, Float mint,
        Float maxt, Float &t, void *temp) const {
    const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
    const Transform &trafo = m_transform->eval(_ray.time);
    Ray ray;
    trafo.inverse()(_ray, ray);
    return kdtree->rayIntersect(ray, mint, maxt, t, temp);
}

bool Instance::rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
    const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
    Ray ray;
    const Transform &trafo = m_transform->eval(_ray.time);
    trafo.inverse()(_ray, ray);
    return kdtree->rayIntersect(ray, mint, maxt);
}

void Instance::adjustTime(Intersection &its, Float time) const {
    Transform trafo = m_transform->eval(its.time).inverse();
    trafo = m_transform->eval(time) * trafo;

    its.dpdu = trafo(its.dpdu);
    its.dpdv = trafo(its.dpdv);
    its.geoFrame = Frame(normalize(trafo(its.geoFrame.n)));
    its.p = trafo(its.p);
    computeShadingFrame(normalize(trafo(its.shFrame.n)), its.dpdu, its.shFrame);
    its.wi = normalize(trafo(its.wi));
    its.instance = this;
    its.time = time;
}

void Instance::fillIntersectionRecord(const Ray &_ray,
    const void *temp, Intersection &its) const {
    const ShapeKDTree *kdtree = m_shapeGroup->getKDTree();
    const Transform &trafo = m_transform->eval(_ray.time);
    Ray ray;
    trafo.inverse()(_ray, ray);
    kdtree->fillIntersectionRecord<false>(ray, temp, its);

    its.shFrame.n = normalize(trafo(its.shFrame.n));
    its.geoFrame = Frame(normalize(trafo(its.geoFrame.n)));
    its.dpdu = trafo(its.dpdu);
    its.dpdv = trafo(its.dpdv);
    its.p = trafo(its.p);
    its.instance = this;
}

void Instance::getNormalDerivative(const Intersection &its,
        Vector &dndu, Vector &dndv, bool shadingFrame) const {
    const Transform &trafo = m_transform->eval(its.time);
    const Transform invTrafo = trafo.inverse();

    /* The following is really super-inefficient, but it's
       needed to be able to deal with general transformations */
    Intersection temp(its);
    temp.p = invTrafo(its.p);
    temp.dpdu = invTrafo(its.dpdu);
    temp.dpdv = invTrafo(its.dpdv);

    /* Determine the length of the transformed normal
       *before* it was re-normalized */
    Normal tn = trafo(normalize(invTrafo(its.shFrame.n)));
    Float invLen = 1 / tn.length();
    tn *= invLen;

    its.shape->getNormalDerivative(temp, dndu, dndv, shadingFrame);

    dndu = trafo(Normal(dndu)) * invLen;
    dndv = trafo(Normal(dndv)) * invLen;

    dndu -= tn * dot(tn, dndu);
    dndv -= tn * dot(tn, dndv);
}

MTS_IMPLEMENT_CLASS_S(Instance, false, Shape)
MTS_EXPORT_PLUGIN(Instance, "Instanced geometry");
MTS_NAMESPACE_END
