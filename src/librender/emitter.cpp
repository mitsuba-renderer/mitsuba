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

#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/track.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

AbstractEmitter::AbstractEmitter(const Properties &props)
 : ConfigurableObject(props), m_shape(NULL), m_type(0) {
    m_worldTransform = props.getAnimatedTransform("toWorld", Transform());
}

AbstractEmitter::AbstractEmitter(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
    m_worldTransform = new AnimatedTransform(stream);
    m_medium = static_cast<Medium *>(manager->getInstance(stream));
    m_shape = static_cast<Shape *>(manager->getInstance(stream));
    m_type = stream->readUInt();
 }

AbstractEmitter::~AbstractEmitter() {
}

void AbstractEmitter::serialize(Stream *stream, InstanceManager *manager) const {
    ConfigurableObject::serialize(stream, manager);
    m_worldTransform->serialize(stream);
    manager->serialize(stream, m_medium.get());
    manager->serialize(stream, m_shape);
    stream->writeUInt(m_type);
}

void AbstractEmitter::addChild(const std::string &name, ConfigurableObject *child) {
    if (child->getClass()->derivesFrom(MTS_CLASS(Medium))) {
        Assert(m_medium == NULL);
        m_medium = static_cast<Medium *>(child);
    } else {
        ConfigurableObject::addChild(name, child);
    }
}

ref<Shape> AbstractEmitter::createShape(const Scene *scene) {
    return NULL;
}

Spectrum AbstractEmitter::samplePosition(PositionSamplingRecord &pRec,
        const Point2 &sample, const Point2 *extra) const {
    NotImplementedError("samplePosition");
}

Spectrum AbstractEmitter::sampleDirection(DirectionSamplingRecord &dRec,
        PositionSamplingRecord &pRec, const Point2 &sample,
        const Point2 *extra) const {
    NotImplementedError("sampleDirection");
}

Spectrum AbstractEmitter::sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
    NotImplementedError("sampleDirect");
}

Spectrum AbstractEmitter::evalPosition(const PositionSamplingRecord &pRec) const {
    NotImplementedError("evalPosition");
}

Spectrum AbstractEmitter::evalDirection(const DirectionSamplingRecord &dRec,
        const PositionSamplingRecord &pRec) const {
    NotImplementedError("evalDirection");
}

Float AbstractEmitter::pdfPosition(const PositionSamplingRecord &pRec) const {
    NotImplementedError("pdfPosition");
}

Float AbstractEmitter::pdfDirection(const DirectionSamplingRecord &dRec,
        const PositionSamplingRecord &pRec) const {
    NotImplementedError("pdfDirection");
}

Float AbstractEmitter::pdfDirect(const DirectSamplingRecord &dRec) const {
    NotImplementedError("pdfDirect");
}

Emitter::Emitter(const Properties &props)
 : AbstractEmitter(props) {
    // Importance sampling weight (used by the luminaire sampling code in \ref Scene)
    m_samplingWeight = props.getFloat("samplingWeight", 1.0f);
}

Emitter::Emitter(Stream *stream, InstanceManager *manager)
 : AbstractEmitter(stream, manager) {
    m_samplingWeight = stream->readFloat();
}

void Emitter::serialize(Stream *stream, InstanceManager *manager) const {
    AbstractEmitter::serialize(stream, manager);

    stream->writeFloat(m_samplingWeight);
}

Spectrum Emitter::sampleRay(Ray &ray,
        const Point2 &spatialSample,
        const Point2 &directionalSample,
        Float time) const {
    NotImplementedError("sampleRay");
}

Spectrum Emitter::eval(const Intersection &its, const Vector &d) const {
    NotImplementedError("eval");
}

Spectrum Emitter::evalEnvironment(const RayDifferential &ray) const {
    NotImplementedError("evalEnvironment");
}

bool Emitter::fillDirectSamplingRecord(DirectSamplingRecord &dRec,
        const Ray &ray) const {
    NotImplementedError("fillDirectSamplingRecord");
}

Emitter::~Emitter() { }

Emitter *Emitter::getElement(size_t index) {
    return NULL;
}

bool Emitter::isCompound() const {
    return false;
}

ref<Bitmap> Emitter::getBitmap(const Vector2i &sizeHint) const {
    NotImplementedError("getBitmap");
}

MTS_IMPLEMENT_CLASS(Emitter, false, AbstractEmitter)
MTS_IMPLEMENT_CLASS(AbstractEmitter, true, ConfigurableObject)
MTS_NAMESPACE_END
