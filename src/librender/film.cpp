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

#include <mitsuba/render/film.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

Film::Film(const Properties &props)
 : ConfigurableObject(props) {
    bool isMFilm = boost::to_lower_copy(props.getPluginName()) == "mfilm";

    /* Horizontal and vertical film resolution in pixels */
    m_size = Vector2i(
        props.getInteger("width", isMFilm ? 1 : 768),
        props.getInteger("height", isMFilm ? 1 : 576)
    );
    /* Crop window specified in pixels - by default, this
       matches the full sensor area */
    m_cropOffset = Point2i(
        props.getInteger("cropOffsetX", 0),
        props.getInteger("cropOffsetY", 0)
    );
    m_cropSize = Vector2i(
        props.getInteger("cropWidth", m_size.x),
        props.getInteger("cropHeight", m_size.y)
    );
    if (m_cropOffset.x < 0 || m_cropOffset.y < 0 ||
        m_cropSize.x <= 0 || m_cropSize.y <= 0 ||
        m_cropOffset.x + m_cropSize.x > m_size.x ||
        m_cropOffset.y + m_cropSize.y > m_size.y )
        Log(EError, "Invalid crop window specification!");

    /* If set to true, regions slightly outside of the film
       plane will also be sampled, which improves the image
       quality at the edges especially with large reconstruction
       filters. */
    m_highQualityEdges = props.getBoolean("highQualityEdges", false);
}

Film::Film(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
    m_size = Vector2i(stream);
    m_cropOffset = Point2i(stream);
    m_cropSize = Vector2i(stream);
    m_highQualityEdges = stream->readBool();
    m_filter = static_cast<ReconstructionFilter *>(manager->getInstance(stream));
}

Film::~Film() { }

void Film::serialize(Stream *stream, InstanceManager *manager) const {
    ConfigurableObject::serialize(stream, manager);
    m_size.serialize(stream);
    m_cropOffset.serialize(stream);
    m_cropSize.serialize(stream);
    stream->writeBool(m_highQualityEdges);
    manager->serialize(stream, m_filter.get());
}

void Film::addChild(const std::string &name, ConfigurableObject *child) {
    const Class *cClass = child->getClass();

    if (cClass->derivesFrom(MTS_CLASS(ReconstructionFilter))) {
        Assert(m_filter == NULL);
        m_filter = static_cast<ReconstructionFilter *>(child);
    } else {
        Log(EError, "Film: Invalid child node! (\"%s\")",
            cClass->getName().c_str());
    }
}

void Film::configure() {
    if (m_filter == NULL) {
        /* No reconstruction filter has been selected. Load a Gaussian filter by default */
        m_filter = static_cast<ReconstructionFilter *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(ReconstructionFilter), Properties("gaussian")));
        m_filter->configure();
    }
}

MTS_IMPLEMENT_CLASS(Film, true, ConfigurableObject)
MTS_NAMESPACE_END
