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

#pragma once
#if !defined(__MITSUBA_RENDER_VOLUME_H_)
#define __MITSUBA_RENDER_VOLUME_H_

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generalized source of volumetric information
 * \ingroup librender
 */
class MTS_EXPORT_RENDER VolumeDataSource : public ConfigurableObject {
public:
	/// Serialize to a binary data stream.
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return the bounding box
	inline const AABB &getAABB() const {
		return m_aabb;
	}

	/// Are float-valued lookups permitted?
	virtual bool supportsFloatLookups() const;

	/// Look up a floating point value by position
	virtual Float lookupFloat(const Point &p) const;

	/// Are spectrum-valued lookups permitted?
	virtual bool supportsSpectrumLookups() const;

	/// Look up a spectrum value by position
	virtual Spectrum lookupSpectrum(const Point &p) const;

	/// Are vector-valued lookups permitted?
	virtual bool supportsVectorLookups() const;

	/// Look up a vector value by position
	virtual Vector lookupVector(const Point &p) const;

	/**
	 * \brief Return the recommended step size for numerical
	 * integration or inifinity if this is not known/applicable
	 */
	virtual Float getStepSize() const = 0;

	/**
	 * \brief Return the maximum floating point value that
	 * could be returned by \ref lookupFloat.
	 *
	 * This is useful when implementing Woodcock-Tracking.
	 */
	virtual Float getMaximumFloatValue() const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~VolumeDataSource();

	/// Protected constructor
	VolumeDataSource(const Properties &props);

	/// Unserialize from a binary data stream
	VolumeDataSource(Stream *stream, InstanceManager *manager);
protected:
	AABB m_aabb;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_VOLUME_H_ */
