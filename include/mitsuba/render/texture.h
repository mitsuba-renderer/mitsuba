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

#if !defined(__TEXTURE_H)
#define __TEXTURE_H

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Base class of all textures. Computes values for an arbitrary surface
 * point. \ref Texture2D is a specialization to UV-based textures.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Texture : public ConfigurableObject, public HWResource {
public:
	/// Return the texture value at \a its
	virtual Spectrum getValue(const Intersection &its) const = 0;

	/// Return the component-wise average value of the texture over its domain
	virtual Spectrum getAverage() const = 0;

	/// Return the component-wise minimum of the texture over its domain
	virtual Spectrum getMinimum() const = 0;

	/// Return the component-wise maximum of the texture over its domain
	virtual Spectrum getMaximum() const = 0;

	/// Return the resolution in pixels, if applicable
	virtual Vector3i getResolution() const;

	/// Return whether the texture takes on a single constant value
	virtual bool isConstant() const = 0;

	/**
	 * \brief Does this texture do pre-filtering when ray 
	 * differentials are available?
	 */
	virtual bool usesRayDifferentials() const = 0;

	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	Texture(const Properties &props);
	Texture(Stream *stream, InstanceManager *manager);

	virtual ~Texture();
};

/**
 * \brief Base class of all 2D textures
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Texture2D : public Texture {
public:
	/// Return the texture value at \a its
	Spectrum getValue(const Intersection &its) const;

	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Texture2D subclass must provide this function
	virtual Spectrum getValue(const Point2 &uv) const = 0;

	/// Texture2D subclass must provide this function
	virtual Spectrum getValue(const Point2 &uv, Float dudx,
			Float dudy, Float dvdx, Float dvdy) const = 0;

	MTS_DECLARE_CLASS()
protected:
	Texture2D(const Properties &props);
	Texture2D(Stream *stream, InstanceManager *manager);

	virtual ~Texture2D();
protected:
	Point2 m_uvOffset;
	Vector2 m_uvScale;
};

MTS_NAMESPACE_END

#endif /* __TEXTURE_H */
