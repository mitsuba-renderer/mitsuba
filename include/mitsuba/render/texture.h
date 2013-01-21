/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
#if !defined(__MITSUBA_RENDER_TEXTURE_H_)
#define __MITSUBA_RENDER_TEXTURE_H_

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
	/**
	 * \brief Return the texture value at \c its
	 * \param filter
	 *    Specifies whether a filtered texture lookup is desired. Note
	 *    that this does not mean that filtering will actually be used.
	 */
	virtual Spectrum eval(const Intersection &its, bool filter = true) const;

	/// Return the component-wise average value of the texture over its domain
	virtual Spectrum getAverage() const;

	/// Return the component-wise minimum of the texture over its domain
	virtual Spectrum getMinimum() const;

	/// Return the component-wise maximum of the texture over its domain
	virtual Spectrum getMaximum() const;

	/// Return the resolution in pixels, if applicable
	virtual Vector3i getResolution() const;

	/// Return whether the texture takes on a single constant value
	virtual bool isConstant() const;

	/**
	 * \brief Does this texture perform any pre-filtering when
	 * ray differentials are available?
	 */
	virtual bool usesRayDifferentials() const;

	/**
	 * \brief Some textures are only proxies for an actual
	 * implementation. This function returns the actual
	 * texture implementation to be used.
	 *
	 * The default implementation returns <tt>this</tt>.
	 */
	virtual ref<Texture> expand();

	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return the underlying bitmap representation (if any)
	virtual ref<Bitmap> getBitmap() const;

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
	/**
	 * \brief Return the texture value at \c its
	 * \param filter
	 *    Specifies whether a filtered texture lookup is desired. Note
	 *    that this does not mean that filtering will actually be used.
	 */
	Spectrum eval(const Intersection &its, bool filter = true) const;

	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Unfiltered texture lookup -- Texture2D subclasses must provide this function
	virtual Spectrum eval(const Point2 &uv) const = 0;

	/// Filtered texture lookup -- Texture2D subclasses must provide this function
	virtual Spectrum eval(const Point2 &uv, const Vector2 &d0,
			const Vector2 &d1) const = 0;

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

#endif /* __MITSUBA_RENDER_TEXTURE_H_ */
