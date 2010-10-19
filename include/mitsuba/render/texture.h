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

#if !defined(__TEXTURE_H)
#define __TEXTURE_H

#include <mitsuba/core/cobject.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_RENDER Texture : public ConfigurableObject, public HWResource {
public:
	virtual Spectrum getValue(const Intersection &its) const = 0;
	virtual Spectrum getAverage() const = 0;
	virtual Spectrum getMaximum() const = 0; /* Component-wise maximum */
	virtual bool usesRayDifferentials() const = 0;
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	Texture(const Properties &props);
	Texture(Stream *stream, InstanceManager *manager);

	virtual ~Texture();
};

MTS_NAMESPACE_END

#endif /* __TEXTURE_H */
