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

#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/**
 * Checkerboard texture
 */
class Checkerboard : public Texture {
public:
	Checkerboard(const Properties &props) : Texture(props) {
		m_reflectance = props.getSpectrum("reflectance", Spectrum(1.0f));
	}

	Checkerboard(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {
		m_reflectance = Spectrum(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		m_reflectance.serialize(stream);
	}

	Spectrum getValue(const Intersection &its) const {
		int x = 2*(((int) its.uv.x) % 2) - 1, y = 2*(((int) its.uv.y) % 2) - 1;
	
		if (x*y == 1)
			return m_reflectance;
		else
			return Spectrum(0.0f);
	}
	
	bool usesRayDifferentials() const {
		return false;
	}
	
	Spectrum getAverage() const {
		return m_reflectance * .5f;
	}
	
	Spectrum getMaximum() const {
		return m_reflectance;
	}

	std::string toString() const {
		return "Checkerboard[]";
	}

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_reflectance;
};

MTS_IMPLEMENT_CLASS_S(Checkerboard, false, Texture)
MTS_EXPORT_PLUGIN(Checkerboard, "Checkerboard texture");
MTS_NAMESPACE_END
