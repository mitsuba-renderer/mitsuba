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

#include <mitsuba/core/fresolver.h>

using namespace mitsuba;

class GeometryConverter {
public:
	inline GeometryConverter() {
		m_srgb = false;
		m_mapSmallerSide = true;
		m_xres = m_yres = -1;
		m_samplesPerPixel = 8;
		m_fov = -1;
		m_filmType = "exrfilm";
	}

	void convert(const fs::path &inputFile, 
		const fs::path &outputDirectory, 
		const fs::path &sceneName,
		const fs::path &adjustmentFile);

	virtual fs::path locateResource(const fs::path &resource) = 0;

	inline void setSRGB(bool srgb) { m_srgb = srgb; }
	inline void setMapSmallerSide(bool mapSmallerSide) { m_mapSmallerSide = mapSmallerSide; }
	inline void setResolution(int xres, int yres) { m_xres = xres; m_yres = yres; }
	inline void setSamplesPerPixel(int samplesPerPixel) { m_samplesPerPixel = samplesPerPixel; }
	inline void setFov(Float fov) { m_fov = fov; }
	inline void setFilmType(const std::string &filmType) { m_filmType = filmType; }
	inline const fs::path &getFilename() const { return m_filename; }
private:
	void convertCollada(const fs::path &inputFile, std::ostream &os,
		const fs::path &textureDirectory,
		const fs::path &meshesDirectory);
	void convertOBJ(const fs::path &inputFile, std::ostream &os,
		const fs::path &textureDirectory,
		const fs::path &meshesDirectory);
public:
	bool m_srgb, m_mapSmallerSide;
	int m_xres, m_yres, m_samplesPerPixel;
	Float m_fov;
	fs::path m_filename;
	std::string m_filmType;
};
