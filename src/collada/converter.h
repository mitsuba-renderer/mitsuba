#include <mitsuba/mitsuba.h>

using namespace mitsuba;

class ColladaConverter {
public:
	inline ColladaConverter() {
		m_srgb = false;
		m_mapSmallerSide = true;
		m_xres = m_yres = -1;
	}

	void convert(const std::string &inputFile, 
		const std::string &outputDirectory, 
		const std::string &sceneName,
		const std::string &adjustmentFile);

	virtual std::string locateResource(const std::string &resource) = 0;

	void setSRGB(bool srgb) { m_srgb = srgb; }
	void setMapSmallerSide(bool mapSmallerSide) { m_mapSmallerSide = mapSmallerSide; }
	void setResolution(int xres, int yres) { m_xres = xres; m_yres = yres; }
public:
	bool m_srgb, m_mapSmallerSide;
	int m_xres, m_yres;
};
