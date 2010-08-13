#include <mitsuba/mitsuba.h>

using namespace mitsuba;

class ColladaConverter {
public:
	inline ColladaConverter() {
		m_srgb = false;
		m_mapSmallerSide = true;
		m_xres = m_yres = -1;
		m_samplesPerPixel = 8;
	}

	void convert(const std::string &inputFile, 
		const std::string &outputDirectory, 
		const std::string &sceneName,
		const std::string &adjustmentFile);

	virtual std::string locateResource(const std::string &resource) = 0;

	inline void setSRGB(bool srgb) { m_srgb = srgb; }
	inline void setMapSmallerSide(bool mapSmallerSide) { m_mapSmallerSide = mapSmallerSide; }
	inline void setResolution(int xres, int yres) { m_xres = xres; m_yres = yres; }
	inline void setSamplesPerPixel(int samplesPerPixel) { m_samplesPerPixel = samplesPerPixel; }
	inline const std::string &getFilename() const { return m_filename; }
public:
	bool m_srgb, m_mapSmallerSide;
	int m_xres, m_yres, m_samplesPerPixel;
	std::string m_filename;
};
