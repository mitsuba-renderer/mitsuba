#include <mitsuba/mitsuba.h>

using namespace mitsuba;

class ColladaConverter {
public:
	inline ColladaConverter() {
		m_srgb = false;
	}

	void convert(const std::string &inputFile, 
		const std::string &outputDirectory, 
		const std::string &sceneName,
		const std::string &adjustmentFile);

	virtual std::string locateResource(const std::string &resource) = 0;

	void setSRGB(bool srgb) { m_srgb = srgb; }
public:
	bool m_srgb;
};
