#include <mitsuba/mitsuba.h>

using namespace mitsuba;

class ColladaConverter {
public:
	inline ColladaConverter() { }

	void convert(const std::string &inputFile, 
		const std::string &outputDirectory, 
		const std::string &sceneName,
		const std::string &adjustmentFile);

	virtual std::string locateResource(const std::string &resource) = 0;
};
