/**
 * Mitsuba COLLADA 1.4 and Wavefront OBJ -> XML converter
 *
 * Takes a DAE, ZAE or OBJ file and turns it into a scene description and separate mesh files
 * using a compact binary format. All associated files are copied into newly created
 * 'textures' and 'meshes' directories
 *
 * Currently supports the following subset of the COLLADA specification:
 * - Arbitrary polygonal meshes
 * - Lambert and Phong materials (allowed to be textured)
 * - Cameras
 * - Spot and Point and Ambient lights
 *
 * When exporting DAE using Maya/FBX, be sure to have it convert all NURBS surfaces into
 * "Software Render Meshes". Triangulation is not required (the code below does this
 * automatically for arbitrary polygonal meshes). The Light and camera export options
 * should be activated, since they are off by default. While modeling the scene, it is
 * advisable to use light sources with an inverse square falloff. Otherwise, the
 * illumination will be competely different when rendering in Mitsuba (the image might
 * be pitch black). Note that most BRDFs in Mitsuba treat surfaces as one-sided, thus they
 * will appear black when seen from the back.
 *
 * The conversion barfs when it gets more than 10MB in one single XML string
 * (error: xmlSAX2Characters: huge text node: out of memory). In this case, split the
 * mesh into smaller pieces or recompile libxml with a higher limit.
 */

// Mitsuba's "Assert" macro conflicts with Xerces' XSerializeEngine::Assert(...).
// This becomes a problem when using a PCH which contains mitsuba/core/logger.h
#if defined(Assert)
# undef Assert
#endif
#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/dom/DOMException.hpp>
#include "converter.h"
#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/version.h>
#if defined(__WINDOWS__)
#include <mitsuba/core/getopt.h>
#endif

class ConsoleGeometryConverter : public GeometryConverter {
public:
	inline ConsoleGeometryConverter() {
	}

	fs::path locateResource(const fs::path &resource) {
		return fs::path();
	}
};

void help() {
	cout << "COLLADA 1.4 & Wavefront OBJ Importer, Copyright (c) " MTS_YEAR " Wenzel Jakob" << endl
		<< "Syntax: mtsimport [options] <DAE/ZAE/OBJ scene> <XML output file> [Adjustment file]" << endl
		<< "Options/Arguments:" << endl
		<<  "   -h          Display this help text" << endl << endl
		<<  "   -a p1;p2;.. Add one or more entries to the resource search path" << endl << endl
		<<  "   -v          Be more verbose" << endl << endl
		<<  "   -s          Assume that colors are in sRGB space." << endl << endl
		<<  "   -m          Map the larger image side to the full field of view" << endl << endl
		<<  "   -z          Import animations" << endl << endl
		<<  "   -y          Don't pack all geometry data into a single file" << endl << endl
		<<  "   -n          Don't import any materials (an adjustments file will be necessary)" << endl << endl
		<<  "   -l <type>   Override the type of film (e.g. 'hdrfilm', 'ldrfilm', ..)" << endl << endl
		<<  "   -r <w>x<h>  Override the image resolution to e.g. 1920x1080" << endl << endl
		<< "Please see the documentation for more information." << endl;
}

int importMain(int argc, char **argv) {
	bool srgb = false, mapSmallerSide = true;
	int optchar;
	char *end_ptr = NULL;
	int xres = -1, yres = -1;
	std::string filmType = "hdrfilm";
	FileResolver *fileResolver = Thread::getThread()->getFileResolver();
	ELogLevel logLevel = EInfo;
	bool packGeometry = true, importMaterials = true,
		 importAnimations = false;

	optind = 1;

	while ((optchar = getopt(argc, argv, "snzvyhmr:a:l:")) != -1) {
		switch (optchar) {
			case 'a': {
					std::vector<std::string> paths = tokenize(optarg, ";");
					for (int i=(int) paths.size()-1; i>=0; --i)
						fileResolver->prependPath(paths[i]);
				}
				break;
			case 's':
				srgb = true;
				break;
			case 'm':
				mapSmallerSide = false;
				break;
			case 'n':
				importMaterials = false;
				break;
			case 'z':
				importAnimations = true;
				break;
			case 'v':
				logLevel = EDebug;
				break;
			case 'l':
				filmType = optarg;
				break;
			case 'y':
				packGeometry = false;
				break;
			case 'r': {
					std::vector<std::string> tokens = tokenize(optarg, "x");
					if (tokens.size() != 2)
						SLog(EError, "Invalid resolution argument supplied!");
					xres = strtol(tokens[0].c_str(), &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Invalid resolution argument supplied!");
					yres = strtol(tokens[1].c_str(), &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Invalid resolution argument supplied!");
				}
				break;
			case 'h':
			default:
				help();
				return -1;
	}
	};

	if (argc-optind < 2) {
		help();
		return -1;
	}

	ref<Logger> log = Thread::getThread()->getLogger();
	log->setLogLevel(logLevel);

	ConsoleGeometryConverter converter;
	converter.setSRGB(srgb);
	converter.setResolution(xres, yres);
	converter.setImportMaterials(importMaterials);
	converter.setImportAnimations(importAnimations);
	converter.setMapSmallerSide(mapSmallerSide);
	converter.setPackGeometry(packGeometry);
	converter.setFilmType(filmType);

	const Logger *logger = Thread::getThread()->getLogger();
	size_t initialWarningCount = logger->getWarningCount();
	converter.convert(argv[optind], "", argv[optind+1], argc > optind+2 ? argv[optind+2] : "");
	size_t warningCount = logger->getWarningCount() - initialWarningCount;

	if (warningCount > 0)
		SLog(EInfo, "Encountered " SIZE_T_FMT " warnings -- please check the "
			"messages above for details.", warningCount);

	return 0;
}

int mts_main(int argc, char **argv) {
	int retval;

	/* Initialize Xerces-C */
	XERCES_CPP_NAMESPACE_USE
	try {
		XMLPlatformUtils::Initialize();
	} catch(const XMLException &toCatch) {
		fprintf(stderr, "Error during Xerces initialization: %s",
			XMLString::transcode(toCatch.getMessage()));
		return -1;
	}

	/* Initialize the core framework */
	Class::staticInitialization();
	Object::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Bitmap::staticInitialization();

	Thread::getThread()->getLogger()->setLogLevel(EInfo);

#if !defined(__WINDOWS__)
	/* Correct number parsing on some locales (e.g. ru_RU) */
	setlocale(LC_NUMERIC, "C");
#endif

	try {
		/* An OpenGL context may be required for the GLU tesselator */
		ref<Session> session = Session::create();
		ref<Device> device = Device::create(session);
		ref<Renderer> renderer = Renderer::create(session);
		renderer->setLogLevel(ETrace);
		renderer->setWarnLogLevel(ETrace);

		session->init();
		device->init();
		renderer->init(device);

		device->makeCurrent(renderer);
		ref<Timer> timer = new Timer();

		retval = importMain(argc, argv);

		if (retval != -1)
			cout << "Finished conversion (took " << timer->getMilliseconds() << " ms)" << endl;

		renderer->shutdown();
		device->shutdown();
		session->shutdown();
	} catch(const XMLException &toCatch) {
		cout << "Caught a Xerces exception: " <<
			XMLString::transcode(toCatch.getMessage()) << endl;
		retval = -1;
	} catch(const DOMException &toCatch) {
		cout << "Caught a Xerces exception: " <<
			XMLString::transcode(toCatch.getMessage()) << endl;
		retval = -1;
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exception: " << e.what() << endl;
		retval = -1;
	} catch (...) {
		std::cerr << "Caught a critical exception of unknown type!" << endl;
		retval = -1;
	}

	XMLPlatformUtils::Terminate();

	/* Shutdown the core framework */
	Bitmap::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	PluginManager::staticShutdown();
	Object::staticShutdown();
	Class::staticShutdown();

	return retval;
}

#if !defined(__OSX__) && !defined(__WINDOWS__)
int main(int argc, char **argv) {
	return mts_main(argc, argv);
}
#endif


