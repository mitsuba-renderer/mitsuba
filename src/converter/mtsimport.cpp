/**
 * Mitsuba COLLADA 1.4 and Wavefront OBJ -> XML converter
 * 
 * Takes a DAE or OBJ file and turns it into a scene description and separate mesh files
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
 *
 * Since Mitsuba does not support per-vertex colors and prefers textures, any vertex colors 
 * part of the input file are not converted and should instead be baked to textures beforehand 
 * (e.g. using Lighting/shading -> Batch bake in Maya).
 */

#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/dom/DOMException.hpp>
#include "converter.h"
#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/core/fresolver.h>
#if defined(WIN32)
#include "../mitsuba/getopt.h"
#endif

XERCES_CPP_NAMESPACE_USE

class ConsoleGeometryConverter : public GeometryConverter {
public:
	inline ConsoleGeometryConverter() {
	}

	std::string locateResource(const std::string &resource) {
		return "";
	}
};

void help() {
	cout << "COLLADA 1.4 & WaveFront OBJ Importer, Version " MTS_VERSION ", Copyright (c) " MTS_YEAR " Wenzel Jakob" << endl
		<< "Syntax: mtsimport [options] <DAE or OBJ scene> <XML destination file> [Adjustment file]" << endl
		<< "Options/Arguments:" << endl
		<<  "   -h          Display this help text" << endl << endl
		<<  "   -p <num>    Use the specified number of samples per pixel." << endl << endl
		<<  "   -s          Assume that colors are in sRGB space." << endl << endl
		<<  "   -m          Map the larger image side to the full field of view" << endl << endl
		<<  "   -r <w>x<h>  Override the image resolution to e.g. 1920×1080" << endl << endl
		<<  "   -f <fov>    Override the field of view to the given value specified in degrees." << endl << endl
		<< "Please see the documentation for more information." << endl;
}

int colladaMain(int argc, char **argv) {
	bool srgb = false, mapSmallerSide = true;
	char optchar, *end_ptr = NULL;
	int xres = -1, yres = -1;
	int samplesPerPixel = 8;
	Float fov = -1;

	optind = 1;

	while ((optchar = getopt(argc, argv, "shmr:p:f:")) != -1) {
		switch (optchar) {
			case 's':
				srgb = true;
				break;
			case 'm':
				mapSmallerSide = false;
				break;
			case 'p':
				samplesPerPixel = strtol(optarg, &end_ptr, 10);
				if (*end_ptr != '\0')
					SLog(EError, "Invalid number of samples per pixel!");
				break;
			case 'f':
				fov = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0')
					SLog(EError, "Invalid field of view value!");
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

	ConsoleGeometryConverter converter;
	converter.setSRGB(srgb);
	converter.setResolution(xres, yres);
	converter.setMapSmallerSide(mapSmallerSide);
	converter.setSamplesPerPixel(samplesPerPixel);
	converter.setFov(fov);
	converter.convert(argv[optind], "", argv[optind+1], argc > optind+2 ? argv[optind+2] : "");

	return 0;
}

int ubi_main(int argc, char **argv) {
	int retval;
	
	/* Initialize Xerces-C */
	try {
		XMLPlatformUtils::Initialize();
	} catch(const XMLException &toCatch) {
		fprintf(stderr, "Error during Xerces initialization: %s",
			XMLString::transcode(toCatch.getMessage()));
		return -1;
	}

	/* Initialize the core framework */
	Class::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();

	Thread::getThread()->getLogger()->setLogLevel(EInfo);

	FileResolver *resolver = FileResolver::getInstance();
#if defined(WIN32)
	char lpFilename[1024];
	if (GetModuleFileNameA(NULL,
		lpFilename, sizeof(lpFilename))) {
		resolver->addPathFromFile(lpFilename);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}
#elif defined(__LINUX__)
	char exePath[PATH_MAX];
	if (getcwd(exePath, PATH_MAX)) {
		resolver->addPathFromFile(exePath);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}
#else
	MTS_AUTORELEASE_BEGIN()
	resolver->addPath(__ubi_bundlepath());
	MTS_AUTORELEASE_END() 
#endif

#if !defined(WIN32)
	/* Correct number parsing on some locales (e.g. ru_RU) */
	setlocale(LC_NUMERIC, "C");
#endif

	try {
		/* An OpenGL context may be required for the GLU tesselator */
		ref<Session> session = Session::create();
		ref<Device> device = Device::create(session);
		ref<Renderer> renderer = Renderer::create(session);

		session->init();
		device->init();
		renderer->init(device);

		device->makeCurrent(renderer);

		retval = colladaMain(argc, argv);

		if (retval != -1)
			cout << "Finished conversion" << endl;

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
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
		retval = -1;
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type!" << endl;
		retval = -1;
	}
	
	XMLPlatformUtils::Terminate();

	/* Shutdown the core framework */
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	Class::staticShutdown();
	
	return retval;
}

#if !defined(__OSX__)
int main(int argc, char **argv) {
	return ubi_main(argc, argv);
}
#endif


