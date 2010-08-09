#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

using namespace mitsuba;

void addImages(Float weight1, const std::string &s1, Float weight2, const std::string &s2, const std::string &s3) {
	ref<FileStream> aFile   = new FileStream(s1, FileStream::EReadOnly);
	ref<FileStream> bFile   = new FileStream(s2, FileStream::EReadOnly);
	ref<FileStream> outFile = new FileStream(s3, FileStream::ETruncReadWrite);

	ref<Bitmap> aBitmap = new Bitmap(Bitmap::EEXR, aFile);
	ref<Bitmap> bBitmap = new Bitmap(Bitmap::EEXR, bFile);
	ref<Bitmap> outBitmap = new Bitmap(aBitmap->getWidth(), aBitmap->getHeight(), 128);

	float *aData = aBitmap->getFloatData();
	float *bData = bBitmap->getFloatData();
	float *outData = outBitmap->getFloatData();
	int width = aBitmap->getWidth();

	for (int y=0; y<aBitmap->getHeight(); ++y) {
		for (int x=0; x<aBitmap->getWidth(); ++x) {
			Float ra = aData[(x + y * width) * 4] * weight1;
			Float ga = aData[(x + y * width) * 4 + 1] * weight1;
			Float ba = aData[(x + y * width) * 4 + 2] * weight1;
			Float rb = bData[(x + y * width) * 4] * weight2;
			Float gb = bData[(x + y * width) * 4 + 1] * weight2;
			Float bb = bData[(x + y * width) * 4 + 2] * weight2;
			outData[(x+y * width) * 4 + 0] = std::max((Float) 0, ra + rb);
			outData[(x+y * width) * 4 + 1] = std::max((Float) 0, ga + gb);
			outData[(x+y * width) * 4 + 2] = std::max((Float) 0, ba + bb);
			outData[(x+y * width) * 4 + 3] = 1;
		}
	}

	outBitmap->save(Bitmap::EEXR, outFile);
}

int main(int argc, char **argv) {
	Class::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	char *end_ptr = NULL;
	try {
		if (argc != 6) {
			cout << "Add the weighted pixel values of two EXR images to produce a new one" << endl;
			cout << "Syntax: addimages <weight 1> <image 1.exr> <weight 2> <image 2.exr> <target.exr>" << endl;
		} else {
			Float weight1 = (Float) strtod(argv[1], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");
			Float weight2 = (Float) strtod(argv[3], &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse floating point value");
			addImages(weight1, argv[2], weight2, argv[4], argv[5]);
		}
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
		exit(-1);
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type! " << std::endl;
		exit(-1);
	}
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Class::staticShutdown();
	return 0;
}
