#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

using namespace mitsuba;

void joinRGB(const std::string &s1, const std::string &s2, const std::string &s3, const std::string &s4) {
	ref<FileStream> rFile   = new FileStream(s1, FileStream::EReadOnly);
	ref<FileStream> gFile   = new FileStream(s2, FileStream::EReadOnly);
	ref<FileStream> bFile   = new FileStream(s3, FileStream::EReadOnly);
	ref<FileStream> outFile = new FileStream(s4, FileStream::ETruncReadWrite);

	ref<Bitmap> rBitmap = new Bitmap(Bitmap::EEXR, rFile);
	ref<Bitmap> gBitmap = new Bitmap(Bitmap::EEXR, gFile);
	ref<Bitmap> bBitmap = new Bitmap(Bitmap::EEXR, bFile);
	ref<Bitmap> outBitmap = new Bitmap(rBitmap->getWidth(), rBitmap->getHeight(), 128);

	float *rData = rBitmap->getFloatData();
	float *gData = gBitmap->getFloatData();
	float *bData = bBitmap->getFloatData();
	float *outData = outBitmap->getFloatData();
	int width = rBitmap->getWidth();

	for (int y=0; y<rBitmap->getWidth(); ++y) {
		for (int x=0; x<rBitmap->getWidth(); ++x) {
			float r = rData[(x + y * width) * 4];
			float g = gData[(x + y * width) * 4];
			float b = bData[(x + y * width) * 4];
			outData[(x+y * width) * 4 + 0] = r;
			outData[(x+y * width) * 4 + 1] = g;
			outData[(x+y * width) * 4 + 2] = b;
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
	try {
		if (argc < 5) {
			cout << "Join three monochromatic EXRs into a colored image" << endl;
			cout << "joinrgb <red.exr> <green.exr> <blue.exr> <combined.exr>" << endl;
		} else {
			joinRGB(argv[1], argv[2], argv[3], argv[4]);
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
