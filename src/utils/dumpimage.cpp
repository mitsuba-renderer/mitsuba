#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <iomanip>

using namespace mitsuba;

Bitmap *downsample(const Bitmap *in) {
	Bitmap *result = new Bitmap(in->getWidth()/2, in->getHeight()/2, 128);
	const float *data = in->getFloatData();
	float *out = result->getFloatData();
	result->clear();
	for (int y=0; y<in->getWidth()/2; ++y) {
		for (int x=0; x<in->getWidth()/2; ++x) {
			Float value1 = data[(2*x    +  2*y     * in->getWidth()) * 4];
			Float value2 = data[(2*x    +  (2*y+1) * in->getWidth()) * 4];
			Float value3 = data[(2*x+1  +  2*y     * in->getWidth()) * 4];
			Float value4 = data[(2*x+1  +  (2*y+1) * in->getWidth()) * 4];
			Float avg = (value1+value2+value3+value4)/4;
			out[(x  +  y * result->getWidth()) * 4 + 0] = avg;
			out[(x  +  y * result->getWidth()) * 4 + 1] = avg;
			out[(x  +  y * result->getWidth()) * 4 + 2] = avg;
			out[(x  +  y * result->getWidth()) * 4 + 3] = 1;
		}
	}
	return result;
}

void dumpImage(const std::string &s1) {
	ref<FileStream> stream = new FileStream(s1, FileStream::EReadOnly);
	ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, stream);
//	bitmap = downsample(bitmap);
//	bitmap = downsample(bitmap);
	stream = new FileStream("downsampled.exr", FileStream::ETruncReadWrite);
	bitmap->save(Bitmap::EEXR, stream);

	float *data = bitmap->getFloatData();

	cout << "A={" << endl;
	cout << std::fixed << endl;
	cout << std::setprecision(12)<< endl;
	for (int y=0; y<bitmap->getWidth(); ++y) {
		cout << "\t{";
		for (int x=0; x<bitmap->getWidth(); ++x) {
			cout << data[(x + y * bitmap->getWidth()) * 4];
			if (x+1 < bitmap->getWidth())
				cout << ", ";

		}
		cout << "}";
		if (y+1 < bitmap->getHeight())
			cout << ",";
		cout << endl;
	}
	cout << "};";
}

int main(int argc, char **argv) {
	Class::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	try {
		if (argc < 2) {
			cout << "dumpimage <image.exr>" << endl;
		} else {
			dumpImage(argv[1]);
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
