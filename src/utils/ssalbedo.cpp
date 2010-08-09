#include <mitsuba/core/aabb.h>
#include <mitsuba/core/fstream.h>

using namespace mitsuba;
	
Float solveSSAlbedo(Float diffAlbedo) {
	if (diffAlbedo == 1 || diffAlbedo == 0)
		return diffAlbedo;

	if (diffAlbedo < 0 || diffAlbedo > 1) {
		cout << "Overflow: "<< diffAlbedo << "!" << endl;
		diffAlbedo = std::max((Float) 0, std::min((Float) 1, diffAlbedo));
	}

	Float l = 0, r = 1;
	while (true) {
		Float m = (l+r)/2, 
				fm = m/2*(1+std::exp(-4.0/3.0*std::sqrt(3*(1-m))))*std::exp(-std::sqrt(3*(1-m)));
		if (fm < diffAlbedo)
			l = m;
		else
			r = m;
		if (std::abs(fm-diffAlbedo) < 1e-3)
			return m;
	}
}

void computeSSAlbedo(const std::string &filename, const std::string target) {
	ref<FileStream> stream = new FileStream(filename, FileStream::EReadOnly);
	stream->setByteOrder(Stream::ELittleEndian);

	int xres = stream->readInt(), yres=stream->readInt(), zres=stream->readInt();
	Vector3i res = Vector3i(xres, yres, zres);
	int channels = stream->readInt();

	size_t nEntries = res.x*res.y*res.z*channels;
	Float xmin = stream->readSingle(), ymin = stream->readSingle(), zmin = stream->readSingle();
	Float xmax = stream->readSingle(), ymax = stream->readSingle(), zmax = stream->readSingle();
	AABB aabb(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));

	SLog(EInfo, "Loading \"%s\": %ix%ix%i (%i channels), %i KiB, %s", filename.c_str(), 
		res.x, res.y, res.z, channels, nEntries*sizeof(float)/1024,
		aabb.toString().c_str());
	float *data = new float[nEntries];
	stream->read(data, nEntries*sizeof(float));
	stream->close();

	SLog(EInfo, "Computing single scattering albedo ..");
	for (size_t i=0; i<nEntries; ++i)
		data[i] = solveSSAlbedo(data[i]);

	SLog(EInfo, "Saving \"%s\": %ix%ix%i (%i channels), %i KiB, %s", target.c_str(), 
		res.x, res.y, res.z, channels, nEntries*sizeof(float)/1024,
		aabb.toString().c_str());

	ref<FileStream> targetStream = new FileStream(target, FileStream::ETruncReadWrite);
	targetStream->setByteOrder(Stream::ELittleEndian);
	res.serialize(targetStream);
	targetStream->writeInt(channels);
	targetStream->writeSingle(xmin); targetStream->writeSingle(ymin); targetStream->writeSingle(zmin);
	targetStream->writeSingle(xmax); targetStream->writeSingle(ymax); targetStream->writeSingle(zmax);
	targetStream->write(data, nEntries*sizeof(float));
	targetStream->close();
}

int main(int argc, char **argv) {
	Class::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();

	try {
		if (argc != 3) {
			cout << "Converts a volume of diffuse color values to " << endl;
			cout << "a single scattering albedo volume by numerically" << endl;
			cout << "inverting (2.4) in the Jensen et al. BSSRDF paper" << endl;
			cout << "Syntax: ssalbedo <source.vol> <target.vol>" << endl;
		} else {
			computeSSAlbedo(argv[1], argv[2]);
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
