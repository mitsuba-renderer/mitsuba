/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#ifndef WIN32
#include <arpa/inet.h>
#endif

using namespace mitsuba;
	
void createVol(const std::string &format, int start, int end, Float dc, std::string outputFile) {
	ref<FileStream> out = new FileStream(outputFile, FileStream::ETruncReadWrite);
	out->setByteOrder(Stream::ELittleEndian);
	ref<FileStream> stream = new FileStream(formatString(format.c_str(), start), FileStream::EReadOnly);
	ref<Bitmap> bitmap = new Bitmap(Bitmap::EPNG, stream);
	Vector3i res(bitmap->getWidth()/2, bitmap->getHeight()/2, end-start+1);
	SAssert(bitmap->getBitsPerPixel() == 16);
	res.serialize(out);
	out->writeInt(1);
	Float zres = (end-start)/(Float)(bitmap->getWidth()); /// XXX verify this computation
	Float xmin = -.5, ymin = -.5, zmin = -zres, xmax = .5, ymax = .5, zmax = zres;
	out->writeSingle(xmin); out->writeSingle(ymin); out->writeSingle(zmin);
	out->writeSingle(xmax); out->writeSingle(ymax); out->writeSingle(zmax);
	float *newData = new float[res.x * res.y];
	for (int i=start; i<=end; ++i) {
		cout << "File " << i << "/" << end-start+1<< endl;
		ref<FileStream> debugStream  = new FileStream("debug.png", FileStream::ETruncReadWrite);
		ref<Bitmap> debugBitmap = new Bitmap(res.x, res.y, 8);
		stream = new FileStream(formatString(format.c_str(), i), FileStream::EReadOnly);
		bitmap = new Bitmap(Bitmap::EPNG, stream);
		uint16_t *data = (uint16_t *) bitmap->getData();
		int pos = 0;
		uint8_t *debugData = (uint8_t *) debugBitmap->getData();
		for (int y=0; y<bitmap->getHeight(); y += 2) {
			for (int x=0; x<bitmap->getWidth(); x+=2) {
				Float value1 = ntohs(data[(x  ) + (y  )*bitmap->getWidth()]) / 63335.0f;
				Float value2 = ntohs(data[(x+1) + (y  )*bitmap->getWidth()]) / 63335.0f;
				Float value3 = ntohs(data[(x+1) + (y+1)*bitmap->getWidth()]) / 63335.0f;
				Float value4 = ntohs(data[(x  ) + (y+1)*bitmap->getWidth()]) / 63335.0f;
				Float value = value1+value2+value3+value4-4*dc;
				value = std::exp(value*4);
				if (value < .2)
					value = 0;
				debugData[pos] = std::min(255, (int) (value * 255));
				newData[pos] = value;
				++pos;
			}
		}
		debugBitmap->save(Bitmap::EPNG, debugStream);
		out->writeSingleArray(newData, res.x * res.y);
	}
	delete[] newData;
	out->close();
}

int main(int argc, char **argv) {
	Class::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();

	try {
		if (argc != 6) {
			cout << "Converts a bunch of 16-bit PNGs into a 3D volume file" << endl;
			cout << "Syntax: createvol <printf-style input filename format> <first index> <last index> <dc> <output.vol>" << endl;
		} else {
			createVol(argv[1], atoi(argv[2]), atoi(argv[3]), (Float) atof(argv[4]), argv[5]);
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
