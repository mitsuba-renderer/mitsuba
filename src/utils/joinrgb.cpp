/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/util.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

MTS_NAMESPACE_BEGIN

class JoinRGB : public Utility {
public:
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

		for (int y=0; y<rBitmap->getHeight(); ++y) {
			for (int x=0; x<rBitmap->getWidth(); ++x) {
				float r = rData[(x + y * width) * 4];
				float g = gData[(x + y * width) * 4 + 1];
				float b = bData[(x + y * width) * 4 + 2];
				outData[(x+y * width) * 4 + 0] = r;
				outData[(x+y * width) * 4 + 1] = g;
				outData[(x+y * width) * 4 + 2] = b;
				outData[(x+y * width) * 4 + 3] = 1;
			}
		}

		outBitmap->save(Bitmap::EEXR, outFile);
	}

	int run(int argc, char **argv) {
		if (argc < 5) {
			cout << "Join three monochromatic EXRs into a colored image" << endl;
			cout << "joinrgb <red.exr> <green.exr> <blue.exr> <combined.exr>" << endl;
		} else {
			joinRGB(argv[1], argv[2], argv[3], argv[4]);
		}
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(JoinRGB, "Join three monochromatic EXRs into a colored image");
MTS_NAMESPACE_END
