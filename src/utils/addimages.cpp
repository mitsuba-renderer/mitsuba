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

class AddImages : public Utility {
public:
	int run(int argc, char **argv) {
		if (argc != 6) {
			cout << "Add the weighted pixel values of two EXR images to produce a new one" << endl;
			cout << "Syntax: mtsutil addimages <weight 1> <image 1.exr> <weight 2> <image 2.exr> <target.exr>" << endl;
			return -1;
		}
		char *end_ptr = NULL;
		Float weight1 = (Float) strtod(argv[1], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");
		Float weight2 = (Float) strtod(argv[3], &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse floating point value");

		ref<FileStream> aFile   = new FileStream(argv[2], FileStream::EReadOnly);
		ref<FileStream> bFile   = new FileStream(argv[4], FileStream::EReadOnly);
		ref<FileStream> outFile = new FileStream(argv[5], FileStream::ETruncReadWrite);

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
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AddImages, "Generate linear combinations of EXR images")
MTS_NAMESPACE_END
