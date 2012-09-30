/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

		ref<Bitmap> rBitmap = new Bitmap(Bitmap::EOpenEXR, rFile);
		ref<Bitmap> gBitmap = new Bitmap(Bitmap::EOpenEXR, gFile);
		ref<Bitmap> bBitmap = new Bitmap(Bitmap::EOpenEXR, bFile);
		rBitmap = rBitmap->separateChannel(0);
		gBitmap = gBitmap->separateChannel(0);
		bBitmap = bBitmap->separateChannel(0);

		std::vector<Bitmap *> sourceBitmaps;
		sourceBitmaps.push_back(rBitmap);
		sourceBitmaps.push_back(gBitmap);
		sourceBitmaps.push_back(bBitmap);

		ref<Bitmap> result = Bitmap::join(Bitmap::ERGBA, sourceBitmaps);
		ref<FileStream> outFile = new FileStream(s4, FileStream::ETruncReadWrite);
		result->write(Bitmap::EOpenEXR, outFile);
	}

	int run(int argc, char **argv) {
		if (argc < 5) {
			cout << "Join three monochromatic images into a RGB-valued EXR file" << endl;
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
