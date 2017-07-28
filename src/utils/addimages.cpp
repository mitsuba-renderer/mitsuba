/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/util.h>

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

        ref<Bitmap> aBitmap = new Bitmap(Bitmap::EOpenEXR, aFile);
        ref<Bitmap> bBitmap = new Bitmap(Bitmap::EOpenEXR, bFile);

        /* A few sanity checks */
        if (aBitmap->getPixelFormat() != bBitmap->getPixelFormat())
            Log(EError, "Error: Input bitmaps have a different pixel format!");
        if (aBitmap->getComponentFormat() != bBitmap->getComponentFormat())
            Log(EError, "Error: Input bitmaps have a different component format!");
        if (aBitmap->getSize() != bBitmap->getSize())
            Log(EError, "Error: Input bitmaps have a different size!");

        ref<Bitmap> outBitmap = new Bitmap(aBitmap->getPixelFormat(),
                aBitmap->getComponentFormat(), aBitmap->getSize());

        size_t nEntries =
            (size_t) aBitmap->getSize().x *
            (size_t) aBitmap->getSize().y *
            aBitmap->getChannelCount();

        switch (aBitmap->getComponentFormat()) {
            case Bitmap::EFloat16: {
                    half *aData = aBitmap->getFloat16Data();
                    half *bData = bBitmap->getFloat16Data();
                    half *outData = outBitmap->getFloat16Data();
                    for (size_t i=0; i<nEntries; ++i)
                        *outData++ = (half) ((float) std::max((Float) 0,
                                weight1 * (Float) (*aData++) +
                                weight2 * (Float) (*bData++)));
                }
                break;

            case Bitmap::EFloat32: {
                    float *aData = aBitmap->getFloat32Data();
                    float *bData = bBitmap->getFloat32Data();
                    float *outData = outBitmap->getFloat32Data();
                    for (size_t i=0; i<nEntries; ++i)
                        *outData++ = (float) std::max((Float) 0,
                                weight1 * (Float) (*aData++) +
                                weight2 * (Float) (*bData++));
                }
                break;

            case Bitmap::EUInt32: {
                    uint32_t *aData = aBitmap->getUInt32Data();
                    uint32_t *bData = bBitmap->getUInt32Data();
                    uint32_t *outData = outBitmap->getUInt32Data();
                    for (size_t i=0; i<nEntries; ++i)
                        *outData++ = (uint32_t) std::max((Float) 0,
                                weight1 * (Float) (*aData++) +
                                weight2 * (Float) (*bData++));
                }
                break;

            default:
                Log(EError, "Unsupported component format!");
        }

        outBitmap->write(Bitmap::EOpenEXR, outFile);
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AddImages, "Generate linear combinations of EXR images")
MTS_NAMESPACE_END
