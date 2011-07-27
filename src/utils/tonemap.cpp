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
#include <mitsuba/core/timer.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/bitmap.h>
#if defined(WIN32)
#include <mitsuba/core/getopt.h>
#endif

MTS_NAMESPACE_BEGIN

class Tonemap : public Utility {
public:
	void help() {
		cout << endl;
		cout << "Synopsis: Loads one or more linear EXR images and writes tonemapped 8-bit PNG/JPGs";
		cout << endl;
		cout << "Usage: mtsutil tonemap [options] <EXR file (s)>" << endl;
		cout << "Options/Arguments:" << endl;
		cout << "   -h             Display this help text" << endl << endl;
		cout << "   -g gamma       Specify the gamma value (The default is -1 => sRGB)" << endl << endl;
		cout << "   -m multiplier  Multiply the pixel values by 'multiplier' (Default = 1)" << endl << endl;
		cout << "   -f fmt         Specifies the output format (png/jpg, default:png)" << endl << endl;
	}

	inline float toSRGB(float value) {
		if (value < 0.0031308f)
			return 12.92f * value;
		return 1.055f * std::pow(value, 0.41666f) - 0.055f;
	}


	int run(int argc, char **argv) {
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		char optchar, *end_ptr = NULL;
		optind = 1;
		Float gamma = -1, multiplier = 1;
		Bitmap::EFileFormat format = Bitmap::EPNG;

		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "hg:m:f:")) != -1) {
			switch (optchar) {
				case 'h': {
						help();
						return 0;
					}
					break;
				case 'g': 
					gamma = (Float) strtod(optarg, &end_ptr);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the gamma value!");
					break;
				case 'f': {
					  std::string fmt = optarg;
					  if (fmt == "png")
						  format = Bitmap::EPNG;
					  else if (fmt == "jpg" || fmt == "jpeg")
						  format = Bitmap::EJPEG;
					  else
						  SLog(EError, "Unknown format! (must be png/jpg)");
					}
					break;
				case 'm': 
					multiplier = (Float) strtod(optarg, &end_ptr);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the multiplier!");
					break;
			};
		}

		if (optind == argc) {
			help();
			return 0;
		}

		Float invGamma = 1.0f/gamma;

		for (int i=optind; i<argc; ++i) {
			fs::path inputFile = fileResolver->resolve(argv[i]);
			Log(EInfo, "Loading EXR image \"%s\" ..", inputFile.file_string().c_str());
			ref<FileStream> is = new FileStream(inputFile, FileStream::EReadOnly);
			ref<Bitmap> input = new Bitmap(Bitmap::EEXR, is);
			ref<Bitmap> output = new Bitmap(input->getWidth(), input->getHeight(), 32);
			float *inputData = input->getFloatData();
			uint8_t *outputData = output->getData();
			for (int y=0; y<input->getHeight(); ++y) {
				for (int x=0; x<input->getWidth(); ++x) {
					size_t idx = y*input->getWidth() + x;
					for (int i=0; i<3; ++i) {
						if (gamma == -1)
							outputData[idx*4+i] = (uint8_t) std::max(std::min((int) (toSRGB(inputData[idx*4+i] * multiplier) * 255), 255), 0);
						else
							outputData[idx*4+i] = (uint8_t) std::max(std::min((int) (std::pow(inputData[idx*4+i] * multiplier, invGamma)*255), 255), 0);
						outputData[idx*4+3] = (uint8_t) std::max(std::min((int) (inputData[idx*4+3] * 255), 255), 0);
					}
				}
			}

			fs::path outputFile = inputFile;
			if (format == Bitmap::EPNG) {
				outputFile.replace_extension(".png");
				Log(EInfo, "Writing tonemapped PNG image \"%s\" ..", outputFile.file_string().c_str());
				ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
			} else if (format == Bitmap::EJPEG) {
				outputFile.replace_extension(".jpg");
				Log(EInfo, "Writing tonemapped JPEG image \"%s\" ..", outputFile.file_string().c_str());
			} else {
				Log(EError, "Unknown format!");
			}

			ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
			output->save(format, os);
		}
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Tonemap, "Simple EXR->PNG tonemapper")
MTS_NAMESPACE_END
