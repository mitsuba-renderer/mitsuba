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
#include <mitsuba/core/timer.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>
#if defined(WIN32)
# include <mitsuba/core/getopt.h>
#endif
#ifdef MTS_OPENMP
# include <omp.h>
#endif

MTS_NAMESPACE_BEGIN

class Tonemap : public Utility {
public:
	void help() {
		cout << endl;
		cout << "Synopsis: Loads one or more EXR/RGBE images and writes tonemapped 8-bit PNG/JPGs";
		cout << endl;
		cout << "Usage: mtsutil tonemap [options] <EXR/RGBE file (s)>" << endl;
		cout << "Options/Arguments:" << endl;
		cout << "   -h             Display this help text" << endl << endl;
		cout << "   -g gamma       Specify the gamma value (The default is -1 => sRGB)" << endl << endl;
		cout << "   -m multiplier  Multiply the pixel values by 'multiplier' (Default = 1)" << endl << endl;
		cout << "   -b r,g,b       Color balance: apply the specified per-channel multipliers" << endl << endl;
		cout << "   -c x,y,w,h     Crop: tonemap a given rectangle instead of the entire image" << endl << endl;
		cout << "   -s w,h         Resize the output image to the specified resolution" << endl << endl;
		cout << "   -r x,y,w,h,i   Add a rectangle at the specified position and intensity, e.g." << endl;
		cout << "                  to make paper figures. The intensity should be in [0, 255]." << endl << endl;
		cout << "   -f fmt         Request a certain output format (png/jpg, default:png)" << endl << endl;
		cout << "   -a             Require the output image to have an alpha channel" << endl << endl;
		cout << "   -p key,burn    Run Reinhard et al.'s photographic tonemapping operator. 'key'" << endl;
		cout << "                  between [0, 1] chooses between low and high-key images and" << endl
			 << "                  'burn' (also [0, 1]) controls how much highlights may burn out" << endl << endl;
		cout << "   -x             Temporal coherence mode: activate this flag when tonemapping " << endl
			 << "                  frames of an animation using the '-p' option to avoid flicker" << endl << endl;
		cout << "   -o file        Save the output with a given filename" << endl << endl;
		cout << "   -t             Multithreaded: process several files in parallel" << endl << endl;
		cout << " The operations are ordered as follows: 1. crop, 2. resize, 3. color-balance, " << endl;
		cout << " 4. tonemap, 5. annotate. To simply process a directory full of EXRs in " << endl;
		cout << " parallel, run the following: 'mtsutil tonemap -t path-to-directory/*.exr'" << endl;
	}

	typedef struct {
		int r[5];
	} Rect;

	int run(int argc, char **argv) {
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		int optchar;
		char *end_ptr = NULL;
		optind = 1;
		Float gamma = -1, multiplier = 1;
		Bitmap::EFileFormat format = Bitmap::EPNG;
		Float cbal[] = {1, 1, 1};
		int crop[] = {0, 0, -1, -1};
		int resize[] = {-1, -1};
		Float tonemapper[] = {-1, -1};
		bool temporalCoherence = false;
		std::vector<Rect> rects;
		std::string outputFilename;
		Bitmap::EPixelFormat pixelFormat = Bitmap::ERGB;
		Float logAvgLuminance = 0, maxLuminance = 0;
		bool runParallel = false;
		ReconstructionFilter *rfilter = NULL;

		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "htxag:m:f:r:b:c:o:p:s:")) != -1) {
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

				case 'x':
					temporalCoherence = true;
					break;

				case 'f': {
					  std::string fmt = boost::to_lower_copy(std::string(optarg));
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

				case 'a':
					pixelFormat = Bitmap::ERGBA;
					break;

				case 'b': {
						std::vector<std::string> tokens = tokenize(optarg, ", ");
						if (tokens.size() != 3)
							Log(EError, "Invalid color balancing parameter!");
						for (int i=0; i<3; ++i) {
							cbal[i] = (Float) std::strtod(tokens[i].c_str(), &end_ptr);
							if (*end_ptr != '\0')
								Log(EError, "Cannot parse floating point number "
									"in color balancing parameter!");
						}
					}
					break;

				case 'c': {
						std::vector<std::string> tokens = tokenize(optarg, ", ");
						if (tokens.size() != 4)
							Log(EError, "Invalid crop parameter!");
						for (int i=0; i<4; ++i) {
							crop[i] = (int) std::strtol(tokens[i].c_str(), &end_ptr, 10);
							if (*end_ptr != '\0')
								Log(EError, "Cannot parse integer in crop parameter!");
						}
					}
					break;

				case 'p': {
						std::vector<std::string> tokens = tokenize(optarg, ", ");
						if (tokens.size() != 2)
							Log(EError, "Invalid tone mapper parameter!");
						for (int i=0; i<2; ++i) {
							tonemapper[i] = (Float) std::strtod(tokens[i].c_str(), &end_ptr);
							if (*end_ptr != '\0')
								Log(EError, "Cannot parse tone mapper parameters!");
						}
					}
					break;

				case 's': {
						std::vector<std::string> tokens = tokenize(optarg, ", ");
						if (tokens.size() != 2)
							Log(EError, "Invalid resize parameter!");
						for (int i=0; i<2; ++i) {
							resize[i] = (int) std::strtol(tokens[i].c_str(), &end_ptr, 10);
							if (*end_ptr != '\0')
								Log(EError, "Cannot parse integer in resize parameter!");
						}
					}
					break;


				case 'r': {
						std::vector<std::string> tokens = tokenize(optarg, ", ");
						if (tokens.size() != 5)
							Log(EError, "Invalid rectangle parameter!");
						Rect r;
						for (int i=0; i<5; ++i) {
							r.r[i] = (int) std::strtol(tokens[i].c_str(), &end_ptr, 10);
							if (*end_ptr != '\0')
								Log(EError, "Cannot parse integer in rectangle parameter!");
						}
						rects.push_back(r);
					}
					break;

				case 'o':
					outputFilename = optarg;
					break;

				case 't':
					runParallel = true;
					break;
			}
		}

		if (runParallel) {
			if (outputFilename != "" || temporalCoherence) {
				Log(EWarn, "Requested multithreaded tonemapping along with incompatible options, disabling threading..");
				runParallel = false;
			} else {
				Log(EInfo, "Performing multithreaded tonemapping ..");
			}
		}

		if (optind == argc) {
			help();
			return 0;
		}

		if (pixelFormat == Bitmap::ERGBA && format == Bitmap::EJPEG)
			Log(EError, "JPEG images do not support an alpha channel!");

		if (resize[0] != -1) {
			/* A resampling operation was requested; use a Lanczos Sinc reconstruction filter by default */
			rfilter = static_cast<ReconstructionFilter *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(ReconstructionFilter), Properties("lanczos")));
			rfilter->configure();
		}

		if (runParallel) {
			ref<Logger> logger = Thread::getThread()->getLogger();

			#if defined(MTS_OPENMP)
				#pragma omp parallel for schedule(static)
			#endif
			for (int i=optind; i<argc; ++i) {
				Thread *thread = Thread::getThread();
				if (!thread) {
					thread = Thread::registerUnmanagedThread("omp");
					thread->setLogger(logger);
				}

				fs::path inputFile = fileResolver->resolve(argv[i]);
				Log(EInfo, "Loading EXR image \"%s\" ..", inputFile.string().c_str());
				ref<FileStream> is = new FileStream(inputFile, FileStream::EReadOnly);
				ref<Bitmap> input = new Bitmap(Bitmap::EOpenEXR, is);

				if (crop[2] != -1 && crop[3] != -1)
					input = input->crop(Point2i(crop[0], crop[1]), Vector2i(crop[2], crop[3]));

				if (resize[0] != -1)
					input = input->resample(rfilter, ReconstructionFilter::EClamp,
							ReconstructionFilter::EClamp, Vector2i(resize[0], resize[1]));

				if (cbal[0] != 1 || cbal[1] != 1 || cbal[2] != 1)
					input->colorBalance(cbal[0], cbal[1], cbal[2]);

				if (tonemapper[0] != -1) {
					Float logAvgLuminance = 0, maxLuminance = 0;
					input->tonemapReinhard(logAvgLuminance, maxLuminance, tonemapper[0], tonemapper[1]);
					Log(EInfo, "Tonemapper reports: log-average luminance = %f, max. luminance = %f",
						logAvgLuminance, maxLuminance);
				}

				ref<Bitmap> output = input->convert(pixelFormat, Bitmap::EUInt8, gamma, multiplier);

				for (size_t i=0; i<rects.size(); ++i) {
					int *r = rects[i].r;
					output->drawRect(Point2i(r[0], r[1]), Vector2i(r[2], r[3]), Spectrum(r[4]/255.0f));
				}

				fs::path outputFile = inputFile;
				if (format == Bitmap::EPNG)
					outputFile.replace_extension(".png");
				else if (format == Bitmap::EJPEG)
					outputFile.replace_extension(".jpg");
				else
					Log(EError, "Unknown target format!");

				Log(EInfo, "Writing tonemapped image to \"%s\" ..", outputFile.string().c_str());

				ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
				output->write(format, os);
			}

		} else {
			for (int i=optind; i<argc; ++i) {
				fs::path inputFile = fileResolver->resolve(argv[i]);
				Log(EInfo, "Loading EXR image \"%s\" ..", inputFile.string().c_str());
				ref<FileStream> is = new FileStream(inputFile, FileStream::EReadOnly);
				ref<Bitmap> input = new Bitmap(Bitmap::EOpenEXR, is);

				if (crop[2] != -1 && crop[3] != -1)
					input = input->crop(Point2i(crop[0], crop[1]), Vector2i(crop[2], crop[3]));

				if (resize[0] != -1)
					input = input->resample(rfilter, ReconstructionFilter::EClamp,
					ReconstructionFilter::EClamp, Vector2i(resize[0], resize[1]));

				if (cbal[0] != 1 || cbal[1] != 1 || cbal[2] != 1)
					input->colorBalance(cbal[0], cbal[1], cbal[2]);

				if (tonemapper[0] != -1) {
					input->tonemapReinhard(logAvgLuminance, maxLuminance, tonemapper[0], tonemapper[1]);
					Log(EInfo, "Tonemapper reports: log-average luminance = %f, max. luminance = %f",
						logAvgLuminance, maxLuminance);
					if (!temporalCoherence) {
						logAvgLuminance = 0;
						maxLuminance = 0;
					}
				}

				ref<Bitmap> output = input->convert(pixelFormat, Bitmap::EUInt8, gamma, multiplier);

				for (size_t i=0; i<rects.size(); ++i) {
					int *r = rects[i].r;
					output->drawRect(Point2i(r[0], r[1]), Vector2i(r[2], r[3]), Spectrum(r[4]/255.0f));
				}

				fs::path outputFile = inputFile;
				if (outputFilename == "") {
					if (format == Bitmap::EPNG)
						outputFile.replace_extension(".png");
					else if (format == Bitmap::EJPEG)
						outputFile.replace_extension(".jpg");
					else
						Log(EError, "Unknown target format!");
				} else {
					outputFile = outputFilename;
				}

				Log(EInfo, "Writing tonemapped image to \"%s\" ..", outputFile.string().c_str());

				ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
				output->write(format, os);
			}
		}
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Tonemap, "Command line batch tonemapper")
MTS_NAMESPACE_END
