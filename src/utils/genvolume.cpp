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
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/util.h>
#if defined(WIN32)
#include <mitsuba/core/getopt.h>
#endif


MTS_NAMESPACE_BEGIN

class GridDataSource {
    public:
    enum EVolumeType {
        EFloat32 = 1,
        EFloat16 = 2,
        EUInt8 = 3,
        EQuantizedDirections = 4
    };
};

class GenerateVolumeDataSource : public Utility {
public:

    void help(){
        cout<<"Generates an input file for the `gridvolume` plugin"<<endl;
        cout<<"mtsutil genvolume [options] <filename>"<<endl;
        cout<<"    -f              Set the output format to Float32 (Default)"<<endl;
        cout<<"    -i              Set the output format to UInt8"<<endl;
        cout<<"    -q              Set the output format to QuantizedDirections (2xUInt8), overwrites channels"<<endl;
        cout<<"    -c <channels>   Set the number of channels (1 or 3), Default: 1"<<endl;
        cout<<"    -v <value>      Set the constant value of the volume, Default: 0.0f"<<endl;
        cout<<endl;
        cout<<"    -b <xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax>"<<endl;
        cout<<"                    Defines the integer bounding box of the filled volume, Default: [0,0,0],[1,1,1]"<<endl;
    }

    int run(int argc, char **argv) {
        ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
        int optchar;
        char *end_ptr = NULL;
        
        //init options with defaults
        GridDataSource::EVolumeType type = GridDataSource::EVolumeType::EFloat32;
        int channels = 1;
        TAABB<Point3i> bounds(Point3i(0),Point3i(1));
        Float value = 0.0f;

        optind = 1;

        /* Parse command-line arguments */
        while ((optchar = getopt(argc, argv, "hfiqc:v:b:")) != -1) {
            switch (optchar) {
                case 'h': {
                        help();
                        return 0;
                    }
                    break;
                case 'f': {
                        type = GridDataSource::EVolumeType::EFloat32;
                    }
                    break;
                case 'i': {
                        type = GridDataSource::EVolumeType::EUInt8;
                    }
                    break;
                case 'q': {
                        type = GridDataSource::EVolumeType::EQuantizedDirections;
                    }
                    break;
                case 'c': {
                        channels = strtol(optarg, &end_ptr,10);
                        if (*end_ptr != '\0')
                            SLog(EError, "Could not parse the channels argument!");
                    }
                    break;
                case 'v': {
                        value = (Float) strtod(optarg, &end_ptr);
                        if (*end_ptr != '\0')
                            SLog(EError, "Could not parse the value argument!");
                    }
                    break;
                case 'b': {
                        std::vector<std::string> tokens = tokenize(optarg, ", ");
                        if (tokens.size() != 6)
                            Log(EError, "Invalid number of bounding box parameters!");
                        for (int i=0; i<3; ++i) {
                            bounds.min[i] = std::strtol(tokens[i].c_str(), &end_ptr, 10);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse integer number "
                                    "in bounding box parameter!");
                        }
                        for (int i=3; i<6; ++i) {
                            bounds.max[i-3] = std::strtol(tokens[i].c_str(), &end_ptr, 10);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse integer number "
                                    "in bounding box parameter!");
                        }
                    }
                    break;
            };
        }

        if (optind == argc || optind+1 < argc) {
            help();
            return 0;
        }
        if(type == GridDataSource::EVolumeType::EQuantizedDirections)
            channels = 2;
        else if(channels != 1 && channels != 3)
            Log(EError, "Invalid number of channels specified for data type. Please specify 1 or 3.");
            

        assert(bounds.isValid());

        cout << bounds.toString().c_str() << endl;
        fs::path outputFile = fileResolver->resolve(argv[optind]);
        ref<FileStream> fs = new FileStream(outputFile, FileStream::ETruncWrite);

        char header[] = {'V','O','L'};
        fs->write(&header, sizeof(header));

        fs->writeChar(3);
        fs->writeInt(type);

        int xres = (bounds.max.x - bounds.min.x)+1;
        int yres = (bounds.max.y - bounds.min.y)+1;
        int zres = (bounds.max.z - bounds.min.z)+1;
        fs->writeInt(xres);
        fs->writeInt(yres);
        fs->writeInt(zres);

        fs->writeInt(channels);

        fs->writeSingle(bounds.min.x);
        fs->writeSingle(bounds.min.y);
        fs->writeSingle(bounds.min.z);
        fs->writeSingle(bounds.max.x);
        fs->writeSingle(bounds.max.y);
        fs->writeSingle(bounds.max.z);

        size_t size = xres * yres * zres * channels;
        void * data;

        switch((GridDataSource::EVolumeType)type) {
            case GridDataSource::EVolumeType::EFloat32: {
                std::vector<float> vector(size, value);
                data = vector.data();
                size *= sizeof(float);
            }
            break;
            
            case GridDataSource::EVolumeType::EQuantizedDirections:
            case GridDataSource::EVolumeType::EUInt8: {
                std::vector<unsigned char> vector(size, value);
                data = vector.data();
                size *= sizeof(unsigned char);
            }
            break;

            default:
                Log(EError, "EVolumeType not supported");
                break;
        };

        fs->write(data, size);
        fs->close();

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(GenerateVolumeDataSource, "Generate a constant volume for `gridvolume` plugin")
MTS_NAMESPACE_END
