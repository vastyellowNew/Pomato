/*************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body
        problem for orbital mechanics applications.

Authors: Wayne Schlei and Xavier Tricoche

Copyright (c) 2013-2018, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/


/// Single Orbit Convolution Pass:
/// Main function for a single convolution pass given NRRD input image and NRRD data
/// Note:  On pass 1, the input image should be colored noise (run noiseGen for generating noise image)

#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <util/wall_timer.hpp>
#include <math/fixed_vector.hpp>
#include <convolution/orbitConvolution.hpp>
#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace nvis;
using namespace OrbitConvolution;

std::string me;
std::string filename;
std::string imageFile;
std::string dataFile;


void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Single Orbit convolution pass computation\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -i  | --image <string>           Input Image File (.nrrd)\n"
            << " -p  | --maxp <unsigned int>      Max number of returns to consider\n"
            << " -d  | --data <string>            Input Data File (ReturnData in .nrrd)\n"
            << " -k  | --pass <int>               Pass number\n"
            << " -np | --numpasses <int>          Total number of passes\n"
            << " -s  | --shortfall <int>          Indicate how to handle shortfall\n"
            << "                                  (0:Drop,1:Penalty,2:Penalty w/Skew,3:UseAnyway)\n"
            << " -o  | --output <string>          Output Image File (.nrrd)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    unsigned int maxp = 50;
    int pass = 0, numPasses = 1;
    int shortfall = 0;
    filename = "none";
    imageFile = "none";
    dataFile = "none";

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-i" || arg == "--image") {
            if (i == argc-1) {
                printUsageAndExit("missing input image");
            }
            imageFile = argv[++i];
        } else if (arg == "-d" || arg == "--data") {
            if (i == argc-1) {
                printUsageAndExit("missing input data");
            }
            dataFile = argv[++i];
        } else if (arg == "-p" || arg == "--maxp") {
            if (i == argc-1) {
                printUsageAndExit("missing max period");
            }
            maxp = atoi(argv[++i]);
        } else if (arg == "-k" || arg == "--pass") {
            if (i == argc-1) {
                printUsageAndExit("missing pass number");
            }
            pass = atoi(argv[++i]);
        } else if (arg == "-np" || arg == "--numpasses") {
            if (i == argc-1) {
                printUsageAndExit("missing total number of passes");
            }
            numPasses = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--shortfall") {
            if (i == argc-1) {
                printUsageAndExit("missing shortfall descriptor");
            }
            shortfall = atoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            filename = argv[++i];
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (filename=="none" || imageFile=="none" || dataFile=="none") {
        printUsageAndExit("");
    }

    //Open Map Data file
    Nrrd* nDataPtr = nrrdNew();
    if( nrrdLoad(nDataPtr,dataFile.c_str(),NULL) ) {
        std::cerr << "Error opening input map file " << dataFile << "\n";
        exit( EXIT_FAILURE );
    }

    //Open current image (initial image made with initializeNoise.cpp)
    Nrrd* nImage = nrrdNew();
    if( nrrdLoad(nImage, imageFile.c_str(), NULL) ) {
        std::cerr << "Error opening image file " << imageFile << "\n";
        exit( EXIT_FAILURE );
    }
    //Load from NRRD file
    const unsigned int colorDim = nImage->axis[0].size;
    const unsigned int resx = nImage->axis[1].size;
    const unsigned int resy = nImage->axis[2].size;
    std::vector<nvis::fvec3> imageData(resx*resy, 0.0);
    //size_t size[3] = {3, resx, resy}; //How it's stored in nrrd
    //Loop through the data and set the values
    const float* nImagePtr = (const float*)nImage->data;
    for (unsigned int n=0; n<resx*resy; ++n) {
        //Set color -> NRRD data goes color, horizontal, vertical
        for (unsigned int k=0; k<colorDim; k++) {
            imageData[n][k] = *(nImagePtr);
            nImagePtr++;
        }
    }
    /*//Sanity check for load noise image into memory
    size_t sd[3] = {3, resx, resy};
    Nrrd *noiseOut = nrrdNew();
    if( nrrdWrap_nva(noiseOut, &imageData.front(), nrrdTypeFloat, 3, sd) ||
        nrrdSave("noiseCheck.nrrd",noiseOut,NULL) )
    {
      std::cerr << "ERROR while checking noise image data load : " << biffGetDone(NRRD) << "\n";
      nrrdNuke(noiseOut);
    }*/


    //The convolution pass
    bool ok = singleOrbitConvolutionPass(pass,numPasses,maxp,nDataPtr,imageData,shortfall);
    if (!ok) {
        std::cerr << "Orbit Convolution Pass failed!\n";
        exit( EXIT_FAILURE );
    }

    //Storing the output
    size_t dims[3] = {3, resx, resy};
    Nrrd* nout = nrrdNew();
    if( nrrdWrap_nva(nout, &imageData.front(), nrrdTypeFloat, 3, dims) ||
            nrrdSave(filename.c_str(), nout, NULL) ) {
        std::cerr << "ERROR while saving image: " << biffGetDone(NRRD) << "\n";
        nrrdNuke(nout);
    }

    return 0;
}
