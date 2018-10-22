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


/// Generating a noise image for orbit convolution

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

std::string me;
std::string filename;


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
            << " -r  | --resolution <int> x 2     Image resolution (eg, 256x256)\n"
            << " -s  | --seed <uint>              Seed\n"
            << " -o  | --output <string>          Output Image File (.nrrd)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    ivec2 res(256);
    unsigned int seed = 41484;

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--seed") {
            if (i == argc-1) {
                printUsageAndExit("missing seed");
            }
            seed = (unsigned int) atoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            filename = argv[++i];
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if ( filename.empty() ) {
        printUsageAndExit("");
    }

    //Random number generator
    boost::mt19937 rng;
    if (seed > 0) {
        rng.seed(seed);
    }
    boost::uniform_01<> dist;

    //Initialize color noise
    const unsigned int resx = res[0];
    const unsigned int resy = res[1];
    std::vector<nvis::fvec3> image(resx*resy);
    for (unsigned int i=0; i<resx*resy; i++) {
        image[i][0] = dist(rng);
        image[i][1] = dist(rng);
        image[i][2] = dist(rng);
    }

    //Storing the output
    size_t dims[3] = {3, resx, resy};
    Nrrd* nout = nrrdNew();
    if( nrrdWrap_nva(nout, &image.front(), nrrdTypeFloat, 3, dims) ||
            nrrdSave(filename.c_str(), nout, NULL) ) {
        std::cerr << "ERROR while creating noise image: " << biffGetDone(NRRD) << "\n";
        nrrdNuke(nout);
    }

    return 0;
}
