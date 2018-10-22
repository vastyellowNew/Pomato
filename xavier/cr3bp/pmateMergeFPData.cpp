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


/// Combine FixedPointData Objects using filters to remove duplicates
/// Author:  Wayne Schlei
/// Date:    7/23/2016

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
// math
#include <util/wall_timer.hpp>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

//API - maps
#include <maps/DP45wrapper.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/poincare_map.hpp>
//API - cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/multipleAngleTracker.hpp>
#include <maps/fixpoints.hpp>
//API - pmate
#include <pmate/FixedPointData.hpp>
#include <pmate/FixedPointDataFilters.hpp>


using namespace nvis;
typedef pmate::FixedPointData<xavier::fixpoint>  FPData;
#if defined(_WIN32) || defined(C_0X)
typedef FPData::FixedPointChain         FixedPointChain;
typedef FPData::FPChainsContainer       FPChainsContainer;
typedef FixedPointChain::iterator       FPChainIter;
typedef FPChainsContainer::iterator     FPChainsContainerIter;
#else
typedef typename FPData::FixedPointChain         FixedPointChain;
typedef typename FPData::FPChainsContainer       FPChainsContainer;
typedef typename FixedPointChain::iterator       FPChainIter;
typedef typename FPChainsContainer::iterator     FPChainsContainerIter;
#endif
typedef xavier::dp5wrapper<double, 42>                                      ode_solver;
typedef orbital::cr3bp                                                      rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>                            section_type;
typedef orbital::MultipleAngleTracker<42, 3>                                tracker_type;
typedef xavier::poincare_map<rhs_type, ode_solver, section_type >           MapType;
typedef xavier::map_analysis_param                                          MapParams;



std::string me("");
std::string fpxFile("none");
std::string fpx2File("none");
std::vector<std::string> inputFiles;
std::string outputFile("none");



void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Merge & Filter fixed point data.  It is assumed that the \n"
            << "fixed points are computed for the same section.  This is built for the \n"
            << "section y=0.  Other sections will require an adjustment to this program.\n"
            << "OPTIONS:\n"
            << " -h  | --help                         Print this information\n"
            << " -i  | --input <string> ... <string>  Input .fpx file[s]\n"
            << " -m  | --mu <double>              Gravity parameter for CR3BP\n"
            << " -d  | --dir  <int>               Direction (>=0:Positive [d],<0:Negative)\n"
            << " -e  | --eps <double>             Tolerance on 'same point' test [d=5e-6]\n"
            << " -p  | --maxp <int>               Max Period cutoff (optional)\n"
            << " -b  | --bounds <float> x 4       Bounds for box filter (optional)\n"
            << " -o  | --output <string>          Output .fpx File\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    double eps = 5.e-6;
    double mu = 1.21505714306e-2;
    int maxp = -1;
    bool maxCutoff = false;
    bool merge = false;
    bool positiveDir = true;
    bool boxFilter = false;
    const double LARGE = std::numeric_limits<double>::max();
    nvis::bbox2 bounds;
    bounds.min()[0] = -LARGE;
    bounds.min()[1] = -LARGE;
    bounds.max()[0] = LARGE;
    bounds.max()[1] = LARGE;
    //Define inputs
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                printUsageAndExit("missing input name");
            }
            fpxFile = argv[++i];
            //Loop remaining files until we find another option:
            int j=i+1;
            std::string tempStr = argv[j]; //Next argument
            while(tempStr.at(0) != '-' && tempStr.at(0) != '\0') {
                //Add file to set of input files
                inputFiles.push_back(tempStr);
                //Set the new string
                j++;
                i++;
                tempStr = argv[j];
            }
            if ((int)inputFiles.size() >= 1) {
                merge = true;
            }
        } else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) {
                printUsageAndExit("missing mu value");
            }
            mu = atof(argv[++i]);
        } else if (arg == "-d" || arg == "--dir") {
            if (i == argc-1) {
                printUsageAndExit("missing positive section direction");
            }
            int d = atoi(argv[++i]);
            positiveDir = (d>=0);
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing tolerance");
            }
            eps = (double) atof(argv[++i]);
        } else if (arg == "-p" || arg == "--maxp") {
            if (i == argc-1) {
                printUsageAndExit("missing max period cutoff");
            }
            maxp = atoi(argv[++i]);
            maxCutoff = (maxp>=0);
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing bounds");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            boxFilter = true;
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            outputFile = argv[++i];
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if ( outputFile=="none" || fpxFile=="none" ) {
        printUsageAndExit("");
    }
    
    
    std::cerr << "---------------------------------------------------------------------\n";
    std::cerr << "PMATE Merge Fixed Point Data \n";
    std::cerr << "---------------------------------------------------------------------\n";
    
    //Open a FixedPointData object from a file
    FPData* theData = new FPData(fpxFile.c_str());
    std::cerr << "Read File : " << fpxFile << " with " << theData->getNumOrbits() << " orbits\n";
    
    //Open a second object from file if input
    if ( merge ) {
        int numFiles = (int) inputFiles.size();
        for(int i=0; i<numFiles; i++) {
            fpx2File = inputFiles[i]; //Secondary files
            FPData fpx2Data( fpx2File.c_str() );
            std::cerr << "Read Secondary File: " << fpx2File << " with " << fpx2Data.getNumOrbits() << " orbits\n";
            //Check that the Hamiltonian values are consistent
            if( theData->getJacobiConstant() != fpx2Data.getJacobiConstant() ) {
                //if not, print error and exit
                std::cerr << "ERROR: Hamiltonian values in input FixedPointData files do NOT match!\n";
                std::cerr << "  Input 1: " << fpxFile << " : H = " << theData->getJacobiConstant() << "\n";
                std::cerr << "  Input 2: " << fpx2File << " : H = " << fpx2Data.getJacobiConstant() << "\n";
                std::cerr << "  Exiting pmateMergeFPData.\n";
                return 1;
            }
        }
    }
    
    //The Hamiltonian value from results:
    double C = theData->getJacobiConstant();
    
    //Load additional data if second input is on
    if (merge) {
        int numFiles = (int) inputFiles.size();
        //Loop through secondary files to add data to main set:
        for(int i=0; i<numFiles; i++) {
            FPData fpx2Data( inputFiles[i].c_str() ); //Read construction
            FPChainsContainer data2;
            fpx2Data.getData( data2 );
            //Add to back of overall data
            theData->insertData( data2 );
        }
    }
    
    //Build Poincare Map object
    //Right-hand Side (CR3BP EOMs)
    rhs_type rhs(C, mu);
    //Section
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    section.isPositive = positiveDir;
    
    //The Poincare Map Engine
    MapType theMap(rhs, section);
    theMap.setPrecision(1.e-12);
    //Only need dummy map parameters object for now...
    MapParams   theMapParams; //Assigns defaults
    theMapParams.max_depth = 1;
    std::vector<double> wTols, wDist;
    for(int i=0; i<3; i++) {
        wTols.push_back( 0.5 );
        wDist.push_back( 1.0 );
    }
    wTols[1] = 1.0;
    wDist[1] = 2000.0; //w_xydot
    theMapParams.winding_convexity_tols = wTols;
    theMapParams.winding_cell_maxDist = wDist;
    theMapParams.the_metric.bounds().min()[0] = -LARGE;
    theMapParams.the_metric.bounds().min()[1] = -2.55;
    theMapParams.the_metric.bounds().max()[0] = LARGE;
    theMapParams.the_metric.bounds().max()[1] = 2.55;
    if(boxFilter) {
        theMapParams.the_metric.bounds() = bounds;
    }
    //Probably the most important part for map params:
    theMapParams.the_metric.periodic()[0] = false; //Not periodic bounds
    theMapParams.the_metric.periodic()[1] = false; //Not periodic bounds
    //If period cutoff is enabled
    if(maxCutoff) {
        theMapParams.max_period = maxp;
    }
    
    //Run the filtering steps by calling pmate function:
    pmate::filterFixedPointData((*theData), theMap, theMapParams, eps, maxCutoff, boxFilter);
    
    
    //Write output of filtering
    theData->write( outputFile.c_str() );
    std::cerr << "---------------------------------------------------------------------\n";
    std::cerr << "Total AFTER filtering = " << theData->getNumOrbits() << " orbits.\n";
    std::cerr << "---------------------------------------------------------------------\n";
    
    //Free info
    delete theData;
    
    return 0;
}
