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


/** Constructing an Adaptive Grid Structure (and output) for PMATE
 *  - Compute a regular grid of Poincare map data in CR3BP
 *    with Adaptive Cell Subdivision and store output data file (AGND)
 *
 *  Author:  Wayne Schlei, Xavier Tricoche (Purdue University)
 *  ----------------------------------------------------------------------
 *  Inputs:
 *  ----------------------------------------------------------------------
 *  1) Map Analysis Parameters (either manual/default or by input file)
 *  2) System Parameters: C, mu, lstar(in km)
 *  3) Grid Bounds (bbox)
 *  4) Grid Resolution
 *  ----------------------------------------------------------------------
 *  Outputs:
 *  ----------------------------------------------------------------------
 *  1) Adaptive Grid Data File
 *  2) Map Analysis Parameters File (if not input)
 *  3) Point cloud file (PSI) for all returns (Optional)
 */
#include <iostream>
#include <vector>
#include <list>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>
#include <limits>
#include <iomanip> //C++11
#include <ctime>

// math
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

//API - xavier
#include <data/grid.hpp>
#include <data/adaptive_grid.hpp>
#include <data/cell.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/definitions.hpp>
#include <maps/adaptiveGridMapInterface.hpp>
#include <maps/map_analysis.hpp>
#include <maps/winding.hpp>
#include <maps/poincare_map.hpp>
//API - cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/cr3bpCellChecker.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/multipleAngleTracker.hpp>
//API - PMATE classes
#include <pmate/AdaptiveSampler.hpp>


//Outputs
#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif
//#include <mpi.h>

using namespace nvis;

// RHS and ODE solver parameters
double eps, C, mu, K;

typedef nvis::fixed_vector<double, 42>                                vec42;
typedef nvis::fixed_vector<double, 6>                                 vec6;
typedef nvis::fixed_matrix<double, 6>                                 mat6;;

//Type Defs
typedef xavier::grid<double, 2>                                              Plane_Grid_Type;
typedef xavier::AdaptiveGrid<std::vector<double>,WindingVecTraits >          AdaptiveGrid_Type;
typedef xavier::map_analysis_param                                           MapParams;
typedef xavier::orbit_data                                                   OrbitDataType;
typedef xavier::AdaptiveGridNodeData<OrbitDataType, AdaptiveGrid_Type,MapParams>        DataSet_Type;
typedef AdaptiveGrid_Type::data_type                                         data_pair_type;
typedef AdaptiveGrid_Type::id_type                                           id_type;
typedef xavier::dp5wrapper<double, 42>                                       ode_solver;
typedef orbital::cr3bp                                                       rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>                             section_type;
typedef orbital::MultipleAngleTracker<42, 3>                                 tracker_type;
typedef xavier::poincare_map<rhs_type, ode_solver, section_type >            map_type;
typedef std::vector<double>                                                  winding_type;
typedef xavier::CR3BP_Convexity<AdaptiveGrid_Type,id_type,winding_type>      CellCheckType;
typedef map_type::return_type                                                return_type;
typedef pmate::AdaptiveSampler<map_type,OrbitDataType,tracker_type,winding_type,WindingVecTraits,CellCheckType,MapParams>  AdaptiveSamplerType;

//------ Analysis ------
typedef nvis::fixed_vector<double,6>                              vec6;
typedef nvis::fixed_vector<double,7>                              vec7;
typedef nvis::fixed_vector<double,42>                             vec42;
typedef nvis::fixed_matrix<double,6,6>                            mat6;
typedef std::vector<nvis::vec2>                                   orbit_type;


//Static declarations
const double invalid_double = std::numeric_limits<double>::max();
//Default Value traits for adaptive grid
const double xavier::DefaultValueTraits::invalid = 10000.;
const double xavier::DefaultValueTraits::unknown = 0.0;
//Winding traits
const std::vector<double> WindingVecTraits::invalid(3,10000.);
const std::vector<double> WindingVecTraits::unknown(3,0.0);
//MaxDepth comparison
int xavier::OnEdgeCompare::maxDepth = 3;
int xavier::ID_Compare::maxDepth = 3;

using namespace xavier;
using namespace pmate;

// *******************************    PARAMETERS     *******************************
std::string me;
std::string filename;
std::string paramsFile;
nvis::ivec2 res(10, 10);

void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Compute Adaptive Grid for PMATE in restricted 3-body problem.\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -ip | --params <string>          Input MapAnalysisParam filename\n"
            << " -wp | --write_params <int>       Option to write params file (0=No,1=Yes)\n"
            << " -md | --maxDepth <int>           Maximum depth level of grid\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -b  | --bounds <float> x 4       Image bounds\n"
            << " -r  | --resolution <int> x 2     Image resolution\n"
            << " -C <double>                      C constant\n"
            << " -m  | --mu <double>              mu constant\n"
            << " -l  | --lstar <double>           Normalization distance in km\n"
            << " -R1 | --radius1 <float>          Exclusion Radius around P1 in CR3BP (nonDim)\n"
            << " -R2 | --radius2 <float>          Exclusion Radius around P2 in CR3BP (nonDim)\n"
            << " -p  | --nump <int>               Number of Map Iterations\n"
            << " -d  | --dir  <int>               Direction (>=0:Positive,<0:Negative)\n"
            << " -t  | --threads <int>            Set number of executed threads\n"
            << " -o  | --output <string>          Output base name (.agnd)\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n"
            << "\n"
            << "NOTE: Input options will overwrite MapAnalysisParam file if stated AFTER '-ip filename' statement\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    bool bounds_set = false;
    bool positiveDir = true;
    int maxp = 100;
    int max_depth = 3;
    filename = "none";
    paramsFile = "none";
    bool writeParams = true;
    
    //Standard Earth-Moon system example
    C = 2.96;
    mu = 1.21505714306e-2;
    double lstar = 384388.174;
    double R1 = 6378.14 / lstar;
    double R2 = 1738.20 / lstar;
    nvis::bbox2 bounds;
    bounds.min()[0] = -0.4;
    bounds.min()[1] = -2.5;
    bounds.max()[0] = 1.1;
    bounds.max()[1] = 2.5;
    
    
    //Map Params
    MapParams   mParams; //Assigns defaults
    //mParams.the_metric = ; //Note, default is non-periodic bounds to +/- infinity
    mParams.nb_iterations = maxp;
    mParams.max_depth = 3;
    xavier::OnEdgeCompare::maxDepth = 3;
    std::vector<double> wTols(3,0.5), wDist(3,1.0);
    wTols[1] = 1.0;
    wDist[1] = 2000.0; //w_xydot
    mParams.winding_convexity_tols = wTols;
    mParams.winding_cell_maxDist = wDist;
    
    
    size_t maxNumThreads = 1;
#if _OPENMP
    maxNumThreads = omp_get_max_threads();
#endif
    size_t nthreads = maxNumThreads;
    
    
    /// Read & Handle Options: --------------------------------------------------
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-C") {
            if (i == argc-1) {
                printUsageAndExit("missing C constant");
            }
            C = strtod(argv[++i],NULL);
        } else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) {
                printUsageAndExit("missing mu constant");
            }
            mu = strtod(argv[++i],NULL);
        } else if (arg == "-l" || arg == "--lstar") {
            if (i == argc-1) {
                printUsageAndExit("mission lstar value");
            }
            lstar = strtod(argv[++i],NULL);
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing bounds");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            mParams.bounds = bounds;
            bounds_set = true;
        } else if (arg == "-ip" || arg == "--params") {
            if (i == argc-1) {
                printUsageAndExit("missing params name");
            }
            paramsFile = argv[++i];
            bool readOK = mParams.read(paramsFile.c_str());
            if (!readOK) {
                throw std::runtime_error( "Param file read failure!" );
            }
            bounds_set = true;
            //Setting some parameters based on MapAnalysisParam input
            maxp = mParams.nb_iterations;
            xavier::OnEdgeCompare::maxDepth = mParams.max_depth;
            writeParams = false;
        } else if (arg == "-wp" || arg == "--write_params") {
            if (i == argc-1) {
                printUsageAndExit( "Missing write params option" );
            }
            if ( atoi(argv[++i]) == 1 ) {
                writeParams = true;
            }
        } else if (arg == "-md" || arg == "--maxDepth") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum depth value");
            }
            max_depth = atoi(argv[++i]);
            if(max_depth != mParams.max_depth) {
                writeParams = true;
            }
            mParams.max_depth = max_depth;
            xavier::OnEdgeCompare::maxDepth = max_depth;
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
            mParams.resolution = res;
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
            mParams.samplingIntegTol = eps;
        } else if (arg == "-p" || arg == "--maxp") {
            if (i == argc-1) {
                printUsageAndExit("missing max period");
            }
            maxp = atoi(argv[++i]);
            if (mParams.nb_iterations != maxp) {
                writeParams = true;
            }
            mParams.nb_iterations = maxp;
        } else if (arg == "-d" || arg == "--dir") {
            if (i == argc-1) {
                printUsageAndExit("missing direction");
            }
            int d = atoi(argv[++i]);
            positiveDir = (d>=0);
        } else if (arg == "-t" || arg == "--threads") {
            if (i == argc-1) {
                printUsageAndExit("missing number of threads");
            }
            nthreads = atoi(argv[++i]);
        }
        
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            filename = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
            mParams.verbose = (verbose<1) ? false : true;
            mParams.debug = (verbose>1) ? true : false;
        }
        
        else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (!bounds_set || filename=="none") {
        printUsageAndExit("");
    }
    
    
    //Threads
    if (nthreads > maxNumThreads) {
        std::cerr << "Warning:  Indicated thread number (" << nthreads
                  << ") is larger than max available (" << maxNumThreads
                  << ").  Setting thread number to " << maxNumThreads << "\n";
        nthreads = maxNumThreads;
    }
#if _OPENMP
    omp_set_num_threads( nthreads );
#endif
    std::cerr << nthreads << " threads employed for computation\n";
    //mParams.nb
    
    
    
    /// Data Objects
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    section.isPositive = positiveDir;
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    
    //Adaptive Sampler Object stores the primary data objects and makes calling simple
    AdaptiveSamplerType* theSampler = new AdaptiveSamplerType(bounds,res,pmap,mParams);
    //Adjust Convexity Checker for CR3BP
    theSampler->cellChecker.setBodyRadii(R1,R2);
    theSampler->cellChecker.setMu(mu);
    theSampler->cellChecker.tols = mParams.winding_convexity_tols;
    theSampler->cellChecker.maxDists = mParams.winding_cell_maxDist;
    
    nvis::timer timer;
    //Call the sampling function
    theSampler->sample();
    
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. \n";
    
    //Store to file
    theSampler->write(filename.c_str());
    
    
    //Write map parameters to file if necessary
    if (writeParams) {
        if (paramsFile == "none") {
            paramsFile = "mapParams.Run";
            time_t rawtime;
            struct tm* timeinfo;
            char buffer[80];
            time (&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer,80,"%d-%m-%Y_%I:%M:%S",timeinfo);
            std::string tstr(buffer);
            paramsFile = tstr + std::string(".param");
        } else {
            paramsFile = "New_" + paramsFile;
        }
        mParams.write( paramsFile.c_str() );
    }
    
    return 0;
}
