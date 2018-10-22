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


/** Creating Manifolds from PMATE FixedPointData in the CR3BP
 *  Authors:  Wayne Schlei, Xavier Tricoche (Purdue University)
 *
 *  Inputs:
 *   1) FixedPointData file
 *   2) MapParameters file (Optional)
 *      (preferably what was used to generate (1) )
 *   3) Manifold Settings (Defaults with modifiable options)
 *      (Can also input ManifoldSettings file (.mset) )
 *   4) Output ManifoldData filename (.im)
 *  Notes:
 *  You MUST input a FixedPointData file (.fpx)!
 *
 *
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
#include <topology/ManifoldClasses.hpp>
#include <pmate/FixedPointData.hpp>
#include <pmate/ManifoldData.hpp>


//Outputs
#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif
//#include <mpi.h>

using namespace nvis;

//Type Definitions
typedef xavier::grid<double, 2>                                             Plane_Grid_Type;
typedef xavier::dp5wrapper<double, 42>                                      ode_solver;
typedef orbital::cr3bp                                                      rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>                            section_type;
typedef orbital::MultipleAngleTracker<42, 3>                                tracker_type;
typedef poincare_map<rhs_type, ode_solver, section_type >                   MapType;
typedef xavier::map_analysis_param                                          MapParams;
#if defined(_WIN32) || defined(C_0X)
typedef MapType::return_type                                                return_type;
typedef MapType::lvec_type                                                  vec_type;
#else
typedef typename MapType::return_type                                       return_type;
typedef typename MapType::lvec_type                                         vec_type;
#endif

//------ Analysis ------
typedef nvis::fixed_vector<double,6>                                        vec6;
typedef nvis::fixed_vector<double,7>                                        vec7;
typedef nvis::fixed_vector<double,42>                                       vec42;
typedef nvis::fixed_matrix<double,6,6>                                      mat6;
typedef std::vector<nvis::vec2>                                             orbit_type;
typedef xavier::orbit_data                                                  OrbitDataType; //Note: 2D
#if defined(_WIN32) || defined(C_0X)
typedef xavier::EdgeRotationFailure<xavier::vec_type>                       EdgeRotFailure;
typedef xavier::EdgeRotationFailure<xavier::vec_type>                       MapDiscont;
#else
typedef xavier::EdgeRotationFailure<vec_type>                               EdgeRotFailure;
typedef xavier::EdgeRotationFailure<vec_type>                               MapDiscont;
#endif
// -- PMATE Objects --
static const int s = rhs_type::numSingularities;
typedef nvis::fixed_vector<double,s+1>                                      ExtendedMapDataVec;
#if defined(_WIN32) || defined(C_0X)
typedef topology::SortableReturnData<xavier::vec_type,ExtendedMapDataVec>   SortableData;
#else
typedef topology::SortableReturnData<vec_type,ExtendedMapDataVec>           SortableData;
#endif
typedef std::vector<xavier::fixpoint>                                       FixedPointChain;
typedef pmate::FixedPointData<xavier::fixpoint>                             FPData;
typedef pmate::ManifoldData<SortableData,xavier::fixpoint>                  ManifoldDataType;


//Static declarations
const double invalid_double = std::numeric_limits<double>::max();
//MaxDepth comparison
int xavier::OnEdgeCompare::maxDepth = 1;
int xavier::ID_Compare::maxDepth = 1;

using namespace xavier;
using namespace topology;
using namespace pmate;

// *******************************    PARAMETERS     *******************************
std::string me;
std::string manSetFileStr = "none";
std::string manDataFileStr = "none";
std::string newFPXFileStr = "none";

///Compose filenames from a base name
void composeFileNames( std::string& baseName )
{
    if (baseName != "none") {
        manSetFileStr = baseName + ".mset";
        newFPXFileStr = baseName + ".fpx";
        manDataFileStr = baseName + ".im";
    } else {
        throw std::runtime_error("Invalid baseName!");
    }
}


void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Computing Invariant Manifolds with PMATE in restricted 3-body problem (y=0 section).\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << "     * Means an element is required to run the program\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE::Manifolds Problem Definition:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -C  | --Hamiltonian <double>     Desired Hamiltonian (or C) constant (Default = 2.96)\n"
            << " -m  | --mu <double>              Gravitational parameter (mu, Default is Earth-Moon)\n"
            << " -e  | --eps <float>              Integration tolerance (default 1e-12)\n"
            << " -l  | --lstar <double>           Normalization distance in km\n"
            << " -d  | --dir  <int>               Direction (>=0:Positive [d],<0:Negative)\n"
            << " -t  | --threads <int>            Set number of executed threads\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE::Manifolds Input Files:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -ip | --params <string>          *Input MapAnalysisParam filename (needed)\n"
            << " -ix | --fpdata <string>          *Input FixedPointData (.fpx) filename\n"
            << " -is | --mSettings <string>       Input ManifoldSettings (.mset) filename\n"
            //<< " -im | --inputMData <string>      Input ManifoldData for continuation option\n"
            << " -ws | --write_mSet <int>         Option to write ManifoldSettings file (0=No,1=Yes[d])\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE::Manifolds Output : (*Need output base name OR ManifoldData filename)\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -o  | --output <string>          Output base name (if declared, identifies all output)\n"
            << " -v  | --verbose <int>            Verbose mode (0: off[d], 1: on)\n"
            << " -dv | --debug_verbose <int>      Debug statements (0: off[d], 1: on)\n"
            << " -om | --imFile <string>          *Computed invariant manifold (ManifoldData) output filename\n"
            << "                                  (always outputs)\n"
            //<< " -od | --write_mapDiscont <int>   Option to write map discont file (0=No[d],1=Yes)\n"
            //<< " -fd | --mapDiscontFile <string>  Map discontinuities output filename\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE::Manifolds  Manifold Settings for Advection with STHManifold algorithm:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -xp | --maxp                     Enable max_period cutoff for manifolds (d: off).  Will skip\n"
            << "                                  any manifold computation that triggers the condition : \n"
            << "                                     |theManifold.thePeriod| > theMapParams.max_period \n"
            << " -b  | --bounds <float> x 4       Computation bounds on section (x,xdot)\n"
            << "                                  Note: These refer to the ManifoldSettings.bounds, which\n"
            << "                                  govern the usable (or display) space for manifolds.\n"
            << "                                  (Default = MapParams.the_metric.bounds )\n"
            << " -ms | --use_manual_step          Utilize a manual step (with parameters) to start manifolds\n"
            << " -me | --meps <float>             Initial step (in 4D)  (Default = 1e-7)\n"
            << " -xs | --max_step <float>         Max allowed initial step (Default = 1e-5)\n"
            << " -si | --max_si <float>           Max allowed stablity index magnitude (Default = 1e10)\n"
            << " -xa | --max_alpha <float>        Max allowed angle (deg) between 3 manifold points\n"
            << "                                  (Default = 17.2deg = 0.3rad)\n"
            //<< " -la | --liberal_alpha <float>    Broader max Alpha (not really used)\n"
            << " -da | --delta_alpha_max <float>  Control parameter on sampling rate (delta*alpha=0 is a line)\n"
            << "                                  (Default = 1e-3)\n"
            << " -dm | --delta_min <float>        Min distance between manifold points (Default = 1e-4 nd)\n"
            << " -dx | --delta_max <float>        Max distance between manifold points (Default = 0.3 nd)\n"
            << " -md | --max_depth <int>          Maximum IM tree depth during manifold advection (stopping condition)\n"
            << "                                  (Default = 8; Warning >10 will take a long time!)\n"
            << " -xg | --enableMaxNumSeg          Enable max number of seeding segments constraint\n"
            << " -ns | --numMaxSegs <int>         Maximum number of segments in manifold\n"
            << "                                  (Default = 8000 segments)\n"
            << " -ng | --numMaxSeeds <int>        Maximum number of seeding segments in manifold\n"
            << "                                  (Default = 500 segments)\n"
            << " -xl | --enableMaxArclength       Enable max arc length stopping criteria\n"
            << " -al | --max_arc_length <float>   Maximum arc length for manifold (stopping condition)\n"
            << "                                  (Default = 30 nd)\n"
            << " -sl | --max_seed_length <float>  Maximum arc length for seeding segments on manifold\n"
            << "                                  (Default = 0.5 nd)\n"
            << " -tm | --dtau_min <float>         Min allowed distance between successive manifold segment \n"
            << "                                  sample points.  Used in creating new manifold arcs and \n"
            << "                                  transversality violation detection. \n"
            << "                                  (Default = 5e-6 nonDim. Note: it may be unwise\n"
            << "                                  to input a value lower than 1e-8.)\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Transversality Violation Detection Parameters:\n"
            << "   => Applied between subsequent points on a line segment existing on the section.\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -R1 | --radius1 <float>          *Exclusion Radius around P1 in CR3BP (nonDim)\n"
            << " -R2 | --radius2 <float>          *Exclusion Radius around P2 in CR3BP (nonDim)\n"
            << "                                  *Both radii should be defined for non-default system (EM)\n"
            << " -st | --sing_tol <float>         Tolerance on checking 'divide by 0' for r^3, so checks\n"
            << "                                  if [r^3 <= sing_tol] which is a singularity intersection\n"
            << "                                  (Default = 1e-20)\n"
            << " -tv | --max_trans_speed <float>  Max transverse speed differential\n"
            << "                                  (Default = 0.5 nonDim)\n"
            << " -td | --max_disp <float> x2      Max map displacement differential\n"
            << "                                  (Default = [0.2,2.0] nondim for [x,xdot])\n"
            << " -tt | --max_abs_time <float>     Max ABSOLUTE differential in time-of-flight\n"
            << "                                  (Default = 2.0 nondim time units)\n"
            << " -tr | --max_rel_time <float>     Max RELATIVE time-of-flight differential factor\n"
            << "                                  (Default = 0.15 or 15% differential)\n"
            << " -tf | --min_fp_tol <float>       Min map displacement magnitude to indicate a potential\n"
            << "                                  fixed point along an edge. (Default = 1e-4)\n"
            //<< " -tn | --ntEdgeDivisions <int>    Max number of subdivisions to consider on non-transverse edge\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE::Manifolds  Fixed Point Refinement Parameters (if needed):\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -sc | --sanityCheck              Run Sanity Check on input fixed points\n"
            << " -re | --refine_eps <float>       Refinement integration precision (Default = 1e-12)\n"
            << " -rt | --refine_tol <float>       Refinement convergence tolerance (Default = 1e-8)\n"
            << " -ri | --refine_max_iters <int>   Refinement maximum number of allowed iterations\n"
            << "                                  (Default = 20)\n"
            << " -rn | --refine_pts <int>         Refinement points per period (Default = 5)\n"
            << " -M  | --linear_stm <int>         Option to utilize linear STM for stability analysis, which\n"
            << "                                  computes the eigen-space basis for each unstable (saddle)\n"
            << "                                  fixed point. (0=Numerical STM, 1=Linear STM [default])\n"
            //<< " -of | --write_failures <int>     Option to write fixed point failures file (0=No[d],1=Yes)\n"
            //<< " -ff | --fpFailFile <string>      Fixed Points that FAIL sanity check output filename\n"
            << "\n"
            << "NOTE: Input options will overwrite ManifoldSettings file if stated AFTER '-is filename' statement\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    double eps = 1.0e-12;
    int verbose = 0;
    bool positiveDir = true;
    bool writeFailures = false;
    bool runSanityCheck = false;
    bool maxCutoff = false;
    std::string fpxFileStr = "none";

    //Standard Earth-Moon system example
    double C = 2.96;
    double mu = 1.21505714306e-2;
    double lstar = 384388.174;
    double R1 = 6378.14 / lstar;
    double R2 = 1738.20 / lstar;
    double sing_tol = 1e-20;


    //Map Analysis Params
    MapParams   theMapParams; //Assigns defaults
    theMapParams.nb_iterations = 50;
    theMapParams.max_depth = 1;
    xavier::OnEdgeCompare::maxDepth = 1;
    xavier::ID_Compare::maxDepth = 1;
    std::vector<double> wTols, wDist;
    for(int i=0; i<3; i++) {
        wTols.push_back( 0.5 );
        wDist.push_back( 1.0 );
    }
    wTols[1] = 1.0;
    wDist[1] = 2000.0; //w_xydot
    theMapParams.winding_convexity_tols = wTols;
    theMapParams.winding_cell_maxDist = wDist;
    theMapParams.samplingIntegTol = eps;
    std::string mpFileStr = "none";
    //Metric bounds influence manifold computation
    const double LARGE = std::numeric_limits<double>::max();
    theMapParams.the_metric.bounds().min()[0] = -LARGE;
    theMapParams.the_metric.bounds().min()[1] = -2.55;
    theMapParams.the_metric.bounds().max()[0] = LARGE;
    theMapParams.the_metric.bounds().max()[1] = 2.55;
    theMapParams.the_metric.periodic()[0] = false; //Not periodic bounds
    theMapParams.the_metric.periodic()[1] = false; //Not periodic bounds

    //Manifold Settings
    ManifoldSettings theManifoldSettings;
    theManifoldSettings.bounds.min() = nvis::vec2(-LARGE, -2.55);
    theManifoldSettings.bounds.max() = nvis::vec2(LARGE, 2.55);
    bool writeSettings = true;
    bool mSetReset = false;

    /// Data Objects
    //Right Hand Side
    rhs_type rhs(C, mu);
    //Section
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    section.isPositive = positiveDir;

    //The Poincare Map Engine
    MapType theMap(rhs, section);
    theMap.setPrecision(eps);


    //Max number of threads available
    size_t maxNumThreads = 1;
#if _OPENMP
    maxNumThreads = omp_get_max_threads();
#endif
    size_t nthreads = maxNumThreads;

    std::cout << "pmateManifoldCreator: input parameters:\n";
    for (int i=0; i<argc; ++i) {
        std::cout << argv[i] << " ";
    }
    std::cout << '\n';


    /// Read & Handle Options: --------------------------------------------------
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        }
        // PMATE:Manifolds Problem Definition ----------------------------------------------------
        else if (arg == "-C" || arg == "--Hamiltonian" ) {
            if (i == argc-1) {
                printUsageAndExit("missing Hamiltonian value");
            }
            C = atof(argv[++i]);
        } else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) {
                printUsageAndExit("missing mu value");
            }
            mu = atof(argv[++i]);
        } else if (arg == "-e" || arg == "--eps" ) {
            if (i == argc-1) {
                printUsageAndExit("missing integration tolerance");
            }
            eps = atof(argv[++i]);
            theMapParams.samplingIntegTol = eps;
        } else if (arg == "-l" || arg == "--lstar") {
            if (i == argc-1) {
                printUsageAndExit("missing lstar value");
            }
            lstar = atof(argv[++i]);
        } else if (arg == "-d" || arg == "--dir") {
            if (i == argc-1) {
                printUsageAndExit("missing positive section direction");
            }
            int d = atoi(argv[++i]);
            positiveDir = (d>=0);
        } else if (arg == "-t" || arg == "--threads") {
            if (i == argc-1) {
                printUsageAndExit("missing number of threads");
            }
            nthreads = atoi(argv[++i]);
        }
        // PMATE:Manifolds Input Files ------------------------------------------------------------
        else if (arg == "-ip" || arg == "--params") {
            if (i == argc-1) {
                printUsageAndExit("missing params name");
            }
            mpFileStr = argv[++i];
            std::cerr << "PMATE: Reading map parameters file...\n";
            bool readOK = theMapParams.read(mpFileStr.c_str());
            if (!readOK) {
                throw std::runtime_error( "Param file read failure!" );
            }
            //Setting some parameters based on MapAnalysisParam input
            xavier::OnEdgeCompare::maxDepth = theMapParams.max_depth;
            xavier::ID_Compare::maxDepth = theMapParams.max_depth;
            //Need to force sampling tolerance to be higher (at least default 1e-12)
            theMapParams.samplingIntegTol = eps;
        } else if (arg == "-ix" || arg == "--fpdata") {
            if (i == argc-1) {
                printUsageAndExit("missing input FixedPointData filename");
            }
            fpxFileStr = argv[++i];
        } else if (arg == "-is" || arg == "--mSettings") {
            if (i == argc-1) {
                printUsageAndExit("missing input ManifoldSettings filename");
            }
            manSetFileStr = argv[++i];
            bool readOK = theManifoldSettings.read(manSetFileStr.c_str());
            if (!readOK) {
                throw std::runtime_error( "ManifoldSettings file read failure!" );
            }
            mSetReset = false;
        } else if (arg == "-ws" || arg == "--write_mSet") {
            if (i == argc-1) {
                printUsageAndExit("missing write grid option");
            }
            writeSettings = (atoi(argv[++i])==1)? true : false;
        }
        // PMATE: Output -----------------------------------------------------------------
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output base name");
            }
            std::string baseStr = argv[++i];
            composeFileNames( baseStr );
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            theMapParams.verbose = (atoi(argv[++i])==1)? true : false;
        } else if (arg == "-dv" || arg == "--debug_verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing debug_verbose flag");
            }
            theMapParams.debug = (atoi(argv[++i])==1)? true : false;
        } else if (arg == "-om" || arg == "--imFile") {
            if (i == argc-1) {
                printUsageAndExit("missing output ManifoldData filename");
            }
            manDataFileStr = argv[++i];
        }
        /*else if (arg == "-fd" || arg == "--mapDiscontFile") {
            if (i == argc-1) printUsageAndExit("missing map discont file name");
            topoExtractor->mapDiscontFileStr = argv[++i];
        }
        else if (arg == "-ff" || arg == "--fpFailFile") {
            if (i == argc-1) printUsageAndExit("missing fixed point failures file name");
            topoExtractor->fpFailFileStr = argv[++i];
        }
        else if (arg == "-od" || arg == "--write_mapDiscont") {
            if (i == argc-1) printUsageAndExit("missing write map discont option");
            writeMapDiscont = (atoi(argv[++i])==1)? true : false;
        }
        else if (arg == "-of" || arg == "--write_failures") {
            if (i == argc-1) printUsageAndExit("missing write fixed point failures option");
            writeFailures = (atoi(argv[++i])==1)? true : false;
        }*/
        // PMATE:Manifolds Sampling Parameters ----------------------------------------------------
        else if (arg == "-xp" || arg == "--maxp") {
            maxCutoff = true;
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing bounds");
            }
            theManifoldSettings.bounds.min()[0] = atof(argv[++i]);
            theManifoldSettings.bounds.min()[1] = atof(argv[++i]);
            theManifoldSettings.bounds.max()[0] = atof(argv[++i]);
            theManifoldSettings.bounds.max()[1] = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-ms" || arg == "--use_manual_step") {
            theManifoldSettings.manualStep = true;
            mSetReset = true;
        } else if (arg == "-me" || arg == "--meps") {
            if (i == argc-1) {
                printUsageAndExit("missing manual inital step");
            }
            theManifoldSettings.eps = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-xs" || arg == "--max_step") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum inital step");
            }
            theManifoldSettings.sdelta_max = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-si" || arg == "--max_si") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum stability index magnitude");
            }
            theManifoldSettings.max_si = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-xa" || arg == "--max_alpha") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum alpha value");
            }
            theManifoldSettings.alpha_max = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-da" || arg == "--delta_alpha_max") {
            if (i == argc-1) {
                printUsageAndExit("missing delta*alpha maximum value");
            }
            theManifoldSettings.delta_alpha_max = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-dm" || arg == "--delta_min") {
            if (i == argc-1) {
                printUsageAndExit("missing delta minimum value");
            }
            theManifoldSettings.delta_min = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-dx" || arg == "--delta_max") {
            if (i == argc-1) {
                printUsageAndExit("missing delta maximum value");
            }
            theManifoldSettings.delta_max = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-xg" || arg == "--enableMaxNumSeg") {
            theManifoldSettings.enableMaxSeg = true;
            mSetReset = true;
        } else if (arg == "-ns" || arg == "--numMaxSegs") {
            if (i == argc-1) {
                printUsageAndExit("missing number of maximum segments");
            }
            theManifoldSettings.maxNumSeg = atoi(argv[++i]);
            mSetReset = true;
        } else if (arg == "-ng" || arg == "--numMaxSeeds") {
            if (i == argc-1) {
                printUsageAndExit("missing number of maximum seeding segments");
            }
            theManifoldSettings.maxNumSeedSeg = atoi(argv[++i]);
            mSetReset = true;
        } else if (arg == "-xl" || arg == "--enableMaxArclength") {
            theManifoldSettings.enableMaxArclength = true;
            mSetReset = true;
        } else if (arg == "-al" || arg == "--max_arc_length") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum arc length");
            }
            theManifoldSettings.max_arc_length = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-sl" || arg == "--max_seed_length") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum seeding arc length");
            }
            theManifoldSettings.maxSeedArclength = atof(argv[++i]);
            mSetReset = true;
        } else if (arg == "-md" || arg == "--max_depth") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum depth for IM tree");
            }
            theManifoldSettings.maxTreeDepth = atoi(argv[++i]);
            mSetReset = true;
        } else if (arg == "-tm" || arg == "--dtau_min") {
            if (i == argc-1) {
                printUsageAndExit("missing delta tau min");
            }
            theManifoldSettings.delta_tau_min = atof(argv[++i]);
            mSetReset = true;
        }
        // PMATE: Transversality Violation Parameters--------------------------------------
        else if (arg == "-R1" || arg == "--radius1") {
            if (i == argc-1) {
                printUsageAndExit("missing exclusion radius for P1");
            }
            R1 = atof(argv[++i]);
        } else if (arg == "-R2" || arg == "--radius2") {
            if (i == argc-1) {
                printUsageAndExit("missing exclusion radius for P2");
            }
            R2 = atof(argv[++i]);
        } else if (arg == "-st" || arg == "--sing_tol") {
            if (i == argc-1) {
                printUsageAndExit("missing singularity tolerance");
            }
            sing_tol = atof(argv[++i]); //Impacts RHS objects
        } else if (arg == "-tv" || arg == "--max_trans_speed") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum transverse speed");
            }
            theMapParams.maxDeltaTransSpeed = atof(argv[++i]);
        } else if (arg == "-td" || arg == "--max_disp") {
            if (i == argc-2) {
                printUsageAndExit("missing maximum displacement");
            }
            theMapParams.maxDeltaMapDisplacement[0] = atof(argv[++i]);
            theMapParams.maxDeltaMapDisplacement[1] = atof(argv[++i]);
        } else if (arg == "-tt" || arg == "--max_abs_time") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum absolute time differntial");
            }
            theMapParams.maxTimeChange = atof(argv[++i]);
        } else if (arg == "-tr" || arg == "--max_rel_time") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum relative time differential");
            }
            theMapParams.maxDeltaTimeFactor = atof(argv[++i]);
        } else if (arg == "-tf" || arg == "--min_fp_tol") {
            if (i == argc-1) {
                printUsageAndExit("missing minimum fixed point detection tolerance");
            }
            theMapParams.min_fp_tol = atof(argv[++i]);
        }
        // PMATE:Manifolds - Fixed Point Refinement Parameters ---------------------------------------
        else if (arg == "-sc" || arg == "--sanityCheck") {
            runSanityCheck = true;
        } else if (arg == "-re" || arg == "--refine_eps") {
            if (i == argc-1) {
                printUsageAndExit("missing refinement integration tolerance");
            }
            theMapParams.refinementIntegTol = atof(argv[++i]);
        } else if (arg == "-rt" || arg == "--refine_tol") {
            if (i == argc-1) {
                printUsageAndExit("missing refinement convergence tolerance");
            }
            theMapParams.refinementConvTol = atof(argv[++i]);
        } else if (arg == "-ri" || arg == "--refine_max_iters") {
            if (i == argc-1) {
                printUsageAndExit("missing refinement max iters");
            }
            theMapParams.refinementMaxIters = atoi(argv[++i]);
        } else if (arg == "-rn" || arg == "--refine_pts") {
            if (i == argc-1) {
                printUsageAndExit("missing refinement points per period");
            }
            theMapParams.refinementPointsPerPeriod = atoi(argv[++i]);
        } else if (arg == "-M" || arg == "--linear_stm") {
            if (i == argc-1) {
                printUsageAndExit("missing refinement STM option");
            }
            theMapParams.linearMonodromy = (atoi(argv[++i])==1)? true : false;
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (fpxFileStr == "none" || mpFileStr =="none" || manDataFileStr == "none") {
        printUsageAndExit("One of fpxFileStr, mpFileStr, manDataFileStr was not defined in input");
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

    //Load FixedPointData
    FPData* fpData = new FPData;
    fpData->read( fpxFileStr.c_str() );
    //Do C and FPData->C match?
    if (C != fpData->getJacobiConstant() ) {
        std::cerr << "Warning::  Hamiltonian mis-match between input (" << C << ") and FixedPointData ("
                  << fpData->getJacobiConstant() << ")!\n";
        std::cerr << "It is on the user to check to make sure the system mup and C values match input data.\n";
        std::cerr << "Adjusting to FixedPointData input C value...\n";
        //return 1;
        C = fpData->getJacobiConstant();
    }

    //Adjust objects to handle new values
    rhs_type adjustedRHS(C,mu);
    adjustedRHS.setSingularityTolerance( sing_tol );
    adjustedRHS.setSingularitySafeDistance(0,R1/lstar);//Applies to pre-filtering manifolds for computation.
    adjustedRHS.setSingularitySafeDistance(1,R2/lstar);
    section.set_rhs( adjustedRHS );
    theMap.setNewComponents( adjustedRHS, section );
    //Start process with sampling integration tolerance
    theMap.setPrecision(theMapParams.samplingIntegTol);



    //Sanity Check on FixedPointData input
    std::vector< FixedPointChain > chains;
#if defined(_WIN32) || defined(C_0X)
    std::vector< FixedPointChain >::iterator chainIT;
#else
    typename std::vector< FixedPointChain >::iterator chainIT;
#endif
    fpData->getData(chains);
    int numOrbits = (int) chains.size(), numSaddles = 0;
    for(chainIT=chains.begin(); chainIT!=chains.end(); ++chainIT) {
        if ((*chainIT)[0].saddle) {
            numSaddles++;
        }
    }
    //Examine full-period propagation to see if it's within tolerances
    if(runSanityCheck) {
        //Run a parallel computation to propagate all orbits (fixed points?)
    }
    //If removing some saddles, we have to store the new FixedPointData object
    newFPXFileStr = fpxFileStr;


    //Build the ManifoldData object and prep for computation:
    ManifoldDataType* theManifoldData =
        new ManifoldDataType((*fpData));
    theManifoldData->setNumThreads( (int) nthreads );
    if (maxCutoff) {
        theManifoldData->setMaxCutoffFlag(true);
    }


    //The Manifold Computation:
    std::cerr << std::setprecision(14) << argv[0] << ": Executing STHManifold algorithm for the system :\n"
              << "---------------------------------------------------------------------------\n"
              << "         mu = " << mu << "\n"
              << "          C = " << C << "\n"
              << "       bbox = [" << theManifoldSettings.bounds.min() << "->"
              << theManifoldSettings.bounds.max() << "]\n"
              << " numSaddles = " << numSaddles << "\n"
              << "---------------------------------------------------------------------------\n";

    nvis::timer timer;

    //Call the advection function
    theManifoldData->compute(theMap,theMapParams,theManifoldSettings);


    double elapsed = timer.elapsed();
    std::cerr << "-------------------------------------------------------------------------\n";
    std::cerr << "\nComputation took " << elapsed/3600.0 << " hrs = " << elapsed << " s. \n";
    std::cerr << "-------------------------------------------------------------------------\n";

    //Store solution
    theManifoldData->write(manDataFileStr.c_str(), newFPXFileStr.c_str());

    //Additional storage if specified
    if (writeSettings && (manSetFileStr != "none") ) {
        theManifoldSettings.write(manSetFileStr.c_str());
    }

    //if (writeMapDiscont && (mapDiscontFileStr != "none") )
    //     theManifoldData->writeMapDiscont();
    //if (writeFailures && (fpFailFileStr != "none") )
    //     writeFailedFixedPoints();


    return 0;
}
