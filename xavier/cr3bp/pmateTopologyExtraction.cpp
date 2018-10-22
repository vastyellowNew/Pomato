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


/** Executing PMATE in the CR3BP over a specified domain
 *  Authors:  Wayne Schlei, Xavier Tricoche (Purdue University)
 *
 *  Execution Steps:
 *  1)Compute a regular grid of Poincare map data in CR3BP
 *    with Adaptive Cell Subdivision and store output data file (AGND)
 *    (like pmateAdaptiveGrid).  This step runs only if no AGND file is
 *    given as input.
 *  2)Evaluate Poincare index on analysis cells (from adaptive grid)
 *    to detect fixed points and evaluate fixed point guesses.
 *  3)Fixed point refinement utilizing a sequence of differential
 *    corrections procedures.
 *  4)Propagate invariant manifolds for unstable fixed points on the
 *    section utilizing the STHmanifold procedure (manifolds via
 *    curve-refinement).
 *  Note:  Invariant manifold propagation for unstable fixed points
 *  may be conducted with a separate function (pmateManifoldCreator).
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
#include <pmate/FixedPointData.hpp>
#include <pmate/FixedPointDataFilters.hpp>
#include <pmate/AdaptiveSampler.hpp>
#include <pmate/TopologyExtractor.hpp>


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
typedef poincare_map<rhs_type, ode_solver, section_type >                   map_type;
typedef xavier::map_analysis_param                                          MapParams;
#if defined(_WIN32) || defined(C_0X)
typedef map_type::return_type                                               return_type;
typedef map_type::lvec_type                                                 lvec_type;
#else
typedef typename map_type::return_type                                      return_type;
typedef typename map_type::lvec_type                                        lvec_type;
#endif

//------ Analysis ------
typedef nvis::fixed_vector<double,6>                                        vec6;
typedef nvis::fixed_vector<double,7>                                        vec7;
typedef nvis::fixed_vector<double,42>                                       vec42;
typedef nvis::fixed_matrix<double,6,6>                                      mat6;
typedef std::vector<nvis::vec2>                                             orbit_type;
typedef std::vector<double>                                                 winding_type;
typedef nvis::ivec3                                                         NodeType;
typedef std::pair<orbit_type, nvis::fvec3>                                  color_orbit_type;
typedef xavier::orbit_data                                                  OrbitDataType; //Note: 2D
typedef xavier::TerminalEdge<NodeType>                                      EdgeType;
typedef xavier::CompositeEdge<NodeType>                                     CompEdgeType;
typedef xavier::EdgeRotationFailure<lvec_type>                              EdgeRotFailure;
#if defined(_WIN32) || defined(C_0X)
typedef std::set<EdgeType>                                                   EdgeSet;
typedef std::set<CompEdgeType>                                               CompositeEdgeSet;
typedef std::set<EdgeType>::const_iterator                                   EdgeSetIter;
typedef std::set<CompEdgeType>::const_iterator                               CompositeEdgeSetIter;
typedef xavier::AdaptiveGrid<winding_type,WindingVecTraits>                  AdaptiveGrid_Type;
typedef AdaptiveGrid_Type::data_type                                         data_pair_type;
typedef AdaptiveGrid_Type::id_type                                           id_type; //Same as NodeType
typedef xavier::AdaptiveGridNodeData<OrbitDataType,AdaptiveGrid_Type,MapParams>  DataSet;
typedef DataSet::DataMapIterator                                             DataMapIterator;
typedef xavier::Cell<DataSet,NodeType>                                       MapCell;
typedef xavier::CR3BP_Convexity<AdaptiveGrid_Type,id_type,winding_type>      CellCheckType;
#else
typedef typename std::set<EdgeType>                                          EdgeSet;
typedef typename std::set<CompEdgeType>                                      CompositeEdgeSet;
typedef typename std::set<EdgeType>::const_iterator                          EdgeSetIter;
typedef typename std::set<CompEdgeType>::const_iterator                      CompositeEdgeSetIter;
typedef xavier::AdaptiveGrid<winding_type,WindingVecTraits>                  AdaptiveGrid_Type;
typedef typename AdaptiveGrid_Type::data_type                                data_pair_type;
typedef typename AdaptiveGrid_Type::id_type                                  id_type; //Same as NodeType
typedef xavier::AdaptiveGridNodeData<OrbitDataType,AdaptiveGrid_Type,MapParams>  DataSet;
typedef typename DataSet::DataMapIterator                                    DataMapIterator;
typedef xavier::Cell<DataSet,NodeType>                                       MapCell;
typedef xavier::CR3BP_Convexity<AdaptiveGrid_Type,id_type,winding_type>      CellCheckType;
#endif
// -- PMATE Objects --
typedef pmate::AdaptiveSampler<map_type,OrbitDataType,tracker_type,winding_type,WindingVecTraits,CellCheckType,MapParams>   AdaptiveSamplerType;
typedef pmate::TopologyExtractor<map_type,OrbitDataType,tracker_type,winding_type,WindingVecTraits,CellCheckType,MapParams> TopologyExtractorType;

//Static declarations
const double invalid_double = std::numeric_limits<double>::max();
//Default Value traits for adaptive grid
const double xavier::DefaultValueTraits::invalid = 10000.;
const double xavier::DefaultValueTraits::unknown = 0.0;
//Winding traits
const std::vector<double> WindingVecTraits::invalid(3,10000.);
const std::vector<double> WindingVecTraits::unknown(3,0.0);
//MaxDepth comparison
int xavier::OnEdgeCompare::maxDepth = 1;
int xavier::ID_Compare::maxDepth = 1;

using namespace xavier;
using namespace pmate;

// *******************************    PARAMETERS     *******************************
std::string me;


void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Extracting Topology with PMATE in restricted 3-body problem (y=0 section).\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << "     * Means an element is required to run the program\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Problem Definition:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -C  | --Hamiltonian <double>     Desired Hamiltonian (or C) constant (Default = 2.96)\n"
            << " -m  | --mu <double>              Gravitational parameter (mu, Default is Earth-Moon)\n"
            << " -b  | --bounds <float> x 4       *Computation bounds on section (x,xdot)\n"
            << " -r  | --resolution <int> x 2     Initial (depth=0) sampling GRID resolution\n"
            << "                                  (Default = [10,10]).  Note that the two values correspond\n"
            << "                                  with the number of vertices in the intial grid,\n"
            << "                                  NOT the number of analysis cells.\n"
            << " -l  | --lstar <double>           Normalization distance in km (set this for not Earth-Moon!)\n"
            << " -d  | --dir  <int>               Direction (>=0:Positive [d],<0:Negative)\n"
            << " -t  | --threads <int>            Set number of executed threads\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Input Files:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -ip | --params <string>          *Input MapAnalysisParam filename (needed if no bounds/res)\n"
            << " -ig | --grid <string>            Input Adaptive Grid Node Data (.agnd) filename\n"
            << " -wp | --write_params <int>       Option to write params file (0=No,1=Yes[d])\n"
            << " -wg | --write_grid <int>         Option to write agnd file (0=No,1=Yes)\n"
            << "                                  Skips writing grid data if input file (-ig) is specified\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Output : (*Need output base name OR fixed point file name)\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -o  | --output <string>          Output base name (if declared, identifies all output)\n"
            << " -v  | --verbose <int>            Verbose mode (0: off[d], 1: on)\n"
            << " -dv | --debug_verbose <int>      Debug statements (0: off[d], 1: on)\n"
            << " -fd | --mapDiscontFile <string>  Map discontinuities output filename\n"
            << " -fg | --fpGuessesFile <string>   Fixed point guesses output filename\n"
            << " -ff | --fpFailFile <string>      Failed fixed points output filename\n"
            << " -fx | --fpxFile <string>         *Computed fixed points output filename (always outputs)\n"
            << " -od | --write_mapDiscont <int>   Option to write map discont file (0=No[d],1=Yes)\n"
            << " -og | --write_guesses <int>      Option to write fixed point guesses file (0=No[d],1=Yes)\n"
            << " -of | --write_failures <int>     Option to write fixed point failures file (0=No[d],1=Yes)\n"
            << " -wc | --watch_cell <int> x3      Watch cell with given ID (i,j,depth)\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Sampling Parameters:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -p  | --sample_iters <int>       Number of Map Iterations for initial sampling and winding\n"
            << "                                  number evaluation (Default = 50 iterates)\n"
            << " -md | --maxDepth <int>           Maximum depth level of adaptive grid (Default = 1)\n"
            << " -e  | --eps <float>              Sampling integration precision (Default = 1e-8)\n"
            << " -wt | --winding_tol <float> x3   Winding number tolerances used in subdivision\n"
            << "                                  (Default = [0.5 1.0 0.5] on the 3 winding numbers)\n"
            << " -wd | --winding_dist <float> x3  Winding number allowed value variance within a cell\n"
            << "                                  (Default = [1.0 2000 1.0] on the 3 winding numbers)\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Poincare Index (Fixed Point Detection) Parameters:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -mp | --min_p <int>              Min analysis period cutoff (Default = 1)\n"
            << " -xp | --max_p <int>              Max analysis period cutoff (Default = 12)\n"
            << " -ma | --max_angle <float>        Max allowed rotation angle (deg) in Poincare index\n"
            << "                                  (Default = 135deg = 3*pi/4)\n"
            << " -lm | --lmin <float>             Min allowed distance between successive cell edge points;\n"
            << "                                  Used in tracking rotation and transversality violation\n"
            << "                                  detection. (Default = 3e-4 nonDim. Note: it may be unwise\n"
            << "                                  to input a value lower than 1e-8.)\n"
            << " --------------------------------------------------------------------------------------------\n"
            << "  PMATE Transversality Violation Detection Parameters:\n"
            << "   => Applied between subsequent points on a line segment existing on the section.\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -R1 | --radius1 <float>          *Exclusion Radius around P1 in CR3BP (km)\n"
            << " -R2 | --radius2 <float>          *Exclusion Radius around P2 in CR3BP (km)\n"
            << "                                  *Both should be input for non-default system (EM)\n"
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
            << "  PMATE Fixed Point Refinement Parameters:\n"
            << " --------------------------------------------------------------------------------------------\n"
            << " -sc | --subcell_res <int> x2     Sub-cell sampling resolution (guess generation)\n"
            << "                                  (Default = [6,6] point resolution)\n"
            << " -re | --refine_eps <float>       Refinement integration precision (Default = 1e-12)\n"
            << " -rt | --refine_tol <float>       Refinement convergence tolerance (Default = 1e-8)\n"
            << " -ri | --refine_max_iters <int>   Refinement maximum number of allowed iterations\n"
            << "                                  (Default = 20)\n"
            << " -rn | --refine_pts <int>         Refinement points per period (Default = 5)\n"
            << " -M  | --linear_stm <int>         Option to utilize linear STM for stability analysis, which\n"
            << "                                  computes the eigen-space basis for each unstable (saddle)\n"
            << "                                  fixed point. (0=Numerical STM, 1=Linear STM [default])\n"
            << "\n"
            << "NOTE: Input options will overwrite MapAnalysisParam file if stated AFTER '-ip filename' statement\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    me = argv[0];
    double eps = 1.0e-8;
    int verbose = 0;
    bool bounds_set = false;
    bool positiveDir = true;
    int niters = 50;
    int max_depth = 1;
    bool writeParams = true;
    bool writeGrid = false; //assume we have input file
    bool writeMapDiscont = false;
    bool writeGuesses = false;
    bool writeFailures = false;

    //Standard Earth-Moon system example
    double C = 2.96;
    double mu = 1.21505714306e-2;
    double lstar = 384388.174; //km
    double R1 = 6378.14 ; //km
    double R2 = 1738.20 ; //km
    double sing_tol = 1e-20;
    nvis::bbox2 bounds;
    bounds.min()[0] = -0.4;
    bounds.min()[1] = -2.5;
    bounds.max()[0] = 1.1;
    bounds.max()[1] = 2.5;
    nvis::ivec2 res(10, 10);


    //Map Params
    MapParams   theMapParams; //Assigns defaults
    theMapParams.bounds = bounds;
    theMapParams.resolution = res;
    theMapParams.nb_iterations = niters;
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
    bool paramsAreReset = false; //Indicates when options may differ than how AGNData was made.
    std::string mpFileStr = "none";
    //Metric bounds influence manifold computation
    const double LARGE = std::numeric_limits<double>::max();
    theMapParams.the_metric.bounds().min()[0] = -LARGE;
    theMapParams.the_metric.bounds().min()[1] = -2.55;
    theMapParams.the_metric.bounds().max()[0] = LARGE;
    theMapParams.the_metric.bounds().max()[1] = 2.55;
    theMapParams.the_metric.periodic()[0] = false; //Not periodic bounds
    theMapParams.the_metric.periodic()[1] = false; //Not periodic bounds



    /// Data Objects
    //Right Hand Side
    rhs_type rhs(C, mu);
    //Section
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    section.isPositive = positiveDir;

    //The Poincare Map Engine
    map_type theMap(rhs, section);
    theMap.setPrecision(eps);


    //Topology Extractor Class
    TopologyExtractorType* topoExtractor = new TopologyExtractorType(theMap, theMapParams);
    topoExtractor->fpFileStr = "FixedPointData.fpx"; //Default, always writes this file
    bool runNameSpecified = false;
    std::string agndFileStr = "none"; //For use if runNameSpecified
    bool readAGND = false;


    //Max number of threads available
    size_t maxNumThreads = 1;
#if _OPENMP
    maxNumThreads = omp_get_max_threads();
#endif
    size_t nthreads = maxNumThreads;

    std::cout << "received command line arguments\n";
    for (int i=0; i<argc; ++i) {
        std::cout << argv[i] << " ";
    }
    std::cout << '\n';


    /// Read & Handle Options: --------------------------------------------------
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("You requested the help output");
        }
        // PMATE: Problem Definition ----------------------------------------------------
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
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing bounds");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            bounds_set = true;
            theMapParams.bounds = bounds;
            theMapParams.the_metric.bounds() = bounds;
            paramsAreReset = true;
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
            theMapParams.resolution = res;
            paramsAreReset = true;
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
        // PMATE: Input Files ------------------------------------------------------------
        else if (arg == "-ip" || arg == "--params") {
            if (i == argc-1) {
                printUsageAndExit("missing params name");
            }
            mpFileStr = argv[++i];
            topoExtractor->mpFileStr = mpFileStr;
            std::cerr << "PMATE: Reading map parameters file...\n";
            bool readOK = theMapParams.read(mpFileStr.c_str());
            if (!readOK) {
                throw std::runtime_error( "Param file read failure!" );
            }
            //Setting some parameters based on MapAnalysisParam input
            niters = theMapParams.nb_iterations;
            bounds = theMapParams.bounds;
            xavier::OnEdgeCompare::maxDepth = theMapParams.max_depth;
            xavier::ID_Compare::maxDepth = theMapParams.max_depth;
            writeParams = false;
            bounds_set = true;
            paramsAreReset = false;
        } else if (arg == "-ig" || arg == "--grid") {
            if (i == argc-1) {
                printUsageAndExit("missing grid data file name");
            }
            agndFileStr = argv[++i];
            topoExtractor->agndFileStr = agndFileStr;
            readAGND = true;
            writeGrid = false;
        } else if (arg == "-wp" || arg == "--write_params") {
            if (i == argc-1) {
                printUsageAndExit("missing write params option");
            }
            writeParams = (atoi(argv[++i])==1)? true : false;
        } else if (arg == "-wg" || arg == "--write_grid") {
            if (i == argc-1) {
                printUsageAndExit("missing write grid option");
            }
            writeGrid = (atoi(argv[++i])==1)? true : false;
        }
        // PMATE: Output -----------------------------------------------------------------
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output base name");
            }
            topoExtractor->composeFileNames( argv[++i] );
            mpFileStr = topoExtractor->mpFileStr;
            runNameSpecified = true;
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
        } else if (arg == "-fd" || arg == "--mapDiscontFile") {
            if (i == argc-1) {
                printUsageAndExit("missing map discont file name");
            }
            topoExtractor->mapDiscontFileStr = argv[++i];
        } else if (arg == "-fg" || arg == "--fpGuessesFile") {
            if (i == argc-1) {
                printUsageAndExit("missing fixed point guesses file name");
            }
            topoExtractor->fpGuessFileStr = argv[++i];
        } else if (arg == "-ff" || arg == "--fpFailFile") {
            if (i == argc-1) {
                printUsageAndExit("missing fixed point failures file name");
            }
            topoExtractor->fpFailFileStr = argv[++i];
        } else if (arg == "-fx" || arg == "--fpxFile") {
            if (i == argc-1) {
                printUsageAndExit("missing fixed point results file name");
            }
            topoExtractor->fpFileStr = argv[++i];
        } else if (arg == "-od" || arg == "--write_mapDiscont") {
            if (i == argc-1) {
                printUsageAndExit("missing write map discont option");
            }
            writeMapDiscont = (atoi(argv[++i])==1)? true : false;
        } else if (arg == "-og" || arg == "--write_guesses") {
            if (i == argc-1) {
                printUsageAndExit("missing write fixed point guesses option");
            }
            writeGuesses = (atoi(argv[++i])==1)? true : false;
        } else if (arg == "-of" || arg == "--write_failures") {
            if (i == argc-1) {
                printUsageAndExit("missing write fixed point failures option");
            }
            writeFailures = (atoi(argv[++i])==1)? true : false;
        } else if (arg == "-wc" || arg == "--watch_cell") {
            if (i == argc-3) {
                printUsageAndExit("missing watch cell id");
            }
            for (int j=0; j<3; j++) {
                topoExtractor->theWatchCell[j] = atoi(argv[++i]);
            }
            topoExtractor->enableWatchCell = true;
        }
        // PMATE: Sampling Parameters ----------------------------------------------------
        else if (arg == "-p" || arg == "--sample_iters") {
            if (i == argc-1) {
                printUsageAndExit("missing sample iterates");
            }
            niters = atoi(argv[++i]);
            if (theMapParams.nb_iterations != niters) {
                writeParams = true;
            }
            theMapParams.nb_iterations = niters;
            paramsAreReset = true;
        } else if (arg == "-md" || arg == "--maxDepth") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum depth value");
            }
            max_depth = atoi(argv[++i]);
            if(max_depth != theMapParams.max_depth) {
                writeParams = true;
            }
            theMapParams.max_depth = max_depth;
            xavier::OnEdgeCompare::maxDepth = max_depth;
            xavier::ID_Compare::maxDepth = max_depth;
            paramsAreReset = true;
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
            theMapParams.samplingIntegTol = eps;
            paramsAreReset = true;
        } else if (arg == "-wt" || arg == "--winding_tol") {
            if (i == argc-3) {
                printUsageAndExit("missing winding number tolerances");
            }
            for(int j=0; j<3; j++) {
                wTols[j] = atof(argv[++i]);
            }
            theMapParams.winding_convexity_tols = wTols;
            paramsAreReset = true;
        } else if (arg == "-wd" || arg == "--winding_dist") {
            if (i == argc-3) {
                printUsageAndExit("missing winding number value variance");
            }
            for(int j=0; j<3; j++) {
                wDist[j] = atof(argv[++i]);
            }
            theMapParams.winding_cell_maxDist = wDist;
            paramsAreReset = true;
        }
        // PMATE: Poincare Index Parameters------------------------------------------------
        else if (arg == "-mp" || arg == "--min_p") {
            if (i == argc-1) {
                printUsageAndExit("missing minimum period cutoff");
            }
            theMapParams.min_period = atoi(argv[++i]);
        } else if (arg == "-xp" || arg == "--max_p") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum period cutoff");
            }
            theMapParams.max_period = atoi(argv[++i]);
        } else if (arg == "-ma" || arg == "--max_angle") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum rotation angle");
            }
            theMapParams.max_angle = M_PI*atof(argv[++i])/180.0;
        } else if (arg == "-lm" || arg == "--lmin") {
            if (i == argc-1) {
                printUsageAndExit("missing min allowed distance lmin");
            }
            theMapParams.lmin = atof(argv[++i]);
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
        // PMATE: Fixed Point Refinement Parameters ---------------------------------------
        else if (arg == "-sc" || arg == "--subcell_res") {
            if (i == argc-2) {
                printUsageAndExit("missing subcell resolution");
            }
            theMapParams.subcellResolution[0] = atof(argv[++i]);
            theMapParams.subcellResolution[1] = atof(argv[++i]);
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
            printUsageAndExit("unrecognized argument: " + arg);
        }
    }
    if (!bounds_set)  {
        printUsageAndExit("Boundary not set");
    }
    else if (mpFileStr =="none") {
        printUsageAndExit("No map file name provided");
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
    topoExtractor->nbthreads = nthreads;

    //Adjust objects to handle new values
    rhs_type adjustedRHS(C,mu);
    adjustedRHS.setSingularityTolerance( sing_tol );
    adjustedRHS.setSingularitySafeDistance(0,R1/lstar);//Applies to filtering!
    adjustedRHS.setSingularitySafeDistance(1,R2/lstar);
    section.set_rhs( adjustedRHS );
    theMap.setNewComponents( adjustedRHS, section );
    //Start process with sampling integration tolerance
    theMap.setPrecision(theMapParams.samplingIntegTol);

    //Utilize the new parameters in the Sampler (read should go here!)
    topoExtractor->setMapParams(theMapParams); //Resets some values within "theSampler"
    topoExtractor->resetGrid();
    //Adjust Convexity Checker for CR3BP
    topoExtractor->theSampler.cellChecker.setBodyRadii(R1/lstar,R2/lstar); //Non-dim inputs!
    topoExtractor->theSampler.cellChecker.setMu(mu);
    topoExtractor->theSampler.cellChecker.tols = wTols;
    topoExtractor->theSampler.cellChecker.maxDists = wDist;

    //Perform the AGND read()
    if (readAGND) {
        if (paramsAreReset) {
            std::cerr << "Warning: Input map parameters are being overwritten with additional inputs!\n";
            std::cerr << "         It is likely that the specified Adaptive Grid Data input file \n";
            std::cerr << "         [ " << agndFileStr << " ]\n";
            std::cerr << "         does NOT correlate with generated data!\n";
        }
        std::cerr << "PMATE: Processing input grid data... \n";
        topoExtractor->theSampler.read( agndFileStr.c_str() );
    } else {
        // We have to sample, but that is handled as part of compute()
    }


    //The Main call:
    std::cerr << setprecision(14) << argv[0] << ": Executing PMATE process for the system :\n"
              << "---------------------------------------------------------------------------\n"
              << "      mu = " << mu << "\n"
              << "       C = " << C << "\n"
              << "    bbox = [" << bounds.min() << "->" << bounds.max() << "]\n"
              << "     res = " << res << "\n"
              << "---------------------------------------------------------------------------\n";

    nvis::timer timer;

    //Call the sampling function
    topoExtractor->compute();


    double elapsed = timer.elapsed();
    std::cerr << "-------------------------------------------------------------------------\n";
    std::cerr << "\nComputation took " << elapsed/3600.0 << " hrs = " << elapsed << " s. \n";
    std::cerr << "-------------------------------------------------------------------------\n";

    //Additional storage if specified
    if (writeParams && runNameSpecified) {
        topoExtractor->writeMapParams();
    } else {
        topoExtractor->mpFileStr = "MapParameters.param";
        topoExtractor->writeMapParams();
    }
    if (writeGrid && runNameSpecified) {
        topoExtractor->writeSamplingData();
    } else if (writeGrid) {
        topoExtractor->agndFileStr = "AdaptiveGridData.agnd";
        topoExtractor->writeSamplingData();
    }
    if (writeMapDiscont && (topoExtractor->mapDiscontFileStr != "none") ) {
        topoExtractor->writeMapDiscont();
    }
    if (writeGuesses && (topoExtractor->fpGuessFileStr != "none") ) {
        topoExtractor->writeFixedPointGuesses();
    }
    if (writeFailures && (topoExtractor->fpFailFileStr != "none") ) {
        topoExtractor->writeFailedFixedPoints();
    }


    return 0;
}
