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


//Testing the two versions of mapDisplacmentUsingCache() side-by-side
// Author: Wayne Schlei
#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>
#include <limits>

// math
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <boost/rational.hpp>
#include <math/rational.hpp>

// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>

// map API
#include <maps/DP45wrapper.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/fixpoints.hpp>
#include <maps/index.hpp>
#include <maps/poincare_map.hpp>
#include <topology/CATtracker.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/EdgeRotationMapCalls.hpp>
#include <topology/AdaptiveEdgeRotation.hpp>

// cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/multipleAngleTracker.hpp>

//Outputs
#include <teem/nrrd.h>
#include <cr3bp/psiWrite.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

// RHS and ODE solver parameters
double eps, C, mu;

typedef orbital::cr3bp                                    rhs_type;
typedef nvis::fixed_vector<double, 42>                                vec42;
typedef nvis::fixed_vector<double, 6>                                 vec6;
typedef nvis::fixed_matrix<double, 6>                                 mat6;
typedef xavier::dp5wrapper<double, 42>                                odesolver_type;
typedef orbital::planar_section<rhs_type, 6, 42>                      section_type;
//typedef orbital::AngleTracker<42, 3>                                  tracker_type;
typedef orbital::MultipleAngleTracker<42, 3>                          tracker_type;
typedef topology::CATtracker<vec42,rhs_type::numSingularities>        TrackerType;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type>  map_type;
typedef map_type::return_type                                         return_type;
typedef xavier::EdgeRotationFailure<nvis::vec2>                     ERotFailure;
static const int s = rhs_type::numSingularities;
typedef nvis::fixed_vector<double,s+1>                                ExtendedMapDataVec;
typedef topology::SortableReturnData<vec2,ExtendedMapDataVec>         SortableData;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

// *******************************    PARAMETERS     *******************************

void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << " [options]\n"
            << "DESCRIPTION: Winding number computation in restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -x0 | --state <float> x2         Testing coordinate on section\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --maxp <int>               Number of map iterations (forward and backward)\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    
    //Earth-Moon system example
    C = 2.96;
    mu = 1.21505714306e-2;
    
    //Initial state for testing
    vec2 x0(0.85375, -0.0220039); //Rotation failure at p=3 (map is good)
    //vec2 x0(0.885641,-0.2085); //Map Fail at beta(p=1)?
    //vec2 x0(0.8759,0.127); //Map Fail at beta(p=2)?
    //vec2 x0(0.875662,0.11175); //Map Fail at beta(p=2)?
    
    
    eps = 1.0e-8;
    int verbose = 0;
    int maxp = 3;
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-C") {
            if (i == argc-1) {
                printUsageAndExit("missing C constant");
            }
            C = atof(argv[++i]);
        } else if (arg == "-m" || arg == "--mu") {
            if (i == argc-1) {
                printUsageAndExit("missing mu constant");
            }
            mu = atof(argv[++i]);
        } else if (arg == "-x0" || arg == "--state") {
            if (i >= argc-2) {
                printUsageAndExit("missing initial state");
            }
            x0[0] = atof(argv[++i]);
            x0[1] = atof(argv[++i]);
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-p" || arg == "--maxp") {
            if (i == argc-1) {
                printUsageAndExit("missing max period");
            }
            maxp = atoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.precision() = eps;
    xavier::map_analysis_param theMapParams;
    theMapParams.min_period = 1;
    theMapParams.max_period = maxp;
    //Default settings for rotation evaluation parameters
    theMapParams.max_angle = M_PI*(135.0/180.0);
    theMapParams.maxDeltaTransSpeed = 0.8;
    theMapParams.maxDeltaMapDisplacement = nvis::vec2(0.2,2.0);
    theMapParams.maxDeltaTimeFactor = 0.25;
    theMapParams.ntEdgeDivisions = 256;
    theMapParams.min_fp_tol = 1.e-6;
    //lmin - don't need for this right now
    //theMapParams.lmin = 0.99*std::min(spacing[0],spacing[1])/256.;
    
    size_t nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cerr << nthreads << " threads available for computation\n";
    
    nvis::timer timer;
    
    //Propagate the map using mapDisplacementUsingCache();
    std::set<SortableData> dCache,bCache,newdCache,newbCache;
    ERotFailure pf;
    
    for (int p=maxp; p>0; --p) {
    
        //Start the Data output
        std::cerr << "===========================================================================\n";
        std::cerr << " Period : " << p << "\n";
        std::cerr << "===========================================================================\n";
        //Original version
        nvis::vec2 beta0, delta0, eta0, state;
        bool originalMapError = false, originalFoundData = true;
        std::set<SortableData>::iterator cit;
        std::cerr << "  Original mapDisplacementUsingCache():\n";
        try {
            beta0 = mapDisplacementUsingCache(x0,pmap,-p,theMapParams,bCache);
            //Look up the state
            cit = bCache.find(SortableData(x0));
            if(cit == bCache.end()) {
                originalFoundData = false;
                throw std::runtime_error("Backward data not found");
            }
            state = cit->returns[p-1]; //Could this be a segfault without being a std::runtime_error?
            std::cerr << "     Iter " << -p << " : " << state << "\n";
            std::cerr << "     beta   = " << beta0 << "\n";
        } catch(...) {
            //Report backward fail
            std::cerr << "     Iter " << -p << " : " << "Failed" << "\n";
            if(!originalFoundData) {
                std::cerr << "     beta   = " << beta0 << "\n";
            } else {
                std::cerr << "     beta   = " << "Failed" << "\n";
            }
            originalMapError = true;
        }
        originalFoundData = true;
        try {
            delta0 = mapDisplacementUsingCache(x0,pmap,p,theMapParams,dCache);
            //Look up the state
            cit = dCache.find(SortableData(x0));
            if(cit == dCache.end()) {
                originalFoundData = false;
                throw std::runtime_error("Forward data not found");
            }
            state = cit->returns[p-1];
            std::cerr << "     Iter " << p << " : " << state << "\n";
            std::cerr << "     delta  = " << delta0 << "\n";
        } catch(...) {
            //Report forward fail
            std::cerr << "     Iter " << p << " : " << "Failed" << "\n";
            if(!originalFoundData) {
                std::cerr << "     delta  = " << delta0 << "\n";
            } else {
                std::cerr << "     delta  = " << "Failed" << "\n";
            }
            originalMapError = true;
        }
        if(originalMapError) {
            std::cerr << "     eta    = " << "Failed" << "\n";
        } else {
            eta0 = (delta0-beta0)/2.0;
            std::cerr << "     eta    = " << eta0 << "\n";
        }
        
        //Modified version
        std::cerr << "---------------------------------------------------------------------------\n";
        std::cerr << "  New mapDisplacementUsingCache() with Interpolation:\n";
        nvis::vec2 beta1,delta1,eta1;
        nvis::vec2 edgeStep(0.0);
        bool mapError = false, dataFound = true;
        try {
            beta1 = mapDisplacementUsingCache(x0,pmap,-p,theMapParams,newbCache,pf,edgeStep);
            //Look up state
            cit = newbCache.find(SortableData(x0));
            if(cit == newbCache.end()) {
                dataFound = false;
                throw std::runtime_error("Backward data not found in new method");
            }
            state = cit->returns[p-1]; //segfault not runtime_error?
            std::cerr << "     Iter " << -p << " : " << state << "\n";
            std::cerr << "     beta   = " << beta1 << "\n";
            //Fixed with interpolation?
        } catch(ERotFailure& erFail) {
            //Report backward fail
            std::cerr << "     ERotFailure: " << erFail.what() << "\n";
            std::cerr << "     Iter " << -p << " : " << "Failed" << "\n";
            if(!dataFound) {
                std::cerr << "     beta   = " << beta1 << "\n";
            } else {
                std::cerr << "     beta   = " << "Failed" << "\n";
            }
            mapError = true;
        } catch(...) {
            std::cerr << "     Error: Not sure what went wrong on BACKWARD step\n";
            std::cerr << "     Iter " << -p << " : " << "Failed" << "\n";
            if(!dataFound) {
                std::cerr << "     beta   = " << beta1 << "\n";
            } else {
                std::cerr << "     beta   = " << "Failed" << "\n";
            }
            mapError = true;
        }
        try {
            delta1 = mapDisplacementUsingCache(x0,pmap,p,theMapParams,newdCache,pf,edgeStep);
            //Look up state
            cit = newdCache.find(SortableData(x0));
            if(cit == newdCache.end()) {
                dataFound = false;
                throw std::runtime_error("Forward data not found in new method");
            }
            state = cit->returns[p-1]; //segfault not runtime_error?
            std::cerr << "     Iter " << p << " : " << state << "\n";
            std::cerr << "     delta   = " << delta1 << "\n";
            //Fixed with interpolation?
        } catch(ERotFailure& erFail) {
            //Report backward fail
            std::cerr << "     ERotFailure: " << erFail.what() << "\n";
            std::cerr << "     Iter " << p << " : " << "Failed" << "\n";
            if(!dataFound) {
                std::cerr << "     delta  = " << delta1 << "\n";
            } else {
                std::cerr << "     delta  = " << "Failed" << "\n";
            }
            mapError = true;
        } catch(...) {
            std::cerr << "     Error: Not sure what went wrong on BACKWARD step\n";
            std::cerr << "     Iter " << p << " : " << "Failed" << "\n";
            if(!dataFound) {
                std::cerr << "     delta  = " << delta1 << "\n";
            } else {
                std::cerr << "     delta  = " << "Failed" << "\n";
            }
            mapError = true;
        }
        if(mapError) {
            std::cerr << "     eta    = " << "Failed" << "\n";
        } else {
            eta1 = (delta1-beta1)/2.0;
            std::cerr << "     eta    = " << eta1 << "\n";
        }
        
    }
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. \n";
    
    return 0;
}
