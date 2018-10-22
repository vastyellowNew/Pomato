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


#include <iostream>
#include <vector>
#include <list>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>
#include <limits>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>

// map API
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/newton.hpp>
#include <maps/fixpoints.hpp>
#include <maps/index.hpp>
#include <maps/poincare_map.hpp>
#include <topology/invariant_manifold.hpp>

// cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

// RHS and ODE solver parameters
double eps, C, mu, K;

typedef orbital::cr3bp                                                rhs_type;
typedef nvis::fixed_vector<double, 42>                                vec42;
typedef nvis::fixed_vector<double, 6>                                 vec6;
typedef nvis::fixed_matrix<double, 6>                                 mat6;
typedef xavier::dp5wrapper<double, 42>                                odesolver_type;
typedef orbital::planar_section<rhs_type, 6, 42>                      section_type;
typedef orbital::AngleTracker<42, 3>                                  tracker_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type>  map_type;
typedef map_type::return_type                                         return_type;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

// *******************************    PARAMETERS     *******************************
std::string me;
std::string filename;
nvis::vec2 L3(-1.00001053340710,0);
nvis::ivec2 res(1000, 1000);

void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Winding number computation in restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -b  | --bounds <float> x 4       Image bounds\n"
            << " -r  | --resolution <int> x 2     Image resolution\n"
            << " -n  | --nseeds <int>             Number of seeds\n"
            << " -x0 | --center <float> x2        Center of rotation\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --maxp <int>               Max considered period\n"
            << " -o  | --output <string>          Output base name\n"
            << " -s  | --seed <int>               Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal)\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    nvis::bbox2 bounds;
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    //std::vector<nvis::vec2> pCenters;
    //pCenters.push_back(L3);
    
    //Earth-Moon system example
    C = 2.96;
    mu = 1.21505714306e-2;
    bounds.min()[0] = -0.3;
    bounds.min()[1] = -0.5;
    bounds.max()[0] = 1.2;
    bounds.max()[1] = 0.5;
    //Possible Centers for Earth-Moon
    std::vector<nvis::vec2> pCenters;
    pCenters.push_back(nvis::vec2(0.83691519554062,0.0)); //L1
    pCenters.push_back(nvis::vec2(1.15568211091105,0.0)); //L2
    pCenters.push_back(nvis::vec2( -1.00506263990269,0.0)); //L3
    pCenters.push_back(nvis::vec2( 0.48784942856940,0.86602540378444)); //L4
    pCenters.push_back(nvis::vec2( 0.48784942856940,-0.86602540378444)); //L5
    pCenters.push_back(nvis::vec2(-mu,0.0)); //P1
    pCenters.push_back(nvis::vec2(1.0-mu,0.0)); //P2
    
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    bool bounds_set = false;
    int maxp = 100;
    filename = "none";
    int nseeds = 100;
    int seeding = 1;
    
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
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing bounds");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            bounds_set = true;
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
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
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            filename = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else if (arg == "-n" || arg == "--nseeds") {
            if (i == argc-1) {
                printUsageAndExit("missing seed number");
            }
            nseeds = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--seeding") {
            if (i == argc-1) {
                printUsageAndExit("missing seeding info");
            }
            seeding = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (!bounds_set || filename=="none") {
        printUsageAndExit("");
    }
    
    nvis::vec2 primary(-mu, 0);
    nvis::vec2 secondary(1-mu, 0);
    nvis::vec2 barycenter(0, 0);
    
    std::cerr << argv[0] << ": computing the winding number of " << nseeds << " orbits in "
              << '[' << bounds.min() << "->" << bounds.max() << "]\n";
              
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.precision() = eps;
    
    float* raster = (float*)calloc(4*res[0]*res[1], sizeof(float));
    float* usedCenter = (float*)calloc(4*res[0]*res[1], sizeof(float));
    double dx = bounds.size()[0] / (double)(nseeds-1);
    int N = res[0]*res[1];
    
    std::vector<nvis::vec2> seeds(nseeds);
    if (seeding == 0) {
        for (int i=0 ; i<nseeds ; ++i) {
            seeds[i] = bounds.min() + nvis::vec2(drand48(), drand48())*bounds.size();
        }
    } else if (seeding == 1) {
        double dx = bounds.size()[0]/(double)(nseeds-1);
        for (int i=0 ; i<nseeds ; ++i) {
            seeds[i][0] = bounds.min()[0] + i*dx;
            seeds[i][1] = 0.5*(bounds.min()[1] + bounds.max()[1]);
        }
    } else if (seeding == 2) {
        double dy = bounds.size()[1]/(double)(nseeds-1);
        for (int i=0 ; i<nseeds ; ++i) {
            seeds[i][0] = 0.5*(bounds.min()[0] + bounds.max()[0]);
            seeds[i][1] = bounds.min()[1] + i*dy;
        }
    } else if (seeding == 3) {
        nvis::vec2 dp = bounds.size() / nvis::vec2(nseeds-1);
        for (int i=0 ; i<nseeds ; ++i) {
            seeds[i] = bounds.min() + i*dp;
        }
    }
    
    size_t nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cerr << nthreads << " threads available for computation\n";
    
    nvis::timer timer;
    int counter=0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<nseeds ; ++n) {
            std::vector<return_type> orbit;
            //Select the tracking center
            //  1)  Locate nearest possible center
            double minDist = 200.0;
            int minID = 0;
            for (int ii=0; ii< (int) pCenters.size(); ii++) {
                //Get seed position (x,y)
                nvis::vec2 pos(seeds[n][0],0.0); //On y=0 plane
                double dist = nvis::norm<double,2>(pos - pCenters[ii]);
                if (dist<minDist) {
                    minID = ii;
                    minDist = dist;
                }
            }
            tracker_type tracker(pCenters[minID]); //JE - L3
            try {
                pmap.map_and_track_complete<tracker_type>(seeds[n], orbit, maxp, tracker);
            } catch(...) {
                continue;
            }
            if (orbit.size() < maxp) {
                continue;
            }
            
            double w = 2.*M_PI*orbit.size()/orbit.back().delta_theta[0];
            if (std::isnan(w) || std::isinf(w)) {
                w = 1000;
            }
            
            orbit.push_back(return_type());
            orbit.back().x = seeds[n];
            for (int i=0 ; i<orbit.size() ; ++i) {
                nvis::vec2 x = orbit[i].x;
                nvis::vec2 c = (x-bounds.min())/bounds.size()*nvis::vec2(res);
                if (c[0]<0 || c[0]>=res[0] || c[1]<0 || c[1]>=res[1]) {
                    continue;
                }
                int m = (int)c[0] + (int)c[1]*res[0];
#pragma openmp atomic
                raster[m] = w;
#pragma openmp atomic
                usedCenter[m] = (float) minID + 1.0;
            }
            
#pragma openmp atomic
            ++counter;
            
            int c = counter;
            std::ostringstream os;
            double dt = timer.elapsed();
            os << "\r" << c << "/" << nseeds << " (" << (float)c*100/nseeds << "%) in "
               << dt << "s. (" << (float)c/dt << "Hz)           \r";
            std::cout << os.str() << std::flush;
        }
    }
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. (" << (float)nseeds/elapsed << "Hz)\n";
    
    size_t dims[] = {res[0], res[1]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, raster, nrrdTypeFloat, 2, dims);
    char* labels[] = { (char*)"X", (char*)"Xdot" };
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
    nrrdSave(filename.c_str(), nout, NULL);
    nrrdNuke(nout);
    
    //Used Center field
    Nrrd* nout2 = nrrdNew();
    nrrdWrap_nva(nout2, usedCenter, nrrdTypeFloat, 2, dims);
    nrrdAxisInfoSet_nva(nout2, nrrdAxisInfoLabel, labels);
    nrrdSave("usedCenters.nrrd", nout2, NULL);
    nrrdNuke(nout2);
    
    return 0;
}
