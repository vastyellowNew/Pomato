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


/** Computing Poincare Map Data for Orbit Convolution
 *  - Compute a regular grid of Poincare map data in CR3BP
 *    and store output such that it's ready for
 *    image processing (NRRD format - 'int' array)
 *  - This is explicitly for the apse_section!
 *
 *  Author:  Wayne Schlei, Xavier Tricoche (Purdue University)
 */
#include <iostream>
#include <vector>
#include <list>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>
#include <limits>

// math
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
#include <maps/poincare_map.hpp>
#include <convolution/orbitConvolution.hpp>

// cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/apse_section.hpp>
#include <cr3bp/tracker.hpp>

//Outputs
#include <teem/nrrd.h>
#include <cr3bp/psiWrite.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif
//#include <mpi.h>

using namespace nvis;


typedef orbital::cr3bp                                    rhs_type;
typedef nvis::fixed_vector<double, 42>                                vec42;
typedef nvis::fixed_vector<double, 6>                                 vec6;
typedef nvis::fixed_matrix<double, 6>                                 mat6;
typedef xavier::dp5wrapper<double, 42>                                odesolver_type;
typedef orbital::apse_section<rhs_type, 6, 42>                        section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type>  map_type;
typedef map_type::return_type                                         return_type;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

// *******************************    PARAMETERS     *******************************
std::string me;
std::string filename;
std::string pointFile;
nvis::ivec2 res(1000, 1000);

void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Compute Map for OC in restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -b  | --bounds <float> x 4       Image bounds\n"
            << " -r  | --resolution <int> x 2     Image resolution\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -s  | --section <int>            Section Type (0:Periapsis, 1:Apoapsis)\n"
            << " -d  | --seedDir <int>            Direction (0:Prograde, 1:Retrograde)\n"
            << " -rb | --refBody <int>            Reference Body (0:Barycenter,1:P1,2:P2)\n"
            << " -p  | --maxp <int>               Max considered period\n"
            << " -t  | --threads <int>            Set number of executed threads\n"
            << " -o  | --output <string>          Output base name\n"
            //<< " -op | --pointFile <string>       Point-clout output (.psi, optional)\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    // RHS and ODE solver parameters
    double eps, C, mu;
    
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    
    //Standard Earth-Moon system example
    C = 2.96;
    mu = 1.21505714306e-2;
    nvis::bbox2 bounds;
    bounds.min()[0] = 0.8;
    bounds.min()[1] = -0.2;
    bounds.max()[0] = 1.2;
    bounds.max()[1] = 0.2;
    
    //Apse section specifics
    int seedDir = 0; //Prograde
    int sectionType = 0; //Periapsis
    int refBody = 2;
    
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    int maxp = 100;
    filename = "none";
    pointFile = "none";
    int nseeds = 100;
    size_t maxNumThreads = 1;
#if _OPENMP
    maxNumThreads = omp_get_max_threads();
#endif
    size_t nthreads = maxNumThreads;
    
    
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
        } else if (arg == "-t" || arg == "--threads") {
            if (i == argc-1) {
                printUsageAndExit("missing number of threads");
            }
            nthreads = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--section") {
            if (i == argc-1) {
                printUsageAndExit("missing section type");
            }
            sectionType = atoi(argv[++i]);
        } else if (arg == "-d" || arg == "--seedDir") {
            if (i == argc-1) {
                printUsageAndExit("missing seed direction");
            }
            seedDir = atoi(argv[++i]);
        }
        
        else if (arg == "-rb" || arg == "--refBody") {
            if (i == argc-1) {
                printUsageAndExit("missing reference body");
            }
            refBody = atoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            filename = argv[++i];
        }
        /*else if (arg == "-op" || arg == "--pointFile") {
            if (i== argc-1) printUsageAndExit("missing point file name");
            pointFile = argv[++i];
        }*/
        else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        }
        
        else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (filename=="none") {
        printUsageAndExit("");
    }
    
    if (nthreads > maxNumThreads) {
        std::cerr << "Warning:  Indicated thread number (" << nthreads
                  << ") is larger than max available (" << maxNumThreads
                  << ").  Setting thread number to " << maxNumThreads << "\n";
        nthreads = maxNumThreads;
    }
#if _OPENMP
    omp_set_num_threads( nthreads );
#endif
    
    //Get number of trajectories to process
    nseeds = res[0] * res[1];
    
    std::cerr << argv[0] << ": Computing Poincare Map for " << nseeds << " orbits in "
              << '[' << bounds.min() << "->" << bounds.max() << "]\n";
              
    const double LARGE = std::numeric_limits<double>::max();
    //CR3BP
    rhs_type rhs(C, mu);
    //Apse Section
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    section.isPeriapsis = (sectionType == 0) ? true : false;
    section.seedPrograde = (seedDir == 0) ? true : false;
    switch (refBody) {
        case 0 :
            section.referenceBody = section_type::BARYCENTER;
            break;
        case 1 :
            section.referenceBody = section_type::P1;
            break;
        case 2 :
            section.referenceBody = section_type::P2;
            break;
        default :
            section.referenceBody = section_type::P2;
            break;
    }
    //The Map
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    
    std::vector<nvis::vec2> seeds(nseeds);
    //Random
    //for (int i=0 ; i<nseeds ; ++i) {
    //  seeds[i] = bounds.min() + nvis::vec2(drand48(), drand48())*bounds.size();
    //}
    //Regular grid
    for (int k=0; k<nseeds; k++) {
        const int ix = k % res[0];
        const int iy = k / res[0];
        seeds[k] = bounds.interpolate( nvis::vec2((ix+0.5)/res[0], (iy+0.5)/res[1]) );
    }
    
    
    std::cerr << nthreads << " threads employed for computation\n";
    
    nvis::timer timer;
    //Cache of data per thread
    std::vector< std::vector<return_type> > orbitDataPerThread[nthreads];
    std::vector< int > seedIDPerThread[nthreads];
    
    int counter=0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<nseeds ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            map_type* theMap = pmap.clone();
            
            std::vector<return_type> orbit;
            //tracker_type tracker(barycenter);
            try {
                //theMap->map_complete(seeds[n], orbit, maxp);
                //In-region version of map that runs 5*maxp iters
                theMap->mapRegion(seeds[n], orbit, maxp, 5*maxp, bounds);
            } catch(...) {
                //Less than maxp is completed
                // Causes:  Singularity, loss of percision, or violates Jacobi Integral,
                // or insufficent returns to region
                // Result:  Just store anyway, if orbit is empty, that's ok here...
                
                //continue;
            }
            
            
            //Store up until available data
            orbitDataPerThread[thread_id].push_back( orbit );
            seedIDPerThread[thread_id].push_back( n );
            
#pragma openmp atomic
            ++counter;
            
            int c = counter;
            std::ostringstream os;
            double dt = timer.elapsed();
            os << "\r" << c << "/" << nseeds << " (" << (float)c*100/nseeds << "%) in "
               << dt << "s. (" << (float)c/dt << "Hz)  Roughly "
               << (((float)nseeds-(float)c)/((float)c/dt))/60.0
               << " min remaining...         \r";
            std::cout << os.str() << std::flush;
        }
    }
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. (" << (float)nseeds/elapsed << "Hz)\n";
    
    //Unpack output per thread
    typedef std::vector<return_type> DataValueType;
    typedef std::map<int,DataValueType> DataMapType;
    typedef std::pair<int,DataValueType> DataPairType;
    DataMapType dataMap;
    for(int tid=0; tid<nthreads; tid++) {
        for(int jj=0; jj<(int) seedIDPerThread[tid].size(); jj++) {
            //Insert each pair (IC_id,returns) into the map container
            dataMap.insert( DataPairType( seedIDPerThread[tid][jj], orbitDataPerThread[tid][jj]) );
        }
    }
    //Create the NRRD file
    OrbitConvolution::convertMapDataToNRRD< DataMapType, DataValueType >
    (res,maxp,bounds,dataMap,filename.c_str() );
    
    //Creating output (if done without convolution functions)
    /*size_t dims[] = {res[0], res[1]};
    Nrrd *nout = nrrdNew();
    nrrdWrap_nva(nout, raster, nrrdTypeFloat, 2, dims);
    nrrdSave(filename.c_str(), nout, NULL);
    nrrdNuke(nout);*/
    
    
    return 0;
}
