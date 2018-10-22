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


// Test function for MultipleShooting.hpp
// Wayne Schlei

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
#include <maps/section.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/fixpoints.hpp>
#include <maps/poincare_map.hpp>
#include <orbital/monodromy.hpp>
#include <orbital/controller.hpp>
#include <orbital/corrections.hpp>

// cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;


typedef orbital::cr3bp                          rhs_type;
typedef nvis::fixed_vector<double, 42>          vec42;
typedef nvis::fixed_vector<double, 6>           vec6;
typedef nvis::fixed_matrix<double, 6>           mat6;
typedef xavier::dp5wrapper<double, 42>          odesolver_type;
typedef orbital::planar_section<rhs_type, 6, 42>                     section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type> map_type;
typedef std::vector<xavier::fixpoint>           FixedPointChain;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;


//bool xavier::record_newton_steps = true;
//std::vector<nvis::vec2> xavier::newton_steps;

// *******************************    PARAMETERS     *******************************
std::string me("");
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Multiple-Shooting method for fixed point search (test) in circular restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -c  | --convergence <float>      Convergence tolerance\n"
            << " -g  | --guess <float> x 2        First guess\n"
            << " -th | --parallel                 Run a parallel test with multiple fixed points\n"
            << " -n  | --nCross <int>             Number of crossings for new guess [default = 1]\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --patchPoints <int>        Points per revolution (init)\n"
            << " -o  | --output <string>          Output file for regular samples\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    // RHS and ODE solver parameters
    
    //Jupiter-Europa system example
    double mu = 2.528017705e-5;
    double C = 3.000;
    //nvis::vec2 x0(-1.205,0); //Assumed IC for testing
    //nvis::vec2 x0(-1.23142, 0); //Saddle with p=4, ||F(X_0)|| = 1.e-7
    //int numCross = 4;
    
    //Earth-Moon Problem cases
    mu = 1.21505714306e-2;
    C = 2.96;
    //Convergence with SINGLE_SHOOT
    nvis::vec2 x0(0.7282608,0.0);  //Saddle p=1, L1 Lyap
    int numCross = 1;
    //Convergence with DISTRIBUTED_ERROR (mixed sampling)
    //nvis::vec2 x0(0.855334,0.0217615);
    //int numCross = 3;
    //Convergence (on wrong orbit) with PERIODICITY_HOMOTOPY
    // but convergence with QUASI_NEWTON with single shooting! :)
    //nvis::vec2 x0(0.302889,0.574074);
    //int numCross = 8;
    
    //Tough cases
    //nvis::vec2 x0(0.718261,0.0); //Bug with HamSecOff & Quasi-Newton? =>FIXED!
    //int numCross = 1;
    //nvis::vec2 x0(0.831884,0.0185185); //Saddle p=3, Period3 DRO
    //int numCross = 3;
    //nvis::vec2 x0(-0.284058,1.38889); //Center p=8
    //int numCross = 8;
    //nvis::vec2 x0(0.288406,0.3888889);
    //int numCross = 9;
    //nvis::vec2 x0(0.36087,0.648148);
    //int numCross = 7;
    
    bool runParallel = false;
    
    
    
    me = argv[0];
    double eps = 1.0e-12;
    double convTol = 1.0e-10;
    int verbose = 1;
    bool seed_set = false;
    int pointsPerRev = 5;
    std::string filename("none");
    bool file_set = false;
    
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
        } else if (arg == "-g" || arg == "--guess") {
            if (i >= argc-2) {
                printUsageAndExit("missing seed");
            }
            x0[0] = atof(argv[++i]);
            x0[1] = atof(argv[++i]);
            seed_set = true;
        } else if (arg == "-th" || arg == "--parallel") {
            runParallel = true;
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-c" || arg == "--convergence") {
            if (i == argc-1) {
                printUsageAndExit("missing convergence tolerance");
            }
            convTol = atof(argv[++i]);
        } else if (arg == "-n" || arg == "--nCross") {
            if (i == argc-1) {
                printUsageAndExit("missing number of crossings");
            }
            numCross = atoi(argv[++i]);
        } else if (arg == "-p" || arg == "--patchPoints") {
            if (i == argc-1) {
                printUsageAndExit("missing number of patch points per rev");
            }
            pointsPerRev = atoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output file");
            }
            filename = argv[++i];
            file_set = true;
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    //if (!seed_set) printUsageAndExit("");
    
    const double LARGE = std::numeric_limits<double>::max();
    if(runParallel) {
        //Force Earth-Moon problem
        mu = 1.21505714306e-2;
        C = 2.96;
        std::cout << "Multiple Shooting Parallel Test : Forcing problem to Earth-Moon examples.\n";
    }
    
    rhs_type rhs(C, mu);
    rhs.setSingularityTolerance( 1e-20 );
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    metric<double, 2> euclidean_metric;
    
    nvis::bbox2 bounds;
    bounds.min() = x0 - nvis::vec2(0.01, 0.01);
    bounds.max() = x0 + nvis::vec2(0.01, 0.01);
    
    //Sample the close proximity around IC
    if (file_set) {
        double* data = (double*)calloc(2*200*200, sizeof(double));
        nvis::vec2 dx = bounds.size() / nvis::vec2(200, 200);
#pragma openmp parallel for
        for (int n=0 ; n<200*200 ; ++n) {
            int i = n%200;
            int j = n/200;
            nvis::vec2 x = bounds.min() + nvis::vec2(i, j)*dx;
            nvis::vec2 f = pmap.map(x, 1) - x;
            data[2*n] = f[0];
            data[2*n+1] = f[1];
            if (!(n%10)) {
                std::ostringstream os;
                os << "\rn=" << n;
                std::cerr << os.str();
            }
        }
        
        size_t s[] = {2, 200, 200};
        double spc[] = {airNaN(), dx[0], dx[1]};
        Nrrd* nout = nrrdNew();
        nrrdWrap_nva(nout, data, nrrdTypeDouble, 3, s);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        nrrdSave(filename.c_str(), nout, NULL);
        delete[] data;
    }
    
    // compute orbit at seed point
    fixpoint fp;
    FixedPointChain fps;
    /*std::vector<nvis::vec2> chain;
    pmap.map(x0, chain, 10*maxp);
    std::cerr << "orbit at " << x0 << ":\n";
    std::copy(chain.begin(), chain.end(), std::ostream_iterator<nvis::vec2>(std::cerr, ", "));
    std::cerr << '\n';
    */
    
    if (!runParallel) {
        //Testing a fixed period
        int p = numCross;
        std::cerr << "Starting Convergence Test Sequence at " << x0 << " for period " << p << "\n";
        bool found = orbital::solveFixedPoints(pmap, euclidean_metric, bounds, x0,
                                               C, 5, p, fp, fps, convTol, true, 20, pointsPerRev, false,
                                               orbital::CorrectionsRegulator::THE_END);
                                               
        if (found) {
            std::cerr << " Resulting Fixed-Point: x = " << fp.pos << " with p = "
                      << fp.K << " and period = " << fp.timePeriod << "\n";
        }
        
        //Test Linear vs Numerical Monodromy approaches
        
    } else {
        // Build some guesses - enough to run on multiple threads (and a bit off to run corrections)
        FixedPointChain guessVec;
        xavier::fixpoint fpGuess;
        fpGuess.pos = nvis::vec2(1.00575725, -0.000000000);
        fpGuess.K = 1;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.81059326, 0.20651224);
        fpGuess.K = 2;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.13224032, 1.54022684);
        fpGuess.K = 3;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.78940592, 0.38076369);
        fpGuess.K = 3;
        guessVec.push_back( fpGuess );
        /*fpGuess.pos = nvis::vec2(0.83956114, -0.27181113);  fpGuess.K = 3;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.45538704, 0.00000000);  fpGuess.K = 3;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.80309362, 0.219482465);  fpGuess.K = 5;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.83844889, -0.18776967);  fpGuess.K = 4;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.80632322, -0.18312229);  fpGuess.K = 4;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.47768620, -0.14861833);  fpGuess.K = 2;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.82150715, -0.15578604);  fpGuess.K = 6;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.62068691, -0.22743231);  fpGuess.K = 7;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.62068691, 0.227432317);  fpGuess.K = 7;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.07331734, 0.000000000);  fpGuess.K = 9;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.77989647, -0.25454964);  fpGuess.K = 9;
        guessVec.push_back( fpGuess );
        fpGuess.pos = nvis::vec2(0.51309213, 0.074916925);  fpGuess.K = 10;
        guessVec.push_back( fpGuess );*/
        int numGuesses = (int) guessVec.size();
        
        size_t nbthreads = 1;
#if _OPENMP
        nbthreads = omp_get_max_threads();
#endif
        
        //Call Multi-Tiered Targeting scheme within each thread to solve for a fixed point
        std::vector<FixedPointChain>* cached_fps = new std::vector<FixedPointChain>[nbthreads]; //Note: store all fixed points to same per thread caches
        FixedPointChain* cached_failed = new FixedPointChain[nbthreads];
        //Temp storage
        FixedPointChain  fixed_points, failed_guesses;
        std::vector<FixedPointChain>  fp_chains;
        
        int counter = 0, foundCount=0;
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for(int n=0; n<numGuesses; ++n) {
                //map_type *amap = pmap.clone(); //new allocator per thread
                int thread_id = 0;
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                //Computation inputs
                nvis::vec2 theGuess = guessVec[n].pos;
                int pp = guessVec[n].K;
                xavier::fixpoint fp_result;
                FixedPointChain resultChain;
                //Parameters in solveFixedPoints():
                nvis::bbox2 bounds2; //Stops points from going too far away from guess
                bounds2.min() = theGuess - nvis::vec2(0.1,0.1);
                bounds2.max() = theGuess + nvis::vec2(0.1,0.1); //Allow this to enter adjacent cells
                //Run multiple shooting with Periodicity, Jacobi constant, and Section constraints
                bool found = orbital::solveFixedPoints<map_type>
                             ( pmap, euclidean_metric, bounds2,
                               theGuess, pmap.rhs().desired_hamiltonian(), 5, pp,
                               fp_result, resultChain, convTol, false, 20, pointsPerRev);
                               
                //Store the result if converged or failed
                if (found) {
                    cached_fps[thread_id].push_back( resultChain );
                    std::cout << "TestMS:  Corrections found the " << (resultChain[0].saddle?"SADDLE ":"CENTER ")
                              << resultChain[0].pos << " for p=" << pp << "\n";
                    #pragma omp atomic
                    foundCount++;
                } else {
                    cached_failed[thread_id].push_back( guessVec[n] );
                }
                
                #pragma omp atomic
                counter++;
                
                if (verbose) {
                    std::ostringstream os("");
                    os << "\rCompleted Corrections process for " << counter << " / " << numGuesses << " ("
                       << 100.*(float)counter/(float)numGuesses << "%), "
                       << foundCount << " fixed points found.     \r" << std::flush;
                    std::cout << os.str();
                }
                
                //Done with map copy
                //delete amap;
                
            }//End parallel for
        }//End parallel statement
        
        //Reset Data
        fixed_points.clear();
        fp_chains.clear();
        failed_guesses.clear();
        //Extract Data per thread
        for (int i=0; i<nbthreads; i++) {
            //Assemble the found fixed points
            for (int j=0; j<(int)cached_fps[i].size(); j++) {
                //Store all fps to external cache
                FixedPointChain& thisChain = cached_fps[i][j];
                for (int chainID =0; chainID < (int) thisChain.size(); ++chainID) {
                    fixed_points.push_back( thisChain[chainID] );
                }
                fp_chains.push_back( thisChain );
            }
            //Assemble the failed fixed points
            for (int k=0; k<(int)cached_failed[i].size(); k++) {
                failed_guesses.push_back( cached_failed[i][k] );
            }
        }
        //Clear thread storage
        for (int i=0; i<nbthreads; i++) {
            cached_fps[i].clear();
            cached_failed[i].clear();
        }
        delete[] cached_fps;
        delete[] cached_failed;
        std::cout << "PMATE : Corrections process found " <<  (int)fp_chains.size() << " orbits.\n";
        
    }
    
    
    return 0;
    
}
