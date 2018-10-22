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


// Test function for Task Scheduling
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

// boost API
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// map API
#include <maps/section.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/map_analysis_param.hpp>
#include <maps/fixpoints.hpp>
#include <maps/poincare_map.hpp>
#include <orbital/monodromy.hpp>
#include <orbital/controller.hpp>
#include <orbital/corrections.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/EdgeRotationMapCalls.hpp>

// cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>

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
typedef xavier::map_analysis_param              MapParams;
typedef std::vector<xavier::fixpoint>           FixedPointChain;

const double invalid_double = std::numeric_limits<double>::max();

//MaxDepth comparison
int xavier::OnEdgeCompare::maxDepth = 1;
int xavier::ID_Compare::maxDepth = 1;


using namespace xavier;


// *******************************    PARAMETERS     *******************************
std::string me("");
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Testing OpenMP 'task' propagation\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    // RHS and ODE solver parameters
    
    //Earth-Moon Problem cases
    double mu = 1.21505714306e-2;
    double C = 2.96;
    //Convergence with SINGLE_SHOOT
    nvis::vec2 x0(0.7282608,0.0);  //Saddle p=1, L1 Lyap
    int numCross = 1;
    
    me = argv[0];
    double eps = 1.0e-12;
    double convTol = 1.0e-10;
    int verbose = 1;
    bool seed_set = false;
    int pointsPerRev = 5;
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
        }
        
        else if (arg == "-v" || arg == "--verbose") {
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
    
    rhs_type rhs(C, mu);
    rhs.setSingularityTolerance( 1e-20 );
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    metric<double, 2> euclidean_metric;
    nvis::bbox2 bounds;
    bounds.min()[0] = -0.4;
    bounds.min()[1] = -2.5;
    bounds.max()[0] = 1.1;
    bounds.max()[1] = 2.5;
    
    //Map parameters
    MapParams   theMapParams; //Assigns defaults
    theMapParams.bounds = bounds;
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
    theMapParams.the_metric.bounds() = bounds;
    theMapParams.the_metric.periodic()[0] = false; //Not periodic bounds
    theMapParams.the_metric.periodic()[1] = false; //Not periodic bounds
    
    
    //Type definitions
    static const int s = rhs_type::numSingularities;
    typedef nvis::fixed_vector<double,s+1>                                ExtendedMapDataVec;
    typedef topology::SortableReturnData<nvis::vec2,ExtendedMapDataVec>   SortableData;
    typedef xavier::EdgeRotationFailure<nvis::vec2>                     MapDiscont;
    typedef topology::MapTaskOutputData<boost::uuids::uuid,SortableData,MapDiscont> TaskStorage;
    
    //Create storage objects
    std::set<SortableData> cache;
    
    
    // Build some guesses - enough to run on multiple threads (and a bit off to run corrections)
    FixedPointChain guessVec;
    std::vector<boost::uuids::uuid> ids;
    boost::uuids::random_generator  gen;
    xavier::fixpoint fpGuess;
    fpGuess.pos = nvis::vec2(1.00575725, -0.000000000);
    fpGuess.K = 1;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.81059326, 0.20651224);
    fpGuess.K = 2;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.13224032, 1.54022684);
    fpGuess.K = 3;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.78940592, 0.38076369);
    fpGuess.K = 3;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.83956114, -0.27181113);
    fpGuess.K = 3;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.45538704, 0.00000000);
    fpGuess.K = 3;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.80309362, 0.219482465);
    fpGuess.K = 5;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.83844889, -0.18776967);
    fpGuess.K = 4;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.80632322, -0.18312229);
    fpGuess.K = 4;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.47768620, -0.14861833);
    fpGuess.K = 2;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.82150715, -0.15578604);
    fpGuess.K = 6;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.62068691, -0.22743231);
    fpGuess.K = 7;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.62068691, 0.227432317);
    fpGuess.K = 7;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.07331734, 0.000000000);
    fpGuess.K = 9;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.77989647, -0.25454964);
    fpGuess.K = 9;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    fpGuess.pos = nvis::vec2(0.51309213, 0.074916925);
    fpGuess.K = 10;
    guessVec.push_back( fpGuess );
    ids.push_back( gen() );
    int numGuesses = (int) guessVec.size();
    
    size_t nbthreads = 1;
#if _OPENMP
    //Disable dynamic scheduling
    omp_set_dynamic(0);
    nbthreads = omp_get_max_threads();
    omp_set_num_threads(nbthreads);
#endif
    
    //Per-thread cache of data
    std::vector<TaskStorage>* perThreadCache = new std::vector<TaskStorage>[nbthreads];
    //Result stored/displayed at the end
    std::vector<double>  displacementMagnitude;
    FixedPointChain      outputFPs;
    
    int counter = 0, foundCount=0;
    bool propDone = false;
    //OpenMP Lock to limit thread access where necessary
#ifdef _OPENMP
    omp_lock_t writeLock;
    omp_init_lock(&writeLock);
#endif
    #pragma omp parallel
    {
        // Job (task) construction on a single thread, but work is proccessed in a pool
        #pragma omp single
        {
            //Extra loop to simulate a manifold
            for(int m=0; m<4; m++) {
            
                //The Propagation loop
                for(int n=0; n<numGuesses; ++n) {
                    //for(int m=0; m<4; m++) {
                    //#pragma omp depend(out:propDone) task
                    #pragma omp task
                    {
#ifdef _OPENMP
                        int thread_id = omp_get_thread_num();
#else
                        int thread_id = 0;
#endif
                        #pragma omp critical
                        std::cout << "(" << thread_id << ") : Propagating...\n";
                        vec_type pm = mapUsingTaskCache(guessVec[n].pos, pmap, guessVec[n].K,
                                                        theMapParams, cache, ids[n],
                                                        perThreadCache[thread_id]);
                        propDone = true;
                        
                        //Indicate done with this step
                        #pragma omp atomic
                        counter++;
                        
                        //Informing user about propagation
                        /*if (verbose) {
                            std::ostringstream os("");
                            os << "\rCompleted propagation for " << counter << " / " << numGuesses << " ("
                            << 100.*(float)counter/(float)numGuesses << "%), "
                            << ".     \r" << std::flush;
                            std::cout << os.str();
                        }*/
                    } //End Propagate task
                    //}  //End group of 4 loop
                    
                    //Method of extraction when a chunk (or one in this case) is done using 'critical'
                    // Note: this doesn't work very well, but it does kinda work
                    /*#pragma omp depend(in:propDone) task
                    {
                      int thread_id = omp_get_thread_num();
                      #pragma omp critical (ptExtract)
                      {
                        std::cout << "(" << thread_id << ") : Storing Task...\n";
                        for(int k=0;k<nbthreads;k++) {
                    typename std::vector<TaskStorage>::iterator mit, tempit;
                    for(mit=perThreadCache[k].begin();mit!=perThreadCache[k].end();++mit) {
                      outputFPs.push_back( guessVec[n] );
                      //Look for relevant data by matching the uuid
                      if (ids[n] == mit->id) {
                        //Check for failure
                        if(mit->failureDetected) {
                          displacementMagnitude.push_back( 1000.0 );
                          std::cerr << "Detected Failure: " << mit->eFail.what() << "\n";
                        } else {
                          //Insert into cache (I don't do anything with it here
                          cache.insert( mit->theData );
                          nvis::vec2 delta = mit->theData.getReturn(guessVec[n].K) - guessVec[n].pos;
                          displacementMagnitude.push_back( nvis::norm(delta) );
                        }
                        //Once done with data, lets pull it out of cache
                        tempit = mit; tempit--;
                        perThreadCache[k].erase(mit);
                        mit = tempit;
                      }
                    }
                        }
                      }//end critical
                    }*/ //end task - store
                }
                
                
                //Method of extraction when all propagations are done (end of loop extraction)
                #pragma omp taskwait
                {
                    //Extract result from per-thread cache  : Could try this AFTER propagation loop...
                    for(int n=0; n<numGuesses; ++n) {
#ifdef _OPENMP
                        int thread_id = omp_get_thread_num();
#else
                        int thread_id = 0;
#endif
                        std::cout << "(" << thread_id << ") : Storing...\n";
                        
                        //Pull out information pertaining to this data run
                        for(int k=0; k<nbthreads; k++) {
#if defined(_WIN32) || defined(C_0X)
                            std::vector<TaskStorage>::iterator mit, tempit;
#else
                            typename std::vector<TaskStorage>::iterator mit, tempit;
#endif
                            for(mit=perThreadCache[k].begin(); mit!=perThreadCache[k].end(); ++mit) {
                                //Look for relevant data by matching the uuid
                                if (ids[n] == mit->id) {
                                    outputFPs.push_back( guessVec[n] );
                                    //Check for failure
                                    if(mit->failureDetected) {
                                        displacementMagnitude.push_back( 1000.0 );
                                        std::cerr << "Detected Failure: " << mit->eFail.what() << "\n";
                                    } else {
                                        //Insert into cache (I don't do anything with it here
                                        cache.insert( mit->theData );
                                        nvis::vec2 delta = mit->theData.getReturn(guessVec[n].K) - guessVec[n].pos;
                                        displacementMagnitude.push_back( nvis::norm(delta) );
                                    }
                                    //Once done with data, lets pull it out of cache
                                    tempit = mit;
                                    tempit--;
                                    perThreadCache[k].erase(mit);
                                    mit = tempit;
                                }
                            }
                        }
                    }
                } //End of task wait region
                //Done with current write
                //omp_unset_lock(&writeLock);
                
            } //end 4 times loop
            
        } //End single directive
    }//End parallel statement
    //Done with lock
#ifdef _OPENMP
    omp_destroy_lock(&writeLock);
#endif
    
    //Show result
    std::cout << " Task-Scheduling Test Result:\n";
    std::cout << " ==================================================================\n";
    std::cout << "   i           x0                                   P(x0)-x0\n";
    std::cout << " -------------------------------------------------------------------\n";
    for (int i=0; i<(int)outputFPs.size(); i++) {
        std::cout << " " << i << "       " << outputFPs[i].pos
                  << "          " << displacementMagnitude[i] << "\n";
    }
    
    
    //Clear thread storage
    for (int i=0; i<nbthreads; i++) {
        perThreadCache[i].clear();
    }
    delete[] perThreadCache;
    
    return 0;
    
}
