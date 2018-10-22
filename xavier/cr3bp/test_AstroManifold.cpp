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


// Test function for computing the manifold of a saddle
// -Astro way:  Subdivide ORBIT and start manifold there.
//  This is the old way, but needed for time and result comparisons
//  for ManBVP and ManCurve.
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
#include <cr3bp/cr3bp.hpp>

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
#include <maps/poincare_map.hpp>
#include <maps/section.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/psiWrite.hpp>
#include <orbital/monodromy.hpp>
#include <orbital/corrections.hpp>
#include <maps/fixpoints.hpp>
#include <topology/ManifoldDataStorage.hpp>
#include <topology/invariant_manifold.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Eigen;
using namespace nvis;



typedef orbital::cr3bp                          rhs_type;
typedef nvis::fixed_vector<double, 42>          vec42;
typedef nvis::fixed_vector<double, 6>           vec6;
typedef nvis::fixed_matrix<double, 6>           mat6;
typedef xavier::dp5wrapper<double, 42>          odesolver_type;
typedef xavier::return_map_info<double, 2>      return_type;
typedef xavier::state_info<double, 6>           return_state;
typedef orbital::planar_section<rhs_type, 6, 42>                     section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type> map_type;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

// *******************************    PARAMETERS     *******************************
std::string me;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Computing invariant manifolds of a saddle-type orbit using the original astrodynamics method\n"
            << " which consists of perturbing equally spaced time increments along the orbital path and computing a \n"
            << " specified number of crossings.\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -mu <float>                      mu constant\n"
            << " -C <float>                       C constant\n"
            << " -lstar                           Distance between primaries (in km)\n"
            << " -g  | --guess <float> x 2        First guess\n"
            << " -p  | --period                   Num returns to section for periodic orbit\n"
            << " -e  | --eps <float>              Fixed point refinement tolerance\n"
            << " -np | --patchpoints              Num patch points per crossing in corrections\n"
            << " -d  | --stepValue <float>        Step Value (in km, default = 20km)\n"
            << " -s  | --seeds <int>              Num Seeds on orbit (default = 1e6)\n"
            << " -sp | --seedLength <int>         Num returns to compute for each seed\n"
            << " -o  | --output <string>          Output file for manifold returns (.psi)\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    // RHS and ODE solver parameters
    double eps, C, mu;
    
    //Jupiter-Europa system example
    //mu = 2.528017705e-5;
    //C = 3.000;
    //double lstar = 6.71100e5;
    //nvis::vec2 x0(-1.23142, 0); //JE Saddle with p=4, ||F(X_0)|| = 1.e-7
    //int numCross = 4;
    //int numPtsPerCross = 4;
    
    //Earth-Moon system example
    C = 2.96;
    mu = 1.21505714306e-2;
    double lstar = 384388.174; //km
    nvis::vec2 x0(0.7282608, 0); //EM L1 Lyapunov (p=1)
    int numCross = 1;
    int numPtsPerCross = 5;
    
    //nvis::vec2 x0(0.84,0); //EM DRO-like saddle (p=3)
    //int numCross = 3;
    //int numPtsPerCross = 5;
    
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    double dValue = 20.0; //km
    bool seed_set = false;
    int numSeeds = 1e6; //seeds from orbit
    int numSeedReturns = 3; //*numCross
    double velocityScale = 0.16;
    std::string filename;
    filename = "none";
    
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-C") {
            if (i == argc-1) {
                printUsageAndExit("missing C constant");
            }
            C = atof(argv[++i]);
        } else if (arg == "-mu" || arg == "--mu") {
            if (i == argc-1) {
                printUsageAndExit("missing mu constant");
            }
            mu = atof(argv[++i]);
        } else if (arg == "-lstar") {
            if (i == argc-1) {
                printUsageAndExit("missing lstar value");
            }
            lstar = atof(argv[++i]);
        } else if (arg == "-g" || arg == "--guess") {
            if (i >= argc-2) {
                printUsageAndExit("missing seed");
            }
            x0[0] = atof(argv[++i]);
            x0[1] = atof(argv[++i]);
            seed_set = true;
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-p" || arg == "--period") {
            if (i == argc-1) {
                printUsageAndExit("missing period (num crossings)");
            }
            numCross = atoi(argv[++i]);
        } else if (arg == "-np" || arg == "--patchpoints") {
            if (i == argc-1) {
                printUsageAndExit("missing patch points per crossing");
            }
            numPtsPerCross = atoi(argv[++i]);
        } else if (arg == "-d" || arg == "--stepValue") {
            if (i == argc-1) {
                printUsageAndExit("missing d-value");
            }
            dValue = atof(argv[++i]);
        } else if (arg == "-s" || arg == "--seeds") {
            if (i == argc-1) {
                printUsageAndExit("missing number of seed points");
            }
            numSeeds = atoi(argv[++i]);
        } else if (arg == "-sp" || arg == "--seedLength") {
            if (i == argc-1) {
                printUsageAndExit("missing number of returns for seed propagation");
            }
            numSeedReturns = atoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output file");
            }
            filename = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (filename=="none") {
        printUsageAndExit("");
    }
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(1.e-12); //High precision for fixed-points and manifolds
    
    nvis::bbox2 bounds;
    bounds.min() = x0 - nvis::vec2(0.1, 0.1);
    bounds.max() = x0 + nvis::vec2(0.1, 0.1);
    metric<double, 2> euclidean_metric; //Non-periodic space
    
    fixpoint fp;
    bool output = (verbose > 0) ? true : false;
    
    
    //Testing a fixed period
    std::vector<fixpoint> _fps;
    int p = numCross;
    std::cerr << "Starting Multiple Shooting at " << x0 << " for period " << p << "\n";
    bool found = orbital::solveFixedPoints<map_type>(pmap, euclidean_metric, bounds, x0,
                 C, 5, p, fp, _fps, eps, ( (verbose > 1) ? true : false ), 20, 5, false);
                 
    std::cerr << " Resulting Fixed-Point Chain: \n";
    for (int i=0; i<(int) _fps.size(); ++i) {
        std::cerr << "    x_" << i << " = " << _fps[i].pos << "\n";
    }
    
    //Compute monodromy matrix
    std::vector<return_state> stateDataPerReturn;
    pmap.map_complete(fp.pos,stateDataPerReturn,p);
    nvis::fixed_matrix<double,6> mMat(0);
    mMat = stateDataPerReturn.back().J; //Full-period STM
    std::cerr << "Monodromy Matrix: " << mMat << "\n";
    
    //Test if output is ok
    if (!found) {
        std::cerr << "Multiple Shooting diverged... Exiting!\n";
        return 0;
    }
    
    
    //Test output
    fp = _fps[0]; //First fixed point
    std::cerr << " Output First Fixed-Point: \n";
    std::cerr << "    x = " << fp.pos << "\n";
    std::cerr << "    p = " << fp.K << "\n";
    std::cerr << "    isSaddle = " << fp.saddle << "\n";
    std::cerr << "    eigVals = [ " << fp.eval[0] << " , " << fp.eval[1] << " ]\n";
    std::cerr << "    log(eigVals) = [ " << std::log(fp.eval[0])
              << " , " << std::log(fp.eval[1]) << " ]\n";
    std::cerr << "    EigVec[0] = " << fp.fullEvec[0] << "\n";
    std::cerr << "    EigVec[1] = " << fp.fullEvec[1] << "\n";
    std::cerr << "    mapEigVec[0] = " << fp.evec[0] << "\n";
    std::cerr << "    mapEigVec[1] = " << fp.evec[1] << "\n";
    
    if (found && !fp.saddle) {
        std::cerr << " Computed fixed point is a center.  Skipping Manifold Test.\n";
        return 0;
    }
    
    
    //Map analysis param
    map_analysis_param _params; //Parameters of map
    _params.nb_iterations = 10*numCross;
    _params.max_period = numCross;
    _params.the_metric = euclidean_metric;
    //_params.verbose = ( (verbose > 1) ? true : false );
    _params.verbose = false;
    _params.bounds = bounds;
    _params.linearMonodromy = true; //Add option later
    
    //Start timer for comparison
    nvis::timer timer;
    
    //Need eigenvectors from monodromy
    vec6 vs0 = fp.fullEvec[0];
    vec6 vu0 = fp.fullEvec[1];
    
    //Compute the points along the orbit & stable/unstable starting conditions
    double dt = fp.timePeriod / (double) numSeeds;
    std::vector< vec6 > startStates; //Set of all initial conditions
    std::vector<int> types; //0 = m_stable, 1 = p_stable, 2 = m_unstable, 3 = p_unstable
    std::vector<int> seedNumber;
    //Store initial info
    vec6 y0 = section.mapToState( fp.pos );
    return_state fInfo;
    fInfo.x = y0;
    fInfo.t = 0.0;
    fInfo.J = mat6::identity();
    mat6 stm = mat6::identity();
    double dStep = dValue / lstar;
    for(int i=0; i<numSeeds; i++) {
        //Mapping should not throw errors here
        fInfo = pmap.integrate_state(y0,dt); //Always starts STM from Identity
        //Compute the STM for this point
        stm = fInfo.J * stm; //Phi(t2,t0) = Phi(t2,t1) * Phi(t1,t0)
        //Transform the eigenvectors
        vec6 vsi = stm * vs0;
        vec6 vui = stm * vu0;
        //Normalize by position
        double vsiPos = sqrt(vsi[0]*vsi[0]+vsi[1]*vsi[1]+vsi[2]*vsi[2]);
        double vuiPos = sqrt(vui[0]*vui[0]+vui[1]*vui[1]+vui[2]*vui[2]);
        vsi /= vsiPos;
        vui /= vuiPos;
        //Setup the 4 initial states
        startStates.push_back( y0 - dStep*vsi);
        types.push_back( 0 ); //m_stable
        seedNumber.push_back( i );
        startStates.push_back( y0 + dStep*vsi);
        types.push_back( 1 ); //p_stable
        seedNumber.push_back( i );
        startStates.push_back( y0 - dStep*vui);
        types.push_back( 2 ); //m_unstable
        seedNumber.push_back( i );
        startStates.push_back( y0 + dStep*vui);
        types.push_back( 3 ); //p_unstable
        seedNumber.push_back( i );
        //The next initial position
        y0 = fInfo.x;
    }
    
    //Setup computation per thread
    size_t nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cerr << nthreads << " threads available for computation\n";
    
    //Point data cloud data caches (for Avizo data file)
    std::vector<PointPSI> ptsPerThread[nthreads], pts;
    std::vector<int> idsPerThread[nthreads], ids;
    std::vector< std::vector<double> > dataPerPointPerThread[nthreads];
    std::vector< std::vector<double> > totalData;
    std::vector< std::string > psi_labels, psi_symbols;
    psi_labels.push_back( std::string("Seed Number") );
    psi_labels.push_back( std::string("Crossing Number") );
    psi_labels.push_back( std::string("Time") );
    psi_labels.push_back( std::string("Thread") );
    psi_labels.push_back( std::string("Type") ); //stable-unstable
    psi_labels.push_back( std::string("Perturbation") ); //Plus-minus
    psi_labels.push_back( std::string("ManifoldType") ); // 0->3
    psi_symbols.push_back( std::string("s"));
    psi_symbols.push_back( std::string("c"));
    psi_symbols.push_back( std::string("t"));
    psi_symbols.push_back( std::string("h"));
    psi_symbols.push_back( std::string("q"));
    psi_symbols.push_back( std::string("w"));
    psi_symbols.push_back( std::string("m"));
    
    
    //Run parallel computation
    int numManifolds = (int) startStates.size();
    int counter = 0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0; n<numManifolds; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            std::vector<return_type> output;
            //Forward or backward?
            int niter = (types[n]<2) ? -1 : 1;
            niter *= numSeedReturns*p;
            try {
                //Run the map given a state
                pmap.map( startStates[n], output, niter);
            } catch (...) {
                //Errors will occur, but store all available points
            }
            
#pragma openmp atomic
            ++counter;
            
            //Prompt user
            std::ostringstream os;
            double elTime = timer.elapsed();
            int c = counter;
            os << "\r" << c << "/" << numManifolds << " ("
               << (float)c*100/numManifolds << "%) in "
               << elTime << "s. (" << (float)c/elTime << "Hz)  Roughly "
               << (((float)numManifolds-(float)c)/((float)c/elTime))/60.0
               << " min remaining...\r";
            std::cout << os.str() << std::flush;
            
            //Storage
            for (int i=0; i<output.size(); i++) {
                nvis::vec2 x = output[i].x;
                ptsPerThread[thread_id].push_back( PointPSI(x[0],0.0,velocityScale*x[1]) );
                idsPerThread[thread_id].push_back( n );
                std::vector<double> dataPerPoint;
                dataPerPoint.push_back( (double) seedNumber[n] );
                dataPerPoint.push_back( (double) i );
                dataPerPoint.push_back( output[i].t );
                dataPerPoint.push_back( (double) thread_id );
                dataPerPoint.push_back( (types[n] == 0 || types[n] == 1) ? 0.0 : 1.0 );
                dataPerPoint.push_back( (types[n] == 0 || types[n] == 2) ? 0.0 : 1.0 ); //0=m,1=p
                dataPerPoint.push_back( (double) types[n] ); //0->3
                dataPerPointPerThread[thread_id].push_back( dataPerPoint );
            }
            
            
        }//End for loop
    }//End parallel statement
    
    //Assemble resulting data
    std::vector<double> sNum, cNum, tVec, tdNum, suMan, pmPert, mType;
    for(int h=0; h<nthreads; h++) {
        for (int i=0; i<(int)ptsPerThread[h].size(); i++) {
            pts.push_back( ptsPerThread[h][i] );
            ids.push_back( idsPerThread[h][i] );
            sNum.push_back( dataPerPointPerThread[h][i][0] );
            cNum.push_back( dataPerPointPerThread[h][i][1] );
            tVec.push_back(   dataPerPointPerThread[h][i][2] );
            tdNum.push_back( dataPerPointPerThread[h][i][3] );
            suMan.push_back( dataPerPointPerThread[h][i][4] );
            pmPert.push_back( dataPerPointPerThread[h][i][5] );
            mType.push_back( dataPerPointPerThread[h][i][6] );
        }
    }
    //Write to file
    std::vector< std::vector<double> > dataVecs;
    dataVecs.push_back(sNum);
    dataVecs.push_back(cNum);
    dataVecs.push_back(tVec);
    dataVecs.push_back(tdNum);
    dataVecs.push_back(suMan);
    dataVecs.push_back(pmPert);
    dataVecs.push_back(mType);
    bool ok = PSI_WriteToFile(pts,ids,dataVecs,psi_labels,psi_symbols,filename.c_str());
    
    //Stop timer
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. (" << (float)numManifolds/elapsed << "Hz)\n";
    
    
    
    return 0;
    
}
