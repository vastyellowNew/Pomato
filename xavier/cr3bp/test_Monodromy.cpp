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


// Test function for testing the monodromy matrix versions
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
#include <maps/poincare_map.hpp>
#include <maps/section.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/fixpoints.hpp>
#include <orbital/monodromy.hpp>
#include <orbital/corrections.hpp>
#include <topology/invariant_manifold.hpp>
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
            << "DESCRIPTION: Multiple-Shooting method for fixed point search (test) in circular restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Solution tolerance\n"
            << " -m  | --method <int>             Method (0:STH, 1:England)\n"
            << " -g  | --guess <float> x 2        First guess\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --numCross <int>           Crossing period (num iterates)\n"
            << " -np | --patchpoints              Number of patch points per crossing in corrections\n"
            << " -o  | --output <string>          Output file for regular samples\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    // RHS and ODE solver parameters
    double eps, C, mu, K;
    
    //Jupiter-Europa system example
    //mu = 2.528017705e-5;
    //C = 3.000;
    //double lstar = 6.71100e5;
    //nvis::vec2 x0(-1.23142, 0); //JE Saddle with p=4, ||F(X_0)|| = 1.e-7
    //int numCross = 4;
    //int numPtsPerCross = 4;
    
    //Earth-Moon system example
    //C = 2.96;
    //mu = 1.21505714306e-2;
    //double lstar = 384388.174; //km
    //nvis::vec2 x0(0.7282608, 0); //EM L1 Lyapunov (p=1)
    //int numCross = 1;
    //int numPtsPerCross = 5;
    
    //nvis::vec2 x0(0.84,0); //EM DRO-like saddle (p=3)
    //int numCross = 3;
    //int numPtsPerCross = 5;
    
    //Saturn-Titan Test case (With faulty eigenvalue pairing)
    C = 3.00;
    mu = 0.00023658080508871;
    double lstar = 1221865.;
    //nvis::vec2 x0(1.005589979268675, -0.044160709534052);
    //nvis::vec2 x0(1.003279305334344, 0.065749205422243);
    //int numCross = 2;
    //nvis::vec2 x0(0.986930325086290, -0.099701691791292);
    //nvis::vec2 x0(0.995201721173873, -0.008304749724716);
    //nvis::vec2 x0(0.968544443136223, 0.123837138162046);
    //int numCross = 3;
    //nvis::vec2 x0(1.011311142520915, -0.003990143379070);
    //int numCross = 4;
    //nvis::vec2 x0(1.004064507827794, 0.000000001774448);
    //nvis::vec2 x0(0.958177537061041, 0.029095017878986);
    nvis::vec2 x0(0.949974162578252, -0.053088266905394);
    int numCross = 5;
    
    
    int numPtsPerCross = 5;
    
    me = argv[0];
    eps = 1.0e-8;
    bool verbose = true;
    bool seed_set = false;
    int maxp = 15;
    std::string filename;
    bool file_set = false;
    int method = 0;
    
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
        } else if (arg == "-g" || arg == "--guess") {
            if (i >= argc-2) {
                printUsageAndExit("missing seed");
            }
            x0[0] = atof(argv[++i]);
            x0[1] = atof(argv[++i]);
            seed_set = true;
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing solution tolerance");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-m" || arg == "--method") {
            if (i == argc-1) {
                printUsageAndExit("missing method");
            }
            method = atoi(argv[++i]);
        } else if (arg == "-p" || arg == "--numCross") {
            if (i == argc-1) {
                printUsageAndExit("missing period");
            }
            numCross = atoi(argv[++i]);
        } else if (arg == "-np" || arg == "--patchpoints") {
            if (i == argc-1) {
                printUsageAndExit("missing patch points per crossing");
            }
            numPtsPerCross = atoi(argv[++i]);
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
            verbose = (atoi(argv[++i])>0)? true:false;
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    //if (!seed_set) printUsageAndExit("");
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(1.e-12); //High precision for fixed-points and manifolds
    
    nvis::bbox2 bounds;
    bounds.min() = x0 - nvis::vec2(0.01, 0.01);
    bounds.max() = x0 + nvis::vec2(0.01, 0.01);
    metric<double, 2> euclidean_metric; //Non-periodic space
    
    fixpoint fp;
    bool output = (verbose > 0) ? true : false;
    
    
    // compute orbit at seed point
    //std::vector<nvis::vec2> chain;
    //pmap.map(x0, chain, 10*maxp);
    //std::cerr << "orbit at " << x0 << ":\n";
    
    //Testing linear STM
    std::vector<fixpoint> _fps;
    int p = numCross;
    std::cerr << "Starting Multiple Shooting at " << x0 << " for period " << p << "\n";
    bool found = orbital::solveFixedPoints<map_type>(pmap, euclidean_metric, bounds, x0,
                 C, 5, p, fp, _fps, eps, verbose, 20, 5, true);
                 
    std::cerr << " Resulting Fixed-Point Chain: \n";
    for (int i=0; i<(int) _fps.size(); ++i) {
        std::cerr << "    x_" << i << " = " << _fps[i].pos << "\n";
    }
    
    //Compute monodromy matrix
    std::vector<return_state> stateDataPerReturn;
    pmap.map_complete(fp.pos,stateDataPerReturn,p);
    nvis::fixed_matrix<double,6> monodromy(0);
    monodromy = stateDataPerReturn.back().J; //Full-period STM
    std::cerr << "Monodromy Matrix: " << monodromy << "\n";
    MatrixXd monoM = orbital::mat2MatrixXd<orbital::mat6>(monodromy);
    std::cerr << "det(linear Monodromy) = " << monoM.determinant() << "\n";
    
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
    std::cerr << "    eigVec[0] = " << fp.evec[0] << "\n";
    std::cerr << "    eigVec[1] = " << fp.evec[1] << "\n";
    
    if (found && !fp.saddle) {
        std::cerr << " Computed fixed point is a center.  Skipping Manifold Test.\n";
        return 0;
    }
    
    //Do it again for Numerical STM
    _fps.clear();
    found = orbital::solveFixedPoints<map_type>(pmap, euclidean_metric, bounds, x0,
            C, 5, p, fp, _fps, eps, verbose, 20, 5, false);
            
    std::cerr << " Resulting Fixed-Point Chain: \n";
    for (int i=0; i<(int) _fps.size(); ++i) {
        std::cerr << "    x_" << i << " = " << _fps[i].pos << "\n";
    }
    
    
    //Test if output is ok
    if (!found) {
        std::cerr << "Multiple Shooting diverged... Exiting!\n";
        return 0;
    }
    
    
    //Test output
    fp = _fps[0]; //First fixed point
    std::cerr << " Output First Fixed-Point (Numerical): \n";
    std::cerr << "    x = " << fp.pos << "\n";
    std::cerr << "    p = " << fp.K << "\n";
    std::cerr << "    isSaddle = " << fp.saddle << "\n";
    std::cerr << "    eigVals = [ " << fp.eval[0] << " , " << fp.eval[1] << " ]\n";
    std::cerr << "    log(eigVals) = [ " << std::log(fp.eval[0])
              << " , " << std::log(fp.eval[1]) << " ]\n";
    std::cerr << "    eigVec[0] = " << fp.evec[0] << "\n";
    std::cerr << "    eigVec[1] = " << fp.evec[1] << "\n";
    
    return 0;
    
}
