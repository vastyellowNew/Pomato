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


// Test function for guesses of fixed points
//  - Employes Corrector class (with CorrectionsRegulator)
// Wayne Schlei

#include <iostream>
#include <vector>
#include <list>
#include <util/wall_timer.hpp>
#include <complex>
#include <sstream>
#include <limits>
#include <exception>


// data structure
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>

// map API
#include <maps/definitions.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/map_analysis.hpp>
#include <maps/fixpoints.hpp>
#include <maps/section.hpp>
#include <maps/poincare_map.hpp>
#include <orbital/controller.hpp>
#include <orbital/corrections.hpp>
// cr3bp
#include <cr3bp/planar_section.hpp>
#include <cr3bp/cr3bp.hpp>

#include <teem/nrrd.h>
#include <cr3bp/psiWrite.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

typedef orbital::cr3bp                          rhs_type;
typedef nvis::fixed_vector<double, 42>          vec42;
typedef nvis::fixed_vector<double, 6>           vec6;
typedef nvis::fixed_matrix<double, 6>           mat6;
typedef xavier::dp5wrapper<double, 42>          odesolver_type;
typedef orbital::CorrectionsRegulator           reg_type;
#if defined(_WIN32) || defined(C_0X)
typedef reg_type::OperatingMode                 mode_type;
#else
typedef typename reg_type::OperatingMode        mode_type;
#endif

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

typedef orbital::planar_section<rhs_type, 6, 42>        section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type> map_type;

template<typename T>
using vec_of_vecs = std::vector< std::vector< T > >;

template<typename T>
void flatten(std::vector<T>& out, const vec_of_vecs<T>& in) {
	out.clear();
	std::for_each(in.begin(), in.end(), [&](const std::vector<T>& v) {
		std::copy(v.begin(), v.end(), std::back_inserter(out));
	});
}


// *******************************    PARAMETERS     *******************************
std::string me;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Testing initial guesses for fixed point search in circular restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -c  | --convergence <float>      Convergence tolerance\n"
            << " -x  | --guess <float> x 2        First guess\n"
            << " -b  | --bounds <float> x 4       Bounds around first guess for sampling box\n"
            << " -g  | --grid <int> x 2           Grid sizes within box in 2D map\n"
            << " -z  | --scale <float>            ScaleFactor on z-plot (xdot) dim\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --patchPoints <int>        Points per revolution (init)\n"
            << " -o  | --output <string>          Output file for regular samples\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

double getMethodValue(mode_type mode)
{
    double value = -1.0; //Fail by default
    switch (mode) {
        case reg_type::SINGLE_SHOOT :
            value = 0.0;
            break;
        case reg_type::RSINGLE_SHOOT :
            value = 1.0;
            break;
        case reg_type::MULTIPLE_SHOOT :
            value = 2.0;
            break;
        case reg_type::RESAMPLE :
            value = 3.0;
            break;
        case reg_type::DISTRIBUTED_ERROR :
            value = 4.0;
            break;
        case reg_type::HAMSEC_OFF :
            value = 5.0;
            break;
        case reg_type::PERIODICITY_HOMOTOPY : //Temporarily skipping homotopy.......
            value = 6.0;
            break;//...........................................................
        case reg_type::QUASI_NEWTON :
            value = 7.0;
            break;
		default:
			throw std::runtime_error("invalid value in double getMethodValue(mod_type)");
    }
    return value;
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    // RHS and ODE solver parameters
    double eps, convTol, C, mu;
    
    //Jupiter-Europa system example
    //mu = 2.528017705e-5;
    //C = 3.000;
    //nvis::vec2 x0(-1.205,0); //Assumed IC for testing
    //nvis::vec2 x0(-1.23142, 0); //Saddle with p=4, ||F(X_0)|| = 1.e-7
    //int numCross = 4;
    
    //Earth-Moon Problem cases
    mu = 1.21505714306e-2;
    C = 2.96;
    //Convergence with SINGLE_SHOOT
    nvis::vec2 guess(0.7282608,0.0);  //Saddle p=1, L1 Lyap
    int numCross = 1;
    //Convergence with DISTRIBUTED_ERROR (mixed sampling)
    //nvis::vec2 x0(0.855334,0.0217615);
    //int numCross = 3;
    //Convergence (on wrong orbit) with PERIODICITY_HOMOTOPY
    //nvis::vec2 x0(0.302889,0.574074);
    //int numCross = 8;
    
    //Others...
    //nvis::vec2 x0(-0.284058,1.38889); //Center p=8
    //int numCross = 8;
    //nvis::vec2 x0(0.288406,0.3888889);
    //int numCross = 9;
    //nvis::vec2 x0(0.36087,0.648148);
    //int numCross = 7;
    
    //Bounding Box
    nvis::bbox2 bounds;
    bounds.min() = guess - nvis::vec2(0.01, 0.06);
    bounds.max() = guess + nvis::vec2(0.01, 0.06);
    //Sampling parameters
    nvis::ivec2 gridRes(2,2);
    double sf = 0.16;  //Scale factor on xdot dimension (for plotting)
    
    
    me = argv[0];
    eps = 1.0e-12;
    convTol = 1.0e-8;
    int verbose = 0;
    bool seed_set = true;
    int maxiter = 20;
    int pointsPerRev = 5;
    std::string filename("ConvTest.psi");
    bool file_set = true;
    
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
        } else if (arg == "-x" || arg == "--guess") {
            if (i >= argc-2) {
                printUsageAndExit("missing seed");
            }
            guess[0] = atof(argv[++i]);
            guess[1] = atof(argv[++i]);
            seed_set = true;
        } else if (arg == "-b" || arg == "--bounds") {
            if (i==argc-1) {
                printUsageAndExit("missing bounds");
            }
            bounds.min() = nvis::vec2(atof(argv[i+1]),atof(argv[i+2]));
            bounds.max() = nvis::vec2(atof(argv[i+3]),atof(argv[i+4]));
            i += 4;
        } else if (arg == "-g" || arg == "--grid") {
            if (i== argc-1) {
                printUsageAndExit("missing number of samples");
            }
            gridRes[0] = atoi(argv[++i]);
            gridRes[1] = atoi(argv[++i]);
        } else if (arg == "-z" || arg == "--scale") {
            if (i== argc-1) {
                printUsageAndExit("missing scale factor");
            }
            sf = atof(argv[++i]);
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
    if (filename=="none") {
        printUsageAndExit("");
    }
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    metric<double, 2> euclidean_metric;
    
    int counter = 0;
    std::vector<fixpoint> fps;
    nvis::bbox2 epsBox;
    epsBox.min() = guess - nvis::vec2(5.e-4,5.e-4);
    epsBox.max() = guess + nvis::vec2(5.e-4,5.e-4);
    
    //PSI file headers
    std::vector< std::string > psi_labels, psi_symbols;
    psi_labels.push_back( std::string("Convergence") );
    psi_labels.push_back( std::string("Method") );
    psi_labels.push_back( std::string("IsSameFP") );
    psi_labels.push_back( std::string("Iters") );
    psi_labels.push_back( std::string("NormF0") );
    psi_labels.push_back( std::string("xEnd") );
    psi_labels.push_back( std::string("xdotEnd") );
    psi_labels.push_back( std::string("ComputeTime_s") );
    psi_symbols.push_back( std::string("c") );
    psi_symbols.push_back( std::string("m") );
    psi_symbols.push_back( std::string("h") );
    psi_symbols.push_back( std::string("i") );
    psi_symbols.push_back( std::string("f0") );
    psi_symbols.push_back( std::string("k") );
    psi_symbols.push_back( std::string("d") );
    psi_symbols.push_back( std::string("t") );
    //PSI data cache
    size_t nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cerr << nthreads << " threads available for computation\n";
    vec_of_vecs<PointPSI> ptsPerThread(nthreads);
	std::vector<PointPSI> pts;
    vec_of_vecs<int> idsPerThread(nthreads);
	std::vector<int> ids;
    std::vector< vec_of_vecs<double> > dataPerPointPerThread(nthreads);
    vec_of_vecs<double> totalData;
    
    //Use the close proximity around IC based on bounds/grid to run corrections
    if (file_set) {
        nvis::timer timer;
        //double *data = (double *)calloc(2*n*m, sizeof(double));
        int numPts = gridRes[0]*gridRes[1];
        nvis::vec2 spacing( bounds.max() - bounds.min() );
        spacing[0] = spacing[0]/((double) gridRes[0]+1.);
        spacing[1] = spacing[1]/((double) gridRes[1]+1.);
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (int kk=0 ; kk<numPts ; ++kk) {
                int i = kk%gridRes[0];
                int j = kk/gridRes[0];
                //guess
                nvis::vec2 delta(spacing[0]*((double)i+1),spacing[1]*((double)j+1));
                nvis::vec2 x0 = bounds.min() + delta;
                nvis::timer threadTimer;
                fixpoint fp;
                std::vector<fixpoint> iterates;
                
                int thread_id = 0;
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                
                //Run the corrections algorithm
                orbital::PeriodicOrbitCorrector<map_type> corrector(
                    pmap, x0, C, numCross, pointsPerRev, verbose);
                corrector.setMaxIters(maxiter);
                corrector.setTolerance(convTol);
                corrector.setStopMode(reg_type::THE_END);
                corrector.setJobIndex(thread_id);
                //Run corrector (with various methods)
                orbital::CorrectionResult theResult = corrector.correct();
                bool found = theResult.converged;
                if (found) {
                    //Compute Map iterates and gather fixed-point information
                    iterates.clear();
                    corrector.setFixedPoints(fp,iterates);
                }
                
                bool rightFP = false;
                if (found && bounds.inside(fp.pos) && epsBox.inside(fp.pos)) {
                    rightFP = true;
                }
                
                //Set data values
                ptsPerThread[thread_id].push_back(PointPSI(x0[0],0,sf*x0[1]));//xdot on zaxis
                idsPerThread[thread_id].push_back( kk );
                std::vector<double> dataPerPoint;
                dataPerPoint.push_back( (found)? 1.0 : 0.0 );
                dataPerPoint.push_back( getMethodValue(corrector.getMode()) );
                dataPerPoint.push_back( (rightFP)? 1.0 : 0.0 );
                dataPerPoint.push_back( (double) theResult.iterations );
                dataPerPoint.push_back( theResult.normF0 );
                dataPerPoint.push_back( fp.pos[0] );//xf - final value after methods
                dataPerPoint.push_back( fp.pos[1] );//xdotf - final value after methods
                double dt = threadTimer.elapsed();
                dataPerPoint.push_back( dt ); //Compute time in seconds
                dataPerPointPerThread[thread_id].push_back(dataPerPoint);
                
                
#pragma openmp atomic
                ++counter;
                dt = timer.elapsed();
                int c = counter;
                std::ostringstream os;
                os << "\r" << c << "/" << numPts << " (" << (float)c*100/((float) numPts) << "%) in "
                   << dt << "s. (" << (float)c/dt << "Hz)           \r";
                std::cout << os.str() << std::flush;
            }
        }
    }
    
    
    //Assemble data caches
    std::vector<double> foundNums, mNums, rfp, itNums, f0Nums, xfs, xdfs, cTimes;
    for (int h=0; h<nthreads; h++) {
        for (int i=0; i<(int)ptsPerThread[h].size(); i++) {
            pts.push_back( ptsPerThread[h][i] );
            ids.push_back( idsPerThread[h][i] );
            foundNums.push_back( dataPerPointPerThread[h][i][0] );
            mNums.push_back( dataPerPointPerThread[h][i][1] );
            rfp.push_back( dataPerPointPerThread[h][i][2] );
            itNums.push_back( dataPerPointPerThread[h][i][3] );
            f0Nums.push_back( dataPerPointPerThread[h][i][4] );
            xfs.push_back( dataPerPointPerThread[h][i][5] );
            xdfs.push_back( dataPerPointPerThread[h][i][6] );
            cTimes.push_back( dataPerPointPerThread[h][i][7] );
        }
    }
    //Write to PSI file
    vec_of_vecs<double> dataVecs;
    dataVecs.push_back( foundNums );
    dataVecs.push_back( mNums );
    dataVecs.push_back( rfp );
    dataVecs.push_back( itNums );
    dataVecs.push_back( f0Nums );
    dataVecs.push_back( xfs );
    dataVecs.push_back( xdfs );
    dataVecs.push_back( cTimes );
    bool ok = PSI_WriteToFile(pts,ids,dataVecs,psi_labels,psi_symbols,filename.c_str());
    
    return 0;
    
}
