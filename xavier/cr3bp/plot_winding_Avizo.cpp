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
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/newton.hpp>
#include <maps/fixpoints.hpp>
#include <maps/invariant_manifold.hpp>
#include <maps/index.hpp>
#include <maps/poincare_map.hpp>

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
double eps, C, mu, K;

typedef orbital::cr3bp                                    rhs_type;
typedef nvis::fixed_vector<double, 42>                                vec42;
typedef nvis::fixed_vector<double, 6>                                 vec6;
typedef nvis::fixed_matrix<double, 6>                                 mat6;
typedef xavier::dp5wrapper<double, 42>                                odesolver_type;
typedef orbital::planar_section<rhs_type, 6, 42>                      section_type;
//typedef orbital::AngleTracker<42, 3>                                  tracker_type;
typedef orbital::MultipleAngleTracker<42, 3>                          tracker_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type>  map_type;
typedef map_type::return_type                                         return_type;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

// *******************************    PARAMETERS     *******************************
std::string me;
std::string filename;
std::string pointFile;
nvis::vec2 L3(-1.00001053340710,0);
nvis::fixed_vector<size_t, 2> res(1000, 1000);

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
            << " -r  | --resolution <uint> x 2    Image resolution\n"
            << " -n  | --nseeds <int>             Number of seeds\n"
            << " -x0 | --center <float> x2        Center of rotation\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --maxp <int>               Max considered period\n"
            << " -o  | --output <string>          Output base name\n"
            << " -op | --pointFile <string>       Point-clout output (.psi)\n"
            << " -s  | --seed <int>               Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal)\n"
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
    nvis::bbox2 bounds;
    bounds.min()[0] = -0.4;
    bounds.min()[1] = -2.5;
    bounds.max()[0] = 1.1;
    bounds.max()[1] = 2.5;
    
    
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    bool bounds_set = false;
    int maxp = 100;
    int p_max = 12;
    filename = "none";
    pointFile = "points.psi";
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
        } else if (arg == "-op" || arg == "--pointFile") {
            if (i== argc-1) {
                printUsageAndExit("missing point file name");
            }
            pointFile = argv[++i];
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
    pmap.setPrecision(eps);
    
    float* raster = (float*)calloc(4*res[0]*res[1], sizeof(float));
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
    
    //Point data cloud data caches (for Avizo data file)
    std::vector<PointPSI> ptsPerThread[nthreads], pts;
    std::vector<int> idsPerThread[nthreads], ids;
    std::vector< std::vector<double> > dataPerPointPerThread[nthreads]; //4xnumPtsPerThread
    std::vector< std::vector<double> > totalData;
    std::vector< std::string > psi_labels, psi_symbols;
    psi_labels.push_back( std::string("Winding Number x-xd") );
    psi_labels.push_back( std::string("Winding Number x-yd") );
    psi_labels.push_back( std::string("Winding Number xd-yd") );
    psi_labels.push_back( std::string("Poloidal Period Approx") );
    psi_labels.push_back( std::string("Crossing Number") );
    psi_labels.push_back( std::string("Time") );
    psi_labels.push_back( std::string("Thread") );
    psi_symbols.push_back( std::string("w"));
    psi_symbols.push_back( std::string("u"));
    psi_symbols.push_back( std::string("v"));
    psi_symbols.push_back( std::string("p"));
    psi_symbols.push_back( std::string("c"));
    psi_symbols.push_back( std::string("t"));
    psi_symbols.push_back( std::string("h"));
    
    
    int counter=0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<nseeds ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            std::vector<return_type> orbit;
            tracker_type tracker(barycenter);
            try {
                pmap.map_and_track_complete<tracker_type>(seeds[n], orbit, maxp, tracker);
            } catch(...) {
                continue;
            }
            if (orbit.size() < maxp) {
                continue;
            }
            
            // toroidal / poloidal (q/p)
            double w = 2.*M_PI*orbit.size()/orbit.back().delta_theta[0];
            if (std::isnan(w) || std::isinf(w)) {
                w = 1000;
            }
            double u = 2.*M_PI*orbit.size()/orbit.back().delta_theta[1];
            if (std::isnan(u) || std::isinf(u)) {
                u = 1000;
            }
            double v = 2.*M_PI*orbit.size()/orbit.back().delta_theta[2];
            if (std::isnan(v) || std::isinf(v)) {
                v = 1000;
            }
            
            //Period Approximation
            boost::rational<int> alpha_w, alpha_u, alpha_v;
            alpha_w = xavier::rational_approx_CF<int>(w,p_max); //Limit max denom as p_max
            alpha_u = xavier::rational_approx_CF<int>(u,p_max);
            alpha_v = xavier::rational_approx_CF<int>(v,p_max);
            //Check for zeros
            //if (abs(alpha_w.numerator()) > 10000) { alpha_w.assign(0,1); }
            //if (abs(alpha_u.numerator()) > 10000) { alpha_u.assign(0,1); }
            //if (abs(alpha_v.numerator()) > 10000) { alpha_v.assign(0,1); }
            
            //LCM of denominators to find highest poloidal period
            int best_p = boost::lcm(alpha_w.denominator(),alpha_u.denominator());
            best_p = boost::lcm(best_p, alpha_v.denominator());
            //GCD of denominators
            //int best_p = boost::gcd(alpha_w.denominator(),alpha_u.denominator());
            //best_p = boost::gcd(best_p, alpha_v.denominator());
            
            //Test the tracking system code
            std::vector<double> wn;
            wn = tracker.getWindingFactors(orbit.size(),orbit.back().delta_theta);
            //Period_range code
            std::vector<double>::iterator it;
            std::vector<boost::rational<int> > alpha_PR;
            std::vector<int> ps;
            for (it = wn.begin(); it!=wn.end(); ++it) {
                //compute rational number via continued fraction
                alpha_PR.push_back( xavier::rational_approx_CF<int>(*it,p_max) );
                ps.push_back( fabs(alpha_PR.back().denominator()) );
            }
            
            //Test output
            if (thread_id == 0) {
                std::cerr << "  Test n = " << n << " : seed = " << seeds[n] << "\n";
                std::cerr << "            Winding Numbers : w = [ " << w << ", " << u << ", " << v << " ] \n";
                std::cerr << "            Tracker (getWF) : w = [ " << wn[0] << ", " << wn[1] << ", " << wn[2] << " ]\n";
                std::cerr << "            Approx Int Ratio: w = [ " << alpha_w << ", " << alpha_u << ", " << alpha_v << " ] \n";
                std::cerr << "            Tracker -> Ratio: w = [ " << alpha_PR[0] << ", " << alpha_PR[1] << ", " << alpha_PR[2] << " ] \n";
                std::cerr << "            Best Periods    : p = [" << alpha_w.denominator() << ", " << alpha_u.denominator()
                          << ", " << alpha_v.denominator() << "]\n";
                std::cerr << "            Tracker -> Per. : p = [" << ps[0] << ", " << ps[1] << ", " << ps[2] << " ]\n";
            }
            
            orbit.push_back(return_type());
            orbit.back().x = seeds[n];
            for (int i=0 ; i<orbit.size() ; ++i) {
                nvis::vec2 x = orbit[i].x;
                //Check if leaves domain bounds
                nvis::vec2 c = (x-bounds.min())/bounds.size()*nvis::vec2(res);
                if (c[0]<0 || c[0]>=res[0] || c[1]<0 || c[1]>=res[1]) {
                    continue;
                }
                //Save data for PSI file
                ptsPerThread[thread_id].push_back(PointPSI(x[0],0.0,x[1]));
                idsPerThread[thread_id].push_back( n );
                std::vector<double> dataPerPoint;
                dataPerPoint.push_back(w);
                dataPerPoint.push_back(u);
                dataPerPoint.push_back(v);
                dataPerPoint.push_back(best_p);
                dataPerPoint.push_back((double)i);
                dataPerPoint.push_back(orbit[i].t);
                dataPerPoint.push_back((double)thread_id);
                dataPerPointPerThread[thread_id].push_back(dataPerPoint);
                
                int m = (int)c[0] + (int)c[1]*res[0];
#pragma openmp atomic
                raster[m] = w;
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
    
    //PSI Data File - always on if run correctly
    //  Unwrap data per thread:
    std::vector<double> wNums, uNums, vNums, pNums, cNums, tVec, threadVec;
    for (int h=0; h<nthreads; h++) {
        for (int i=0; i<(int)ptsPerThread[h].size(); i++) {
            pts.push_back( ptsPerThread[h][i] );
            ids.push_back( idsPerThread[h][i] );
            wNums.push_back( dataPerPointPerThread[h][i][0] );
            uNums.push_back( dataPerPointPerThread[h][i][1] );
            vNums.push_back( dataPerPointPerThread[h][i][2] );
            pNums.push_back( dataPerPointPerThread[h][i][3] );
            cNums.push_back( dataPerPointPerThread[h][i][4] );
            tVec.push_back( dataPerPointPerThread[h][i][5] );
            threadVec.push_back( dataPerPointPerThread[h][i][6] );
        }
    }
    //Write to file
    std::vector< std::vector<double> > dataVecs;
    dataVecs.push_back(wNums);
    dataVecs.push_back(uNums);
    dataVecs.push_back(vNums);
    dataVecs.push_back(pNums);
    dataVecs.push_back(cNums);
    dataVecs.push_back(tVec);
    dataVecs.push_back(threadVec);
    bool ok = PSI_WriteToFile(pts,ids,dataVecs,psi_labels,psi_symbols,pointFile.c_str());
    return 0;
}
