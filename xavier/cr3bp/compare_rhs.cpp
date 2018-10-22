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

// maps API
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/poincare_map.hpp>
#include <graphics/colors.hpp>

// cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/pcr3bp.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>


#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

nvis::ivec2 res(1000,1000);
std::string filename;
double eps, C, mu, K;
int seeding;
nvis::vec2 L3(-1.00001053340710, 0);
int maxp = 100;

typedef orbital::cr3bp                                                  rhs_type1;
typedef orbital::pcr3bp                                                 rhs_type2;
typedef orbital::pcr3bp_reduced                                         rhs_type3;
typedef xavier::dp5wrapper<double, 42>                                  odesolver_type1;
typedef xavier::dp5wrapper<double, 20>                                  odesolver_type2;
typedef xavier::dp5wrapper<double, 4>                                   odesolver_type3;
typedef orbital::planar_section<rhs_type1, 6, 42>                       section_type1;
typedef orbital::planar_section<rhs_type2, 4, 20>                       section_type2;
typedef orbital::planar_section_reduced<rhs_type3, 4>                   section_type3;
typedef xavier::poincare_map<rhs_type1, odesolver_type1, section_type1> map_type1;
typedef xavier::poincare_map<rhs_type2, odesolver_type2, section_type2> map_type2;
typedef xavier::poincare_map<rhs_type3, odesolver_type3, section_type3> map_type3;
typedef orbital::AngleTracker<42, 3>                                    tracker_type1;
typedef orbital::AngleTracker<20, 2>                                    tracker_type2;
typedef orbital::AngleTracker<4,  2>                                    tracker_type3;
typedef xavier::planar_rotation_winding_number                          winding_type;

using namespace xavier;

nvis::bbox2 bounds;

inline int pos(const nvis::vec2& x)
{
    const double& minx = bounds.min()[0];
    const double& miny = bounds.min()[1];
    const double& maxx = bounds.max()[0];
    const double& maxy = bounds.max()[1];
    
    int i = floor(res[0] * (x[0] - minx) / (maxx - minx));
    int j = floor(res[1] * (x[1] - miny) / (maxy - miny));
    if (i < 0 || i >= res[0] || j < 0 || j >= res[1]) {
        return -1;
    }
    return i + j*res[0];
}

inline void print_color(const nvis::vec3& col, int i, float* data)
{
    if (i < 0) {
        return;
    }
    data[3*i  ] += col[0];
    data[3*i+1] += col[1];
    data[3*i+2] += col[2];
}

void generate_seeds(std::vector< nvis::vec2 >& p, int seeding, int nsamples)
{
    const double& minx = bounds.min()[0];
    const double& miny = bounds.min()[1];
    const double& maxx = bounds.max()[0];
    const double& maxy = bounds.max()[1];
    
    p.clear();
    switch (seeding) {
        case 0: {
            for (int n = 0 ; n < nsamples ; ++n)
                p.push_back(nvis::vec2(minx + drand48()*(maxx - minx),
                                       miny + drand48()*(maxy - miny)));
            break;
        }
        case 1: {
            double dx = (maxx - minx) / (double)nsamples;
            for (int n = 0 ; n < nsamples ; ++n) {
                p.push_back(nvis::vec2(minx + (double)n*dx, 0.5*(miny + maxy)));
            }
            break;
        }
        case 2: {
            double dy = (maxy - miny) / (double)nsamples;
            for (int n = 0 ; n < nsamples ; ++n) {
                p.push_back(nvis::vec2(0.5*(minx + maxx), miny + (double)n*dy));
            }
            break;
        }
        case 3: {
            double dx = (maxx - minx) / (double)nsamples;
            double dy = (maxy - miny) / (double)nsamples;
            for (int n = 0 ; n < nsamples ; ++n) {
                p.push_back(nvis::vec2(minx + (double)n*dx, miny + (double)n*dy));
            }
            break;
        }
        default: {
            std::cout << "unknown seeding pattern. using random sampling" << std::endl;
            seeding = 0;
            generate_seeds(p, seeding, nsamples);
        }
    }
}

std::string me;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Compute Poincare plot of circular restricted 3-body problem\n"
            << "             and color-code the result by winding number around L3.\n"
            << "             The computation is carried out for the full circular restricted\n"
            << "             3-body problem, the full planar CR3BP, and the planar CR3BP\n"
            << "             without integration of the state transformation matrix\n"
            << "OPTIONS:\n"
            << " -h | --help                     Print this information\n"
            << " -o | --output <string>          Output file name\n"
            << " -r | --res <int> (x2)           Resolution\n"
            << " -p | --maxp <int>               Number of iterations\n"
            << " -b | --bounds <float> (x4)      Computation bounds\n"
            << " -e | --eps <float>              Integration precision\n"
            << " -n | --samples <int>            Number of orbits to compute\n"
            << " -C <float>                      C constant\n"
            << " -m | --mu <float>               mu constant\n"
            << " -s | --seed <int> [<float>x2]   Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal, 4: single, 5: bounds)\n"
            << " -v | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

template<typename Map, typename Tracker>
void do_compute(const std::vector<nvis::vec2>& seeds, const Map& pmap,
                const std::string& filename)
{
    typedef typename Map::return_type return_type;
    typedef typename Map::rhs_type    rhs_type;
    typedef Tracker                   tracker_type;
    
    int nseeds = seeds.size();
    
    float* raster = (float*)calloc(res[0]*res[1], sizeof(float));
    
    nvis::timer timer;
    int started=0, completed=0;
    int cmp, st;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<nseeds ; ++n) {
        
#pragma openmp atomic
            ++started;
            
            st = started;
            cmp = completed;
            
            std::ostringstream os;
            double dt = timer.elapsed();
            os << "\r" << cmp << "/" << nseeds << " orbits completed (" << (float)cmp*100/nseeds << "%) in "
               << dt << "s. (" << (float)(n+1)/dt << "Hz) / " << st << " orbits in progress        \r";
            std::cout << os.str() << std::flush;
            
            std::vector<return_type> orbit;
            tracker_type t(L3);
            try {
                pmap.map_and_track_complete(seeds[n], orbit, maxp, t);
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
                raster[m] = w;
            }
            
#pragma openmp atomic
            ++completed;
#pragma openmp atomic
            --started;
            
            st = started;
            cmp = completed;
            
            os.clear();
            os.str("");
            dt = timer.elapsed();
            os << "\r" << cmp << "/" << nseeds << " orbits completed (" << (float)cmp*100/nseeds << "%) in "
               << dt << "s. (" << (float)(n+1)/dt << "Hz) / " << st << " orbits pending        \r";
            std::cout << os.str() << std::flush;
            
        }
    }
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. (" << (float)nseeds/elapsed << " Hz)\n";
    
    size_t dims[] = {res[0], res[1]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, raster, nrrdTypeFloat, 2, dims);
    char* labels[] = { (char*)"X", (char*)"Xdot" };
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
    nrrdSave(filename.c_str(), nout, NULL);
    nrrdNuke(nout);
    std::cout << filename << " saved\n";
}

int main(int argc, char* argv[])
{
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    bool bounds_set = false;
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
    
    const double LARGE = std::numeric_limits<double>::max();
    
    std::cout << "Running computation with cr3bp...\n";
    rhs_type1 rhs1(C, mu);
    section_type1 section1(rhs1);
    section1.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section1.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type1 pmap1(rhs1, section1);
    pmap1.setPrecision(eps);
    std::ostringstream os;
    os << filename << "_cr3bp.nrrd";
    do_compute<map_type1, tracker_type1>(seeds, pmap1, os.str());
    
    std::cout << "Running computation with pcr3bp...\n";
    rhs_type2 rhs2(C, mu);
    section_type2 section2(rhs2);
    section2.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section2.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type2 pmap2(rhs2, section2);
    pmap2.setPrecision(eps);
    os.clear();
    os.str("");
    os << filename << "_pcr3bp.nrrd";
    do_compute<map_type2, tracker_type2>(seeds, pmap2, os.str());
    
    std::cout << "Running computation with pcr3bp_reduced...\n";
    rhs_type3 rhs3(C, mu);
    section_type3 section3(rhs3);
    section3.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section3.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type3 pmap3(rhs3, section3);
    pmap3.setPrecision(eps);
    os.clear();
    os.str("");
    os << filename << "_pcr3bp_reduced.nrrd";
    do_compute<map_type3, tracker_type3>(seeds, pmap3, os.str());
    
    return 0;
}
