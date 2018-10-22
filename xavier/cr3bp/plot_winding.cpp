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

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

// map API
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
using namespace xavier;
using namespace orbital;

nvis::ivec2 res(1000,1000);
std::string filename;
double eps, C, mu, K;
int seeding;
int maxp = 100;
int which_rhs=2;
int verbose = 0;
size_t nthreads;
nvis::vec2 L3(-1.00001053340710, 0);
nvis::vec2 center=L3;
nvis::bbox2 bounds, sbounds;

typedef cr3bp                                                       rhs_type1;
typedef pcr3bp                                                      rhs_type2;
typedef pcr3bp_reduced                                              rhs_type3;
typedef dp5wrapper<double, 42>                                      odesolver_type1;
typedef dp5wrapper<double, 20>                                      odesolver_type2;
typedef dp5wrapper<double, 4>                                       odesolver_type3;
typedef planar_section<rhs_type1, 6, 42>                            section_type1;
typedef planar_section<rhs_type2, 4, 20>                            section_type2;
typedef planar_section_reduced<rhs_type3, 4>                        section_type3;
typedef poincare_map<rhs_type1, odesolver_type1, section_type1>     map_type1;
typedef poincare_map<rhs_type2, odesolver_type2, section_type2>     map_type2;
typedef poincare_map<rhs_type3, odesolver_type3, section_type3>     map_type3;
typedef AngleTracker<42, 3>                                         tracker_type1;
typedef AngleTracker<20, 2>                                         tracker_type2;
typedef AngleTracker<4,  2>                                         tracker_type3;
typedef planar_rotation_winding_number                              winding_type;

#define PLOT_QUIET 1

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

void generate_seeds(std::vector< nvis::vec2 >& p, int seeding, int nsamples)
{
    const double& minx = sbounds.min()[0];
    const double& miny = sbounds.min()[1];
    const double& maxx = sbounds.max()[0];
    const double& maxy = sbounds.max()[1];
    
    p.resize(nsamples);
    switch (seeding) {
        case 0: {
            orbital::cr3bp rhs(C, mu);
            for (int n = 0 ; n < nsamples ;) {
                p[n] = nvis::vec2(minx + drand48()*(maxx - minx),
                                  miny + drand48()*(maxy - miny));
                try {
                    rhs.yd(p[n][0], p[n][1]);
                } catch(...) {
                    continue;
                }
                ++n;
            }
            
            break;
        }
        case 1: {
            double dx = (maxx - minx) / (double)nsamples;
            for (int n = 0 ; n < nsamples ; ++n) {
                p[n] = nvis::vec2(minx + (double)n*dx, 0.5*(miny + maxy));
            }
            break;
        }
        case 2: {
            double dy = (maxy - miny) / (double)nsamples;
            for (int n = 0 ; n < nsamples ; ++n) {
                p[n] = nvis::vec2(0.5*(minx + maxx), miny + (double)n*dy);
            }
            break;
        }
        case 3: {
            double dx = (maxx - minx) / (double)nsamples;
            double dy = (maxy - miny) / (double)nsamples;
            for (int n = 0 ; n < nsamples ; ++n) {
                p[n] = nvis::vec2(minx + (double)n*dx, miny + (double)n*dy);
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
bool save_coord = false;
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
            << " -h  | --help                     Print this information\n"
            << " -o  | --output <string>          Output file name\n"
            << " -sc | --savecoord                Save orbit coordinates to file\n"
            << " -r  | --res <int> (x2)           Resolution\n"
            << " -p  | --maxp <int>               Number of iterations\n"
            << " -b  | --bounds <float> (x4)      Plotting bounds\n"
            << " -sb | --seedbounds <float> (x4)  Seeding bounds\n"
            << " -c  | --center <float> (x2)      Poloidal rotation center\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -rhs <int>                       Right hand side (0: cr3bp, 1: pcr3bp, 2: pcr3bp_reduced)\n"
            << " -n  | --samples <int>            Number of orbits to compute\n"
            << " -C  <float>                      C constant\n"
            << " -m  | --mu <float>               mu constant\n"
            << " -s  | --seed <int>               Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal)\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

template<typename T>
void plot_on_raster(const std::vector<T>& orbit, float* raster, float value)
{
    for (int i=0 ; i<orbit.size() ; ++i) {
        nvis::vec2 x = orbit[i].x;
        nvis::vec2 c = (x-bounds.min())/bounds.size()*nvis::vec2(res);
        if (c[0]<0 || c[0]>=res[0] || c[1]<0 || c[1]>=res[1]) {
            continue;
        }
        int m = (int)c[0] + (int)c[1]*res[0];
        raster[m] = value;
    }
}

template<typename Map, typename Tracker>
void do_compute(const std::vector<nvis::vec2>& seeds, const Map& pmap,
                const std::string& filename)
{
    typedef typename Map::return_type     return_type;
    typedef typename Map::rhs_type        rhs_type;
    typedef Tracker                       tracker_type;
    typedef std::vector<nvis::vec2>       curve_type;
    typedef std::pair<double, curve_type> orbit_type;
    
    int nseeds = seeds.size();
    
    float* raster;
    if (!save_coord) {
        raster = (float*)calloc(res[0]*res[1], sizeof(float));
    }
    float* failed_raster = (float*)calloc(res[0]*res[1], sizeof(float));
    
    std::vector<orbit_type> all_orbits[nthreads];
    
    nvis::timer timer;
    int started=0, completed=0, success=0;
    int cmp, st, suc;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<nseeds ; ++n) {
            int th_id = 0;
#ifdef _OPENMP
            th_id = omp_get_thread_num();
#endif
            
#pragma openmp atomic
            ++started;
            
            st = started;
            cmp = completed;
            suc = success;
            
            std::ostringstream os;
            double dt = timer.elapsed();
            os << "\r" << cmp << "/" << nseeds << " orbits completed (" << (float)cmp*100/nseeds << "%) in "
               << dt << "s. (" << (float)cmp/dt << "Hz) / " << st << " orbits in progress / "
               << suc << " (" << 100.*(float)suc/(float)cmp << "%) successful integrations         \r";
            std::cout << os.str() << std::flush;
            
            std::vector<return_type> orbit;
            tracker_type t(L3);
            try {
                pmap.map_and_track_complete(seeds[n], orbit, maxp, t);
            } catch(...) {}
#pragma openmp atomic
            --started;
#pragma openmp atomic
            ++completed;
            
            if (orbit.empty()) {
                continue;
            }
            
            double w = 2.*M_PI*orbit.size()/orbit.back().delta_theta[0];
            if (std::isnan(w) || std::isinf(w)) {
                w = 1000;
            }
            orbit.push_back(return_type());
            orbit.back().x = seeds[n];
            if (orbit.size() < maxp+1) {
                plot_on_raster<return_type>(orbit, failed_raster, orbit.size());
            } else if (!save_coord) {
#pragma openmp atomic
                ++success;
                
                plot_on_raster<return_type>(orbit, raster, w);
            } else {
#pragma openmp atomic
                ++success;
                
                all_orbits[th_id].push_back(orbit_type());
                orbit_type& orb = all_orbits[th_id].back();
                orb.first = w;
                orb.second.resize(orbit.size());
                for (int i=0 ; i<orbit.size() ; ++i) {
                    orb.second[i] = orbit[i].x;
                }
            }
            
            st = started;
            cmp = completed;
            suc = success;
            
            os.clear();
            os.str("");
            dt = timer.elapsed();
            os << "\r" << cmp << "/" << nseeds << " orbits completed (" << (float)cmp*100/nseeds << "%) in "
               << dt << "s. (" << (float)cmp/dt << "Hz) / " << st << " orbits in progress / "
               << suc << " (" << 100.*(float)suc/(float)cmp << "%) successful integrations         \r";
            std::cout << os.str() << std::flush;
        }
    }
    double elapsed = timer.elapsed();
    std::cout << "\nComputation took " << elapsed << " s. (" << (float)nseeds/elapsed << " Hz)\n";
    std::cout << success << " orbits (" << 100.*(float)success/(float)nseeds << "%) were successfully integrated\n";
    
    {
        size_t dims[] = {res[0], res[1]};
        Nrrd* nout = nrrdNew();
        nrrdWrap_nva(nout, failed_raster, nrrdTypeFloat, 2, dims);
        char* labels[] = { (char*)"X", (char*)"Xdot" };
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
        if (nrrdSave("failed.nrrd", nout, NULL)) {
            std::cerr << "Error exporting failed.nrrd: " << biffGetDone(NRRD) << '\n';
        }
        nrrdNuke(nout);
    }
    
    if (!save_coord) {
        size_t dims[] = {res[0], res[1]};
        Nrrd* nout = nrrdNew();
        nrrdWrap_nva(nout, raster, nrrdTypeFloat, 2, dims);
        char* labels[] = { (char*)"X", (char*)"Xdot" };
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
        if (nrrdSave(filename.c_str(), nout, NULL)) {
            std::cerr << "Error exporting " << filename << ": " << biffGetDone(NRRD) << '\n';
        }
        nrrdNuke(nout);
    } else {
        size_t ncurves = 0;
        for (int i=0 ; i<nthreads ; ++i) {
            ncurves += all_orbits[i].size();
        }
        
        float* data = (float*)calloc(2*ncurves*(maxp+2), sizeof(float));
        int n=0;
        int stride = 2*(maxp+2);
        for (int i=0 ; i<nthreads ; ++i) {
            for (int j=0 ; j<all_orbits[i].size() ; ++j, ++n) {
                const orbit_type& orb = all_orbits[i][j];
                data[n*stride] = orb.first;
                for (int k=0 ; k<orb.second.size() ; ++k) {
                    data[n*stride+2*(k+1)  ] = orb.second[k][0];
                    data[n*stride+2*(k+1)+1] = orb.second[k][1];
                }
            }
        }
        
        size_t dims[] = {2*(maxp+2), ncurves};
        Nrrd* nout = nrrdNew();
        nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, dims);
        char* labels[] = { (char*)"curveID", (char*)"W,0,x0,y0,x1,y1,..." };
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
        if (nrrdSave(filename.c_str(), nout, NULL)) {
            std::cerr << "Error exporting " << filename << ": " << biffGetDone(NRRD) << '\n';
        }
        nrrdNuke(nout);
    }
    
    std::cout << filename << " saved\n";
}

int main(int argc, char* argv[])
{
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    
    me = argv[0];
    eps = 1.0e-8;
    bool bounds_set = false;
    filename = "none";
    int nseeds = 100;
    int seeding = 1;
    sbounds.min() = sbounds.max() = nvis::vec2(0,0);
    
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
        } else if (arg == "-sb" || arg == "--seedbounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing seeding bounds");
            }
            sbounds.min()[0] = atof(argv[++i]);
            sbounds.min()[1] = atof(argv[++i]);
            sbounds.max()[0] = atof(argv[++i]);
            sbounds.max()[1] = atof(argv[++i]);
        } else if (arg == "-c" || arg == "--center") {
            if (i >= argc-2) {
                printUsageAndExit("missing poloidal center");
            }
            center[0] = atof(argv[++i]);
            center[1] = atof(argv[++i]);
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
        } else if (arg == "-rhs") {
            if (i == argc-1) {
                printUsageAndExit("missing right hand side id");
            }
            which_rhs = atoi(argv[++i]);
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
        } else if (arg == "-sc" || arg == "--savecoord") {
            save_coord = true;
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
    if (nvis::all(sbounds.size() == nvis::vec2(0,0))) {
        sbounds = bounds;
    }
    
    nvis::vec2 primary(-mu, 0);
    nvis::vec2 secondary(1-mu, 0);
    nvis::vec2 barycenter(0, 0);
    
    std::cerr << argv[0] << ": computing the winding number of " << nseeds << " orbits in "
              << '[' << bounds.min() << "->" << bounds.max() << "] over " << maxp << " iterations\n";
              
    std::vector<nvis::vec2> seeds(nseeds);
    generate_seeds(seeds, seeding, nseeds);
    
    nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cerr << nthreads << " threads available for computation\n";
    
    const double LARGE = std::numeric_limits<double>::max();
    
    if (which_rhs == 0) {
        std::cout << "Running computation with cr3bp...\n";
        rhs_type1 rhs1(C, mu);
        section_type1 section1(rhs1);
        section1.bounds().min() = nvis::vec2(-LARGE, -LARGE);
        section1.bounds().max() = nvis::vec2(LARGE, LARGE);
        map_type1 pmap1(rhs1, section1);
        pmap1.precision() = eps;
        do_compute<map_type1, tracker_type1>(seeds, pmap1, filename);
    } else if (which_rhs == 1) {
        std::cout << "Running computation with pcr3bp...\n";
        rhs_type2 rhs2(C, mu);
        section_type2 section2(rhs2);
        section2.bounds().min() = nvis::vec2(-LARGE, -LARGE);
        section2.bounds().max() = nvis::vec2(LARGE, LARGE);
        map_type2 pmap2(rhs2, section2);
        pmap2.precision() = eps;
        do_compute<map_type2, tracker_type2>(seeds, pmap2, filename);
    } else {
        std::cout << "Running computation with pcr3bp_reduced...\n";
        rhs_type3 rhs3(C, mu);
        section_type3 section3(rhs3);
        section3.bounds().min() = nvis::vec2(-LARGE, -LARGE);
        section3.bounds().max() = nvis::vec2(LARGE, LARGE);
        map_type3 pmap3(rhs3, section3);
        pmap3.precision() = eps;
        do_compute<map_type3, tracker_type3>(seeds, pmap3, filename);
    }
    return 0;
}
