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

//Data
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
//Map API
#include <maps/DP45wrapper.hpp>
#include <maps/section.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/poincare_map.hpp>
//cr3bp
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/multipleAngleTracker.hpp>
//Display
#include <graphics/colors.hpp>
#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

nvis::fixed_vector<size_t, 2> res;
int maxi, nsamples;
std::string filename;
double eps, C, mu, K;
int seeding;
nvis::vec2 x0;


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

using namespace xavier;

nvis::bbox2 bounds, sbounds;

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
    const double& minx = sbounds.min()[0];
    const double& miny = sbounds.min()[1];
    const double& maxx = sbounds.max()[0];
    const double& maxy = sbounds.max()[1];

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
            << "DESCRIPTION: compute Poincare plot of circular restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h | --help                     Print this information\n"
            << " -o | --output <string>          Output file name\n"
            << " -r | --resolution <int> (x2)    Resolution\n"
            << " -i | --iterations <int>         Number of iterations\n"
            << " -b | --bounds <float> (x4)      Computation bounds\n"
            << " -e | --eps <float>              Integration precision\n"
            << " -n | --samples <int>            Number of sample trajectories\n"
            << " -C <float>                      C constant\n"
            << " -m | --mu <float>               mu constant\n"
            << " -s | --seed <int> [<float>x2]   Seeding pattern (0: random, 1: x-axis, 2: y-axis, 3: diagonal, 4: single, 5: bounds)\n"
            << " -v | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

nvis::fvec3 b2y(double u)
{
    const nvis::fvec3 blue(0,0,1);
    const nvis::fvec3 yellow(1,1,0);
    return (1.-u)*blue + u*yellow;
}

int main(int argc, char* argv[])
{
    me = argv[0];
    bounds.min() = bounds.max() = nvis::vec2(0,0);
    C = 3.000;
    mu = 2.528017705e-5;
    res = nvis::vec2(512, 512);
    eps = 1.0e-7;
    nsamples = 100;
    maxi = 1000;
    filename = "none";
    int verbose = 0;

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                printUsageAndExit("missing output name");
            }
            filename = argv[++i];
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
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
        } else if (arg == "-i" || arg == "--iterations") {
            if (i == argc-1) {
                printUsageAndExit("missing number of iterations");
            }
            maxi = atof(argv[++i]);
        } else if (arg == "-n" || arg == "--samples") {
            if (i == argc-1) {
                printUsageAndExit("missing number of samples");
            }
            nsamples = atoi(argv[++i]);
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing precision");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-s" || arg == "--seeding") {
            if (i == argc-1) {
                printUsageAndExit("missing seeding info");
            }
            seeding = atoi(argv[++i]);
            if (seeding == 4) {
                if (i >= argc-2) {
                    printUsageAndExit("missing seeding location");
                }
                x0[0] = atof(argv[++i]);
                x0[1] = atof(argv[++i]);
            } else if (seeding >= 5) {
                if (i >= argc-4) {
                    printUsageAndExit("missing seeding bounds");
                }
                sbounds.min()[0] = atof(argv[++i]);
                sbounds.min()[1] = atof(argv[++i]);
                sbounds.max()[0] = atof(argv[++i]);
                sbounds.max()[1] = atof(argv[++i]);
            }
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument" + arg);
        }
    }
    if (nvis::all(bounds.min() == bounds.max()) || filename == "none") {
        printUsageAndExit("");
    }

    if (seeding < 5) {
        sbounds = bounds;
    } else {
        seeding -= 4;
    }

    std::cout << "parameters: resolution = " << res << '\n';
    std::cout << bounds << std::endl;

    if (nsamples == -1) {
        nsamples = res[0];
    }
    if (maxi == -1) {
        maxi = res[0];
    }

    srand48(987654321);
    std::vector< nvis::vec2 > seeds;
    if (seeding == 4) {
        nsamples = 1;
        seeds.resize(1);
        seeds[0] = x0;
    } else {
        generate_seeds(seeds, seeding, nsamples);
    }

    float* hits = (float*)calloc(3 * res[0] * res[1], sizeof(float));
    if (hits == 0) {
        std::cerr << "Unable to allocate target image\n";
        exit(-1);
    }
    unsigned int count = 0;
    float pct, last_pct = 0;

    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds() = bounds;

    nvis::timer _timer;

    size_t nb_threads = 1;
#ifdef _OPENMP
    nb_threads = omp_get_max_threads();
#endif
    std::cerr << "there are " << nb_threads << " threads available\n";

    std::vector<std::string> thread_status(nb_threads);
    for (int i=0 ; i<nb_threads ; ++i) {
        thread_status[i] = "idle";
    }
    //Storage - this assumes you have 1 winding number per orbit (hard to compute in CR3BP)
    std::vector<std::list<std::pair<double, std::list<nvis::vec2> > > > th_orbits(nb_threads);

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < nsamples ; ++n) {

            nvis::vec2 x = seeds[n];
            nvis::vec3 color = b2y((x[0]-sbounds.min()[0])/sbounds.size()[0]);
            print_color(color, pos(x), hits);

            section_type s(section);
            s.seed = x;

            map_type _map(rhs, s);
            _map.setPrecision(eps);
            if (seeding == 4) {
                _map.verbose() = (verbose==2);
            }

            int th_id = 0;
#ifdef _OPENMP
            th_id = omp_get_thread_num();
#endif
            nvis::timer _t;
            if (verbose>=1) {
                std::ostringstream os;
                os << "\t" << th_id << ": seeding at " << x << '\n';
                std::cerr << os.str();
            }

            std::ostringstream os;
            os << "integrating from " << x;
            thread_status[th_id] = os.str();

            th_orbits[th_id].push_back(std::pair<double, std::list<nvis::vec2> >());
            std::vector< vec2 > __hits;
            try {
                _map.map(x, __hits, maxi);
            } catch (...) {
            }

            std::vector<nvis::vec2> vec;
            vec.push_back(x);
            std::copy(__hits.begin(), __hits.end(), std::back_inserter(vec));
            xavier::metric_type m;
            th_orbits[th_id].back().first = xavier::best_period(vec, maxi/4, m);//**PROBLEM!**
            th_orbits[th_id].back().second.push_back(x);
            std::copy(__hits.begin(), __hits.end(), std::back_inserter(th_orbits[th_id].back().second));

            os.clear();
            os.str("");
            os << "idle";
            thread_status[th_id] = os.str();

            if (verbose>=1) {
                std::ostringstream os;
                os << "\t" << th_id << ": completed in " << _t.elapsed() << '\n';
                std::cerr << os.str();
            }

            if (verbose>=3) {
                std::ostringstream os;
                os << "\nMap iterates from " << x << ":\n";
                for (int i=0 ; i<__hits.size() ; ++i) {
                    os << '\t' << i << ": " << __hits[i] << " (" << nvis::norm(__hits[i]-x) << ")\n";
                }
                std::cerr << os.str();
            }

            #pragma omp atomic
            ++count;

            int c = count;
            pct = (float)c / (float)nsamples * 100.;
            if (true) {
                std::ostringstream os;
                float dt = _timer.elapsed();
                os << '\r' << c << " curves computed in " << dt << " s. ("
                   << pct << "%, " << (float)n/dt << " Hz)              \r";
                std::cout << os.str() << std::flush;
                last_pct = pct;
            }
            if (verbose==2) {
                std::ostringstream os;
                os << "progress:\n";
                for (int i=0 ; i<nb_threads ; ++i) {
                    os << "#" << i << ": " << thread_status[i] << "\n";
                }
                os << '\n';
                std::cerr << os.str();
            }
        }
    }
    double elapsed = _timer.elapsed();
    std::cout << "\nintegration of " << nsamples << " orbits completed in "
              << elapsed << " seconds (" << (double)nsamples/elapsed << "Hz)\n";

    std::vector<double> all_periods;
    typedef std::list<std::pair<double, std::list<nvis::vec2> > > orbit_list;
    for (int i=0 ; i<nb_threads ; ++i) {
        const orbit_list& orb = th_orbits[i];
        for (orbit_list::const_iterator it=orb.begin() ; it!=orb.end() ; ++it) {
            all_periods.push_back(it->first);
        }
    }
    std::sort(all_periods.begin(), all_periods.end());

    std::cerr << "min period = " << all_periods.front() << '\n';
    std::cerr << "max period = " << all_periods.back() << '\n';
    // control points
    std::vector<double> cps;
#if 0
    size_t stride = all_periods.size()/8;
    cps.push_back(all_periods.front());
    cps.push_back(all_periods[1*stride]);
    cps.push_back(all_periods[2*stride]);
    cps.push_back(all_periods[3*stride]);
    cps.push_back(all_periods[4*stride]);
    cps.push_back(all_periods[5*stride]);
    cps.push_back(all_periods[6*stride]);
    cps.push_back(all_periods[7*stride]);
    cps.push_back(all_periods.back());
#else
    double dx = sbounds.size()[0]/8.;
    double minx = sbounds.min()[0];
    cps.push_back(minx);
    cps.push_back(minx+dx);
    cps.push_back(minx+2*dx);
    cps.push_back(minx+3*dx);
    cps.push_back(minx+4*dx);
    cps.push_back(minx+5*dx);
    cps.push_back(minx+6*dx);
    cps.push_back(minx+7*dx);
    cps.push_back(minx+8*dx);
#endif

    // scolor scale
    std::vector<nvis::fvec3> scale;
    scale.push_back(midnight);
    scale.push_back(blue);
    scale.push_back(cyan);
    scale.push_back(green);
    scale.push_back(yellow);
    scale.push_back(orange);
    scale.push_back(red);
    scale.push_back(magenta);
    scale.push_back(white);
    xavier::fixed_color_map<double> cmap(cps, scale);

    nvis::bbox2 plot_bounds;

    for (int i=0 ; i<nb_threads ; ++i) {
        const orbit_list& orb = th_orbits[i];
        for (orbit_list::const_iterator it=orb.begin() ; it!=orb.end() ; ++it) {
            std::ostringstream os;
            os << "orbit_" << it->second.front()[0] << ".csv";
            std::fstream f(os.str().c_str(), std::ios::out);
            nvis::vec3 col = cmap(it->second.front()[0]);
            f << "# curve seeded at " << it->second.front() << '\n';
            for (std::list<nvis::vec2>::const_iterator jit=it->second.begin() ; jit!=it->second.end() ; ++jit) {
                print_color(col, pos(*jit), hits);
                plot_bounds.add(*jit);
                f << (*jit)[0] << "," << (*jit)[1] << '\n';
            }
            f.close();
        }
    }

    std::cerr << "plot bounds = " << plot_bounds << '\n';

    Nrrd* nout = nrrdNew();
    size_t size[3] = {3, res[0], res[1]};
    if (nrrdWrap_nva(nout, hits, nrrdTypeFloat, 3, size) ||
            nrrdSave(filename.c_str(), nout, NULL)) {
        std::cout << "ERROR while exporting FTLE file: " << biffGetDone(NRRD)
                  << std::endl;
    }
    nrrdNuke(nout);

    return 0;
}
