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

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <cr3bp/cr3bp.hpp>

// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>

// map API
#include <maps/DP45wrapper.hpp>
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/index.hpp>
#include <maps/poincare_map.hpp>

// cr3bp
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

// RHS and ODE solver parameters
double eps, C, mu, K;


template<typename Map>
struct step_wrapper {
    step_wrapper(const Map& pmap, size_t p) : _map(pmap), _p(p) {}
    nvis::vec2 operator()(const nvis::vec2& x) const
    {
        return _map.map(x, _p) - x;
    }
    
    const Map& _map;
    size_t _p;
};

typedef pcr3bp_reduced                                         rhs_type;
typedef fixed_matrix<double, 4>                                vec_type;
typedef dp5wrapper<double, 4>                                  odesolver_type;
typedef planar_section_reduced<rhs_type, 4>                    section_type;
typedef poincare_map<rhs_type, odesolver_type, section_type>   map_type;
typedef step_wrapper<map_type>                                 step_type;

const double invalid_double = std::numeric_limits<double>::max();

bool xavier::record_newton_steps = true;
std::vector<nvis::vec2> xavier::newton_steps;


template<typename Step>
double rotation_angle(const vec_type& x0, const vec_type& v0,
                      const vec_type& x1, const vec_type& v1,
                      const Step& step, double dx,
                      double max_theta=M_PI/8.)
{

    if (nvis::norm(v0) < dx || nvis::norm(v1) < dx) {
        index_step_size_underflow e("unable to track continous vector rotation");
        e.where = 0.5*(x0 + x1);
        throw e;
    }
    
    double theta = xavier::signed_angle(v0, v1);
    
    if (fabs(theta) <= max_theta) {
        return theta;
    } else if (nvis::norm(x1 - x0) <= dx) {
        // only bad things will come out of trying to proceed in such cases...
        xavier::index_step_size_underflow e("unable to track continous vector rotation");
        e.where = 0.5*(x0 + x1);
        e.step = nvis::norm(x1 - x0);
        e.theta = fabs(theta);
        
        std::cerr << "found an ambiguous point at " << e.where << ", angle = " << e.theta << '\n';
        throw e;
    }
    double u = secant(v0, v1);
    if (std::isnan(u) || std::isinf(u)) {
        std::cerr << "secant method returned NaN" << std::endl;
        xavier::invalid_secant_value e("invalid secant value");
        e.x0 = x0;
        e.x1 = x1;
        e.v0 = v0;
        e.v1 = v1;
        e.theta = theta;
        throw e;
    }
    vec_type x = (1-u)*x0 + u*x1;
    vec_type v = step(x);
    std::cerr.precision(12);
    std::cerr << x[0] << "," << x[1] << "," << v[0] << "," << v[1] << '\n';
    
    return rotation_angle<Step>(x0, v0, x, v, step, dx) +
           rotation_angle<Step>(x, v, x1, v1, step, dx);
}


template<typename Step>
double edge_rotation(const Step& step, const nvis::vec2& x0, const nvis::vec2& x1,
                     double lmin, int n)
{

    std::cerr << "\n\nedge rotation...\n";
    
    double theta = 0;
    if (lmin > 0) {
        double du = 1./(double)n;
        std::vector<nvis::vec2> x(n+1), v(n+1);
        double u = du;
        
        for (int i=0; i<=n ; ++i, u+=du) {
            x[i] = (1.-u)*x0 + u*x1;
            v[i] = step(x[i]);
        }
        for (int i=0 ; i<n ; ++i) {
            theta += rotation_angle<Step>(x[i], v[i], x[i+1], v[i+1], step, lmin);
        }
    } else {
        nvis::vec2 v0 = step(x0);
        nvis::vec2 v1 = step(x1);
        theta = signed_angle(v0, v1);
    }
    return theta;
}

double dx = 1.0e-12;
template<typename Step>
int cell_index(const Step& step, const nvis::vec2& x0, const nvis::vec2& x2)
{
    nvis::vec2 x1(x2[0], x0[1]);
    nvis::vec2 x3(x0[0], x2[1]);
    std::vector<nvis::vec2> y;
    nvis::vec2 x[4] = { x0, x1, x2, x3 };
    double theta = 0;
    for (int i=0 ; i<4 ; ++i) {
        std::cerr << "edge #" << i << '\n';
        try {
            theta += edge_rotation(step, x[i], x[(i+1)%4], dx, 2);
        } catch(xavier::index_step_size_underflow& e) {
            y.push_back(e.where);
        }
    }
    if (y.size()) {
        std::cerr << y.size() << " ambiguous locations detected\n";
        std::copy(y.begin(), y.end(), std::ostream_iterator<nvis::vec2>(std::cerr, ", "));
        for (int i=0 ; i<y.size() ; ++i) {
            std::cerr << "norm at " << y[i] << " = " << nvis::norm(step(y[i])) << '\n';
        }
        if (y.size() == 2) {
            xavier::ambiguous_index e("ambiguous index");
            e.x0 = y[0];
            e.x1 = y[1];
            throw e;
        } else {
            xavier::index_step_size_underflow e("index step size underflow");
            e.where = y[0];
            throw e;
        }
    }
    
    return 0.5*theta/M_PI;
}

double gauss(double x, double s, double x0)
{
    double a = 1/(s*sqrt(2.*M_PI));
    return a*exp(-(x-x0)*(x-x0)/(2*s*s));
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
            << "DESCRIPTION: Index computation in restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -b  | --bounds <float> x 4       Bounds of square region to analyze\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --maxp                     Max considered period\n"
            << " -dx <float>                      Lower bound on sampling steps\n"
            << " -o  | --output <string>          Output file for regular samples\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    
    me = argv[0];
    eps = 1.0e-8;
    int verbose = 0;
    nvis::bbox2 bounds;
    bool seed_set = false;
    int maxp = 15;
    std::string filename;
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
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit("missing bounds");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            bounds_set = true;
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
    if (!seed_set) {
        printUsageAndExit("");
    }
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.precision() = eps;
    metric<double, 2> euclidean_metric;
    
    fixpoint fp;
    nvis::bbox2 bounds;
    bounds.min() = x0 - nvis::vec2(0.01, 0.01);
    bounds.max() = x0 + nvis::vec2(0.01, 0.01);
    
    // compute orbit at seed point
    std::vector<nvis::vec2> chain;
    pmap.map(x0, chain, 10*maxp);
    std::cout << "orbit at " << x0 << ":\n";
    std::copy(chain.begin(), chain.end(), std::ostream_iterator<nvis::vec2>(std::cout, ", "));
    std::cout << '\n';
    
    // determine period at seed location
    std::vector<double> p2d(maxp);
    for (int i=0 ; i<maxp ; ++i) {
        p2d[i] = average_distance(chain, i+1, euclidean_metric);
    }
    {
        double* data = (double*)calloc(500, sizeof(double));
        for (int i=0 ; i<maxp ; ++i) {
            int p = i+1;
            for (int j=0 ; j<500 ; ++j) {
                double x = 1. + (float)j/50.;
                for (int k=1 ; k<=p ; ++k) {
                    double x0 = (double)p/(double)k;
                    data[j] += gauss(x, p2d[i], x0);
                }
            }
        }
        size_t sz[] = { 500 };
        Nrrd* nout = nrrdNew();
        nrrdWrap_nva(nout, data, nrrdTypeDouble, 1, sz);
        nrrdSave("winding_pdf.nrrd", nout, NULL);
    }
    std::vector<int> periods;
    best_periods(periods, chain, maxp, maxp, euclidean_metric);
    std::cout << "best periods at " << x0 << " are: ";
    std::copy(periods.begin(), periods.end(), std::ostream_iterator<int>(std::cout, "  "));
    std::cout << '\n';
    
    size_t nthreads = 1;
#if _OPENMP
    nthreads = omp_get_num_threads();
#endif
    
    std::cerr << "there are " << nthreads << " threads available\n";
    
    std::sort(periods.begin(), periods.end());
    for (int i=0 ; i<periods.size() ; ++i) {
        const int p = periods[i];
        std::cout << "computing index for period " << p << '\n';
        step_type step(pmap, p);
        
        if (file_set) {
            double* data = (double*)calloc(2*200*200, sizeof(double));
            nvis::vec2 dx = bounds.size() / nvis::vec2(200, 200);
            #pragma omp parallel
            {
                #pragma omp for schedule(dynamic,1)
                for (int n=0 ; n<200*200 ; ++n) {
                    int i = n%200;
                    int j = n/200;
                    nvis::vec2 x = bounds.min() + nvis::vec2(i, j)*dx;
                    nvis::vec2 f = step(x);
                    data[2*n] = f[0];
                    data[2*n+1] = f[1];
                    if (!(n%10)) {
                        std::ostringstream os;
                        os << "\rn=" << n;
                        std::cerr << os.str();
                    }
                }
            }
            std::cerr << '\n';
            
            std::ostringstream os;
            os << filename << "_p=" << p << ".nrrd";
            size_t s[] = {2, 200, 200};
            double spc[] = {airNaN(), dx[0], dx[1]};
            Nrrd* nout = nrrdNew();
            nrrdWrap_nva(nout, data, nrrdTypeDouble, 3, s);
            nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
            nrrdSave(os.str().c_str(), nout, NULL);
        }
        
        try {
            int idx = cell_index<step_type>(step, bounds.min(), bounds.max());
            std::cout << "\tindex = " << idx << '\n';
        } catch(xavier::ambiguous_index& e) {
            std::cerr << "\t\t**ambiguous index** caught during index computation\n";
            std::cerr << "\t\t  two ambiguous locations are " << e.x0 << " and " << e.x1 << '\n';
            std::cerr << "\t\t  associated norms are " << nvis::norm(step(e.x0)) << " and " << nvis::norm(step(e.x1)) << '\n';
        } catch(xavier::index_step_size_underflow& e) {
            std::cerr << "\t\t**step size underflow** caught during index computation\n";
            std::cerr << "\t\t  error occured at " << e.where << '\n';
            std::cerr << "\t\t  norm at that location = " << e.step << '\n';
            std::cerr << "\t\t  angle at that location = " << e.theta << "\n\n";
        } catch(xavier::invalid_secant_value& e) {
            std::cerr << "\t\t**invalid secant value** caught during index computation\n";
            std::cerr << "\t\t  error occurred between " << e.x0 << " and " << e.x1 << '\n';
            std::cerr << "\t\t  corresponding vectors: " << e.v0 << " and " << e.v1 << '\n';
            std::cerr << "\t\t  angle at that location = " << e.theta << "\n\n";
        }
    }
    
    return 0;
}
