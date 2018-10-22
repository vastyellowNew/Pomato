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
#include <maps/section.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/poincare_map.hpp>
#include <maps/newton.hpp>
#include <maps/fixpoints.hpp>
#include <topology/invariant_manifold.hpp>

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

// RHS and ODE solver parameters
double eps, C, mu, K;

typedef orbital::cr3bp                          rhs_type;
typedef nvis::fixed_vector<double, 42>          vec42;
typedef nvis::fixed_vector<double, 6>           vec6;
typedef nvis::fixed_matrix<double, 6>           mat6;
typedef xavier::dp5wrapper<double, 42>          odesolver_type;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

class plane : public orbital::section<double, 42, 2> {
public:
    bool verbose;
    lvec_type seed;
    
    plane(const rhs_type& rhs) : _rhs(rhs), verbose(false), _counter(0), _dist(0) {}
    plane(const plane& other) : _rhs(other._rhs), verbose(false), _counter(0), _dist(0) {}
    
    const lbox_type& bounds() const
    {
        return _bounds;
    }
    
    lbox_type& bounds()
    {
        return _bounds;
    }
    
    std::pair<lvec_type, lmat_type> project(const gvec_type& x) const
    {
        std::pair<lvec_type, lmat_type> r;
        r.first[0] = x[0];
        r.first[1] = x[3];
        r.second(0,0) = x[6];
        r.second(0,1) = x[9];
        r.second(1,0) = x[24];
        r.second(1,1) = x[27];
        return r;
    }
    
    gvec_type unproject(const lvec_type& x) const
    {
        gvec_type gv(0);
        gv[0] = x[0];
        gv[3] = x[1];
        gv[4] = _rhs.yd(x[0], x[1]);
        // initialize flow map Jacobian to Identity
        gv[6] = 1;
        gv[13] = 1;
        gv[20] = 1;
        gv[27] = 1;
        gv[34] = 1;
        gv[41] = 1;
        return gv;
    }
    
    value_type distance(const gvec_type& x) const
    {
        return x[1];
    }
    
    gvec_type distance_first_partials(const gvec_type& x) const
    {
        gvec_type gv(0);
        gv[1] = 1;
        return gv;
    }
    
private:
    rhs_type    _rhs;
    lbox_type   _bounds;
    
    mutable size_t      _counter;
    mutable double      _dist;
    mutable gvec_type   _last;
};

typedef plane section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type> map_type;

bool xavier::record_newton_steps = true;
std::vector<nvis::vec2> xavier::newton_steps;

// *******************************    PARAMETERS     *******************************
std::string me;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Newton fixed point search test in circular restricted 3-body problem\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -g  | --guess <float> x 2        First guess\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --maxp                     Max considered period\n"
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
    nvis::vec2 x0;
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
    pmap.setPrecision(eps);
    metric<double, 2> euclidean_metric;
    
    fixpoint fp;
    nvis::bbox2 bounds;
    bounds.min() = x0 - nvis::vec2(0.01, 0.01);
    bounds.max() = x0 + nvis::vec2(0.01, 0.01);
    
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
    }
    
    // compute orbit at seed point
    std::vector<nvis::vec2> chain;
    pmap.map(x0, chain, 10*maxp);
    std::cerr << "orbit at " << x0 << ":\n";
    std::copy(chain.begin(), chain.end(), std::ostream_iterator<nvis::vec2>(std::cerr, ", "));
    std::cerr << '\n';
    
    // determine period at seed location
    std::vector<int> periods;
    best_periods(periods, chain, maxp, 5, euclidean_metric);
    std::cerr << "best periods at " << x0 << " are: ";
    std::copy(periods.begin(), periods.end(), std::ostream_iterator<int>(std::cerr, "  "));
    std::cerr << '\n';
    
    int p = *std::min_element(periods.begin(), periods.end());
    std::cerr << "Starting Newton at " << x0 << " for period " << p << '\n';
    
    //Run a newton line search for the fixed point (Minimizing f(x) = ||P^p(x)-x||)
    bool found = simple_newton<map_type>(pmap, euclidean_metric, bounds, x0, 5, p,
                                         fp, chain, 1.0e-5, 1.0e-5, true, 50, true);
                                         
    return 0;
    
}
