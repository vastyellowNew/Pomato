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
#include <maps/poincare_map.hpp>
#include <maps/section.hpp>
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
// sampling parameters
nvis::ivec2 res;
int maxi, maxp;
nvis::bbox2 bounds;
// graphics parameters
bool print_friendly;

typedef orbital::cr3bp              rhs_type;
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
    ///Mirror theorem
    lvec_type mirror(const lvec_type& x) const
    {
        return lvec_type(x[0],-x[1]);
    }
    
    /// Symmetry test
    bool isSymmetric() const
    {
        return true;
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

// *******************************    PARAMETERS     *******************************
std::string me;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Test flow and Jacobian integration\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -t  | --tmax <float>             Integration time\n"
            << " -n  | --nseeds <int>             Number of seeds\n"
            << " -d <float>                       Stencil size\n"
            << " -b  | --bounds <float> (x4)      Computation bounds\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
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
    bounds.min() = bounds.max() = nvis::vec2(0,0);
    eps = 1.0e-8;
    maxi = 100;
    int verbose = 0;
    int nseeds = 1;
    double tmax = 10.;
    double h = 0.005;
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-n" || arg == "--nseeds") {
            if (i == argc-1) {
                printUsageAndExit("missing number of seeds");
            }
            nseeds = atoi(argv[++i]);
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
        } else if (arg == "-t" || arg == "--tmax") {
            if (i == argc-1) {
                printUsageAndExit("missing integration time");
            }
            tmax = atof(argv[++i]);
        } else if (arg == "-d") {
            if (i == argc-1) {
                printUsageAndExit("missing stencil size");
            }
            h = atof(argv[++i]);
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (nvis::all(bounds.min() == bounds.max()) || tmax==0 ) {
        printUsageAndExit("");
    }
    
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds() = bounds;
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    
    for (int i=0 ; i<nseeds ; ++i) {
        std::cout << "Testing integration from seed #" << i+1 << "/" << nseeds << '\n';
        nvis::vec2 x0 = bounds.min() + nvis::vec2(drand48(), drand48())*bounds.size();
        std::cout << "starting point = " << x0 << '\n';
        map_type::return_state r0(section.unproject(x0),0.0);
        std::cout << "Full starting point = " << r0.x << '\n';
        map_type::return_state r; //Full state & Jacobian Data structure
        r = pmap.flow_map(x0, tmax);
        //Test print full Jacobian
        std::cout << " dt = " << r.t << '\n';
        std::cout << "Full end position = " << r.x << '\n';
        std::cout << "Full end Jacobian = " << r.J << '\n';
        std::cerr << "Determinant = " << nvis::det(r.J) << '\n';
        
        vec42 y = r.getState();
        std::pair<vec2,mat2> pr = section.project(y);
        nvis::mat2 J = pr.second;
        std::cout << "end position = " << pr.first << '\n';
        std::cout << "end Jacobian = " << J << '\n';
        std::cerr << "Determinant = " << nvis::det(J) << '\n';
        
        
        r = pmap.flow_map(x0+nvis::vec2(h, 0), tmax);
        pr = section.project(r.getState());
        nvis::vec2 fr = pr.first;
        r = pmap.flow_map(x0-nvis::vec2(h, 0), tmax);
        pr = section.project(r.getState());
        nvis::vec2 fl = pr.first;
        r = pmap.flow_map(x0+nvis::vec2(0, h), tmax);
        pr = section.project(r.getState());
        nvis::vec2 ft = pr.first;
        r = pmap.flow_map(x0-nvis::vec2(0, h), tmax);
        pr = section.project(r.getState());
        nvis::vec2 fb = pr.first;
        nvis::vec2 dfdx = (fr - fl)/(2.*h);
        nvis::vec2 dfdy = (ft - fb)/(2.*h);
        nvis::mat2 Jcd;
        Jcd(0,0) = dfdx[0];
        Jcd(0,1) = dfdy[0];
        Jcd(1,0) = dfdx[1];
        Jcd(1,1) = dfdy[1];
        std::cerr << "Central difference Jacobian = " << Jcd << '\n';
        std::cerr << "Determinant = " << nvis::det(Jcd) << '\n';
        std::cerr << "Difference = " << nvis::norm(J-Jcd)/nvis::norm(Jcd) << '\n';
    }
    
    return 0;
    
}
