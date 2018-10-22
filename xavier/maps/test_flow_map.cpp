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

// display
#include <graphics/colors.hpp>
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>

// map API
#include <cr3bp/poincare_map.hpp>
#include <cr3bp/section.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/newton.hpp>
#include <maps/fixpoints.hpp>
#include <maps/invariant_manifold.hpp>

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
int nogfx, line_width, point_size, display_size[2], width, height;
bool print_friendly;

typedef orbital::cr3bp                          rhs_type;
typedef nvis::fixed_vector<double, 42>          vec42;
typedef nvis::fixed_vector<double, 6>           vec6;
typedef nvis::fixed_matrix<double, 6>           mat6;
typedef xavier::dp5wrapper<double, 42>          odesolver_type;

const double invalid_double = std::numeric_limits<double>::max();

using namespace xavier;

class plane : public section<double, 42, 2> {
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
            << "DESCRIPTION: Topological analysis of circular restricted 3-body problem\n"
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
    pmap.precision() = eps;
    
    for (int i=0 ; i<nseeds ; ++i) {
        std::cout << "Testing integration from seed #" << i+1 << "/" << nseeds << '\n';
        nvis::vec2 x0 = bounds.min() + nvis::vec2(drand48(), drand48())*bounds.size();
        std::cout << "starting point = " << x0 << '\n';
        map_type::return_type r;
        r = pmap.flow_map(x0, tmax);
        std::cout << "end position = " << r.x << '\n';
        std::cout << "end Jacobian = " << r.J << '\n';
        nvis::mat2 J = r.J;
        std::cerr << "Determinant = " << nvis::det(J) << '\n';
        
        r = pmap.flow_map(x0+nvis::vec2(h, 0), tmax);
        nvis::vec2 fr = r.x;
        r = pmap.flow_map(x0-nvis::vec2(h, 0), tmax);
        nvis::vec2 fl = r.x;
        r = pmap.flow_map(x0+nvis::vec2(0, h), tmax);
        nvis::vec2 ft = r.x;
        r = pmap.flow_map(x0-nvis::vec2(0, h), tmax);
        nvis::vec2 fb = r.x;
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
