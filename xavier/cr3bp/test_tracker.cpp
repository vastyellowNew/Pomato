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
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/poincare_map.hpp>
#include <maps/newton.hpp>
#include <maps/fixpoints.hpp>
#include <maps/index.hpp>
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

struct MyTracker {
    MyTracker() : _center(-1.00001053340710,0), _path(), _theta(0) {}
    
    void initialize(const vec42& v)
    {
        _path.clear();
        _path.push_back(nvis::vec4(v[0], v[1], v[3], v[4]));
        _last = nvis::vec2(v[0], v[3]) - _center;
        _theta = 0;
        _angles.clear();
        _angles.push_back(0);
    }
    
    double operator()(const vec42& v)
    {
        nvis::vec2 cur = nvis::vec2(v[0], v[3]) - _center;
        _theta += signed_angle(_last, cur);
        _path.push_back(nvis::vec4(v[0], v[1], v[3], v[4]));
        _angles.push_back(_theta);
        _last = cur;
        return _theta;
    }
    
    void mark_crossing() {}
    
    nvis::vec2 _center, _last;
    std::list<nvis::vec4> _path;
    std::list<double> _angles;
    double _theta;
};

typedef plane                                                         section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type>  map_type;
typedef map_type::return_type                                         return_type;

// *******************************    PARAMETERS     *******************************
std::string me;
std::string filename;
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
            << " -x  | --seed <float> x 2         Seed point\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -p  | --maxp <int>               Max considered period\n"
            << " -o  | --output <string>          Output base name\n"
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
    filename = "none";
    
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
        } else if (arg == "-x" || arg == "--seed") {
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
                printUsageAndExit("missing output base name");
            }
            filename = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            if (i == argc-1) {
                printUsageAndExit("missing verbose flag");
            }
            verbose = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (!seed_set || filename=="none") {
        printUsageAndExit("");
    }
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    
    // compute orbit at seed point
    std::vector<return_type> chain;
    MyTracker tracker;
    pmap.map_and_track_complete<MyTracker>(x0, chain, maxp, tracker);
    
    std::string name1, name2;
    name1 = filename + "_all.nrrd";
    name2 = filename + "_plot.nrrd";
    
    float* data = (float*)calloc(5*tracker._path.size(), sizeof(float));
    std::list<nvis::vec4>::const_iterator itv = tracker._path.begin();
    std::list<double>::const_iterator ita = tracker._angles.begin();
    for (int n=0; itv!=tracker._path.end() && ita!=tracker._angles.end() ; ++itv, ++ita) {
        data[n++] = (*itv)[0];
        data[n++] = (*itv)[1];
        data[n++] = (*itv)[3];
        data[n++] = (*itv)[4];
        data[n++] = *ita;
    }
    size_t dims[] = {5, tracker._path.size()};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, dims);
    char* labels[] = { (char*)"x,y,xdot,ydot,theta", (char*)"t" };
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
    nrrdSave(name1.c_str(), nout, NULL);
    nrrdNuke(nout);
    
    dims[0] = 3;
    dims[1] = chain.size();
    data = (float*)calloc(3*chain.size(), sizeof(float));
    for (int i=0 ; i<chain.size() ; ++i) {
        data[3*i  ] = chain[i].x[0];
        data[3*i+1] = chain[i].x[1];
        data[3*i+2] = chain[i].delta_theta;
    }
    nout = nrrdNew();
    nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, dims);
    labels[0] = (char*)"x,xdot,theta";
    labels[1] = (char*)"iteration";
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoLabel, labels);
    nrrdSave(name2.c_str(), nout, NULL);
    nrrdNuke(nout);
    
    return 0;
}
