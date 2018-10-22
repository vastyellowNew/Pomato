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
#include <vis/integral_curve.hpp>

// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>

// map API
#include <maps/section.hpp>
#include <maps/DP45wrapper.hpp>
#include <cr3bp/cr3bp.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

// RHS and ODE solver parameters
double eps, C, mu, K, tmax, xmin, xmax;
int    nsamples;

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

typedef plane                section_type;
typedef integral_curve<42>   int_type;

// *******************************    PARAMETERS     *******************************
std::string me;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Visualize PCR3BP trajectories\n"
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -b  | --bounds <float> x 2       Seeding bounds (along x axis)\n"
            << " -n  | --nsamples <int>           Number of orbits\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -t  | --tmax                     Max integration time\n"
            << " -o  | --output <string>          Output file\n";
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
    
    std::string filename = "none";
    xmin = -1.5;
    xmax = -1;
    nsamples = 100;
    tmax = 100.;
    
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
            if (i >= argc-2) {
                printUsageAndExit("missing seed");
            }
            xmin = atof(argv[++i]);
            xmax = atof(argv[++i]);
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
        } else if (arg == "-n" || arg == "--nsamples") {
            if (i == argc-1) {
                printUsageAndExit("missing number of samples");
            }
            nsamples = atoi(argv[++i]);
        } else if (arg == "-t" || arg=="--tmax") {
            if (i == argc-1) {
                printUsageAndExit("missing integration time");
            }
            tmax = atof(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (filename == "none") {
        printUsageAndExit("");
    }
    
    const double LARGE = std::numeric_limits<double>::max();
    rhs_type rhs(C, mu);
    section_type section(rhs);
    
    typedef std::vector<nvis::vec3> orbit;
    std::vector<orbit> all_orbits(nsamples);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int i=0 ; i<nsamples ; ++i) {
            std::ostringstream os;
            os << '\r' << i*100/nsamples << '\r';
            std::cout << os.str();
            double u = xmin + (float)i/(float)(nsamples-1)*(xmax-xmin);
            nvis::vec2 x(u, 0);
            vec42 y = section.unproject(x);
            int_type solver(y);
            try {
                solver.advance(rhs_type, tmax);
            } catch(...) {}
            orbit_type o = all_orbits[i];
            double T = solver.t_max();
            for (double t=0 ; t<T ; t+=tmax/500.) {
                o.push_back(solver(t));
            }
        }
    }
    int npts = 0;
    for (int i=0 ; i<all_orbits.size() ; ++i) {
        npts += all_orbits[i].size();
    }
    
    std::cout << '\n';
    std::fstream out(filename.c_str(), std::ios::out);
    out << "# vtk DataFile Version 2.0\n"
        << "PCR3BP computed by " << argv[0] << '\n'
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << npts << " float\n";
    for (int i=0 ; i<all_orbits.size() ; ++i) {
    
    }
    
    
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
    std::vector<int> periods;
    best_periods(periods, chain, maxp, 5, euclidean_metric);
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
        } catch(xavier::index_step_size_underflow& e) {
            std::cerr << "\t\t**step size underflow** caught during index computation\n";
            std::cerr << "\t\t  error occured at " << e.where << '\n';
            std::cerr << "\t\t  norm at that location = " << e.step << '\n';
            std::cerr << "\t\t  angle at that location = " << e.theta << '\n';
        } catch(xavier::invalid_secant_value& e) {
            std::cerr << "\t\t**invalid secant value** caught during index computation\n";
            std::cerr << "\t\t  error occurred between " << e.x0 << " and " << e.x1 << '\n';
            std::cerr << "\t\t  corresponding vectors: " << e.v0 << " and " << e.v1 << '\n';
            std::cerr << "\t\t  angle at that location = " << e.theta << '\n';
        }
    }
    
    return 0;
}
