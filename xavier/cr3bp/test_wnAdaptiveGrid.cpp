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
#include <data/adaptive_grid.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/poincare_map.hpp>
#include <graphics/colors.hpp>
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

nvis::ivec2 res(10,10);
nvis::ivec2 sres(10, 10);
std::string filename;
double eps, C, mu, K, tol=0.1;
int maxdepth=4; //Default
int maxp = 100;
int which_rhs=2;
int verbose = 0;
size_t nthreads;
//JE - L3 point
nvis::vec2 L3(-1.00001053340710, 0);
//nvis::vec2 center=L3;
nvis::bbox2 bounds, sbounds;

struct value_traits {
    static const double invalid = 1000.;
    static const double unknown = 0.;
};

const double invalid_value = value_traits::invalid;

//Model Defs
//typedef pcr3bp_reduced                                          rhs_type;
//typedef dp5wrapper<double, 4>                                 odesolver_type;
//typedef planar_section_reduced<rhs_type, 4>                     section_type;
typedef nvis::fixed_vector<double, 42>              vec42;
typedef nvis::fixed_vector<double, 6>               vec6;
typedef xavier::dp5wrapper<double, 42>              odesolver_type;
typedef orbital::cr3bp                      rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>        section_type;
typedef orbital::AngleTracker<42, 1>                tracker_type;
typedef poincare_map<rhs_type, odesolver_type, section_type>    map_type;

//Grid Defs
typedef xavier::AdaptiveGrid<double, value_traits>              adaptive_grid_type;
typedef adaptive_grid_type::Node                                node_type;
typedef adaptive_grid_type::id_type                             id_type;
typedef adaptive_grid_type::cell_data_type                      cell_data_type;
typedef adaptive_grid_type::vertex_container_type               vertex_container_type;
typedef adaptive_grid_type::cell_container_type                 cell_container_type;
//Orbit Defs
typedef std::vector<nvis::vec2>                                 curve_type;
typedef std::pair<double, curve_type>                           orbit_type;

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

std::string me;
bool save_coord = false;
void printUsageAndExit(const std::string& what)
{
    if (what.size()) {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
            << "USAGE: " << me << " [options]\n"
            << "DESCRIPTION: Compute continuous approximation of winding number in circular\n"
            << "             restricted 3-body problem using convexity criterion with\n"
            << "             winding number and adaptive refinement.\n"
            //<< "             The computation is carried out for the planar CR3BP\n"
            //<< "             without integration of the state transformation matrix for\n"
            //<< "             increased efficiency."
            << "OPTIONS:\n"
            << " -h  | --help                     Print this information\n"
            << " -o  | --output <string>          Output file name\n"
            << " -sc | --savecoord                Save orbit coordinates to file\n"
            << " -r  | --res <int> (x2)           Resolution\n"
            << " -sr | --sres <int> (x2)          Sampling resolution\n"
            << " -d  | --maxd <int>               Maximum subdivision depth\n"
            << " -p  | --maxp <int>               Number of iterations\n"
            << " -b  | --bounds <float> (x4)      Plotting bounds\n"
            << " -c  | --center <float> (x2)      Poloidal rotation center\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -t  | --tolerance <float>        Relative tolerance on convexity criterion\n"
            << " -C  <float>                      C constant\n"
            << " -m  | --mu <float>               mu constant\n"
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

/// Singular winding number - test only x-xdot value right now - full version in maps/map_analysis.hpp
double winding_number(curve_type& curve, const nvis::vec2& x, const map_type& pmap,
                      const nvis::vec2& center)
{
    typedef map_type::return_type return_type;
    
    curve.clear();
    std::vector<return_type> orb;
    tracker_type tracker(center);
    try {
        pmap.map_and_track_complete<tracker_type>(x, orb, maxp, tracker);
    } catch(...) {}
    
    if (orb.empty() || (int)orb.size() < maxp) {
        return invalid_value;
    }
    
    //Compute winding number manually (could also call tracker.getWindingFactors()
    double w = 2.*M_PI*orb.size()/orb.back().delta_theta[0];
    if (std::isnan(w) || std::isinf(w)) {
        w = invalid_value;
    }
    
    curve.resize(orb.size());
    for (int i=0 ; i<orb.size() ; ++i) {
        curve[i] = orb[i].x;
    }
    return w;
}

struct OrbitInfo {
    OrbitInfo() : seed_id(0)
    {
        orbit.first = invalid_value;
        orbit.second.clear();
    }
    OrbitInfo(const id_type& id) : seed_id(id)
    {
        orbit.first = invalid_value;
        orbit.second.clear();
    }
    id_type     seed_id;
    orbit_type  orbit;
};

/// Sample the seeds that need to be filled from
void sample_grid(adaptive_grid_type& grid, const std::vector<id_type>& seeds,
                 const map_type& pmap, bool subres = false, nvis::vec2 center = nvis::vec2(0,0))
{
    typedef map_type::return_type  return_type;
    
    int nseeds = seeds.size();
    
    std::vector<OrbitInfo> orbits[nthreads];
    
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
            
            orbits[th_id].push_back(OrbitInfo(seeds[n]));
            OrbitInfo&  the_info  = orbits[th_id].back();
            orbit_type& the_orbit = the_info.orbit;
            curve_type& the_curve = the_orbit.second;
            
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
            
            //This runs the map for maxp iters
            nvis::vec2 seed = grid.getVertex(seeds[n]);
            the_orbit.first = winding_number(the_curve, seed, pmap, center);
            
            #pragma omp atomic
            --started;
            #pragma omp atomic
            ++completed;
            
            if (the_orbit.first != invalid_value)
#pragma openmp atomic
                ++success;
                
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
    
    timer.restart();
    for (int n=0 ; n<nthreads ; ++n) {
        for (int i=0 ; i<orbits[n].size() ; ++i) {
            const OrbitInfo&  info   = orbits[n][i];
            const id_type&    id     = info.seed_id;
            const double&     w      = info.orbit.first;
            const curve_type& curve  = info.orbit.second;
            
            grid.setVertexValue(id, w);
            if (subres) {
                grid.insertCellValue(grid.getVertex(id), w);
            }
            for (int k=0 ; k<curve.size() ; ++k) {
                try {
                    grid.insertCellValue(curve[k], w);
                } catch(...) {}
            }
        }
    }
    elapsed = timer.elapsed();
    std::cout << "Data structure update took " << elapsed << " s. (" << (float)nseeds/elapsed << " Hz)\n";
}

//Check a cell needs refinement as it is not satisfying convexity in winding number
bool check_cell(const adaptive_grid_type& grid, const id_type& cell_id)
{
    const cell_data_type& data = grid.getCellData(cell_id);
    //Note this assumes a single winding number value - Not true in general
    std::vector<double> corner_values(4);
    corner_values[0] = grid.getVertexValue(cell_id);
    corner_values[1] = grid.getVertexValue(cell_id + id_type(1,0,0));
    corner_values[2] = grid.getVertexValue(cell_id + id_type(1,1,0));
    corner_values[3] = grid.getVertexValue(cell_id + id_type(0,1,0));
    double min = *std::min_element(corner_values.begin(), corner_values.end());
    double max = *std::max_element(corner_values.begin(), corner_values.end());
    double eps = tol*(max-min);
    min -= eps;
    max += eps;
    typedef cell_data_type::const_iterator iterator_type;
    for (iterator_type it=data.begin() ; it!=data.end() ; ++it) {
        if (it->second<min || it->second>max) {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    //Jupiter-Europa system example
    //mu = 2.528017705e-5;
    //C = 3.000;
    
    //Earth-Moon Problem cases
    mu = 1.21505714306e-2;
    C = 2.96;
    nvis::vec2 center(0,0);
    
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
        } else if (arg == "-c" || arg == "--center") {
            if (i >= argc-2) {
                printUsageAndExit("missing poloidal center");
            }
            center[0] = atof(argv[++i]);
            center[1] = atof(argv[++i]);
        } else if (arg == "-r" || arg == "--res") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
        } else if (arg == "-sr" || arg == "--sres") {
            if (i >= argc-2) {
                printUsageAndExit("missing sampling resolution");
            }
            sres[0] = atoi(argv[++i]);
            sres[1] = atoi(argv[++i]);
        } else if (arg == "-d" || arg == "--maxd") {
            if (i == argc-1) {
                printUsageAndExit("missing maximum refinement depth");
            }
            maxdepth = atoi(argv[++i]);
        } else if (arg == "-e" || arg == "--eps") {
            if (i == argc-1) {
                printUsageAndExit("missing integration precision");
            }
            eps = atof(argv[++i]);
        } else if (arg == "-t" || arg == "--tolerance") {
            if (i == argc-1) {
                printUsageAndExit("missing convexity tolerance");
            }
            tol = atof(argv[++i]);
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
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (!bounds_set || filename=="none") {
        printUsageAndExit("");
    }
    
    //nvis::vec2 p1(-mu, 0);
    //nvis::vec2 p2(1-mu, 0);
    
    adaptive_grid_type grid(bounds, sres[0], sres[1]);
    std::vector<id_type> seed_ids;
    
    nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::cerr << nthreads << " threads available for computation\n";
    
    const double LARGE = std::numeric_limits<double>::max();
    //std::cout << "Running computation with pcr3bp_reduced...\n";
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds().min() = nvis::vec2(-LARGE, -LARGE);
    section.bounds().max() = nvis::vec2(LARGE, LARGE);
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    
    for (int depth=0 ; depth<=maxdepth ; ++depth) {
        std::cout << "sampling map at depth " << depth << '\n';
        
        seed_ids.clear();
        grid.getUndefinedVertices(seed_ids);
        
        sample_grid(grid, seed_ids, pmap, false);
        std::vector<id_type> leaves_id;
        grid.getLeaves(leaves_id);
        
        // determine empty leaves
        std::vector<id_type> empty;
        for (int i=0 ; i<leaves_id.size() ; ++i) {
            if (grid.getCellData(leaves_id[i]).empty()) {
                empty.push_back(leaves_id[i]);
            }
        }
        
        // create new seeds at their center
        std::vector<id_type> new_seeds(empty.size());
        for (int i=0 ; i<empty.size() ; ++i) {
            new_seeds[i] = grid.sub_id(empty[i]) + id_type(1,1,0);
        }
        
        // compute additional orbits
        sample_grid(grid, new_seeds, pmap, true, center);
        
        // determine which cells must be refined
        if (depth == maxdepth) {
            break;
        }
        
        std::vector<id_type> _refine[nthreads];
        int nleaves = leaves_id.size();
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic, 1)
            for (int n=0 ; n<nleaves ; ++n) {
                int th_id = 0;
#ifdef _OPENMP
                th_id = omp_get_thread_num();
#endif
                
                if (!check_cell(grid, leaves_id[n])) {
                    _refine[th_id].push_back(leaves_id[n]);
                }
            }
        }
        std::vector<id_type> to_refine;
        for (int n=0 ; n<nthreads ; ++n) {
            std::copy(_refine[n].begin(), _refine[n].end(), std::back_inserter(to_refine));
        }
        
        // refine identified cells
        for (int i=0 ; i<to_refine.size() ; ++i) {
            grid.refine(to_refine[i]);
        }
        
        std::cout << to_refine.size() << " cells have been refined.\n";
    }
    
    //Work with one winding number to debug
    typedef std::pair<nvis::vec2, double> point_type;
    
    // collect all the data points that have been computed
    std::vector<point_type> all_points;
    const vertex_container_type& vertices = grid.getVertices();
    vertex_container_type::const_iterator vcit;
    for (vcit=vertices.begin() ; vcit!=vertices.end() ; ++vcit) {
        all_points.push_back(point_type(grid.getVertex(vcit->first), vcit->second));
    }
    std::vector<id_type> leaves;
    grid.getLeaves(leaves);
    for (int i=0 ; i<leaves.size() ; ++i) {
        const cell_data_type& data = grid.getCellData(leaves[i]);
        cell_data_type::const_iterator cdit;
        for (cdit=data.begin() ; cdit!=data.end() ; ++cdit) {
            all_points.push_back(point_type(cdit->first, cdit->second));
        }
    }
    std::vector<nvis::bbox2> all_cells(leaves.size());
    for (int i=0 ; i<leaves.size() ; ++i) {
        all_cells[i] = grid.cellBounds(leaves[i]);
    }
    
    std::ostringstream os;
    os << filename << "-points.nrrd";
    float* _data = (float*)calloc(3*all_points.size(), sizeof(float));
    for (int i=0 ; i<all_points.size() ; ++i) {
        _data[3*i  ] = all_points[i].first[0];
        _data[3*i+1] = all_points[i].first[1];
        _data[3*i+2] = all_points[i].second;
    }
    Nrrd* nout = nrrdNew();
    size_t sz[] = {3, all_points.size()};
    nrrdWrap_nva(nout, _data, nrrdTypeFloat, 2, sz);
    nrrdSave(os.str().c_str(), nout, NULL);
    nrrdNuke(nout);
    
    os.clear();
    os.str("");
    os << filename << "-cells.nrrd";
    _data = (float*)calloc(4*all_cells.size(), sizeof(float));
    for (int i=0 ; i<all_cells.size() ; ++i) {
        _data[4*i  ] = all_cells[i].min()[0];
        _data[4*i+1] = all_cells[i].min()[1];
        _data[4*i+2] = all_cells[i].max()[0];
        _data[4*i+3] = all_cells[i].max()[1];
    }
    nout = nrrdNew();
    sz[0] = 4;
    sz[1] = all_cells.size();
    nrrdWrap_nva(nout, _data, nrrdTypeFloat, 2, sz);
    nrrdSave(os.str().c_str(), nout, NULL);
    nrrdNuke(nout);
    
    return 0;
}
