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
#include <map>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET 1

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <cr3bp/cr3bp.hpp>
#include <math/rational.hpp>

// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>

// display
#include <graphics/colors.hpp>
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>

// map API
#include <maps/section.hpp>
#include <orbital/corrections.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/winding.hpp>
#include <maps/map_analysis.hpp>
#include <maps/definitions.hpp>
#include <maps/newton.hpp>
#include <maps/fixpoints.hpp>
#include <maps/poincare_map.hpp>
#include <topology/invariant_manifold.hpp>

// cr3bp
#include <cr3bp/planar_section.hpp>
#include <cr3bp/cr3bp.hpp>

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

// RHS and ODE solver parameters
double eps, C, mu, K;
nvis::vec2 center;
// sampling parameters
nvis::ivec2 res;
int maxi, maxp;
nvis::bbox2 bounds;
// graphics parameters
int nogfx, line_width, point_size, display_size[2], width, height;
bool print_friendly;

typedef orbital::cr3bp                                                rhs_type;
typedef nvis::fixed_vector<double, 42>          vec42;
typedef nvis::fixed_vector<double, 6>                vec6;
typedef nvis::fixed_matrix<double, 6>                mat6;
typedef xavier::dp5wrapper<double, 42>                odesolver_type;

const double invalid_double = std::numeric_limits<double>::max();

bool xavier::record_newton_steps = false;
std::vector<nvis::vec2> xavier::newton_steps;

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
    
    mutable size_t                _counter;
    mutable double                _dist;
    mutable gvec_type         _last;
};

typedef plane section_type;
typedef xavier::poincare_map<rhs_type, odesolver_type, section_type> map_type;

struct BasicTracker {
    BasicTracker(const nvis::vec2& center) : _center(center), _theta(0), _values() {}
    
    void initialize(const vec42& v)
    {
        _last = nvis::vec2(v[0], v[3]) - _center;
        _theta = 0;
        _values.clear();
    }
    
    std::vector<double> operator()(const vec42& v)
    {
        nvis::vec2 cur = nvis::vec2(v[0], v[3]) - _center;
        _theta += signed_angle(_last, cur);
        _last = cur;
        return std::vector<double>(1,_theta);
    }
    
    void mark_crossing()
    {
        _values.push_back(_theta);
    }
    
    nvis::vec2 _center, _last;
    std::vector<double> _values;
    double _theta;
};


// -------------------------
//
//                         Analysis
//
// -------------------------
typedef grid<double, 2>                plane_type;

plane_type*                                _plane;
nvis::ivec2                                _plane_res;
nvis::bbox2                                 _plane_bounds;
metric_type                                _plane_metric;
nvis::vec2                                _plane_spacing;
dataset_type*                                _dataset;

typedef rational_surface_found::edge_type                                segment_type;
typedef rational_surface_found::edge_point                                 edge_point;
typedef std::vector<nvis::vec2>                                                orbit_type;
typedef std::pair<orbit_type, nvis::fvec3>                                color_orbit_type;

// Less Than operator for 2D bounding boxes
struct Lt_bbox {
    nvis::lexicographical_order Lt;
    
    bool operator()(const nvis::bbox2& b0, const nvis::bbox2& b1) const
    {
        if (Lt(b0.min(), b1.min())) {
            return true;
        } else if (Lt(b1.min(), b0.min())) {
            return false;
        }
        return Lt(b0.max(), b1.max());
    }
};

typedef tagged<nvis::bbox2, Lt_bbox>                        p_box_type;

std::map<p_edge_type, double>                                        _edge_angles;
std::set<p_edge_type>                                                _unresolved_edges;
std::vector<tagged<p_cell_type>        >                                _cells_indices;
std::vector<tagged<p_box_type> >                                _boxes_indices;
std::vector<p_edge_type>                                        _failed_edges;
std::vector<std::vector<fixpoint> >                                _chains;
std::vector<std::vector<separatrix> >                                _separatrices;
std::vector<fixpoint>                                                _fixpoints;
std::map<double, rational_type>                                        _valid_rationals;

// -------------------------
//
//                         Display
//
// -------------------------
std::vector<color_orbit_type>                                        _orbits;
xavier::discrete_color_map<int>*                                _cmap;
xavier::map_analysis_param                                        _params, _debug_params;
std::vector<nvis::ivec2>                                        _saddle_cells, _center_cells;
nvis::vec2                                                        _last;
std::vector<nvis::vec2>                                                _problematic_seeds, _failed_newton;
std::vector<color_orbit_type>                                        _rational_surfaces;
std::vector<nvis::vec2>                                                _saddles;
std::vector<nvis::vec2>                                                _centers;
std::vector<fixpoint>                                                _unfiltered_fixed_points;
std::vector<nvis::vec2>                                                _all_seeds;
bool show_failed_edges, show_failed_newton;

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
            << " -r  | --res <int> (x2)           Sampling resolution\n"
            << " -i  | --iterations <int>         Number of iterations\n"
            << " -b  | --bounds <float> (x4)      Computation bounds\n"
            << " -e  | --eps <float>              Integration precision\n"
            << " -c  | --center <float> (x2)      Poloidal rotation center\n"
            << " -C <float>                       C constant\n"
            << " -mu <float>                      mu constant\n"
            << " -mp | --maxperiod <int>          Max considered period\n"
            << " -d  | --display <int> (x2)       Display resolution\n"
            << " -ps | --point <int>              Point size\n"
            << " -ls | --line <int>               Line width\n"
            << " -v  | --verbose <int>            Verbose mode (0: no, 1: basic, 2: full)\n";
    exit(1);
}

// *******************************    GRAPHICS     *******************************
int main_window;
nvis::vec2 wc;
void keyboard(unsigned char key, int x, int y);
void display();

void idle(void)
{
    // switch context back to main window after GLUI activity
    glutSetWindow(main_window);
}

inline void draw_circle(const nvis::vec2& c, double r, double dz=0)
{
    float dt = M_PI/20;
    glBegin( GL_TRIANGLE_FAN );
    glVertex3f(c[0], c[1], dz);
    for( float t = 0; t <=2*M_PI+dt; t += dt ) {
        glVertex3f( c[0] + r*cos(t), c[1] + r*sin(t), dz );
    }
    glEnd();
}

double current_size()
{
    nvis::vec3 min = GLUT_helper::world_coordinates(0,0);
    nvis::vec3 d = GLUT_helper::world_coordinates(width, 0) - min;
    d[0] = fabs(d[0]);
    d[1] = fabs(d[1]);
    d[2] = fabs(d[2]);
    return *std::max_element(&d[0], &d[3]);
}

void display_matrix(const GLfloat m[16])
{
    std::cerr << "(" << m[0] << ", " << m[1] << ", " << m[2] << ", " << m[3] << ")\n"
              << "(" << m[4] << ", " << m[5] << ", " << m[6] << ", " << m[7] << ")\n"
              << "(" << m[8] << ", " << m[9] << ", " << m[10] << ", " << m[11] << ")\n"
              << "(" << m[12] << ", " << m[13] << ", " << m[14] << ", " << m[15] << ")\n";
}

std::set<int> periods_to_show;
void draw(void)
{
    static int count = 0;
    
    glDisable(GL_DEPTH_TEST);
    
    if (!print_friendly) {
        glClearColor(0, 0, 0, 1);
    } else {
        glClearColor(1, 1, 1, 1);
    }
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw separatrices
    glEnable(GL_BLEND);
    if (!print_friendly) {
        glEnable(GL_LINE_SMOOTH);
    } else {
        glDisable(GL_LINE_SMOOTH);
    }
    glLineWidth(line_width);
    glBegin(GL_LINES);
    for (int i=0 ; i<_separatrices.size() ; ++i) {
        for (int j=0 ; j<_separatrices[i].size() ; ++j) {
            int s = j%4;
            if (!print_friendly) {
                if (s<2) {
                    glColor3f(0, 0, 1);
                } else {
                    glColor3f(1, 0, 0);
                }
            } else {
                glColor3f(0, 0, 0);
            }
            const separatrix& sep = _separatrices[i][j];
            const std::vector<nvis::vec2>& pts = sep.manifold;
            for (int n=0 ; n<pts.size()-1 ; ++n) {
                glVertex2f(pts[n][0], pts[n][1]);
                nvis::vec2 x = pts[n] + _plane_metric.displacement(pts[n], pts[n+1]);
                glVertex2f(x[0], x[1]);
            }
        }
    }
    glEnd();
    
    // draw fixpoints
    for (int i=0 ; i<_fixpoints.size() ; ++i) {
        fixpoint& fp = _fixpoints[i];
        if (fp.saddle) {
            glColor3f(1, 0, 0);
        } else {
            glColor3f(0, 1, 0);
        }
        //fp.pos = _plane_metric.modulo(fp.pos);
        double cur_sz = current_size();
        draw_circle(fp.pos, 0.002*cur_sz, 0.5);
    }
}

void display(void)
{
    GLUT_helper::setup_display(draw);
    glutSwapBuffers();
}

void guiCallback(int)
{
    display();
}

bool save_to_file;
void mouse(int button, int state, int x, int y)
{
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    
    std::ostringstream os;
    os << "Cursor at " << wc;
    glutSetWindowTitle(os.str().c_str());
    GLUT_helper::glut_helper_mouse(button, state, x, y);
    display();
}

bool arrow_pressed = false;
int __x0__, __y0__;
void keySpecial(int key, int x, int y)
{

    if (!arrow_pressed) {
        arrow_pressed = true;
        __x0__ = x;
        __y0__ = y;
        mouse(GLUT_MIDDLE_BUTTON, GLUT_DOWN, __x0__, __y0__);
    } else {
        switch (key) {
            case GLUT_KEY_LEFT: {
                __x0__ -= 1;
                break;
            }
            case GLUT_KEY_RIGHT: {
                __x0__ += 1;
                break;
            }
            case GLUT_KEY_UP: {
                __y0__ += 1;
                break;
            }
            case GLUT_KEY_DOWN: {
                __y0__ -= 1;
            }
            default: {
            }
        }
        nvis::vec3 _wc = GLUT_helper::world_coordinates(__x0__, __y0__);
        wc = nvis::vec2(_wc[0], _wc[1]);
        if (nvis::all(wc == _last)) {
            return;
        } else {
            _last = wc;
        }
        std::ostringstream os;
        os << "Cursor at " << wc;
        glutSetWindowTitle(os.str().c_str());
        
        GLUT_helper::glut_helper_motion(__x0__, __y0__);
    }
}

void keySpecialUp(int key, int x, int y)
{
    arrow_pressed = false;
    mouse(GLUT_MIDDLE_BUTTON, GLUT_UP, __x0__, __y0__);
}

void keyboard(unsigned char key, int x, int y)
{
    static double sensitivity = 1000;
    if(key == 'r') {
        GLUT_helper::resetCamera();
    }
    
    glutPostRedisplay();
}

// *******************************    INITIALIZATION    *******************************
size_t nbthreads;
static void init()
{
    nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    _plane_res = res;
    _plane_bounds = bounds;
    _plane_metric.bounds() = _plane_bounds;
    _plane_metric.periodic()[0] = false;
    _plane_metric.periodic()[1] = false;
    int npoints = res[0] * res[1];
    
    std::cerr << "plane resolution = " << res << std::endl;
    
    GLUT_helper::box = _plane_bounds;
    
    // allowable winding factors
    std::cerr << "valid rationals for period range [1, " << maxp << "]=";
    for (int pol_r=1 ; pol_r<=maxp ; ++pol_r) {
        for (int tor_r=1 ; tor_r<=3*maxp ; ++tor_r) {
            rational_type q(pol_r, tor_r);
            _valid_rationals[value<int, double>(q)] = q;
        }
    }
    for (std::map<double, rational_type>::const_iterator it=_valid_rationals.begin() ; it!=_valid_rationals.end() ; ++it) {
        std::cout << it->first << "(" << it->second << "), ";
    }
    std::cout << '\n';
    
    _plane = new plane_type(_plane_res, _plane_bounds);
    _dataset = new dataset_type(*_plane, orbit_data());
    _params.nb_iterations = maxi;
    _params.max_period = maxp;
    _params.the_metric = _plane_metric;
    _plane_spacing = _plane->spacing();
    
    _debug_params.verbose = true;
    _debug_params.the_metric = _plane_metric;
    
    std::cerr << "plane spacing = " << _plane_spacing << std::endl;
}

// ***************************    CELL ANALYSIS    ***************************

// per-cell period range
void compute_cell_periods(std::vector<p_cell_type>& cells)
{
    nvis::ivec2 cellres(_plane_res[0]-1, _plane_res[1]-1);
    int ncells = cellres[0]*cellres[1];
    
    std::vector<p_cell_type> cached_cells[nbthreads];
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<ncells ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            int i = n % cellres[0];
            int j = n / cellres[0];
            nvis::ivec2 cell_id(i,j);
            std::vector<int> periods;
            period_range(periods, *_dataset, cell_id, _valid_rationals);
            for (int k=0 ; k<periods.size() ; ++k) {
                cached_cells[thread_id].push_back(p_cell_type(cell_id, periods[k]));
            }
            std::ostringstream os;
            os << "\rcomputed period of " << n << " / " << ncells
               << " (" << 100*n/ncells << "%), " << cached_cells[thread_id].size() << " cached cells      \r"
               << std::flush;
            std::cout << os.str();
        }
    }
    std::cerr << "\nperiod range computation took " << _timer.elapsed() << '\n';
    
    for (int i=0 ; i<nbthreads ; ++i) {
        std::copy(cached_cells[i].begin(), cached_cells[i].end(),
                  std::back_inserter(cells));
    }
}

template<typename Map>
inline nvis::vec2 p_step(const nvis::vec2& x, const Map& map, int period)
{
    return map.map(x, period)-x;
}

template<typename Map>
void __edge_rotation(const Map& pmap, std::vector<p_edge_type>& failed,
                     const std::vector<p_edge_type>& edges, double lmin,
                     int nb_sub=4)
{

    typedef std::pair<p_edge_type, double> pair_type;
    std::vector<pair_type> cached_angles[nbthreads];
    std::vector<p_edge_type> cached_edges[nbthreads];
    nvis::timer _timer;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<edges.size() ; ++n) {
            int thread_id = 0;
            _debug_params.verbose = false;
            std::ostringstream os, os1;
            
            os1 << "\rcomputed " << n << " edge angles / "
                << edges.size() << " (" << 100*n/edges.size()
                << "%)          \r" << std::flush;
            std::cout << os1.str();
            
            const p_edge_type& e = edges[n];
            nvis::ivec2 i0 = e.object()[0];
            nvis::ivec2 i1 = e.object()[1];
            
            if (!(*_dataset)(i0).valid() || !(*_dataset)(i1).valid()) {
                continue;
            }
            int p = e.tag();
            nvis::vec2 v0 = p_vector((*_dataset)(i0), p, _plane_metric);
            nvis::vec2 v1 = p_vector((*_dataset)(i1), p, _plane_metric);
            if (_debug_params.verbose) {
                os << "v0 = " << v0 << ", v1 = " << v1 << std::endl;
            }
            
            map_type* amap = pmap.clone();
            
            double theta = 0;
            if (lmin > 0) {
                nvis::vec2 x0 = (*_plane)(i0);
                nvis::vec2 x1 = (*_plane)(i1);
                double du = 1./(double)(nb_sub+1);
                std::vector<nvis::vec2> x(nb_sub+2), v(nb_sub+2);
                double u = du;
                try {
                    for (int i=1; i<=nb_sub ; ++i, u+=du) {
                        x[i] = (1.-u)*x0 + u*x1;
                        v[i] = p_step(x[i], pmap, p);
                    }
                } catch(...) {
                    continue;
                }
                x[0] = x0;
                x[nb_sub+1] = x1;
                v[0] = v0;
                v[nb_sub+1] = v1;
                try {
                    for (int i=0 ; i<=nb_sub ; ++i) {
                        theta += adaptive_rotation_angle(x[i], v[i], x[i+1], v[i+1], *amap, p, lmin, _debug_params);
                    }
                    if (_debug_params.verbose) {
                        os << "total dtheta for edge = " << theta << std::endl;
                    }
                    cached_angles[thread_id].push_back(pair_type(e, theta));
                } catch (index_step_size_underflow& err) {
                    if (_debug_params.verbose) {
                        os << "step size underflow caught" << std::endl;
                    }
                    cached_edges[thread_id].push_back(e);
                } catch(...) {
                }
            } else {
                theta = signed_angle(v0, v1);
                if (_debug_params.verbose) {
                    os << "direct angle = " << theta << std::endl;
                }
                if (fabs(theta) < MAX_ANGLE_VARIATION) {
                    cached_angles[thread_id].push_back(std::pair<p_edge_type, double>(edges[n], theta));
                    if (_debug_params.verbose) {
                        os << "valid angle. done" << std::endl;
                    }
                } else {
                    if (_debug_params.verbose) {
                        os << "angle is too large." << std::endl;
                    }
                    cached_edges[thread_id].push_back(edges[n]);
                }
            }
            if (_debug_params.verbose) {
                std::cerr << os.str();
            }
        }
    }
    std::cerr << "processing of " << edges.size() << " edges took "
              << _timer.elapsed() << " s.                          \n";
              
    failed.clear();
    for (int i=0 ; i<nbthreads ; ++i) {
        for (int j=0 ; j<cached_angles[i].size() ; ++j) {
            const pair_type& p = cached_angles[i][j];
            _edge_angles[p.first] = p.second;
        }
        std::copy(cached_edges[i].begin(), cached_edges[i].end(),
                  std::back_inserter(failed));
    }
    std::cerr << "We were unable to process " << failed.size() << " edges out of "
              << edges.size() << " input edges ("
              << 100*failed.size()/edges.size() << "%)\n";
}

int __index(const nvis::ivec2& cell_id, int p,
            std::vector<tagged<p_box_type> >& idx_boxes,
            bool verbose=false)
{
    nvis::ivec2 pt[5];
    pt[0] = cell_id;
    pt[1] = cell_id + nvis::ivec2(1,0);
    pt[2] = cell_id + nvis::ivec2(1,1);
    pt[3] = cell_id + nvis::ivec2(0,1);
    pt[4] = cell_id;
    
    double theta = 0;
    bool valid = true;
    std::ostringstream os;
    for (int i=0 ; i<4 ; ++i) {
        edge_type e(pt[i], pt[i+1]);
        p_edge_type et(e, p);
        if (_unresolved_edges.find(et) != _unresolved_edges.end()) {
            //__multi_cell_index(cell_id, p, idx_boxes, verbose);
            return 0;
        }
        std::map<p_edge_type, double>::iterator it = _edge_angles.find(et);
        if (it == _edge_angles.end()) {
            // should not happen: all edges are initialized with invalid value
            if (verbose) {
                os << "edge " << pt[i] << "-" << pt[i+1]
                   << " was not found in database\n";
            }
            valid = false;
            break;
        }
        double dtheta = it->second;
        if (dtheta == invalid_double) {
            if (verbose) {
                os << "edge " << pt[i] << "-" << pt[i+1]
                   << " has an invalid angle\n";
            }
            valid = false;
            break;
        }
        if (verbose) {
            os << "edge " << pt[i] << "-" << pt[i+1]  << " has angle " << dtheta << '\n';
        }
        if (i>=2) {
            if (verbose) {
                os << "edge " << pt[i] << "-" << pt[i+1]
                   << " was computed the other way around\n";
            }
            dtheta *= -1;
        }
        theta += dtheta;
    }
    if (verbose) {
        os << "total angle is " << theta << std::endl;
        std::cerr << os.str();
    }
    
    if (valid) {
        theta /= 2.*M_PI;
        return lround(theta);
    } else {
        return 0;
    }
}


// per edge vector rotation
template<typename Map>
void compute_edge_rotation(const Map& pmap, const std::vector<p_cell_type>& cells)
{
    // assume cells sorted by period
    
    std::cerr << "there are " << cells.size() << " cells in input\n";
    
    std::set<p_edge_type> unique_edges;
    std::cerr << "identifying edges... " << std::flush;
    edge_type e;
    for (int n=0 ; n<cells.size() ; ++n) {
        nvis::ivec2 id = cells[n].object();
        int p = cells[n].tag();
        e = edge_type(id, id + nvis::ivec2(1,0));
        unique_edges.insert(p_edge_type(e, p));
        e = edge_type(id + nvis::ivec2(1,0), id + nvis::ivec2(1,1));
        unique_edges.insert(p_edge_type(e, p));
        e = edge_type(id + nvis::ivec2(1,1), id + nvis::ivec2(0,1));
        unique_edges.insert(p_edge_type(e, p));
        e = edge_type(id + nvis::ivec2(0,1), id);
        unique_edges.insert(p_edge_type(e, p));
    }
    std::cerr << ": " << unique_edges.size() << " edges found\n";
    
    for (std::set<p_edge_type>::const_iterator it=unique_edges.begin() ;
            it!=unique_edges.end() ; ++it) {
        // initial value is invalid angle
        _edge_angles[*it] = invalid_double;
    }
    std::vector<p_edge_type> all_edges(unique_edges.begin(), unique_edges.end());
    std::cerr << "edge map created\n";
    
    std::vector<p_edge_type> difficult_edges, problematic_edges;
    // easy
    // __edge_rotation(pmap, difficult_edges, all_edges, 0);
    // difficult
    double lmin = 0.99*std::min(_plane_spacing[0], _plane_spacing[1])/64.;
    __edge_rotation(pmap, problematic_edges, all_edges, lmin);
    // problematic
    lmin = 0.99*std::min(_plane_spacing[0], _plane_spacing[1])/2048;
    __edge_rotation(pmap, difficult_edges, problematic_edges, lmin);
    _unresolved_edges.insert(difficult_edges.begin(), difficult_edges.end());
}

void compute_cell_index(const std::vector<p_cell_type>& cells)
{
    std::vector<tagged<p_cell_type> > cached_indices[nbthreads];
    std::vector<tagged<p_box_type> >  cached_indexed_boxes[nbthreads];
    
    nvis::timer _timer;
    std::ostringstream os;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<cells.size() ; ++n) {
            int thread_id = 0;
            bool verbose = false;
            bool __verbose=false;
            
            nvis::ivec2 cell_id = cells[n].object();
            
            int p = cells[n].tag();
            int index = __index(cell_id, p, cached_indexed_boxes[thread_id]);
            if (index != 0) {
                cached_indices[thread_id].push_back(tagged<p_cell_type>(cells[n], index));
            }
            
            std::ostringstream os;
            os << "\rcomputed the index of " << n << " cells / "
               << cells.size() << " (" << 100*n/cells.size() << "%)         \r"
               << std::flush;
            std::cout << os.str();
        }
    }
    std::cerr << "\ncomputation of cell indices took " << _timer.elapsed() << " s.                          \n";
    
    int counter = 0;
    for (int i=0 ; i<nbthreads ; ++i) {
        std::copy(cached_indices[i].begin(), cached_indices[i].end(),
                  std::back_inserter(_cells_indices));
        std::copy(cached_indexed_boxes[i].begin(), cached_indexed_boxes[i].end(),
                  std::back_inserter(_boxes_indices));
    }
    counter = _cells_indices.size();
    std::cerr << "There are " << counter << " cells containing fixed points "
              << "(" << 100*counter/cells.size() << "%)\n";
    std::cerr << "these cells are: " << std::endl;
    for (int i=0 ; i<_cells_indices.size() ; ++i) {
        std::cerr << "id = " << _cells_indices[i].object().object() << ", period = " << _cells_indices[i].object().tag()
                  << ", index = " << _cells_indices[i].tag() << std::endl;
    }
    std::vector<int> nbsaddles(maxp), nbcenters(maxp);
    for (int i=0 ; i<_cells_indices.size() ; ++i) {
        int idx = _cells_indices[i].tag();
        int p = _cells_indices[i].object().tag();
        if (idx < 0) {
            ++nbsaddles[p-1];
            nvis::ivec2 id = _cells_indices[i].object().object();
            nvis::vec2 x = (*_plane)(id) + 0.5*_plane_spacing;
            _saddles.push_back(x);
        } else if (idx > 0) {
            ++nbcenters[p-1];
            nvis::ivec2 id = _cells_indices[i].object().object();
            nvis::vec2 x = (*_plane)(id) + 0.5*_plane_spacing;
            _centers.push_back(x);
        }
    }
    std::cerr << "statistics by period:\n";
    for (int i=0 ; i<maxp ; ++i) {
        std::cerr << "period " << i+1 << ": " << nbsaddles[i] << " saddles, "
                  << nbcenters[i] << " centers\n";
    }
}

// ***************************    FIXED POINTS    ***************************

template<typename InputIterator>
InputIterator epsilon_find(const InputIterator first, const InputIterator last,
                           const nvis::vec2& x, const metric_type& _metric, double eps)
{
    // return element corresponding to min distance under epsilon to
    // reference position
    std::map<double, InputIterator> norm_to_it;
    for (InputIterator i=first ; i!=last ; ++i) {
        const nvis::vec2& y = i->pos;
        norm_to_it[_metric.distance(x, y)] = i;
    }
    if (norm_to_it.begin()->first < eps) {
        return norm_to_it.begin()->second;
    }
    
    return last;
}

template<typename Map>
void __fixed_points(const Map& pmap, std::vector<fixpoint>& fixpoints,
                    std::map<p_cell_type, nvis::vec2>& cell_to_seed,
                    std::map<p_cell_type, nvis::vec2>& completed_cells,
                    bool prefound = false)
{

    std::map<p_cell_type, nvis::vec2> detected_cells;
    
    std::vector<p_cell_type> cells;
    typedef std::map<p_cell_type, nvis::vec2>::iterator        iterator;
    for (iterator it=cell_to_seed.begin() ; it!=cell_to_seed.end() ; ++it) {
        cells.push_back(it->first);
    }
    
    size_t nbcells = cells.size();
    std::cerr << "processing " << nbcells << " cells in input\n";
    
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<cells.size() ; ++n) {
        
            const p_cell_type& pcell = cells[n];
            nvis::ivec2 cell_id = pcell.object();
            int period = pcell.tag();
            
            std::cerr << "looking for fixed point in cell " << cell_id << std::endl;
            
            //First guess is the center of the cell (see compute_fixed_points())
            nvis::vec2 x = cell_to_seed[pcell];
            std::vector<nvis::vec2> chain;
            map_type* amap = pmap.clone();
            double jc = pmap.rhs().desired_hamiltonian();
            fixpoint fp;
            nvis::bbox2 bounds;
            bounds.min() = (*_plane)(cell_id);
            bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
            //Run a MultipleDim Newton search for the fixed point
            amap->setPrecision(1.0e-12);
            bool found = orbital::solveFixedPoints<map_type>(*amap, _plane_metric,
                         bounds, x, jc, 5, period, fp, fixpoints, 1.e-8, true, 10);
                         
            if (found) {
#pragma openmp atomic
                completed_cells[pcell] = fp.pos;
#pragma openmp atomic
                //fixpoints.push_back(fp);
                for (int i=1 ; i<(int)fixpoints.size() ; ++i) {
                    p_cell_type id(_plane->local_coordinates(fixpoints[i].pos).first, period);
                    iterator it = cell_to_seed.find(id);
                    if (it!=cell_to_seed.end())
#pragma openmp atomic
                        it->second = fixpoints[i].pos; // this cell was already identified
                        
#pragma openmp atomic
                    detected_cells[id] = fixpoints[i].pos; // this cell is new
                }
            }
            
#pragma openmp atomic
            int nb_found = fixpoints.size();
            
            std::ostringstream os;
            os << "\rprocessed " << n << " candidate cells from " << cells.size()
               << " (" << 100*n/cells.size() << "%), found "
               << nb_found << " fixed points     \r" << std::flush;
            std::cout << os.str();
        }
    }
    std::cerr << "\nthis round of fixpoint extraction took " << _timer.elapsed() << '\n';
    
    cell_to_seed.clear();
    int nb_found_already = 0;
    for (iterator it=detected_cells.begin() ; it != detected_cells.end() ; ++it) {
        iterator jt = completed_cells.find(it->first);
        if (jt != completed_cells.end()) {
            ++nb_found_already;
        } else {
            cell_to_seed[it->first] = it->second;
        }
    }
    
    std::cerr << completed_cells.size() << " fixed points were found\n";
    std::cerr << cell_to_seed.size() << " additional fixed points were discovered\n";
}

std::vector<std::vector<nvis::vec2> > _tmp_seeds;
template<typename Map>
void __fixed_points_from_boxes(const Map& pmap, std::vector<fixpoint>& fixpoints,
                               std::map<p_cell_type, nvis::vec2>& cell_to_seed,
                               std::map<p_cell_type, nvis::vec2>& completed_cells)
{

    // per-thread cache
    typedef std::pair<p_cell_type, nvis::vec2> solution_type;
    std::vector< std::vector<solution_type> >         completed(nbthreads);
    std::vector< std::vector<solution_type> >          discovered(nbthreads);
    std::vector< std::vector<fixpoint> >                 fps(nbthreads);
    std::vector< std::vector<nvis::vec2> >                seeds(nbthreads);
    
    std::map<p_cell_type, nvis::vec2> detected_cells;
    
    std::cerr << "processing " << _boxes_indices.size() << " boxes in input\n";
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<_boxes_indices.size() ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            const p_box_type b = _boxes_indices[n].object();
            int period = b.tag();
            const nvis::bbox2& bounds = b.object();
            bool verbose = false;
            
            //First guess is center of the cell
            nvis::vec2 x = bounds.center();
            seeds[thread_id].push_back(x);
            
            std::vector<fixpoint> chain;
            map_type* amap = pmap.clone();
            double jc = pmap.rhs().desired_hamiltonian();
            fixpoint fp;
            //Run a Multiple Dim Newton search for the fixed point
            amap->setPrecision(1.0e-12);
            bool found = orbital::solveFixedPoints<map_type>(*amap, _plane_metric,
                         bounds, x, jc, 5, period, fp, chain, 1.e-8, true, 10);
            if (verbose) {
                std::cerr << "Newton search was " << (found ? "successful" : "unsuccessful") << '\n';
            }
            
            if (found) {
                nvis::ivec2 cell_id = _plane->local_coordinates(chain[0].pos).first;
                p_cell_type pcell(cell_id, period);
                
                completed[thread_id].push_back(solution_type(pcell, fp.pos));
                fps[thread_id].push_back(fp);
                
                for (int i=1 ; i<chain.size() ; ++i) {
                    p_cell_type id(_plane->local_coordinates(chain[i].pos).first, period);
                    discovered[thread_id].push_back(solution_type(id, chain[i].pos));
                }
                std::ostringstream os;
                os << fp << " found" << std::endl;
                std::cerr << os.str();
            }
        }
    }
    std::cerr << "this round of fixpoint extraction took " << _timer.elapsed() << '\n';
    
    int nb_found = 0, nb_found_already = 0;
    for (int thid=0 ; thid<nbthreads ; ++thid) {
        std::copy(seeds[thid].begin(), seeds[thid].end(), std::back_inserter(_all_seeds));
        nb_found += completed[thid].size();
        for (int i=0 ; i<completed[thid].size() ; ++i) {
            completed_cells[completed[thid][i].first] = completed[thid][i].second;
        }
        for (int i=0 ; i<discovered[thid].size() ; ++i) {
            std::map<p_cell_type, nvis::vec2>::iterator it = completed_cells.find(discovered[thid][i].first);
            if (it != completed_cells.end()) {
                ++nb_found_already;
            } else {
                cell_to_seed[discovered[thid][i].first] = discovered[thid][i].second;
            }
        }
        for (int i=0 ; i<fps[thid].size() ; ++i) {
            fixpoints.push_back(fps[thid][i]);
        }
    }
    
    std::cerr << nb_found << " fixed points were found\n";
    std::cerr << cell_to_seed.size() << " additional fixed points were discovered\n";
}

void uniquify_fixed_points(std::vector<fixpoint>& fps)
{
    std::vector<std::list<fixpoint> > p2fp(maxp+1);
    for (int i=0 ; i<fps.size() ; ++i) {
        int p = fps[i].K;
        const nvis::vec2& x = fps[i].pos;
        std::list<fixpoint>& _list = p2fp[p];
        if (_list.empty() ||
                epsilon_find(_list.begin(), _list.end(), x, _plane_metric, 1.0e-5) == _list.end()) {
            _list.push_back(fps[i]);
        }
    }
    fps.clear();
    for (int i=1 ; i<=maxp ; ++i) {
        std::copy(p2fp[i].begin(), p2fp[i].end(), std::back_inserter(fps));
    }
}

template<typename Map>
void compute_fixed_points(const Map& pmap)
{
    std::sort(_cells_indices.begin(), _cells_indices.end());
    
    std::map<p_cell_type, nvis::vec2>        cell_to_seed;
    for (int i=0 ; i<_cells_indices.size() ; ++i) {
        const nvis::ivec2& cell_id = _cells_indices[i].object().object();
        cell_to_seed[_cells_indices[i].object()] = (*_plane)(cell_id) + 0.5*_plane_spacing;
        _all_seeds.push_back((*_plane)(cell_id) + 0.5*_plane_spacing);
    }
    
    std::map<p_cell_type, nvis::vec2>                                        completed_cells;
    std::vector<fixpoint>&                                                                 fixpoints = _fixpoints;
    
    for (int i=0 ; cell_to_seed.size() ; ++i) {
        std::cerr << "\n\n\n\nSTARTING ROUND " << i << " of fixed point extraction\n";
        __fixed_points(pmap, fixpoints, cell_to_seed, completed_cells, (i>0));
        if (!i) {
            __fixed_points_from_boxes(pmap, fixpoints, cell_to_seed, completed_cells);
        }
        std::cerr << "COMPLETED ROUND " << i << " of fixed point extraction\n";
        std::cerr << "there are currently " << fixpoints.size() << " fixed points identified\n";
        uniquify_fixed_points(fixpoints);
        std::cerr << "of those " << fixpoints.size() << " are unique\n";
        for (int p=1 ; p<=maxp ; ++p) {
            std::cerr << "PERIOD " << p << '\n';
            for (int n=0 ; n<fixpoints.size() ; ++n) {
                const fixpoint& fp = fixpoints[n];
                if (fp.K == p) {
                    std::cerr << fp << '\n';
                }
            }
        }
        std::cerr << "\n\n\n\n";
    }
    std::cerr << fixpoints.size() << " total fixpoints found\n";
    _unfiltered_fixed_points.clear();
    std::copy(fixpoints.begin(), fixpoints.end(), std::back_inserter(_unfiltered_fixed_points));
}


// ******************************    MAIN    ******************************
int main(int argc, char* argv[])
{
    //Jupiter-Europa system example
    mu = 2.528017705e-5;
    C = 3.000;
    center = nvis::vec2(-1.00001053340710,0); //L3 in JE
    
    me = argv[0];
    bounds.min() = bounds.max() = nvis::vec2(0,0);
    res = nvis::ivec2(-1, -1);
    eps = 1.0e-8;
    maxi = 100;
    maxp = 15;
    int verbose = 0;
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            printUsageAndExit("");
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                printUsageAndExit("missing resolution");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
        } else if (arg == "-c" || arg == "--center") {
            if (i >= argc-2) {
                printUsageAndExit("missing poloidal rotation center");
            }
            center[0] = atof(argv[++i]);
            center[1] = atof(argv[++i]);
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
        } else if (arg == "-mp" || arg == "--maxperiod") {
            if (i == argc-1) {
                printUsageAndExit("missing max period");
            }
            maxp = atof(argv[++i]);
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
        } else if (arg == "-d" || arg == "--display") {
            if (i >= argc-2) {
                printUsageAndExit("missing display size");
            }
            display_size[0] = atoi(argv[++i]);
            display_size[0] = atoi(argv[++i]);
        } else if (arg == "-ps" || arg == "--point") {
            if (i == argc-1) {
                printUsageAndExit("missing point size");
            }
            point_size = atoi(argv[++i]);
        } else if (arg == "-ls" || arg == "--line") {
            if (i == argc-1) {
                printUsageAndExit("missing line width");
            }
            line_width = atoi(argv[++i]);
        } else {
            printUsageAndExit("unrecognized argument");
        }
    }
    if (nvis::all(bounds.min() == bounds.max()) || res[0]<0 || res[1]<0 ) {
        printUsageAndExit("");
    }
    
    width = display_size[0];
    height = display_size[1];
    
    // initialize GLUT
    if (!nogfx) {
        glutInit(&argc, argv);
        glutInitDisplayString("samples rgba double alpha");
        glutInitWindowSize(width, height);
        glutInitWindowPosition(20, 20);
        main_window = glutCreateWindow(argv[0]);
        
        // configure OpenGL for aesthetic results
        if (!print_friendly) {
            glEnable(GL_LINE_SMOOTH);
        }
        if (!print_friendly) {
            glHint(GL_LINE_SMOOTH, GL_NICEST);
        }
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_POLYGON_SMOOTH);
        glHint(GL_POINT_SMOOTH, GL_NICEST);
        glHint(GL_POLYGON_SMOOTH, GL_NICEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    
    std::vector<nvis::fvec3> colors;
    xavier::spiral_scale(colors, maxp, 1);
    std::vector<int> reference_values(maxp);
    for (int i=1 ; i<=maxp ; ++i) {
        reference_values[i-1] = i;
    }
    _cmap = new xavier::discrete_color_map<int>(reference_values, colors);
    
    if (!nogfx) {
        // draw mesh
        double x = _plane_bounds.min()[0];
        for (int i=0 ; i<_plane_res[0] ; ++i, x+=_plane_spacing[0]) {
            _params.edges.push_back(nvis::vec2(x, _plane_bounds.min()[1]));
            _params.edges.push_back(nvis::vec2(x, _plane_bounds.max()[1]));
        }
        double y = _plane_bounds.min()[1];
        for (int j=0 ; j<_plane_res[1] ; ++j, y+=_plane_spacing[1]) {
            _params.edges.push_back(nvis::vec2(_plane_bounds.min()[0], y));
            _params.edges.push_back(nvis::vec2(_plane_bounds.max()[0], y));
        }
    }
    init();
    
    int npoints = _plane_res[0] * _plane_res[1];
    
    rhs_type rhs(C, mu);
    section_type section(rhs);
    section.bounds() = bounds;
    map_type pmap(rhs, section);
    pmap.setPrecision(eps);
    
    // define angular tracker for winding number computation
    BasicTracker tracker(center);
    
    nvis::timer _timer;
    std::cerr << "Sampling map at " << npoints << " grid vertices...\n";
    sample_raster<map_type, BasicTracker>(*_dataset, *_plane, pmap, _params, tracker, true);
    double dt = _timer.elapsed();
    std::cerr << "map sampling at grid vertices took " << dt << " seconds "
              << "(" << _plane_res[0]*_plane_res[1]/dt << "Hz)\n";
              
    // check resulting values
    for (int n=0 ; n<npoints ; ++n) {
        int i = n % _plane_res[0];
        int j = n / _plane_res[0];
        nvis::ivec2 id(i,j);
        if (!(*_dataset)(id).valid()) {
            _problematic_seeds.push_back((*_plane)(id));
        }
    }
    
    std::vector<p_cell_type> cells;
    compute_cell_periods(cells);
    compute_edge_rotation(pmap, cells);
    compute_cell_index(cells);
    compute_fixed_points(pmap);
    
    return 0;
    
}
