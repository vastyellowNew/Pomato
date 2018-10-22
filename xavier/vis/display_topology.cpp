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


/*************************************************************************
Copyright (c) 2018, Xavier Tricoche, Purdue University

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
**********************************************************************/


//Standard Libs
#include <boost/format.hpp>
#include <boost/limits.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <vector>

// cxxopts: command line parsing
#include <misc/cxxopts.hpp>
#include <misc/filename.hpp>
// nvis
#include <math/rational.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/wall_timer.hpp>
// xavier
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/multipleAngleTracker.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/fixpoints.hpp>
#include <maps/poincare_map.hpp>
#include <maps/DP45wrapper.hpp>
#include <orbital/corrections.hpp>
#include <pmate/ManifoldData.hpp>
#include <pmate/FixedPointData.hpp>
#include <topology/invariant_manifold.hpp>
#include <vis/vtk_utils.hpp>

#include "display_topology.hpp"

#include "vtkCoordinate.h"
#include "vtkOutlineSource.h"
#include "vtkInteractorStyleImage.h"
#include "vtkSphereSource.h"
#include "vtkTextureMapToSphere.h"
#include "vtkProperty.h"
#include "vtkTexture.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkColor.h"
#include "vtkPropPicker.h"

#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;
using namespace topology;

//Type Definitions
typedef nvis::fixed_vector<double, 2>                                  vec2;
typedef nvis::fixed_vector<double, 3>                                  vec3;
typedef nvis::fixed_vector<double, 6>                                  vec6;
typedef nvis::fixed_vector<double, 42>                                 vec42;
typedef xavier::grid<double, 2>                                        grid_type;
typedef xavier::dp5wrapper<double, 42>                                 ode_solver;
typedef orbital::cr3bp                                                 rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>                       section_type;
typedef orbital::MultipleAngleTracker<42, 3>                           tracker_type;
typedef poincare_map<rhs_type, ode_solver, section_type >              map_type;
typedef xavier::state_info<rhs_type::value_type, rhs_type::dimension>  return_state;
typedef map_type::return_type                                          return_type;
typedef map_type::lvec_type                                            vec_type;
typedef xavier::map_analysis_param                                     map_parameters;
typedef std::vector<xavier::fixpoint>                                  chain_type;
typedef pmate::FixedPointData<xavier::fixpoint>                        fixpt_data;
typedef std::vector<nvis::vec2>                                        orbit_type;
typedef nvis::bounding_box<vec2>                                       bounds_type;

static const int s = rhs_type::numSingularities;
typedef nvis::fixed_vector<double,s+1>                                 xmap_data_type;
typedef topology::SortableReturnData<vec_type, xmap_data_type>         sortable_data_type;
typedef pmate::ManifoldData<sortable_data_type, xavier::fixpoint>      manifold_data;
typedef xavier::EdgeRotationFailure<vec_type>                          rotation_failure_type;
typedef xavier::EdgeRotationFailure<vec_type>                          map_discont_type;
typedef manifold_data::ManifoldSegType                                 manifold_seg_type;

//Static declarations
const double invalid_double = std::numeric_limits<double>::max();

typedef std::array<vec2, 2> segment_type;
typedef std::vector<vec2> curve_type;
typedef std::pair<int, int> segment_index_type;

int orbit_id, fp_id;
bool display_eigenvectors, display_stable_manifolds,
     display_unstable_manifolds, display_fixpoints, display_discontinuities,
     display_topology, display_all, verbose, fit, save_and_exit;
std::array<int, 2> res;
std::string fp_data_name, man_data_name, system_name, camera_name, title;
double eigen_length, vscale;
int pt_size;
std::array<float, 3> bg_color;
std::array<double, 4> box({{invalid_double, invalid_double, invalid_double, invalid_double}});

size_t image_counter=0;
std::string image_basename;

class cr3bp_style : public vtkInteractorStyleImage {
public:
    static cr3bp_style* New();
    vtkTypeMacro(cr3bp_style, vtkInteractorStyleImage);

    void set_renderer(vtkRenderer* renderer) { m_renderer = renderer; image_counter = 0; }

    virtual void OnKeyPress() override
    {
      // Get the pressed key by symbol
      vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();

      vtkRenderWindow* window = rwi->GetRenderWindow();

      if (key == "c") {
          vtk_utils::print_camera_settings(std::cout, m_renderer);
      }
      else if (key == "s") {
          std::ostringstream oss;
          oss << image_basename << "_" << std::setw(4) << std::setfill('0') << image_counter++ << ".jpg";
          vtk_utils::save_frame(window, oss.str());
      }
      else if (key == "at") {
          int x, y;
          rwi->GetMousePosition(&x, &y);
          VTK_CREATE(vtkCoordinate, coord);
          coord->SetCoordinateSystemToDisplay();
          coord->SetValue(x,y,0);
          double* p = coord->GetComputedWorldValue(m_renderer);
          std::cout << "Display coordinates: (" << x << ", " << y << "), world coordinates: (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
      }
      // Forward events
      vtkInteractorStyleImage::OnKeyPress();
    }

    virtual void OnLeftButtonDown() override
    {
        int* clickPos = this->GetInteractor()->GetEventPosition();
        VTK_CREATE(vtkPropPicker, picker);
        picker->Pick(clickPos[0], clickPos[1], 0, m_renderer);
        vtkActor* actor = picker->GetActor();
        if (actor) {
            if (modified_actors.find((long int)actor) != modified_actors.end()) {
                if (verbose) std::cout << "unselecting actor\n";
                actor->GetProperty()->SetPointSize(pt_size);
                modified_actors.erase((long int)actor);
            }
            else {
                if (verbose) {
                    std::cout << "selected actor\n";
                    actor->PrintSelf(std::cout, vtkIndent());
                }
                actor->GetProperty()->SetPointSize(5*pt_size);
                modified_actors.insert((long int)actor);
            }
        }
        else {
            if (verbose) std::cout << "no actor selected\n";
        }

        vtkInteractorStyleImage::OnLeftButtonDown();
    }

private:
    vtkSmartPointer<vtkRenderer> m_renderer;
    size_t image_counter;
    std::set<long int> modified_actors;

};
vtkStandardNewMacro(cr3bp_style);

typedef enum { STABLE, UNSTABLE, NONE, UNKNOWN } filtering_type;

filtering_type manifold_filter(NONE);

inline vtkSmartPointer<vtkActor> make_actor(vtkPolyData* data) {
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(data);
    vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    return actor;
}

template<typename Point_>
inline Point_ scaled(const Point_& pt) {
    Point_ _scaled = pt;
    _scaled[1] *= vscale;
    return _scaled;
}

template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T, N>& a) {
    os << "[";
    for (size_t i=0; i<N-1; ++i) os << a[i] << ", ";
    os << a.back() << ']';
    return os;
}

bool parse_keyvalue(std::string& key, std::string& value, std::fstream& fs) {
    std::string buffer;
    std::getline(fs, buffer);
    std::cout << "line = " << buffer << '\n';
    if (buffer[0] == '#') return false;
    size_t i = buffer.find_first_of(':');
    if (i == std::string::npos) return false;
    key = buffer.substr(0, i);
    std::string valuestr = buffer.substr(i+1);
    if (valuestr.empty()) return false;
    value = valuestr.substr(valuestr.find_first_not_of(" "));
    return true;
}

VTK_SMART(vtkActor) make_planet(const std::string& filename, double radius, double x) {
    VTK_CREATE(vtkSphereSource, sphere);
    sphere->SetRadius(radius);
    sphere->SetThetaResolution(45);
    sphere->SetPhiResolution(90);
    sphere->SetCenter(x, 0, 0);
    sphere->Update();
    VTK_CREATE(vtkTextureMapToSphere, tosphere);
    tosphere->SetInputConnection(sphere->GetOutputPort());
    tosphere->SetPreventSeam(0);
    VTK_CREATE(vtkPolyDataMapper, mapper);
    mapper->SetInputConnection(tosphere->GetOutputPort());
    mapper->ScalarVisibilityOff();
    VTK_CREATE(vtkJPEGReader, reader);
    reader->SetFileName(filename.c_str());
    reader->Update();
    VTK_CREATE(vtkTexture, texture);
    texture->SetInputConnection(reader->GetOutputPort());
    texture->EdgeClampOn();
    VTK_CREATE(vtkProperty, property);
    property->BackfaceCullingOn();
    VTK_CREATE(vtkActor, planet);
    planet->SetMapper(mapper);
    planet->SetTexture(texture);
    planet->SetProperty(property);
    planet->RotateX(90);
    return planet;
}

int main(int argc, const char* argv[]) {
    std::string me = argv[0];

    try {
        cxxopts::Options options(argv[0], "Visualize CR3BP topology extracted by PMATE");
        options
            .positional_help("[optional args]")
                .show_positional_help();

        options.add_options("Required")
            ("f,fps", "Fixed points data input filename", cxxopts::value<std::string>(fp_data_name))
            ("m,manis", "Invariant manifolds data input filename", cxxopts::value<std::string>(man_data_name));
        options.add_options("Display control")
            ("fpid", "Selected fixed point index", cxxopts::value<int>(fp_id)->default_value("-1"))
            ("orbid", "Selected orbit index", cxxopts::value<int>(orbit_id)->default_value("-1"))
            ("show_fps", "Display fixed points", cxxopts::value<bool>(display_fixpoints)->default_value("false"))
            ("show_evs", "Display eigenvectors", cxxopts::value<bool>(display_eigenvectors)->default_value("false"))
            ("show_stable", "Display stable manifolds", cxxopts::value<bool>(display_stable_manifolds)->default_value("false"))
            ("show_unstable", "Display unstable manifolds", cxxopts::value<bool>(display_unstable_manifolds)->default_value("false"))
            ("show_discont", "Display map discontinuities", cxxopts::value<bool>(display_discontinuities)->default_value("false"))
            ("show_topo", "Display fixed points and manifolds", cxxopts::value<bool>(display_topology)->default_value("false"))
            ("show_all", "Display complete topology", cxxopts::value<bool>(display_all)->default_value("false")->implicit_value("true"));
        options.add_options("Graphics settings")
            ("p,ptsize", "Fixed point size (in pixels)", cxxopts::value<int>(pt_size)->default_value("4"))
            ("l,length", "Length of eigenvectors", cxxopts::value<double>(eigen_length)->default_value("0.001"))
            ("r,res", "Window resolution", cxxopts::value< std::array<int, 2> >(res)->default_value("1024 768"), "\"<width> <height>\"")
            ("bgcolor", "Background color", cxxopts::value< std::array<float, 3> >(bg_color)->default_value("0 0 0"), "\"<r> <g> <b>\"")
            ("s,vscale", "Vertical scaling", cxxopts::value<double>(vscale)->default_value("0.16"))
            ("b,bounds", "Bounding box", cxxopts::value< std::array<double, 4> >(box), "\"xmin xdotmin xmax xdotmax\"")
            ("fit", "Fit bounds in window", cxxopts::value<bool>(fit)->default_value("false"))
            ("system", "Planet/moon system data file", cxxopts::value<std::string>(system_name))
            ("snapshot", "Snapshot basename", cxxopts::value<std::string>(image_basename)->default_value("snapshot"))
            ("camera", "Camera setup", cxxopts::value<std::string>(camera_name))
            ("title", "Legend title", cxxopts::value<std::string>(title))
            ("exit", "Save snapshot and exit", cxxopts::value<bool>(save_and_exit)->default_value("false"));
        options.add_options("Miscellaneous")
            ("v,verbose", "Verbose output", cxxopts::value<bool>(verbose)->default_value("false"))
            ("h,help", "Display this message");

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help(std::vector<std::string>()/*{""}*/) << std::endl;
            exit(0);
        }

        if (!result.count("fps") && !result.count("manis")) {
            std::cerr << "ERROR: Required parameters are missing\n";
            std::cout << options.help(std::vector<std::string>()) << std::endl;
            exit(1);
        }
    }
    catch (const cxxopts::OptionException& e) {
        std::cerr << "error parsing command line options: " << e.what() << std::endl;
        exit(1);
    }

    if (display_all) { // overrule individual selections
        display_stable_manifolds = display_unstable_manifolds =
            display_fixpoints = display_eigenvectors =
            display_discontinuities = true;
    }
    else if (display_topology) {
        display_stable_manifolds = display_unstable_manifolds = display_fixpoints = true;
    }

    fixpt_data fp_data(fp_data_name.c_str());

    bool no_manifold = man_data_name.empty();

    if (display_stable_manifolds || display_unstable_manifolds || display_discontinuities) {
        if (no_manifold) {
            display_stable_manifolds = display_unstable_manifolds = display_discontinuities = false;
            if (verbose) {
                std::cout << "No manifold data filename provided. Skipping manifold visualization.\n";
            }
        }
    }

    if (verbose) {
        std::cout << "resolution is " << res << '\n';
        std::cout << "display stable manifolds:   " << (display_stable_manifolds ? "true" : "false") << '\n';
        std::cout << "display unstable manifolds: " << (display_unstable_manifolds ? "true" : "false") << '\n';
        std::cout << "display fixed points:       " << (display_fixpoints ? "true" : "false") << '\n';
        std::cout << "display eigenvectors:       " << (display_eigenvectors ? "true" : "false") << '\n';
        std::cout << "display discontinuities:    " << (display_discontinuities ? "true" : "false") << '\n';
        std::cout << "bounds:                     " << box << '\n';
        std::cout << "background color:           " << bg_color << '\n';
        std::cout << "Fixed point size:           " << pt_size << '\n';
    }

    manifold_data mdata(fp_data);
    if (!no_manifold) mdata.read(man_data_name.c_str());

    VTK_SMART(vtkActor) P1_actor, P2_actor;
    if (!system_name.empty()) {
        std::fstream sysf(system_name.c_str(), std::ios::in);
        double mu=-1, r1=-1, r2=-1;
        std::string p1_name, p2_name, path;
        while (!sysf.eof()) {
            std::string kw, val;
            bool ok = parse_keyvalue(kw, val, sysf);
            if (!ok) {
                //std::cout << "skipped line\n";
                continue;
            }
            else if (kw == "system" || kw == "System") {
                if (verbose) std::cout << "Considered system is " << val << '\n';
            }
            else if (kw == "mu") {
                mu = std::stod(val);
                if (verbose) std::cout << "mu value=" << mu << '\n';
            }
            else if (kw == "P1_radius" || kw == "R1" || kw== "r1") {
                r1 = std::stod(val);
                if (verbose) std::cout << "R1=" << r1 << '\n';
            }
            else if (kw == "P2_radius" || kw == "R2" || kw== "r2") {
                r2 = std::stod(val);
                if (verbose) std::cout << "R2=" << r2 << '\n';
            }
            else if (kw == "P1_image") {
                p1_name = val;
                if (verbose) std::cout << "P1 image is " << p1_name << '\n';
            }
            else if (kw == "P2_image") {
                p2_name = val;
                if (verbose) std::cout << "P2 image is " << p2_name << '\n';
            }
            else if (kw == "path") {
                path = val;
            }
            else {
                std::cout << "unrecognized keyword: " << kw << '\n';
            }
        }
        sysf.close();
        if (mu > 0 && r1 > 0 && r2 > 0 && !p1_name.empty() && !p2_name.empty()) {
            P1_actor = make_planet(path + '/' + p1_name, r1, -mu);
            P2_actor = make_planet(path + '/' + p2_name, r2, 1-mu);
        }
    }

    std::vector< std::vector< segment_type > > stable_eigenvectors, unstable_eigenvectors;
    std::vector< std::vector< curve_type > > stable_manifolds, unstable_manifolds;
    std::vector< std::vector<vec2> > saddles, centers;
    std::vector< std::vector< std::pair<vec2, int> > > fails;
    double eigen_length = 0.001;
    // std::map<segment_index_type, size_t> segment_to_manifold;

    std::vector< vtkSmartPointer<vtkActor> > unstable_evec_actors, stable_evec_actors,
    saddles_actors, centers_actors, fails_actors, unstable_manifolds_actors,
    stable_manifolds_actors;


    // Create a rendering queue based on display/filter options
    //-----------------------------------------------------------------------------------------
    std::queue<int> manifold_render_queue;
    int numManifolds = 0;
    if (!no_manifold)
        numManifolds = (int) mdata.mapManifolds.size();
    if (orbit_id >= 0) { // Valid orbit ID
        // Only one orbit, check if we are doing only one fixed point or not
        if (fp_id >= 0) { // Valid fixed point ID
            // Only one fixed point of orbit chain
            // Select which manifolds from this orbit and particular fixed point
            if (display_stable_manifolds && !display_unstable_manifolds) {
                // Only stable manifolds
                for (int mID=0;mID<numManifolds;mID++) {
                    int baseorbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
                    int fpIndex = mdata.mapManifolds[mID].fpdPointIdx;
                    bool isUnstable = mdata.mapManifolds[mID].isForward();
                    if (!isUnstable && (baseorbit_id==orbit_id) && (fpIndex==fp_id))
                        manifold_render_queue.push(mID);
                }
            }
            else if (display_unstable_manifolds && !display_stable_manifolds) {
                // Only unstable manifolds
                for (int mID=0;mID<numManifolds;mID++) {
                    int baseorbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
                    int fpIndex = mdata.mapManifolds[mID].fpdPointIdx;
                    bool isUnstable = mdata.mapManifolds[mID].isForward();
                    if (isUnstable && (baseorbit_id==orbit_id) && (fpIndex==fp_id))
                        manifold_render_queue.push(mID);
                }
            }
            else {
                // All manifolds from orbit
                for (int mID=0;mID<numManifolds;mID++) {
                    int baseorbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
                    int fpIndex = mdata.mapManifolds[mID].fpdPointIdx;
                    if((baseorbit_id == orbit_id) && (fpIndex==fp_id)) manifold_render_queue.push(mID);
                }
            }
        }
        else {
            // Select which manifolds from this orbit_id
            if (display_stable_manifolds && !display_unstable_manifolds) {
                // Only stable manifolds
                for (int mID=0;mID<numManifolds;mID++) {
                    int baseorbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
                    bool isUnstable = mdata.mapManifolds[mID].isForward();
                    if (!isUnstable && (baseorbit_id==orbit_id)) manifold_render_queue.push(mID);
                }
            }
            else if (display_unstable_manifolds && !display_stable_manifolds) {
                // Only unstable manifolds
                for (int mID=0;mID<numManifolds;mID++) {
                    int baseorbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
                    bool isUnstable = mdata.mapManifolds[mID].isForward();
                    if (isUnstable && (baseorbit_id==orbit_id)) manifold_render_queue.push(mID);
                }
            }
            else {
                // All manifolds from orbit
                for (int mID=0;mID<numManifolds;mID++) {
                    int baseorbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
                    if(baseorbit_id == orbit_id) manifold_render_queue.push(mID);
                }
            }
        }
    }
    else {
        // Gather manifolds based on type
        if (display_stable_manifolds && !display_unstable_manifolds) {
            // Only stable manifolds
            for (int mID=0;mID<numManifolds;mID++) {
                bool isUnstable = mdata.mapManifolds[mID].isForward();
                if (!isUnstable) manifold_render_queue.push(mID);
            }
        }
        else if (display_unstable_manifolds && !display_stable_manifolds) {
            // Only unstable manifolds
            for (int mID=0;mID<numManifolds;mID++) {
                bool isUnstable = mdata.mapManifolds[mID].isForward();
                if (isUnstable) manifold_render_queue.push(mID);
            }
        }
        else {
            if (verbose) {
                std::cout << "Displaying both manifold types\n";
            }
            // All manifolds
            for (int mID=0;mID<numManifolds;mID++) manifold_render_queue.push(mID);
        }
    }

    //Render eigenVectors as line segments => NOTE: Renders ALL orbits/fixedpoints for now.
    //-----------------------------------------------------------------------------------------
    if (display_eigenvectors) {
        //Loop through all available fixed points in FixedPointData to define points
        int numfp =  fp_data.getNumFixedPoints();
        int numOrbits = fp_data.getNumOrbits();
        if (verbose) {
            std::cout << "Creating eigenvectors for " << numfp
                << " fixed points for " << numOrbits << " orbits\n";
        }
        // check that there is indeed something to do
        if (numfp > 0) {
            int eigUPoints = 0, eigSPoints = 0;
            int eigULineNum = 0, eigSLineNum = 0;
            //For each orbit
            for (int k=0; k<numOrbits; k++) {
                //Only look at unstable orbits (saddles)
                if (!(fp_data.isSaddle(k))) continue;

                unstable_eigenvectors.push_back(std::vector< segment_type >() );
                stable_eigenvectors.push_back(std::vector< segment_type >());
                for (int i=0; i < fp_data.getOrbitPeriod(k); i++) {
                    //First chain only here - just one orbit
                    const fixpoint& fp = fp_data.getFixedPoint(k,i);
                    if (verbose) {
                        std::cout << me  << ": fp " << i << " : x0 = [" << fp.pos[0] << ", " << fp.pos[1]
                            << "], eigenvalues = [" << fp.eval[0] << ", " << fp.eval[1] << "], v_stable = ["
                            << fp.evec[0][0] << " , " << fp.evec[0][1] << "] unstable eigenvector = ["
                            << fp.evec[1][0] << " , " << fp.evec[1][1] << "]\n";
                    }

                    double s = eigen_length;
                    vec2 pt0 = fp.pos - s*fp.evec[1];
                    vec2 pt1 = fp.pos + s*fp.evec[1];
                    unstable_eigenvectors.back().push_back(segment_type{{scaled(pt0), scaled(pt1)}});
                    eigULineNum++;

                    pt0 = fp.pos - s*fp.evec[0];
                    pt1 = fp.pos + s*fp.evec[0];
                    stable_eigenvectors.back().push_back(segment_type{{scaled(pt0), scaled(pt1)}});
                    eigSLineNum++;
                }
                unstable_evec_actors.push_back(make_actor(vtk_utils::make_polylines(unstable_eigenvectors.back())));
                unstable_evec_actors.back()->GetProperty()->SetColor(1, 0, 0);
                stable_evec_actors.push_back(make_actor(vtk_utils::make_polylines(stable_eigenvectors.back())));
                stable_evec_actors.back()->GetProperty()->SetColor(0, 0, 1);
            }
        }
    } //End Eigenvector Render

    //-----------------------------------------------------------------------------------------
    //Render manifolds
    //-----------------------------------------------------------------------------------------

    //Determine if we will show the manifold objects:
    //Fill Coordinate Arrays for each Manifold object (from rendering queue)
    typedef std::pair<int,int> IntPair;
    int unstableSegTotal = 0, unstablePoints = 0;
    int stableSegTotal = 0, stablePoints = 0;

    //Make the main render queue as a priority_queue that puts long, highly unstable orbits first
    // (putting them on the bottom of the visualization).
    ManifoldDataRenderPriority render_priority(mdata);
    std::priority_queue<int,std::vector<int>,ManifoldDataRenderPriority> main_render_queue(render_priority);
    //Emplace queue objects
    while(!manifold_render_queue.empty()) {
        main_render_queue.push(manifold_render_queue.front());
        manifold_render_queue.pop();
    }
    //For each manifold object in queue
    while (!main_render_queue.empty()) {
        int mID = main_render_queue.top(); //Priority rendering
        int ptsPerLine = 0;
        bool isUnstable = mdata.mapManifolds[mID].isForward();
        if ((isUnstable && !display_unstable_manifolds) ||
            (!isUnstable && !display_stable_manifolds)) {
            main_render_queue.pop();
            continue;
        }

        vec2 prev(0.0,0.0);
        bool firstSeg = true, connected = true;
        //int forbit_id = fp_data.getorbit_id(mdata.mapManifolds[mID].fpdorbit_idx,
        //                                  mdata.mapManifolds[mID].fpdPointIdx);
        int forbit_id = mdata.mapManifolds[mID].fpdOrbitIdx;
        //For each segment

        curve_type* current_curve = NULL;
        std::vector<manifold_seg_type>::iterator segit;
        for (segit =  mdata.mapManifolds[mID].segments.begin();
             segit != mdata.mapManifolds[mID].segments.end(); segit++) {
            //The two points of this segment
            vec2& x0 = (*segit)[0];
            vec2& x1 = (*segit)[1];
            IntPair manIntPair(mID,segit->segID);

            //Check if this point will be within bounds
            // bool inBounds = true;
            // nvis::bbox2 bbox(nvis::vec2(portXbounds.getValue(0),portYbounds.getValue(0)),
            //                  nvis::vec2(portXbounds.getValue(1),portYbounds.getValue(1)) );
            // if ( !bbox.inside(x0) || !bbox.inside(x1) ) {
            //   connected = false;
            //   inBounds = false;
            // }
            // if (!inBounds) continue; //Go to next segment

            //Is this connected to the last segment
            if (!firstSeg) {
                if (nvis::all( x0 == prev )) {
                    connected = true;
                }
                else {
                    connected = false;
                }
            }
            //If first segment or not connected, add both points
            if (firstSeg || !connected) {
                //If not connected, we must end the last line
                if (isUnstable) {
                    // FIX ME: need reference from segment to position in corresponding manifold
                    unstable_manifolds.push_back(std::vector<curve_type>());
                    unstable_manifolds.back().push_back(curve_type());
                    unstable_manifolds.back().back().push_back(scaled(x0));
                    unstable_manifolds.back().back().push_back(scaled(x1));
                    unstableSegTotal++;
                }
                else {
                    stable_manifolds.push_back(std::vector<curve_type>());
                    stable_manifolds.back().push_back(curve_type());
                    stable_manifolds.back().back().push_back(scaled(x0));
                    stable_manifolds.back().back().push_back(scaled(x1));
                    stableSegTotal++;
                }

                //Update to no longer be first segment
                firstSeg = false;
                //Restart the pointsPerLine Counter
                ptsPerLine = 2;
            }
            else { //Already connected so just add a single point
                if(isUnstable) {
                    unstable_manifolds.back().back().push_back(scaled(x1));
                    unstableSegTotal++;
                }
                else {
                    stable_manifolds.back().back().push_back(scaled(x1));
                    stableSegTotal++;
                }
                ptsPerLine++;
            }
            //Store the last point to see if we need it again
            prev = x1;
        }
        unstable_manifolds_actors.push_back(make_actor(vtk_utils::make_polylines(unstable_manifolds.back())));
        stable_manifolds_actors.push_back(make_actor(vtk_utils::make_polylines(stable_manifolds.back())));
        unstable_manifolds_actors.back()->GetProperty()->SetColor(1,0,0);
        stable_manifolds_actors.back()->GetProperty()->SetColor(0,0,1);


        //Progress queue
        main_render_queue.pop();

    } //End while loop for main render

    if (verbose) {
        // Output
        std::cout << __FILE__ << ": Stable Lines = " << stable_manifolds.size() << " with " << stableSegTotal << " segments\n";
        std::cout << __FILE__ << ": Unstable Lines = " << unstable_manifolds.size() << " with " << unstableSegTotal << " segments\n";
    }

    //Render fixed points - Only currently renders ALL
    //---------------------------------------------------------------------------------------
    if (display_fixpoints) {
        //After manifolds to make them appear "on top"
        //Fill Coordinate arrays
        int numOrbits = fp_data.getNumOrbits();
        std::vector<vtkColor3ub> colors;
        size_t n_saddles=0;
        size_t n_centers=0;
        srand48(1234);
        for(int k=0;k<numOrbits;k++) {
            colors.clear();
            unsigned short saturation = 256*drand48();
            if (fp_data.isSaddle(k)) {
                saddles.push_back(std::vector< vec2 >());
            }
            else {
                centers.push_back(std::vector< vec2 >());
            }
            for(int i=0;i<fp_data.getOrbitPeriod(k); i++) {
                const xavier::fixpoint& fp = fp_data.getFixedPoint(k,i);
                if (fp_data.isSaddle(k)) {
                    saddles.back().push_back(scaled(fp.pos));
                    colors.push_back(vtkColor3ub(255, 255-saturation, 255-saturation));
                }
                else {
                    centers.back().push_back(scaled(fp.pos));
                    colors.push_back(vtkColor3ub(255-saturation, 255-saturation, 255));
                }
            }
            vtkSmartPointer<vtkPolyData> pd;
            if (fp_data.isSaddle(k))
                pd = vtk_utils::make_points(saddles.back());
            else
                pd = vtk_utils::make_points(centers.back());
            vtk_utils::add_vertices(pd);
            vtk_utils::add_colors(pd, colors);
            if (fp_data.isSaddle(k)){
                saddles_actors.push_back(make_actor(pd));
            }
            else {
                centers_actors.push_back(make_actor(pd));
            }

        }
        if (verbose) {
            std::cout << saddles.size() << " saddle chains and " << centers.size() << " center chains to display\n";
        }
    }


    //---------------------------------------------------------------------------------------
    //Render the Map Discontinuities
    //---------------------------------------------------------------------------------------
    if (display_discontinuities) {
        //Iterate through the sepList and grab all points from manifold
        // std::vector<vec2> xFail;
        // std::vector<SbColor> failColors;
        std::list< map_discont_type >::const_iterator failIT;
        int failPtIdx = 0;
        std::queue<int> renderQueue(manifold_render_queue);
        //For each Manifold in render queue
        while (!renderQueue.empty()) {
            int k = renderQueue.front();
            std::list<map_discont_type> sepList;
            mdata.mapManifolds[k].getMapDiscontinuityList(sepList);
            //For each MapDiscont
            if (!sepList.empty())
                fails.push_back(std::vector< std::pair<vec2, int> >());
            for(failIT = sepList.begin(); failIT!=sepList.end(); ++failIT) {
                std::pair<vec2, int> fail_pt(failIT->where(), failIT->type);
                fails.back().push_back(fail_pt);
                //Increment counter
                failPtIdx++;
            } //End SepList Loop

            //Mark the last seed point - End of processed manifold
            vec2 lastSeedPt = (mdata.mapManifolds[k].getWorkingSegment())[1];
            fails.back().push_back(std::make_pair(lastSeedPt, -1));
            failPtIdx++;

            std::vector<int> cases(fails.back().size());
            std::vector<vec2> pts(fails.back().size());
            for (size_t i=0; i<pts.size(); ++i) {
                pts[i] = fails.back()[i].first;
                cases[i] = fails.back()[i].second;
            }
            vtkSmartPointer<vtkPolyData> pd = vtk_utils::make_points(pts);
            vtk_utils::add_scalars(pd, cases, true, "discontinuity/failure case");
            vtk_utils::add_vertices(pd);
            fails_actors.push_back(make_actor(pd));
            fails_actors.back()->GetProperty()->SetColor(1,0,1);

            //Move to next manifold
            renderQueue.pop();
        } //End While Loop
    }

    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
    window->SetSize(res[0], res[1]);

    // renWin->SetAlphaBitPlanes(1);
    // renWin->SetMultiSamples(0);
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    window->AddRenderer(renderer);
    if (display_eigenvectors) {
        std::for_each(unstable_evec_actors.begin(), unstable_evec_actors.end(), [&](VTK_SMART(vtkActor) a)
            { renderer->AddActor(a); });
        std::for_each(stable_evec_actors.begin(), stable_evec_actors.end(), [&](VTK_SMART(vtkActor) a)
            { renderer->AddActor(a); });
    }
    if (P1_actor) renderer->AddActor(P1_actor);
    if (P2_actor) renderer->AddActor(P2_actor);
    if (display_fixpoints) {
        std::for_each(saddles_actors.begin(), saddles_actors.end(), [&](VTK_SMART(vtkActor) a)
            {
                a->GetProperty()->SetPointSize(pt_size);
                a->GetProperty()->SetRenderPointsAsSpheres(1);
                renderer->AddActor(a);
            });

        std::for_each(centers_actors.begin(), centers_actors.end(), [&](VTK_SMART(vtkActor) a)
            {
                a->GetProperty()->SetPointSize(pt_size);
                a->GetProperty()->SetRenderPointsAsSpheres(1);
                renderer->AddActor(a);
            });
    }
    if (display_discontinuities) {
        std::for_each(fails_actors.begin(), fails_actors.end(), [&](VTK_SMART(vtkActor) a)
            {
                a->GetProperty()->SetPointSize(pt_size);
                a->GetProperty()->SetRenderPointsAsSpheres(1);
                renderer->AddActor(a);

            });
    }
    if (display_unstable_manifolds) {
        std::for_each(unstable_manifolds_actors.begin(), unstable_manifolds_actors.end(), [&](VTK_SMART(vtkActor) a)
            {
                renderer->AddActor(a);
            });
    }
    if (display_stable_manifolds) {
        std::for_each(stable_manifolds_actors.begin(), stable_manifolds_actors.end(), [&](VTK_SMART(vtkActor) a)
            {
                renderer->AddActor(a);
            });
    }
    renderer->SetBackground(bg_color[0], bg_color[1], bg_color[2]);

    if (std::all_of(box.begin(), box.end(), [](double x){return x!=invalid_double;})) {
        if (fit) {
            if (verbose) std::cout << "setting camera scope to valid bounds\n";
            vtk_utils::fill_window(renderer, bounds_type(vec2(box[0], box[1]), vec2(box[2], box[3])));
        }
        VTK_CREATE(vtkOutlineSource, frame);
        frame->SetBounds(box[0], box[2], box[1], box[3], 0, 0);
        VTK_CREATE(vtkPolyDataMapper, fmapper);
        VTK_PLUG(fmapper, frame);
        VTK_CREATE(vtkActor, factor);
        factor->SetMapper(fmapper);
        factor->GetProperty()->SetLineStipplePattern(0xf0f0);
        factor->GetProperty()->SetLineStippleRepeatFactor(1);
        factor->GetProperty()->SetPointSize(1);
        factor->GetProperty()->SetLineWidth(1);
        factor->GetProperty()->SetColor(0.5, 0.5, 0.5);
        renderer->AddActor(factor);
    }
    else if (!camera_name.empty()) {
        vtk_utils::import_camera_settings(camera_name, renderer);
    }
    else {
        renderer->ResetCamera();
    }

    if (!title.empty()) {
        VTK_CREATE(vtkTextActor, txt);
        txt->SetInput(title.c_str());
        txt->GetTextProperty()->SetFontFamilyToArial();
        txt->GetTextProperty()->SetFontSize(22);
        txt->GetTextProperty()->SetColor(1,1,1);
        txt->SetDisplayPosition(20, 20);
        txt->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedViewport();
        txt->GetPosition2Coordinate()->SetValue(0.1, 0.1);
        renderer->AddActor(txt);
    }

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(window);
    VTK_CREATE(cr3bp_style, style);
    style->set_renderer(renderer);
    interactor->SetInteractorStyle(style);

    window->Render();

    if (save_and_exit) {
        std::string basename = xavier::filename::remove_extension(image_basename);
        vtk_utils::save_frame(window, basename + ".jpg");
        exit(0);
    }

    interactor->Start();

    return 0;
}
