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


#ifndef __XAVIER_MAPS_INDEX_HPP__
#define __XAVIER_MAPS_INDEX_HPP__

#include <stdexcept>
#include <exception>
#include <iostream>
#include <maps/definitions.hpp>
#include <maps/mapExceptions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/misc.hpp>
//Avizo
#ifdef HX_HAS_STD
#include <QString>
#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#endif

namespace xavier {


/// Compute the Displacement vector between a starting point and it's pth map iterate
template <class DATA, class METRIC>
inline vec_type p_vector(const vec_type& x0, const DATA& o, int p, const METRIC& _metric)
{
    //Double-checking for an error (likely not to be thrown at this point...
    if (p>(int) o.steps.size()) {
#ifdef HX_HAS_STD
        theMsg->printf("ERROR: Period Exceeds available data!!!");
        theMsg->printf("       Requested p = %d:  Available Steps = %d.", p, o.steps.size());
        theWorkArea->setProgressInfo("ERROR: Period Exceeds available data!!!");
        //theWorkArea->stopWorking();
#else
        std::cerr << "requested period " << p << " while steps contains " << o.steps.size() << " positions\n";
        //Throw an error and kill it
        throw std::runtime_error("period exceeds available data");
#endif
    }
    //return _metric.displacement(o.steps[0], o.steps[p]); //RasterData
    return _metric.displacement(x0, o.steps[p-1]); //AdaptiveGridNodeData
}

/// Compute the Displacement vector between a starting point and it's pth map iterate
template <class DATA, class METRIC>
inline vec_type p_vector(const DATA& o, int p, const METRIC& _metric)
{
    //Double-checking for an error (likely not to be thrown at this point...)
    if (p>=(int) o.steps.size()) {
#ifdef HX_HAS_STD
        theMsg->printf("ERROR: Period Exceeds available data!!!");
        theMsg->printf("       Requested p = %d:  Available Steps = %d.", p, o.steps.size());
        theWorkArea->setProgressInfo("ERROR: Period Exceeds available data!!!");
        //theWorkArea->stopWorking();
#else
        std::cerr << "requested period " << p << " while steps contains " << o.steps.size() << " positions\n";
        //Throw an error and kill it
        throw std::runtime_error("period exceeds available data");
#endif
    }
    return _metric.displacement(o.steps[0], o.steps[p]); //RasterData
}


/** Adaptive rotation angle computation using bisection
 * Assumes you have only one period to compute per edge!
 *  NOTE:  This is insanely inefficient for the
 * CR3BP as well as any function that is costly to integrate.
 * Created some new classes to handle these functions in data/edge.hpp,
 * data/cell.hpp, topology/Edge_Rotation.hpp.
 */
template<typename MAP>
double adaptive_rotation_angle(const vec_type& x0, const vec_type& v0,
                               const vec_type& x1, const vec_type& v1,
                               const MAP& pmap, unsigned int period,
                               double dx, map_analysis_param& params);
                               
template<typename MAP>
double robust_rotation_angle(const vec_type& x0, const vec_type& v0,
                             const vec_type& x1, const vec_type& v1,
                             const MAP& pmap, unsigned int period,
                             double dx, map_analysis_param& params);
                             
// this function may throw one of two exceptions:
// 1) rational_surface_suspected: indicating that the magnitude of the
//    p-map is vanishing smoothly along one of the cell's edges
// 2) higher_precision_required: meaning that we are unable to
//    track the smooth rotation of the p-map with the current
//    precision of the map. Of course, if the map is analytical, this
//    diagnostic is moot and one has to give up.
template<typename MAP>
int poincare_index(const grid_type& mesh, const ivec_type& cell_id,
                   const dataset_type& dataset,
                   const MAP& pmap, int period, double lmin,
                   map_analysis_param& params);
                   
template<typename MAP>
int __poincare_index(const grid_type& mesh,
                     const ivec_type& cell_id,
                     const dataset_type& dataset,
                     const MAP& pmap, int period, double lmin,
                     bool scrutinize,
                     std::vector<std::pair<int, vec_type> >& sing_edges,
                     map_analysis_param& params);
                     
template<typename MAP>
int __fast_poincare_index(const grid_type& mesh,
                          const ivec_type& cell_id,
                          const dataset_type& dataset,
                          const MAP& pmap, int period, double lmin,
                          bool scrutinize,
                          std::vector<std::pair<int, vec_type> >& sing_edges,
                          map_analysis_param& params);
                          
} //end xavier


//---------------------------------------------------------------------------------
//Template source implementation (Compiler needs access to this)
/** Compute rotation angle between two points (x0,x1) on an edge
 *  v0,v1 are the map vectors given the period
 *  Throws an index_step_size_underflow error if unable to resolve rotation
 */
template<typename MAP>
double xavier::adaptive_rotation_angle(const vec_type& x0, const vec_type& v0,
                                       const vec_type& x1, const vec_type& v1,
                                       const MAP& pmap, unsigned int period,
                                       double dx, map_analysis_param& params)
{

    if (nvis::norm(v0) < 1.0e-6 || nvis::norm(v1) < 1.0e-6) {
        index_step_size_underflow e("unable to track continuous vector rotation");
        e.where = 0.5*(x0 + x1);
        throw e;
    }
    
    double theta = xavier::signed_angle(v0, v1);
    
    bool verbose = params.verbose;
    
    if (fabs(theta) <= params.max_angle) {
        return theta;
    } else if (nvis::norm(x1 - x0) <= dx) {
        if (verbose) {
            std::cerr << "angle " << theta << " between " << v0 << " at " << x0 << " and " << v1 << " at " << x1
                      << " remains too large. failed\n";
        }
        
        // only bad things will come out of trying to proceed in such cases...
        index_step_size_underflow e("unable to track continuous vector rotation");
        e.where = 0.5*(x0 + x1);
        throw e;
    }
#ifndef SECANT_METHOD
    vec_type x = 0.5 * (x0 + x1);
#else
    double u = secant(v0, v1);
    if (std::isnan(u) || std::isinf(u)) {
        std::cerr << "secant method returned NaN" << std::endl;
    }
    vec_type x = (1-u)*x0 + u*x1;
#endif
    vec_type v = rhs(x, pmap, period, params);
    
    if (verbose) {
        std::cerr << x0[0] << ", " << x0[1] << ", " << v0[0] << ", " << v0[1] << '\n';
        std::cerr << x1[0] << ", " << x1[1] << ", " << v1[0] << ", " << v1[1] << '\n';
    }
    
    return adaptive_rotation_angle(x0, v0, x, v, pmap, period, dx, params) +
           adaptive_rotation_angle(x, v, x1, v1, pmap, period, dx, params);
}


template<typename MAP>
double xavier::robust_rotation_angle(const vec_type& x0, const vec_type& v0,
                                     const vec_type& x1, const vec_type& v1,
                                     const MAP& pmap, unsigned int period,
                                     double dx, map_analysis_param& params)
{
    double theta = xavier::signed_angle(v0, v1);
    bool verbose = params.verbose;
    if (fabs(theta) <= params.max_angle) {
        return theta;
    }
    
    vec_type x = 0.5 * (x0 + x1);
    vec_type v = rhs(x, pmap, period, params);
    
    // we stop when one of two things happens:
    // 1) the norm of the vector at the halfway point is less than a
    //    prescribed epsilon (dx)
    // 2) the rotation from left and right is not monotonic
    double norm = nvis::norm(v);
    if (norm < dx) {
        index_step_size_underflow e("unable to track continuous vector rotation");
        e.where = x;
        throw e;
    }
    
    return adaptive_rotation_angle(x0, v0, x, v, pmap, period, dx, params, verbose) +
           adaptive_rotation_angle(x, v, x1, v1, pmap, period, dx, params, verbose);
}

template<typename MAP>
int xavier::__poincare_index(const grid_type& mesh,
                             const ivec_type& cell_id,
                             const dataset_type& dataset,
                             const MAP& pmap, int period, double lmin,
                             bool scrutinize,
                             std::vector<std::pair<int, vec_type> >& sing_edges,
                             map_analysis_param& params)
{

    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
    
    vec_type v[5];
    vec_type x[5];
    
    sing_edges.clear();
    
    // corner values
    for (int i=0 ; i<4 ; ++i) {
        ivec_type delta(ids[i][0], ids[i][1]);
        x[i] = mesh(cell_id + delta);
        v[i] = p_vector(dataset(cell_id + delta), period, params.the_metric);
    }
    // close the loop for convenience
    x[4] = x[0];
    v[4] = v[0];
    
    double dtheta[4];
    for (int i=0 ; i<4 ; ++i) {
        dtheta [i] = signed_angle(v[i], v[i+1]);
    }
    
    // compute Poincare index of the map for current period
    double theta = 0;
    for (int edgeID=0 ; edgeID<4 ; ++edgeID) {
        if (dtheta[edgeID] < params.max_angle) {
            theta += dtheta[edgeID];
            continue;
        }
        
        try {
            theta += adaptive_rotation_angle(x[edgeID], v[edgeID], x[edgeID+1], v[edgeID+1],
                                             pmap, period, lmin, params);
        } catch (index_step_size_underflow& e) {
            if (scrutinize) {
                sing_edges.push_back(std::pair<int, vec_type>(edgeID, e.where));
                if (sing_edges.size() == 2) {
                    throw rational_surface_suspected("none");
                }
                continue; // no need to further look at this edge. it is broken anyway
            } else {
                return 0;    // refined tracking failed, too.
            }
        } catch(...) {
            return 0;
        }
    }
    if (sing_edges.size()) {
        throw higher_precision_required("none");
    }
    
    // theta is a valid angle. turn it into an index
    long int idx = lrint(0.5*theta/M_PI);
    if (params.verbose) {
        std::cerr << "integral index = " << idx << '\n';
    }
    return idx;
}

template<typename MAP>
int xavier::__fast_poincare_index(const grid_type& mesh,
                                  const ivec_type& cell_id,
                                  const dataset_type& dataset,
                                  const MAP& pmap, int period, double lmin,
                                  bool scrutinize,
                                  std::vector<std::pair<int, vec_type> >& sing_edges,
                                  map_analysis_param& params)
{

    typedef std::pair<vec_type, int>    data_type;
    
    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
    
    vec_type v[5];
    vec_type x[5];
    
    sing_edges.clear();
    
    for (int i=0 ; i<4 ; ++i) {
        ivec_type delta(ids[i][0], ids[i][1]);
        x[i] = mesh(cell_id + delta);
        v[i] = p_vector(dataset(cell_id + delta), period, params.the_metric);
    }
    x[5] = x[0];
    v[5] = v[0];
    
    // compute Poincare index of the map for current period
    double theta = 0;
    for (int edgeID=0 ; edgeID<4 ; ++edgeID) {
        try {
            theta += robust_rotation_angle(x[edgeID], v[edgeID], x[edgeID+1], v[edgeID+1],
                                           pmap, period, lmin, params);
        } catch (index_step_size_underflow& e) {
            if (scrutinize) {
                sing_edges.push_back(std::pair<int, vec_type>(edgeID, e.where));
                if (sing_edges.size() == 2) {
                    throw rational_surface_suspected("none");
                }
                continue; // no need to further look at this edge. it is broken anyway
            } else {
                return 0;    // refined tracking failed, too.
            }
        } catch(...) {
            return 0;
        }
    }
    if (sing_edges.size()) {
        throw higher_precision_required("none");
    }
    
    // theta is a valid angle. turn it into an index
    long int idx = lrint(0.5*theta/M_PI);
    if (params.verbose) {
        std::cerr << "integral index = " << idx << '\n';
    }
    return idx;
}

template<typename MAP>
int xavier::poincare_index(const grid_type& mesh, const ivec_type& cell_id,
                           const dataset_type& dataset,
                           const MAP& pmap, int period, double lmin,
                           map_analysis_param& params)
{

    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1}, {0,0} };
    
    std::vector<std::pair<int, vec_type> > sing_edges;
    
    try {
        return __poincare_index(mesh, cell_id, dataset, pmap, period, lmin,
                                true, sing_edges, params);
    } catch (rational_surface_suspected& e) {
        if (params.verbose) {
            std::cerr << "rational surface suspected in cell " << cell_id << std::endl;
            std::cerr << "first increase the precision\n";
        }
        sing_edges.clear();
        try {
            return __poincare_index(mesh, cell_id, dataset, pmap, period, lmin/10.,
                                    false, sing_edges, params);
        } catch(rational_surface_suspected& f) {
            if (params.verbose) {
                std::cerr << "rational surface suspicion remains at higher precision\n";
            }
            throw f;
        } catch (...) {
            if (params.verbose) {
                std::cerr << "Failed to track the orientation of the vector field. giving up.\n";
            }
            return 0;
        }
    } catch (higher_precision_required& e) {
        throw e;
    } catch(...) {
        return 0;
    }
}


#endif
