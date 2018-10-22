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


#ifndef __mapNd_hpp
#define __mapNd_hpp

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <stdexcept>
#include "metric.hpp"
#include <vector>

namespace xavier {

///Map Returns on section
template<typename T, int N>
struct return_map_info {
    typedef T                            value_type;
    typedef nvis::fixed_vector<T, N>     vec_type;
    typedef nvis::fixed_matrix<T, N>     mat_type;
    
    return_map_info()
        : x(0), J(0), t(0), nsteps(0), total_nsteps(0) {}
        
    return_map_info(const vec_type& _x, const mat_type& _J)
        : x(_x), J(_J), t(0), nsteps(0), total_nsteps(0) {}
        
    vec_type x;
    mat_type J;
    double t;
//    double delta_theta;
    ///Vector of angles being tracked
    std::vector<double> delta_theta;
    unsigned int nsteps;
    unsigned int total_nsteps;
};

/// ND State structure (for easy reference)
template<typename T, int N>
struct state_info {
    typedef T                               value_type;
    typedef nvis::fixed_vector<T, N>        vec_type;
    typedef nvis::fixed_vector<T, N* (1+N)>  xstate_type;
    typedef nvis::fixed_matrix<T, N>        mat_type;
    
    state_info()
        : x(0), J(0), t(0), nsteps(0), total_nsteps(0) {}
        
    state_info(const vec_type& _x, const mat_type& _J)
        : x(_x), J(_J), t(0), nsteps(0), total_nsteps(0) {}
        
    state_info(const xstate_type& _y, const double& _t)
        : t(_t), nsteps(0), total_nsteps(0)
    {
        for(int i=0; i<N; i++) {
            x[i] = _y[i];
        }
        //Jacobian - Assuming integrated in tandum as STM
        for(int ii=0; ii<N; ii++) {
            for(int jj=0; jj<N; jj++) {
                J[ii][jj] = _y[N + N*ii + jj];
            }
        }
    }
    
    void setState(const xstate_type& _y)
    {
        for(int i=0; i<N; i++) {
            x[i] = _y[i];
        }
        //Jacobian - Assuming integrated in tandum
        for(int ii=0; ii<N; ii++) {
            for(int jj=0; jj<N; jj++) {
                J[ii][jj] = _y[N + N*ii + jj];
            }
        }
    }
    
    xstate_type getState()
    {
        xstate_type y(0);
        for(int i=0; i<N; i++) {
            y[i] = x[i];
        }
        for(int ii=0; ii<N; ii++) {
            for(int jj=0; jj<N; jj++) {
                y[N + N*ii + jj] = J[ii][jj];
            }
        }
        return y;
    }
    
    vec_type x;
    mat_type J;
    double t;
    unsigned int nsteps;
    unsigned int total_nsteps;
};

/// Abstract Map class in ND
template<typename T, int N>
struct mapNd {
    typedef T                            value_type;
    typedef nvis::fixed_vector<T, N>     vec_type;
    typedef nvis::fixed_vector<bool, N>  bvec_type;
    typedef nvis::fixed_matrix<T, N>     mat_type;
    typedef nvis::bounding_box<vec_type> box_type;
    typedef return_map_info<T, N>        return_type;
    typedef mapNd<T, N>                  self_type;
    typedef metric<T, N>                 metric_type;
    
    static const int dim = N;
    
    virtual vec_type map(const vec_type& in, int niter) const = 0;
    virtual void map(const vec_type& in, std::vector<vec_type>& out, int niter) const = 0;
    
    virtual return_type map_complete(const vec_type& in, int niter) const = 0;
    virtual void map_complete(const vec_type& in, std::vector<return_type>& out, int niter) const = 0;
    
    virtual const box_type& bounds() const = 0;
    virtual bvec_type periodic() const = 0;
    
    /// Computing winding number if it is available explicitly through map information
    virtual double winding_number(const vec_type& x0, const std::vector<vec_type>& orbit) const = 0;
    
    /// Computing winding number if it is available explicitly through map information
    virtual double winding_number(const vec_type& x0, int n=100) const
    {
        std::vector<vec_type> orbit;
        map(x0, orbit, n);
        return winding_number(x0, orbit);
    }
    
    virtual self_type* clone() const = 0;
};

typedef mapNd<double, 2>    map2d;

}

#endif
