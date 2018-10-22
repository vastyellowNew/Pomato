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


/////////////////////////////////////////
// A class for the Standard Map        //
/////////////////////////////////////////
#ifndef __standard_map_hpp
#define __standard_map_hpp

#include <vector>
#include <exception>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdexcept>
#include <map>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include "mapNd.hpp"
#include "metric.hpp"
#include <math.h>

namespace xavier {

class standard_map : public map2d {
    //static const double 2.0*M_PI = 6.28318530717958647688;
    
    inline void _fwd(double& x, double& y) const
    {
        y -= _k / 2.0*M_PI * sin(2.0*M_PI * x);
        x += y;
    }
    
    inline nvis::mat2 fwd_jac(double x, double y) const
    {
        nvis::mat2 J;
        J(0,0) = 1 - _k*cos(2.0*M_PI * x);
        J(0,1) = 1;
        J(1,0) = -_k*cos(2.0*M_PI*x);
        J(1,1) = 1;
        return J;
    }
    
    inline void _bwd(double& x, double& y) const
    {
        x -= y;
        y += _k / 2.0*M_PI * sin(2.0*M_PI * x);
    }
    
    inline nvis::mat2 bwd_jac(double x, double y) const
    {
        nvis::mat2 J;
        J(0,0) = 1;
        J(0,1) = -1;
        J(1,0) = _k*cos(2.0*M_PI*(x-y));
        J(1,1) = 1 - _k*cos(2.0*M_PI*(x-y));
        return J;
    }
    
    inline nvis::vec2 _map(const nvis::vec2& x, bool fwd = true) const
    {
        nvis::vec2 y = x;
        if (fwd) {
            _fwd(y[0], y[1]);
        } else {
            _bwd(y[0], y[1]);
        }
        return y;
    }
    
    inline std::pair<nvis::vec2,nvis::mat2> _map_and_jacobian(const nvis::vec2& x, bool fwd = true) const
    {
        std::pair<nvis::vec2, nvis::mat2> r;
        r.first = x;
        if (fwd) {
            r.second = fwd_jac(x[0], x[1]);
            _fwd(r.first[0], r.first[1]);
        } else {
            r.second = bwd_jac(x[0], x[1]);
            _bwd(r.first[0], r.first[1]);
        }
        return r;
    }
    
public:

    standard_map(double k) : _k(k), _bounds(vec_type(0,0),vec_type(1,1)) {}
    
    const double& parameter() const
    {
        return _k;
    }
    
    nvis::vec2 map(const nvis::vec2& x, int niter) const
    {
        nvis::vec2 y = x;
        for (unsigned int i = 0 ; i < abs(niter) ; ++i) {
            y = _map(y, (niter > 0));
        }
        return y;
    }
    
    void map(const nvis::vec2& x, std::vector< nvis::vec2 >& hits, int niter) const
    {
        hits.resize(abs(niter));
        nvis::vec2 y = x;
        for (unsigned int i = 0; i < (unsigned int) abs(niter) ; ++i) {
            y = _map(y, (niter > 0));
            hits[i] = y;
        }
    }
    
    return_type map_complete(const nvis::vec2& x, int niter) const
    {
        nvis::mat2 last, cur;
        nvis::vec2 y = x;
        for (unsigned int i = 0; i < (unsigned int) abs(niter) ; ++i) {
            std::pair<nvis::vec2, nvis::mat2> res = _map_and_jacobian(y, (niter > 0));
            cur = res.second;
            if (i>0) {
                cur = cur * last;
            }
            last = cur;
            y = res.first;
        }
        return return_type(y, last);
    }
    
    void map_complete(const nvis::vec2& x, std::vector<return_type>& out, int niter) const
    {
        out.resize(abs(niter));
        nvis::vec2 y = x;
        for (unsigned int i = 0; i < (unsigned int) abs(niter) ; ++i) {
            std::pair<nvis::vec2, nvis::mat2> res = _map_and_jacobian(y, (niter > 0));
            out[i].x = res.first;
            out[i].J = res.second;
            if (i>0) {
                out[i].J = out[i].J * out[i-1].J;
            }
            y = out[i].x;
        }
    }
    
    ///Return the bounds for (x,y) => ([0,1],[0,1])
    //nvis::bbox2 bounds() const {
    //    return nvis::bbox2(nvis::vec2(0,0), nvis::vec2(1,1));
    //}
    const box_type& bounds() const
    {
        return _bounds;
    }
    
    nvis::fixed_vector<bool, 2> periodic() const
    {
        return nvis::fixed_vector<bool, 2>(true, true);
    }
    
    standard_map* clone() const
    {
        return new standard_map(*this);
    }
    
    double winding_number(const nvis::vec2& x0, const std::vector<nvis::vec2>& orbit) const
    {
        double w = 0;
        double number = 0;
        for (int i=0 ; i<(int)orbit.size() ; ++i) {
            double d = fabs(metric_type::periodic_subtraction(orbit[i][0], x0[0], 1));
            if (d == 0) {
                return (orbit[i][0] - x0[0])/(float)(i+1);
            } else {
                w += 1./d; // weight is inverse propertional to periodic distance to x0
                number += (orbit[i][0] - x0[0])/(float)(i+1)/d;
            }
        }
        return number/w;
    }
    
private:
    //Chaos Parameter
    double _k;
    //Bounding Box
    box_type _bounds;
};

}

#endif //__standard_map_hpp