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


#ifndef __XAVIER_METRIC_HPP__
#define __XAVIER_METRIC_HPP__

#include <vector>
#include <limits>
#include <map>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <utility>
#include <assert.h>

namespace xavier {
template<typename T, int N>
class metric {
    static T __modulo(T a, T b)
    {
        T r = fmod(a, b);
        if (!r) {
            return 0;
        } else {
            return a >= 0 ? r : b + r;
        }
    }
    
    static T __min_norm(T a, T b)
    {
        if (2*a > b) {
            return a - b;
        } else if (2*a + b < 0) {
            return a + b;
        } else {
            return a;
        }
    }
    
    static T sign(T a)
    {
        return (a < 0) ? -1 : 1;
    }
    
public:
    typedef T                                   scalar_type;
    typedef nvis::fixed_vector<scalar_type, N>  vec_type;
    typedef nvis::fixed_vector<bool, N>         bvec_type;
    typedef nvis::bounding_box<vec_type>        bounds_type;
    
    metric() :
        __bounds(vec_type(0), vec_type(0)), __periodic(false) {}
    metric(const bounds_type& bounds) :
        __bounds(bounds), __periodic(false) {}
    metric(const bounds_type& bounds, const bvec_type& periodic) :
        __bounds(bounds), __periodic(periodic) {}
        
    const bounds_type& bounds() const
    {
        return __bounds;
    }
    
    bounds_type& bounds()
    {
        return __bounds;
    }
    
    const bvec_type& periodic() const
    {
        return __periodic;
    }
    
    bvec_type& periodic()
    {
        return __periodic;
    }
    
    vec_type size() const
    {
        return __bounds.size();
    }
    
    scalar_type diameter() const
    {
        return nvis::norm(size());
    }
    
    static scalar_type periodic_subtraction(scalar_type x0, scalar_type x1, scalar_type width)
    {
        scalar_type dx = x1 - x0;
        scalar_type d = sign(dx) * fmod(fabs(dx), width);
        return __min_norm(d, width);
    }
    
    vec_type modulo(const vec_type& x) const
    {
        vec_type y = x - __bounds.min();
        const vec_type& sz = size();
        for (int i = 0 ; i < N ; ++i) {
            if (__periodic[i]) {
                y[i] = __modulo(y[i], sz[i]);
            }
        }
        return y + __bounds.min();
    }
    
    vec_type displacement(const vec_type& a, const vec_type& b) const
    {
        vec_type d = modulo(b) - modulo(a);
        const vec_type& sz = size();
        for (int i = 0 ; i < N ; ++i) {
            if (__periodic[i]) {
                d[i] = __min_norm(d[i], sz[i]);
            }
        }
        return d;
    }
    
    scalar_type distance(const vec_type& a, const vec_type& b) const
    {
        return nvis::norm(displacement(a, b));
    }
    
    scalar_type angle(const vec_type& x0, const vec_type& x1, const vec_type& x2) const
    {
        vec_type v0 = displacement(x0, x1);
        vec_type v1 = displacement(x1, x2);
        v0 /= nvis::norm(v0);
        v1 /= nvis::norm(v1);
        scalar_type alpha = acos(nvis::inner(v0, v1));
        if (N==2) {
            // if we are in 2D, angle sign makes sense, otherwise not
            scalar_type sin_alpha = v0[0] * v1[1] - v0[1] * v1[0];
            if (sin_alpha < 0) {
                alpha *= -1;
            }
        }
        return alpha;
    }
    
    void clip_segment(std::vector<std::pair<vec_type, vec_type> >& segments,
                      const vec_type& a, const vec_type& b) const
    {
        segments.clear();
        vec_type d = displacement(a, b);
        vec_type __a = modulo(a);
        if (!__bounds.inside(__a)) {
            return;
        }
        if (__bounds.inside(__a + d)) {
            segments.resize(1);
            segments[0] = std::make_pair(__a, __a + d);
        } else {
            std::vector<scalar_type> t(N, std::numeric_limits<scalar_type>::max());
            vec_type __b = modulo(b);
            vec_type __x, __y;
            vec_type dir = d / nvis::norm(d);
            // compute intersections with boundary
            for (int i = 0 ; i < N ; ++i) {
                scalar_type w = (dir[i] > 0 ? __bounds.max()[i] : __bounds.min()[i]);
                t[i] = (w - __a[i]) / dir[i];
            }
            // which intersection is closest?
            int imin = std::distance(t.begin(), std::min_element(t.begin, t.end()));
            __x = __a + t[imin] * dir;
            // split segment accordingly
            if (__periodic[imin]) {
                __y = __x;
                __y[imin] = (dir[imin] > 0 ? __bounds.min()[imin] : __bounds.max()[imin]);
                segments.resize(2);
                segments[0] = std::make_pair(__a, __x);
                segments[1] = std::make_pair(__y, __b);
            }
            // if not periodic along that direction, simply discard second half
            else {
                segments.resize(1);
                segments[0] = std::make_pair(__a, __x);
            }
        }
    }
    
private:
    bounds_type __bounds;
    bvec_type   __periodic;
    
};

static metric<double, 2> __default_metric;

} // namespac xavier

#endif