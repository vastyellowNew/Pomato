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


#ifndef __CR3BP_TRACKER_HPP__
#define __CR3BP_TRACKER_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <math/angle.hpp>

namespace orbital {

/// A single angle tracking system - tracking the x-xdot rotation
/// N = coordinate dim (xstate_type, n=42), M = Num winding numbers (1 here)
template<int N, int M=1>
struct AngleTracker {
    typedef nvis::fixed_vector<double, N> gvec_type;
    static const int space_dim = M;
    
    AngleTracker(const nvis::vec2& center)
        : _center(center), _values(), _theta(0) {}
        
    void initialize(const gvec_type& v)
    {
        _last = nvis::vec2(v[0], v[3]) - _center;//x-xdot
        //_last = nvis::vec2(v[0], v[4]) - _center;//x-ydot
        //_last = nvis::vec2(v[3], v[4]) - _center;//xdot-ydot
        _theta = 0;
        _values.clear();
    }
    
    std::vector<double> operator()(const double& t, const gvec_type& v)
    {
        nvis::vec2 cur = nvis::vec2(v[0], v[3]) - _center;//x-xdot
        //nvis::vec2 cur = nvis::vec2(v[0],v[4]) - _center;//x-ydot
        //nvis::vec2 cur = nvis::vec2(v[3], v[4]) - _center;//xdot-ydot
        _theta += xavier::signed_angle(_last, cur);
        _last = cur;
        return std::vector<double>(1,_theta);
    }
    
    void mark_crossing()
    {
        _values.push_back(_theta);
    }
    
    /// Get a set of default winding factors (when integration hits a singularity)
    std::vector<double> getDefaultWindingFactors()
    {
        std::vector<double> factors;
        factors.push_back( -1.0 );
        return factors;
    }
    
    /// Output the winding factors (see sample_raster() in map_analysis.hpp)
    std::vector<double> getWindingFactors(uint vectorSize, std::vector<double>& delta_theta)
    {
        std::vector<double> factors;
        //q/p - poloidal/toroidal factors
        double w = 2.*M_PI*vectorSize/delta_theta[0];
        if (std::isnan(w) || std::isinf(w)) {
            w = 1000;
        }
        factors.push_back(w);
        
        return factors;
    }
    
    nvis::vec2 _center, _last;
    std::vector<double> _values;
    double _theta;
};

}

#endif
