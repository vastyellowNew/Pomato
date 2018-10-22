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


#ifndef __MULTIPLE_ANGLE_TRACKER_HPP__
#define __MULTIPLE_ANGLE_TRACKER_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <math/angle.hpp>

namespace orbital {

template<int N, int M>
struct MultipleAngleTracker {
    typedef nvis::fixed_vector<double, N> gvec_type;
    typedef nvis::fixed_vector<double, M> angle_type;
    typedef std::vector<double> angle_output;
    static const int space_dim = M;
    
    MultipleAngleTracker() : _values(), _theta(0)
    {
        for (int i=0; i<M; i++) {
            _center[i] = nvis::vec2(0,0);
        }
    }
    
    MultipleAngleTracker(const nvis::vec2& center)
        : _values(), _theta(0)
    {
        for (int i=0; i<M; i++) {
            _center[i] = center;
        }
    }
    
    void setCenter(const int idx, const nvis::vec2& center)
    {
        _center[idx] = center;
    }
    void setCenter(const nvis::vec2& center)
    {
        for (int i=0; i<M; i++) {
            _center[i] = center;
        }
    }
    void setCenter(const nvis::vec2* centerPtr)
    {
        //if you have a set of centers
        for (int i=0; i<M; i++) {
            _center[i] = centerPtr[i];
        }
    }
    
    void initialize(const gvec_type& v)
    {
        //x-xdot
        _last[0] = nvis::vec2(v[0], v[3]) - _center[0];
        _theta[0] = 0;
        //x-ydot
        _last[1] = nvis::vec2(v[0], v[4]) - _center[1];
        _theta[1] = 0;
        //xdot-ydot
        _last[2] = nvis::vec2(v[3], v[4]) - _center[2];
        _theta[2] = 0;
        
        //Note:  What state dimensions to check may be tied to the Section
        //In planarCR3BP y=0 map, the relevant angles to check for winding
        //numbers are the dimensions without y.
        
        //Maybe template this with
        //a SECTION and add functions to get parallel and perpendicular
        //dimensions to the section.
        
        _values.clear();
    }
    
    /// Output the vector of tracked angles -> feeds to return_map_info.delta_theta
    angle_output operator()(const double& t, const gvec_type& v)
    {
        //Convert _theta to angle_output
        angle_output out;
        
        for (int i=0; i<M; i++) {
            nvis::vec2 cur;
            if (i==0) {
                //x-xdot
                cur = nvis::vec2(v[0], v[3]) - _center[i];
            } else if (i==1) {
                //x-ydot
                cur = nvis::vec2(v[0], v[4]) - _center[i];
            } else {
                //xdot-ydot
                cur = nvis::vec2(v[3], v[4]) - _center[i];
            }
            _theta[i] += xavier::signed_angle(_last[i], cur);
            out.push_back( _theta[i] );
            _last[i] = cur;
        }
        
        return out;
    }
    
    void mark_crossing()
    {
        _values.push_back(_theta);
    }
    
    /// Get a set of default winding factors (when integration hits a singularity)
    std::vector<double> getDefaultWindingFactors()
    {
        std::vector<double> factors;
        for(int i=0; i<M; i++) {
            factors.push_back( -1.0 );
        }
        return factors;
    }
    
    
    /// Output the winding factors (see sample_raster() in map_analysis.hpp)
    std::vector<double> getWindingFactors(unsigned int vectorSize, std::vector<double>& delta_theta)
    {
        std::vector<double> factors;
        //p/q - toroidal/poloidal factors
        for (int i=0; i<M; i++) {
            double w = 2.*M_PI*vectorSize/delta_theta[i]; //Winding factor (i)
            bool isnanResult = ISNANCALL(w);
            bool isinfResult = ISINFCALL(w);
            if (isnanResult || isinfResult) {
                w = 1000;
            }
            factors.push_back(w);
        }
        return factors;
    }
    
    
    nvis::vec2 _center[space_dim];
    nvis::vec2 _last[space_dim];
    std::vector<angle_type> _values;
    angle_type _theta;
};

}


/** WindingVecTraits mechanism to indicate unknown or invalid values
 */
struct WindingVecTraits {
public:
    WindingVecTraits() {}
    //Static members signifying "invalid" values and "unknown" values
    static const std::vector<double> invalid;
    static const std::vector<double> unknown;
};


#endif
