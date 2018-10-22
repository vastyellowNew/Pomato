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


/** Close Approach distance and Time tracking system
 *  Author: Wayne Schlei (Purdue University)
 *  Use:  For tracking minimum close approach distances to singularities
 *  and propagation times between iterations during a call to poincare_map::map().
 *  This is typically only employed within the evaulation of the Poincare Index
 *  (i.e., namely within Edge_Rotation functions).
 *
 *  Outputs in CR3BP, for example:
 *  CA1: Minimum close approach distance to P1 between returns
 *  CA2: Minimum close approach distance to P2 between returns
 *  T  : Propagation time between returns (if backward propagation, value is negative)
 *
 *  Future work:  This needs some implementation work such that it is generalized
 *  to problems other than CR3BP!
 *
 */
#ifndef __CAT_TRACKER_HPP__
#define __CAT_TRACKER_HPP__

#include <vector>
#include <math/fixed_vector.hpp>

namespace topology {

///Close Approach distance and Time tracking system [CA1 CA2 ... CAM T]
template<typename V, int NUMSINGULARITIES=2>
struct CATtracker {
    typedef V                                              gvec_type;
    typedef nvis::fixed_vector<double, NUMSINGULARITIES+1> DataVec;
    typedef std::vector<double>                            DataOutput;
    static const int numDataValues = NUMSINGULARITIES+1;
    
    CATtracker() : _values(), _data(0)
    {
        for (int i=0; i<NUMSINGULARITIES; i++) {
            _center[i] = nvis::vec2(0,0);
        }
    }
    
    CATtracker(const nvis::vec2& center)
        : _values(), _data(0)
    {
        for (int i=0; i<NUMSINGULARITIES; i++) {
            _center[i] = center;
        }
    }
    
    void setSingularity(const int idx, const nvis::vec2& center)
    {
        _center[idx] = center;
    }
    
    
    void initialize(const gvec_type& v)
    {
        for(int i=0; i<NUMSINGULARITIES; i++) {
            //CA - Pi
            _lastR[i] = nvis::vec2(v[0], v[1]) - _center[i];
            _data[i] = nvis::norm(_lastR[i]);
        }
        //Time
        _lastTime = 0.0;
        _data[NUMSINGULARITIES]  = 0; //Summed time from last crossing
        
        //Note:  HARD-CODED for the Planar CR3BP (2singularities, 2D positions)
        _values.clear();
    }
    
    /// Output the vector of tracked quantities
    DataOutput operator()(const double& t, const gvec_type& v)
    {
        //Convert _theta to angle_output
        DataOutput out;
        for (int i=0; i<NUMSINGULARITIES; i++) {
            //CA Radius & distance
            _lastR[i] = nvis::vec2(v[0], v[1]) - _center[i];
            if ( nvis::norm(_lastR[i]) < _data[i]) {
                _data[i] = nvis::norm(_lastR[i]);
            }
            out.push_back( _data[i] );
        }
        //Time - elapsed since last crossing
        _data[NUMSINGULARITIES] = t - _lastTime;
        
        return out;
    }
    
    /// Update values when map crossing occurs
    void mark_crossing()
    {
        //Store value vector per crossing
        _values.push_back( _data );
        //Reset values with current information
        for(int i=0; i<NUMSINGULARITIES; i++) {
            _data[i] = nvis::norm(_lastR[i]);
        }
        _lastTime += _data[NUMSINGULARITIES]; //sum total
    }
    
    /// Gather the output
    std::vector<DataVec>& getResult()
    {
        return _values;
    }
    
    nvis::vec2 _center[NUMSINGULARITIES];
    nvis::vec2 _lastR[NUMSINGULARITIES];
    double _lastTime;
    std::vector<DataVec> _values;
    DataVec _data;
};

}


#endif
