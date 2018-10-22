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


/* PeriodicOrbit class object
 * Author: Wayne Schlei (Purdue University)
 */
#ifndef PERIODIC_ORBIT_HPP
#define PERIODIC_ORBIT_HPP


#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
//CR3BP objects
//#include <State.h>
#include <design/Trajectory.hpp>

using namespace xavier;

namespace pmateDesign {

/** A specific Trajectory object representing a periodic orbit.
 *  This adds a simple implementation of a linear parameter (alpha)
 *  ranging from [0,1].
 *
 *  Note:  The last state must be a duplicate of the first state
 *  but at the total period of the orbit such that the 'deltaT'
 *  member variable represents the orbital time period.
 */
template <typename STATE>
class PeriodicOrbit : public Trajectory<STATE> {
public:
    typedef STATE                  value_type;
    typedef PeriodicOrbit<STATE>   SelfType;
    typedef Trajectory<STATE>      BaseType;
    
    /// Constructor
    PeriodicOrbit() : BaseType(), fpIndex(-1) {}
    
    /// Constructor with states and times
    PeriodicOrbit(const std::vector<STATE>& in, const std::vector<double>& t) :
        BaseType(in,t), fpIndex(-1) {}
        
    /// Copy constructor (a bit tricky)
    PeriodicOrbit(const SelfType& other) :
        BaseType(), fpIndex(other.fpIndex)
    {
        BaseType::setValidity( other.isValid() );
        if(other.getNumStates() == 0) {
            return;
        }
        BaseType::setStates( other.states, other.times );
    }
    
    /// Modulate linear parameter to be between [0,1]
    double moduloAlpha(const double& a) const
    {
        double alpha = a;
        if(alpha>=0.0 && alpha < 1.0) {
            return alpha;
        } else if (alpha >= 1.0) {
            //Recursive call
            return moduloAlpha( alpha - 1.0 );
        } else {
            //Below 0, recursive call
            return moduloAlpha( alpha + 1.0 );
        }
    }
    
    /// Get the time given a linear parameter on [0,1] representing time
    double getTime(const double& alpha) const
    {
        return (moduloAlpha(alpha)) * BaseType::getTimeOfFlight();
    }
    
    /// Get a state given a linear parameter on [0,1] representing time
    void getStateFromAlpha(const double& alpha, STATE& out) const
    {
        double t = getTime(alpha);
        BaseType::getState( t, out);
    }
    
    /// Get the periodic offset (||xf - x0||)
    double getPeriodicError(const int stateDim=6) const
    {
        double e = 0.0;
        for(int i=0; i<stateDim; i++) {
            double dx = 0.0;
            dx = BaseType::states.back()[i] - BaseType::states[0][i];
            e += dx*dx;
        }
        return sqrt(e);
    }
    
    /// Index in Fixed Point Data
    int fpIndex;
};



} //end pmateDesign


#endif // PERIODIC_ORBIT_HPP