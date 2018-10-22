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


#ifndef __DP45wrapper_hpp__
#define __DP45wrapper_hpp__

#include <math/dopri5.hpp>
#include <math/fixed_vector.hpp>
#include <iostream>

// extremely lightweight wrapper of nvis::dopri5<T>
// to illustrate the implicit API of ODESolver in poincare_map

namespace xavier {

template<typename T, int N>
class dp5wrapper : public nvis::dopri5<nvis::fixed_vector<T, N> > {
public:
    typedef T                           value_type;
    typedef nvis::fixed_vector<T, N>    vector_type;
    typedef nvis::dopri5<vector_type>   int_type;
    
    dp5wrapper() : int_type() {}
    
    void set_init_cond(const vector_type& x, double t0)
    {
        this->y = x;
        this->t = t0;
    }
    
    void set_t_max(double _t)
    {
        this->t_max = _t;
    }
    
    void set_init_step(double _h)
    {
        this->h = _h;
    }
    
    void set_precision(double e)
    {
        this->reltol = this->abstol = e;
    }
};

}

#endif // __DP45wrapper_hpp__