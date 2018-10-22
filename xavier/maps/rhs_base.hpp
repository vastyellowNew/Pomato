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


//Right-hand side abstract base class
//Author: Wayne Schlei & Xavier Tricoche
#ifndef __rhs_base_hpp
#define __rhs_base_hpp

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <stdexcept>

/** Abstract base class for Right-Hand Side objects
 *
 *  Sub-classes of this type represent the differential equations of a dynamical system.
 *  They can be analytical or numerical (vector-field) EOMs
 *
 *  Please provide:
 *  1)The EOMs with STM under operator()
 *  2)The EOMs for states under value()
 *  3)The Jacobian (A) matrix for the EOMs under derivative()
 *  4)A way to evaluate the Hamiltonian
 *  5)The first partial derivatives of the hamiltonian wrt the state
 *  6)A constant value for time-invariant systems under desired_hamiltonian()
 *  7)Any singularities of the problem (template param S)
 */
template<typename T, int N, int S=0>
class rhs_base {
public:
    typedef T                               value_type;
    typedef nvis::fixed_vector<T, N>        state_type;
    typedef nvis::fixed_vector<T, N* (1+N)>  xstate_type;
    typedef nvis::fixed_matrix<T, N>        mat_type;
    typedef nvis::bounding_box<state_type>  box_type;
    typedef rhs_base<T, N, S>               self_type;
    
    static const int dimension = N;
    static const int numSingularities = S;
    
    struct undefined_point : public std::runtime_error {
        undefined_point() : std::runtime_error( "undefined point" )
        {
        }
    };
    
    virtual xstate_type operator()(const value_type& t, const xstate_type& y) const = 0;
    virtual state_type value(const value_type& t, const state_type& y) const = 0;
    virtual mat_type derivative(const value_type& t, const state_type& y) const = 0;
    
    /// Energy constant (i.e., Hamiltonian for conservative, time-invariant systems)
    virtual value_type hamiltonian(const value_type& t, const xstate_type& y) const = 0;
    /// Energy constant value (i.e., Hamiltonian value for analysis)
    virtual value_type desired_hamiltonian() const = 0;
    /// First order partials of Energy Constant evaluated at a state
    virtual state_type hamiltonian_first_partials(const value_type& t, const xstate_type& y) const = 0;
    /// Array of singularities unique to the problem (if none, return NULL)
    virtual const xstate_type* singularities() const = 0;
    /// Number of singularities in EOMs
    virtual int getNumSingularities() const = 0;
    /// Get a minimum numerically safe distance from a specific singularity (on map vector displacement)
    virtual value_type getSingularitySafeDistance(const int& idx) const = 0;
    
    virtual self_type* clone() const = 0;
    
    virtual const box_type& bounds() const = 0;
};

#endif // __rhs_base_hpp
