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


/** General Error Handling classes/objects for various Poincare Map
    Topology Extraction procedures
    Author:  Wayne Schlei
 */

#ifndef __POINCARE_MAP_EXCEPTIONS_HPP
#define __POINCARE_MAP_EXCEPTIONS_HPP

#include <string>
#include <exception>
#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <maps/definitions.hpp>
#include <maps/mapNd.hpp>

namespace xavier {

//Note: vec_type is pre-defined as an nvis::vec2!!!!

template<typename VEC>
struct MapSingularityError : public std::runtime_error {
    MapSingularityError() :
        std::runtime_error("Singularity in EOMs"),
        where(0.0)
    {}
    MapSingularityError(const std::string& s) :
        std::runtime_error("Singularity" + s),
        where(0.0)
    {}
    MapSingularityError(const VEC& v) :
        std::runtime_error("Singularity in EOMs"),
        where(v)
    {}
    //Location of error
    VEC where;
};

//Poincare Map Engine---------------------------------------------
/// A Standard Error with a call to a Poincare Map object
struct MapUndefined : public std::runtime_error {
    MapUndefined() : std::runtime_error("Map Undefined") {}
    MapUndefined(const std::string& s) : std::runtime_error("Map Undefined" + s) {}
    ///Location of error (in section space)
    vec_type where;
};

///Failed integration within the poincare_map class
template <typename VEC, typename GVEC>
struct FailedIntegration : public std::runtime_error {
    ///Constructor
    explicit FailedIntegration(const std::string& __arg)
        : std::runtime_error(__arg) {}
    ///Error code
    int error_code;
    ///Independent variable
    double t;
    ///Local State vector ('Initial Position') of error
    VEC where;
    ///The Dependent variable at the error
    GVEC y;
};

/// Poincare Map Engine was unable to reach the desired number of iterations
struct MapUnderflow : public std::underflow_error {
    explicit MapUnderflow(const std::string& __arg)
        : std::underflow_error(__arg) {}
    vec_type where;
};

/// Unknown map error that returns location
struct MapUnknownError : public std::runtime_error {
    MapUnknownError() : std::runtime_error("Unknown Map Error") {}
    MapUnknownError(const std::string& s) : std::runtime_error("Unknown Map Error: " + s)
    {}
    ///Location of error (in section space)
    vec_type where;
};

// IndexComputation ------------------------------------------------
// -> Note:  These may need to be templatized!
struct index_step_size_underflow : public std::underflow_error {
    explicit index_step_size_underflow(const std::string& __arg)
        : std::underflow_error(__arg) {}
        
    vec_type where;
    double   step;
    double   theta;
};


struct ambiguous_index : public std::runtime_error {
    explicit ambiguous_index(const std::string& __arg)
        : std::runtime_error(__arg) {}
        
    vec_type x0, x1;
};


struct invalid_secant_value: public std::runtime_error {
    explicit invalid_secant_value(const std::string& __arg)
        : std::runtime_error(__arg) {}
        
    vec_type x0, x1;
    vec_type v0, v1;
    double   theta;
};

/// Error handling on edge if a fixed point is suspected
class FixedPoint_Suspected : public std::exception {
public :
    virtual const char* what() const throw()
    {
        return "Fixed Point Suspected";
    }
    vec_type where;
    int period;
    double displacementError;
    double theta;
};


} //xavier

#endif
