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


#ifndef __FIXED_POINT_TESTS_HPP__
#define __FIXED_POINT_TESTS_HPP__

#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <algorithm>
//xavierAPI
#include <math/fixed_vector.hpp>
#include <maps/definitions.hpp>
#include <maps/mapExceptions.hpp>
#include <maps/misc.hpp>
#include <maps/index.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/EdgeRotationFailure.hpp>


using namespace xavier;

namespace topology {



/// Test if the provided point is close to a fixed point for the given period
template <class VEC, class PARAM>
bool isFixedPointSuspected(
    const VEC& x0, const VEC& delta, const int& p,
    const PARAM& params, EdgeRotationFailure<VEC>& fobj )
{
    typedef EdgeRotationFailure<VEC> ERotFailure;
    //Fixed point test is if the norms are vanishing
    const double& MIN_FP_TOL = params.min_fp_tol;
    if ( nvis::norm(delta) <= MIN_FP_TOL ) {
        fobj.setType(ERotFailure::FIXED_POINT_SUSPECTED);
        fobj.period = p;
        fobj.failurePos = x0;
        return true;
    }
    return false;
}
/// Test if there's a fixed point between two adjacent points (Forward Displacement)
template <class VEC, class PARAM>
bool isFixedPointSuspected(
    const VEC& x0, const VEC& delta0,
    const VEC& x1, const VEC& delta1,
    const int& p, const PARAM& params, EdgeRotationFailure<VEC>& fobj )
{
    typedef EdgeRotationFailure<VEC> ERotFailure;
    //Compute signed angle, fixed points will usually only be at a flip
    double theta = xavier::signed_angle(delta0,delta1);
    fobj.theta = theta;
    //Likely not a fixed point if we don't have a flip in map displacement along edge
    if ( fabs(fabs(theta) - M_PI) > 15.*M_PI/180.) {
        return false;
    }
    
    //Fixed point test is if the norms are vanishing
    const double& MIN_FP_TOL = params.min_fp_tol;
    if ( nvis::norm(delta0) <= MIN_FP_TOL && nvis::norm(delta1)<= MIN_FP_TOL ) {
        fobj.setType(ERotFailure::FIXED_POINT_SUSPECTED);
        fobj.period = p;
        fobj.failurePos = (x0+x1)/2.0; //midpoint
        return true;
    }
    return false;
}


/// Test if the provided point is close to a fixed point for the given period (Map Tangent Version)
template <class VEC, class PARAM>
bool isFixedPointSuspectedTangent(
    const VEC& x0, const VEC& delta, const VEC& beta, const VEC& eta,
    const int& p, const PARAM& params, EdgeRotationFailure<VEC>& fobj )
{
    typedef EdgeRotationFailure<VEC> ERotFailure;
    //Fixed point test is if the norms are vanishing
    const double& MIN_FP_TOL = params.min_fp_tol;
    if ( nvis::norm(eta) <= MIN_FP_TOL ) {
        //Figure out if it's a double period
        if ( nvis::norm(delta) >= 100.0*MIN_FP_TOL && nvis::norm(beta) >= 100.0*MIN_FP_TOL ) {
            fobj.setType(ERotFailure::DOUBLE_PERIOD_FIXED_POINT);
            fobj.period = 2*p;
        } else {
            fobj.setType(ERotFailure::FIXED_POINT_SUSPECTED);
            fobj.period = p;
        }
        fobj.failurePos = x0;
        return true;
    }
    return false;
}
/// Test if there's a fixed point between two adjacent points (Map Tangent Version)
template <class VEC, class PARAM>
bool isFixedPointSuspectedTangent(
    const VEC& x0, const VEC& delta0, const VEC& beta0, const VEC& eta0,
    const VEC& x1, const VEC& delta1, const VEC& beta1, const VEC& eta1,
    const int& p, const PARAM& params, EdgeRotationFailure<VEC>& fobj )
{
    typedef EdgeRotationFailure<VEC> ERotFailure;
    //Compute signed angle, fixed points will usually only be at a flip
    double theta = xavier::signed_angle(eta0,eta1);
    fobj.theta = theta;
    //Likely not a fixed point if we don't have a flip in map displacement along edge
    if ( fabs(fabs(theta) - M_PI) > 15.*M_PI/180.) {
        return false;
    }
    
    //Fixed point test is if the norms are vanishing
    const double& MIN_FP_TOL = params.min_fp_tol;
    if ( nvis::norm(eta0) <= MIN_FP_TOL && nvis::norm(eta1)<= MIN_FP_TOL ) {
        //Figure out if it's a double period
        if ( ( nvis::norm(delta0) >= 1.e2*MIN_FP_TOL && nvis::norm(beta0) >= 1.e2*MIN_FP_TOL ) ||
                ( nvis::norm(delta1) >= 1.e2*MIN_FP_TOL && nvis::norm(beta1) >= 1.e2*MIN_FP_TOL ) ) {
            fobj.setType(ERotFailure::DOUBLE_PERIOD_FIXED_POINT);
            fobj.period = 2*p;
        } else {
            fobj.setType(ERotFailure::FIXED_POINT_SUSPECTED);
            fobj.period = p;
        }
        fobj.failurePos = (x0+x1)/2.0; //midpoint
        return true;
    }
    return false;
}

} //end topology

#endif