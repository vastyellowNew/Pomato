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


// AdaptiveEdgeRotation : Functions for computing Poincare index components per edge
// Author:  Wayne Schlei

#ifndef __ADAPTIVE_EDGE_ROTATION_HPP__
#define __ADAPTIVE_EDGE_ROTATION_HPP__

#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <algorithm>
//xavierAPI
#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
#include <data/edge.hpp>
#include <maps/definitions.hpp>
#include <maps/mapExceptions.hpp>
#include <maps/misc.hpp>
#include <maps/index.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/AdaptiveEdge.hpp>
#include <topology/FixedPointTests.hpp>
#include <topology/SectionTransversality.hpp>
#include <topology/EdgeRotationMapCalls.hpp>

using namespace xavier;

namespace topology {

/** Tracking Conditions for Edge Rotation
 *  Returns true if the supposed step satisfies edge rotation constraints:
 *    - Smoothness, No fixed points, No Section Separation for non-transverse sections
 */
template <class MAP, class PARAM, class SRDATA>
bool getEdgeRotationTrackingCondition(
    const vec_type& x0, const vec_type& v0,
    const vec_type& x1, const vec_type& v1,
    const MAP& pmap, const int& period, const PARAM& params,
    std::set<SRDATA>& cache
)
{
    typedef typename MAP::lvec_type vec_type;
    // Note:  Inputs v0,v1 correspond to the map-displacement vectors
    // Thresholds in params
    const double& MAX_ANGLE_VARIATION = params.max_angle;
    
    //The main problem : Rotation angle is below threshold (smoothness test)
    double theta = xavier::signed_angle(v0,v1);
    if ( !(fabs(theta) <= MAX_ANGLE_VARIATION) ) {
        return false;
    }
    
    //SECTION TRANSVERSALITY:  We also must check if the section is transverse to the flow
    bool forwardFail  = detectSectionSeparation<MAP,PARAM,SRDATA>(x0,x1,cache,pmap,period,params);
    if(forwardFail) {
        return false;
    }
    
    //The segment passes the smoothness & transversality tests.
    return true;
}

/// Test angled rotation value between two vectors
template <class MAP, class PARAM>
bool runAngleTest( const vec_type& v0, const vec_type& v1, const MAP& pmap, const PARAM& params)
{
    typedef typename MAP::lvec_type vec_type;
    // Note:  Inputs v0,v1 correspond to the map-displacement
    // Thresholds in params
    const double& MAX_ANGLE_VARIATION = params.max_angle;
    
    //The main problem : Rotation angle is below threshold (smoothness test)
    double theta = xavier::signed_angle(v0,v1);
    if ( !(fabs(theta) <= MAX_ANGLE_VARIATION) ) {
        return false;
    }
    return true;
}



/** Tracking Conditions for Edge Rotation
 *  Returns true if the supposed step satisfies edge rotation constraints:
 *    - Smoothness, No fixed points, No Section Separation for non-transverse sections
 *    - Tests the Map Tangent vector
 */
template <class MAP, class PARAM, class SRDATA>
bool getEdgeRotationTrackingConditionTangent(
    const vec_type& x0, const vec_type& v0,
    const vec_type& x1, const vec_type& v1,
    const MAP& pmap, const int& period, const PARAM& params,
    std::set<SRDATA>& cache, std::set<SRDATA>& backCache
)
{
    typedef typename MAP::lvec_type vec_type;
    // Note:  Inputs v0,v1 correspond to the map-tangent vectors
    // Thresholds in params
    const double& MAX_ANGLE_VARIATION = params.max_angle;
    
    //The main problem : Rotation angle is below threshold (smoothness test)
    double theta = xavier::signed_angle(v0,v1);
    if ( !(fabs(theta) <= MAX_ANGLE_VARIATION) ) {
        return false;
    }
    
    //SECTION TRANSVERSALITY:  We also must check if the section is transverse to the flow
    bool forwardFail  = detectSectionSeparation<MAP,PARAM,SRDATA>(x0,x1,cache,pmap,period,params);
    bool backwardFail = detectSectionSeparation<MAP,PARAM,SRDATA>(x0,x1,backCache,pmap,-period,params);
    if(forwardFail || backwardFail) {
        return false;
    }
    
    //The segment passes the smoothness & transversality tests.
    return true;
}


/** Compute rotation angle between two points (x0,x1) on an edge
   v0,v1 are the map vectors given the period

   -Throws an EdgeRotationFailure exception when the rotation can't be resolved

   -Differs from previous version in that it stores and utilizes map iterates
    for other p values (p is descending in order).

   -There is a possibility that an "Unresolved" rotation can still occur without
    detecting a separation condition.  The main instance of this is a Stable
        Manifold crossing (and these are marked as such with the specific
        EdgeRotationFailure type, but no error is thrown).  Other cases may exist.

   -Note: IF using v as the "tangent vector" which is an approximation
    of the map tangent if map curves were continuous.  So
      v = (P^p(x) - P^-p(x))/2;
      ->Kinda like low-pass filter or Central-Difference Approx
      ->This works great when you don't have section tangency violations
 */
template<class MAP,class PARAM,class SRDATA>
double adaptiveRotationAngle(
    const vec_type& x0, const vec_type& v0, const vec_type& x1, const vec_type& v1,
    const MAP& pmap, const int period, const PARAM& params,
    std::set<SRDATA>& cache, std::set<SRDATA>& backCache,
    std::vector<EdgeRotationFailure<typename MAP::lvec_type> >& smcFailVec
)
{
    typedef EdgeRotationFailure<typename MAP::lvec_type> ERotFailure;
    
    //Evaluate the Angle displacement
    double theta = xavier::signed_angle(v0, v1);
    bool verbose = params.verbose;
    const double& lmin = params.lmin;
    
    //Compute the forward and backward mappings for the endpoints
    // ->they should exist except when singularity at node, which throws discontinuity
    ERotFailure eFail;
    vec_type delta0, delta1;
    try {
        delta0 = mapDisplacementUsingCache(x0, pmap, period, params, cache, eFail);
    } catch(ERotFailure& ef) {
        delta0 = vec_type(50);
        throw ef;
    }
    try {
        delta1 = mapDisplacementUsingCache(x1, pmap, period, params, cache, eFail);
    } catch(ERotFailure& ef) {
        delta1 = vec_type(50);
        throw ef;
    }
    /*vec_type beta0, beta1;
    try {
        beta0  = mapDisplacementUsingCache(x0, pmap, -period, params, backCache, eFail);
    } catch(ERotFailure& ef) {
        beta0 = vec_type(50);
        throw ef;
    }
    try {
        beta1  = mapDisplacementUsingCache(x1, pmap, -period, params, backCache, eFail);
    } catch(ERotFailure& ef) {
        beta1 = vec_type(50);
        throw ef;
    }*/
    
    
    /// Check if we have success (i.e., the tracking conditions are all satisfied)
    if ( // --  Tracking Conditions = OK ? ---
        getEdgeRotationTrackingCondition<MAP,PARAM>(x0,v0,x1,v1,pmap,period,params,cache)
        //getEdgeRotationTrackingCondition<MAP,PARAM>(x0,v0,x1,v1,pmap,period,params,cache,backCache)
    ) {
        return theta;
    }
    /// If not successful, check if it's a section separation condition (have to test both forward and backward)
    /*else if (detectSectionSeparation(x0,x1,backCache,pmap,period,params) ) {
      // -> Tell the level above to switch to the non-transverse version
      ERotFailure err(ERotFailure::BACKWARD_SECTION_SEP, 0.5*(x0+x1), period);
      err.theta = theta;
      if (verbose) {
            std::cerr << "Period (" << period << "): Detected BACKWARD-time section separation between "
              << v0 << " at " << x0 << " and " << v1 << " at " << x1
                     << ".\n Executing Non-transverse version of adaptive rotation...\n";
      }
      throw err;
    }*/
    else if (detectSectionSeparation(x0,x1,cache,pmap,period,params) ) {
        // -> Tell the level above to switch to the non-transverse version
        ERotFailure err(ERotFailure::FORWARD_SECTION_SEP, 0.5*(x0+x1), period);
        err.theta = theta;
        /*if (verbose) {
              std::cerr << "Period (" << period << "): Detected FORWARD-time section separation between "
                << v0 << " at " << x0 << " and " << v1 << " at " << x1
                       << ".\n Executing Non-transverse version of adaptive rotation...\n";
        }*/
        throw err;
    }
    /// If no success and at the smallest allowed step, find out what went wrong...
    else if (nvis::norm(x1 - x0) <= lmin) {
        //We have a failure, so we must figure out what it is and fix it if we can.
        // -> These should NOT be map failures as those are checked elsewhere.
        if (verbose) {
            std::cerr << "Period (" << period << "): angle " << theta << " between " << v0 << " at " << x0
                      << " and " << v1 << " at " << x1 << " remains too large. failed\n";
        }
        
        //Check if we have fixed point conditions
        ERotFailure fperr;
        if( isFixedPointSuspected(x0,delta0,x1,delta1,period,params,fperr) ) {
            //if( isFixedPointSuspectedTangent(x0,delta0,beta0,v0,x1,delta1,beta1,v1,period,params,fperr) ) {
            if(verbose) {
                std::cout << "   Edge_Rot: Fixed Point Detected = " << eFail.where() << " at p=" << period << "\n";
            }
            fperr.theta = theta;
            throw fperr;
        }
        
        // IF we are at lmin, have not resolved rotation, AND NO separation is detected
        // THEN, this is a special UNRESOLVED case where we are crossing a stable invariant manifold
        // near the fixed point (magnitude of displacement is vanishing)
        // -> Just accept all points where the midpoint works as stable manifold crossings
        bool midpointWorked = true;
        vec_type xMid = 0.5*(x0+x1);
        vec_type deltaMid;
        try {
            deltaMid = mapDisplacementUsingCache(xMid, pmap, period, params, cache, eFail);
        } catch(ERotFailure& ef) {
            deltaMid = vec_type(50.0);
            midpointWorked = false;
        }
        if(midpointWorked) {
            unsigned int lowMagCount = 0;
            double deltaMidNorm = nvis::norm(deltaMid);
            if (deltaMidNorm < 0.1) {
                lowMagCount++;
            }
            double delta0Norm = nvis::norm(delta0);
            if (delta0Norm < 0.1) {
                lowMagCount++;
            }
            double delta1Norm = nvis::norm(delta1);
            if (delta1Norm < 0.1) {
                lowMagCount++;
            }
            //Only a stable manifold crossing if magnitude of delta is dropping within
            if(lowMagCount >= 2) {
                //Likely a stable manifold crossing
                ERotFailure smcFail;
                smcFail.setType(ERotFailure::STABLE_MANIFOLD_CROSS);
                smcFail.failurePos = xMid;
                smcFail.period = period;
                double thetaA = xavier::signed_angle(delta0,deltaMid);
                double thetaB = xavier::signed_angle(deltaMid,delta1);
                smcFail.theta = thetaA + thetaB;
                bool useAsGuess = (deltaMidNorm <= 10.0*params.min_fp_tol);
                if (params.debug) {
                    std::cout << "Period (" << period << "): However, a stable-manifold crossing is suspected\n";
                    std::cout << "  StableMan Cross: Partial rotations : " << thetaA << " + "
                              << thetaB << " = " << thetaA+thetaB << "\n";
                    if(useAsGuess) {
                        std::cout << "  StableMan Cross is close enough (" << deltaMidNorm << ") for guess\n";
                    }
                }
                //Throw the error to store
                //throw smcFail;
                smcFailVec.push_back(smcFail);
                //Instead of throwing the error, just accept rotation value as smooth
                return smcFail.theta;
            }
        }
        
        //  We know that there is a problem, but we were unable to fix it or know what it is...
        ERotFailure uf; //Unknown Error by default
        uf.period = period;
        uf.failurePos = 0.5*(x0+x1);
        if (verbose) {
            std::cerr << "Period (" << period << "): Unable to track rotation on edge, reason is 'UNKNOWN'\n";
        }
        throw uf;
    }
    
    //Without failure or not at minimum distance, we compute a subdivision and it's rotation
#ifndef SECANT_METHOD
    vec_type x = 0.5 * (x0 + x1);
#else
    double u = secant<vec_type>(v0, v1);
    if (std::isnan(u) || std::isinf(u)) {
        std::cerr << "secant method returned NaN" << std::endl;
    }
    vec_type x = (1-u)*x0 + u*x1;
#endif
    
    //Displacement vectors
    vec_type delta;
    ///Forward and Backward Map Calls
    ERotFailure pf;
    //Forward map
    delta = mapDisplacementUsingCache(x,pmap,period,params,cache,pf);
    //Backward map
    //vec_type beta;
    //beta = mapDisplacementUsingCache(x,pmap,-period,params,backCache,pf);
    
    //Map Tangent vector:  Like a low-pass filter on the map, makes displacement field more like canonical saddle
    //vec_type eta = (delta - beta) / 2.0;
    /*if (verbose) {
        std::cerr << x0[0] << ", " << x0[1] << ", " << v0[0] << ", " << v0[1] << '\n';
        std::cerr << x1[0] << ", " << x1[1] << ", " << v1[0] << ", " << v1[1] << '\n';
    }*/
    
    return adaptiveRotationAngle(x0, v0, x, delta, pmap, period, params, cache, backCache, smcFailVec) +
           adaptiveRotationAngle(x, delta, x1, v1, pmap, period, params, cache, backCache, smcFailVec);
    //return adaptiveRotationAngle(x0, v0, x, eta, pmap, period, params, cache, backCache) +
    //       adaptiveRotationAngle(x, eta, x1, v1, pmap, period, params, cache, backCache);
}


/// Function to recursively subdivide along an edge for use in
/// Non-Transverse version of adaptiveRotationAngle()
template<class MAP,class PARAM,class SRDATA,class EROTFAILURE>
void ntAdaptiveRotation(
    const vec_type& x0, const vec_type& v0, const vec_type& x1, const vec_type& v1,
    const MAP& pmap, const int period, const PARAM& params, const double edgeLength,
    std::set<SRDATA>& cache, std::set<SRDATA>& backCache,
    std::set<SRDATA>& dispSet,
    std::list<EROTFAILURE>& sepFailureList
)
{
    typedef typename MAP::lvec_type vec_type;
    typedef typename MAP::lmat_type mat_type;
    typedef typename MAP::xstate_type xstate_type;
    typedef std::pair<vec_type,mat_type>  ProjPair;
    //typedef EdgeRotationFailure<vec_type> ERotFailure;
    typedef EROTFAILURE ERotFailure;
    //Adaptive Edge data structures
    typedef AdaptiveEdge<SRDATA,vec_type>      EdgeType;
    typedef typename EdgeType::IdType          IdType;
    typedef std::pair<IdType,IdType>           Segment;
    
    //static const int numSing = MAP::rhs_type::numSingularities;
    
    //Ideally, we would like to continue to subdivide segments to max dx (=lmin)
    //for every section separation location, but this is time-prohibitive.
    
    //Evaluate the Angle displacement
    //bool verbose = params.verbose;
    bool debug = params.debug;
    //const int& ntEdgeDivisions = params.ntEdgeDivisions;
    //const double lmin = edgeLength / ((double)ntEdgeDivisions-1.0); //Subdivide based on current edge
    const double lmin = params.lmin; //Fixed distance
    const double lmax = edgeLength / 32.0; //Force at least 5 subdivisions
    //In theory, we don't need as many subdivisions knowing we have Discontinuities in Line integral
    
    //Adaptive Edge Data structure
    EdgeType adaptiveEdge(x0,x1);
    //std::cout << __FILE__ << ": NTRotation Testing Edge [" << x0 << " , " << x1 << "]\n";
    //List of segments to try
    std::vector<Segment>  segsToCheck;
    segsToCheck.push_back( Segment(IdType(0,0),IdType(1,0)) );
    
    //Compute the forward and backward mappings for the endpoints (they should exist)
    ERotFailure eFail;
    vec_type delta0, delta1;
    //vec_type eta0 = v0;
    //vec_type eta1 = v1;
    try {
        delta0 = mapDisplacementUsingCache(x0, pmap, period, params, cache, eFail);
    } catch(ERotFailure& ef) {
        ef.edgeID = IdType(0,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(0,0),ef);
        delta0 = vec_type(50);
    }
    try {
        delta1 = mapDisplacementUsingCache(x1, pmap, period, params, cache, eFail);
    } catch(ERotFailure& ef) {
        ef.edgeID = IdType(1,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(1,0),ef);;
        delta1 = vec_type(50);
    }
    /*vec_type beta0, beta1;
    try {
        beta0  = mapDisplacementUsingCache(x0, pmap, -period, params, backCache, eFail);
    } catch(ERotFailure& ef) {
        ef.edgeID = IdType(0,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(0,0),ef);;
        beta0 = vec_type(50);
    }
    try {
        beta1  = mapDisplacementUsingCache(x1, pmap, -period, params, backCache, eFail);
    } catch(ERotFailure& ef) {
        ef.edgeID = IdType(1,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(1,0),ef);;
        beta1 = vec_type(50);
    }*/
    
    //Should check initial nodes for fixed points
    if( isFixedPointSuspected(x0,delta0,period,params,eFail) ) {
        //if( isFixedPointSuspectedTangent(x,delta,beta,eta,period,params,eFail) ) {
        eFail.edgeID = IdType(0,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(0,0),eFail);
        if(debug) {
            std::cout << "   Edge_Rot: Fixed Point Detected = " << eFail.where() << " at p=" << period << "\n";
        }
    }
    if( isFixedPointSuspected(x1,delta1,period,params,eFail) ) {
        //if( isFixedPointSuspectedTangent(x,delta,beta,eta,period,params,eFail) ) {
        eFail.edgeID = IdType(1,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(1,0),eFail);
        if(debug) {
            std::cout << "   Edge_Rot: Fixed Point Detected = " << eFail.where() << " at p=" << period << "\n";
        }
    }
    
    //Add the displacement information into AdaptiveEdge
    std::vector<vec_type> dispValues;
    dispValues.push_back( delta0 );
    //dispValues.push_back( beta0 );
    //dispValues.push_back( eta0 );
    adaptiveEdge.setValue( IdType(0,0), SRDATA(x0,dispValues) );
    dispValues.clear();
    dispValues.push_back( delta1 );
    //dispValues.push_back( beta1 );
    //dispValues.push_back( eta1 );
    adaptiveEdge.setValue( IdType(1,0), SRDATA(x1,dispValues) );
    
    
    //Edge subdivision algorithm
    int numSegsToCheck = 1;
    while ( numSegsToCheck > 0 ) {
        std::vector<Segment> newSegs;
        typename std::vector<Segment>::iterator segIT;
        
        //Check each segment for subdivision    ----------------------------
        for(segIT=segsToCheck.begin(); segIT!=segsToCheck.end(); ++segIT) {
            //std::cout << "Testing segment : " << segIT->first << " , " << segIT->second << "\n";
            bool subdivide = false;
            //if segment length is bigger than max segment length, subdivide
            vec_type xA = adaptiveEdge.getVertex(segIT->first);
            vec_type xB = adaptiveEdge.getVertex(segIT->second);
            double length = nvis::norm(xB-xA);
            if(length > lmax) {
                subdivide = true;
            }
            //Check for a discontinuity point, subdivide if still above lmin
            if(!subdivide && (length > lmin)) {
                //Are Node 0 and Node 1 both Singularities? stop
                if(adaptiveEdge.isSingularityPoint(segIT->first) &&
                        adaptiveEdge.isSingularityPoint(segIT->second) ) {
                    continue;
                }
                
                //Is Node 0 or Node 1 a discontinuity point? subdivide?
                if(adaptiveEdge.isDiscontinuity(segIT->first) ||
                        adaptiveEdge.isDiscontinuity(segIT->second)) {
                    //subdivide the node
                    subdivide = true;
                }
            }
            //Checking tracking conditions
            if(!subdivide) {
                //Get relevant data
                SRDATA dataA = adaptiveEdge.getValue(segIT->first);
                SRDATA dataB = adaptiveEdge.getValue(segIT->second);
                vec_type deltaA = dataA.returns[0], deltaB = dataB.returns[0];
                //vec_type betaA = dataA.returns[1],  betaB = dataB.returns[1];
                //vec_type etaA = dataA.returns[2],   etaB = dataB.returns[2];
                
                //Edge Rotation tracking conditions met?
                bool rotConditionOK = getEdgeRotationTrackingCondition(
                                          xA,deltaA,xB,deltaB,pmap,period,params,cache);
                //MapTangent Vector Version
                //bool rotConditionOK = getEdgeRotationTrackingConditionTangent(
                //                            xA,etaA,xB,etaB,pmap,period,params,cache,backCache);
                
                //Decide what to do based on tracking conditions:
                if (rotConditionOK) {
                    continue;    //Done with segment
                } else if ( length <= lmin ) {
                    //std::cout << "  Segment at lmin: " << length << " < lmin=" << lmin << "\n";
                    ERotFailure rotErr; //Empty discontinuity for checking errors
                    //First, check if nodes are discontinuities (any type),
                    //  if yes, no further errors to report
                    if ( adaptiveEdge.isDiscontinuity(segIT->first) ||
                            adaptiveEdge.isDiscontinuity(segIT->second) ) {
                        continue; //Done
                    }
                    //Check for fixed point conditions
                    else if(
                        isFixedPointSuspected(xA,deltaA,xB,deltaB,period,params,rotErr)
                        //isFixedPointSuspectedTangent(xA,deltaA,betaA,etaA,xB,deltaB,betaB,etaB,period,params,rotErr)
                    ) {
                        if(debug) {
                            std::cout << "   Edge_Rot: Fixed Point Detected = " << eFail.where() << " at p=" << period << "\n";
                        }
                        rotErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                        adaptiveEdge.insertDiscontinuity(rotErr);
                        continue; //Done
                    }
                    /// If at lmin without successful rotation conditions nor a fixed point suspected,
                    /// check if it's a section-separation condition (test both forward and backward):
                    bool discontinuityFound = false;
                    /*if (detectSectionSeparation(xA,xB,backCache,pmap,-period,params) ) {
                      // If we are at lmin, store the backward separation point
                      vec_type intersection = xA+xB;
                      intersection /= 2.0;
                      ERotFailure sepErr(ERotFailure::BACKWARD_SECTION_SEP,intersection,-period);
                      if (detectSingularityIntersection(xA,xB,backCache,pmap,-period,params) ) {
                        sepErr.setType(ERotFailure::BACKWARD_SINGULARITY);
                        //Data output for debugging
                        typename std::set<SRDATA>::iterator cit0,cit1;
                        cit0 = cache.find( SRDATA(x0) ); //Note:  Should be found
                        cit1 = cache.find( SRDATA(x1) ); //Note:  Should be found
                        vec_type p0 = (*cit0).returns[std::abs(period)-1]; //P^p(x0)
                        vec_type p1 = (*cit1).returns[std::abs(period)-1]; //P^p(x1)
                        double dt0 = (*cit0).data[std::abs(period)-1][numSing];
                        double dt1 = (*cit1).data[std::abs(period)-1][numSing];
                        const xstate_type* singPtr = pmap.rhs().singularities();
                        std::cout << " Backward Singularity detected :\n";
                        std::cout << "  period = " << -period << "\n";
                        std::cout << "  x0 = " << xA << " x1 = " << xB << "\n";
                        std::cout << "  p0 = " << p0 << " p1 = " << p1 << "\n";
                        std::cout << "  dt0 = " << dt0 << " dt1 = " << dt1 << "\n";
                        for(int k=0;k<numSing;k++) {
                          std::cout << "  Singularity " << k << ":\n";
                          double singularitySafeDistance = pmap.rhs().getSingularitySafeDistance(k);
                          double minCloseApproachDist0 = (*cit0).data[std::abs(period)-1][k];
                          double minCloseApproachDist1 = (*cit1).data[std::abs(period)-1][k];
                          //Check for sign change in position coordinate from singularity
                          ProjPair singPair = pmap.section().project( singPtr[k] );
                          vec_type s = singPair.first;
                          vec_type r0 = p0 - s;
                          vec_type r1 = p1 - s;
                          std::cout << "    r0 = " << r0 << "  r1 = " << r1 << "\n";
                          std::cout << "    (r0[0]*r1[0])<0) = " << ((r0[0]*r1[0])<0) << "\n";
                          std::cout << "    (p0[1]*p1[1])<0) = " << ((p0[1]*p1[1])<0) << "\n";
                          std::cout << "    SafeDist = " << singularitySafeDistance
                                    << " MinCA0 = " << minCloseApproachDist0
                                    << " MinCA1 = " << minCloseApproachDist1 << "\n";
                        }
                      }
                      adaptiveEdge.insertDiscontinuity( sepErr );
                      discontinuityFound = true;
                      std::cout << " Discontinuity Found: " << sepErr.what() << "\n";
                    }*/
                    //Forward
                    if (detectSectionSeparation(xA,xB,cache,pmap,period,params) ) {
                        // If we are at lmin, store the forward separation point
                        vec_type intersection = xA+xB;
                        intersection /= 2.0;
                        ERotFailure sepErr(ERotFailure::FORWARD_SECTION_SEP,intersection,period);
                        sepErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1); //For sorting
                        if (detectSingularityIntersection(xA,xB,cache,pmap,period,params) ) {
                            sepErr.setType(ERotFailure::FORWARD_SINGULARITY);
                            //Data output for debugging
                            /*typename std::set<SRDATA>::iterator cit0,cit1;
                            cit0 = cache.find( SRDATA(x0) ); //Note:  Should be found
                            cit1 = cache.find( SRDATA(x1) ); //Note:  Should be found
                            vec_type p0 = (*cit0).returns[std::abs(period)-1]; //P^p(x0)
                            vec_type p1 = (*cit1).returns[std::abs(period)-1]; //P^p(x1)
                            double dt0 = (*cit0).data[std::abs(period)-1][numSing];
                            double dt1 = (*cit1).data[std::abs(period)-1][numSing];
                            const xstate_type* singPtr = pmap.rhs().singularities();
                            std::cout << " Forward Singularity detected :\n";
                            std::cout << "  period = " << period << "\n";
                            std::cout << "  x0 = " << xA << " x1 = " << xB << "\n";
                            std::cout << "  p0 = " << p0 << " p1 = " << p1 << "\n";
                            std::cout << "  dt0 = " << dt0 << " dt1 = " << dt1 << "\n";
                            for(int k=0;k<numSing;k++) {
                              std::cout << "  Singularity " << k << ":\n";
                              double singularitySafeDistance = pmap.rhs().getSingularitySafeDistance(k);
                              double minCloseApproachDist0 = (*cit0).data[std::abs(period)-1][k];
                              double minCloseApproachDist1 = (*cit1).data[std::abs(period)-1][k];
                              //Check for sign change in position coordinate from singularity
                              ProjPair singPair = pmap.section().project( singPtr[k] );
                              vec_type s = singPair.first;
                              vec_type r0 = p0 - s;
                              vec_type r1 = p1 - s;
                              std::cout << "    r0 = " << r0 << "  r1 = " << r1 << "\n";
                              std::cout << "    (r0[0]*r1[0])<0) = " << ((r0[0]*r1[0])<0) << "\n";
                              std::cout << "    (p0[1]*p1[1])<0) = " << ((p0[1]*p1[1])<0) << "\n";
                              std::cout << "    SafeDist = " << singularitySafeDistance
                                        << " MinCA0 = " << minCloseApproachDist0
                                        << " MinCA1 = " << minCloseApproachDist1 << "\n";
                            }*/
                        }
                        adaptiveEdge.insertDiscontinuity( sepErr );
                        discontinuityFound = true;
                        //std::cout << " Discontinuity Found: " << sepErr.what() << "\n";
                    }
                    //Check for Stable manifold crossings - usually only close to fixed points
                    if (!discontinuityFound) {
                        //Try one last level of subdivision to see if we crossed a stable manifold
                        vec_type xMid = 0.5*(xA+xB);
                        vec_type deltaMid(50.0);
                        bool midpointWorked = true;
                        ERotFailure tempFail;
                        try {
                            deltaMid = mapDisplacementUsingCache(xMid,pmap,period,params,cache,tempFail);
                        } catch(ERotFailure& ef) {
                            deltaMid = vec_type(50.0);
                            midpointWorked=false;
                        }
                        if(midpointWorked) {
                            unsigned int lowMagCount = 0;
                            double deltaMidNorm = nvis::norm(deltaMid);
                            if (deltaMidNorm < 0.1) {
                                lowMagCount++;
                            }
                            double delta0Norm = nvis::norm(deltaA);
                            if (delta0Norm < 0.1) {
                                lowMagCount++;
                            }
                            double delta1Norm = nvis::norm(deltaB);
                            if (delta1Norm < 0.1) {
                                lowMagCount++;
                            }
                            if(lowMagCount >= 2) {
                                //Likely stable manifold crossing
                                tempFail.setType(ERotFailure::STABLE_MANIFOLD_CROSS);
                                tempFail.failurePos = xMid;
                                tempFail.period = period;
                                double thetaA = xavier::signed_angle(deltaA,deltaMid);
                                double thetaB = xavier::signed_angle(deltaMid,deltaB);
                                tempFail.theta = thetaA + thetaB;
                                adaptiveEdge.insertDiscontinuity( tempFail );
                                discontinuityFound = true;
                            }
                        }
                        //If midpoint failed, move on...
                    }
                    
                    //At this point, we don't know what went wrong, so should throw 'unknown'
                    // HOWEVER, it is likely that these are still separation conditions (likely singularity
                    // intersections) that are not detected by the heuristics because the subdivision isn't
                    // refined enough.  Thus, lets pass these as UNRESOLVED cases
                    if (!discontinuityFound) {
                        vec_type intersection = xA+xB;
                        intersection /= 2.0;
                        ERotFailure singErr(ERotFailure::UNRESOLVED,intersection,period);
                        singErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1); //For sorting
                        /*bool isDeltaSmooth = runAngleTest(deltaA,deltaB,pmap,params);
                        if (isDeltaSmooth) {
                          //singErr.setType(ERotFailure::BACKWARD_SINGULARITY);
                          singErr.period = -period;
                        }*/
                        adaptiveEdge.insertDiscontinuity( singErr );
                        if(debug) {
                            std::cout << " Unknown Discontinuity Found: " << singErr.what() << "\n";
                        }
                        
                    }
                    continue; //Force a completion of this segment if at lmin
                } else {
                    subdivide = true;    //rotation not ok, but not at lmin so subdivide
                }
            }
            if(subdivide) {
                //Call the refine command on the left/bottom node
                adaptiveEdge.refine(segIT->first);
                //Set the new vertex (could be midpoint or something else like Secant method)
                vec_type xMidPoint = xA+xB;
                xMidPoint /= 2.0;
                IdType newID = adaptiveEdge.getChildID(segIT->first,1);
                adaptiveEdge.setVertex(newID, xMidPoint);
                //std::cout << " Subdivision values: \n";
                //std::cout << "   New Point : " << newID << " -> " << xMidPoint << "\n";
            }
        } //End Segment check loop ---------------------------------------------------------------
        //Add only new segments to list for next step
        adaptiveEdge.getNewSegments(newSegs);
        
        //Run map and set new information (mapDisplacement,ERotFailure)
        std::vector<IdType> seedIDs;
        adaptiveEdge.getEmptyLeaves(seedIDs);
        typename std::vector<IdType>::iterator seedIT;
        for(seedIT=seedIDs.begin(); seedIT!=seedIDs.end(); ++seedIT) {
            //std::cout << " Propagation of New Nodes:\n";
            //std::cout << "   SeedID = " << (*seedIT) << "\n";
            //Get position node from structure
            vec_type x = adaptiveEdge.getVertex((*seedIT));
            //std::cout << "   Seed x = " << x << "\n";
            vec_type delta,beta,eta;
            bool detectedSeparation = false;
            
            //Run forward map
            try {
                delta = mapDisplacementUsingCache(x,pmap,period,params,cache,eFail);
            } catch(ERotFailure& ef) {
                ef.edgeID = (*seedIT);
                adaptiveEdge.insertDiscontinuityNode((*seedIT),ef);
                delta = vec_type(50);
                detectedSeparation = true;
                if(debug) {
                    std::cout << "   Delta(x) Encountered Error = " << ef.what() << "\n";
                }
            }
            
            //Run backward map
            /*try {
                beta = mapDisplacementUsingCache(x,pmap,-period,params,backCache,eFail);
            } catch(ERotFailure& ef) {
                ef.edgeID = (*seedIT);
                adaptiveEdge.insertDiscontinuityNode((*seedIT),ef);
                beta = vec_type(50);
                detectedSeparation = true;
                std::cout << "   Beta(x) Encountered Error = " << ef.what() << "\n";
            }
            eta = (delta-beta)/2.0;*/
            
            //Check for a fixed point at this new point
            if( !detectedSeparation && isFixedPointSuspected(x,delta,period,params,eFail) ) {
                //if( !detectedSeparation && isFixedPointSuspectedTangent(x,delta,beta,eta,period,params,eFail) ) {
                eFail.edgeID = (*seedIT);
                adaptiveEdge.insertDiscontinuityNode((*seedIT),eFail);
                detectedSeparation = true;
                if(debug) {
                    std::cout << "   Edge_Rot: Fixed Point Detected = " << eFail.where() << " at p=" << period << "\n";
                }
            }
            //Run test for transversality condition (dydot(P^p(x)) != 0)
            else if( !detectedSeparation ) {
                //Forward map
                if(detectSectionSeparation(x,cache,pmap,period,params)) {
                    ERotFailure nodeTransFail(ERotFailure::FORWARD_SECTION_SEP_NODE,x,period);
                    nodeTransFail.edgeID = (*seedIT);
                    adaptiveEdge.insertDiscontinuityNode((*seedIT),nodeTransFail);
                    detectedSeparation = true;
                    if(debug) {
                        std::cout << "   Forward Separation Node = " << nodeTransFail.where() << "\n";
                    }
                }
                /*//Backward map
                if(detectSectionSeparation(x,backCache,pmap,-period,params)) {
                  ERotFailure nodeTransFail(ERotFailure::BACKWARD_SECTION_SEP_NODE,x,-period);
                  adaptiveEdge.insertDiscontinuityNode((*seedIT),nodeTransFail);
                  detectedSeparation = true;
                  std::cout << "   Backward Separation Node = " << nodeTransFail.where() << "\n";
                }*/
            }
            //Generate map displacement data
            dispValues.clear();
            dispValues.push_back( delta );
            //dispValues.push_back( beta );
            //dispValues.push_back( eta );
            //Insert into adaptiveEdge structure
            adaptiveEdge.setValue( (*seedIT), SRDATA(x,dispValues) );
        } //End mapping loop for new points
        
        //Fill segsToCheck with new information
        segsToCheck.clear();
        segsToCheck = newSegs;
        numSegsToCheck = (int) segsToCheck.size();
        //std::cout << " After Propagation: " << numSegsToCheck << " Segments to test on next loop\n";
        
    } //End subdivision algorithm
    
    /// Fill output data:
    adaptiveEdge.getData(dispSet);
    adaptiveEdge.getDiscontinuityList(sepFailureList);
    //Uniquify sep failures?
}

/** Compute rotation angle between two points (x0,x1) on an edge
   * This is specific to a Section Transversality violation
   * Stores all steps, restart angle summation at each section separation point
     through the tracking process so this segments the summation into a
     piecewise computation.
   * Throws an EdgeRotationFailure when encountering an unfixable mapping error
     or a suspected fixed-point
   * Also requires the edge-set information so it can store all section separation locations
   -Note:  Using v as the "tangent vector" which is an approximation
    of the map tangent if map curves were continuous.  So
      v = (P^p(x) - P^-p(x))/2;
      //Kinda like low-pass filter or Central-Difference Approx
 */
template<class EDGE,class MAP,class PARAM, class SRDATA, class FP>
double adaptiveRotationAngleNonTransverse( EDGE& e, const double& edgeLength,
        const std::vector<typename MAP::lvec_type>& xEdge, const std::vector<typename MAP::lvec_type>& vEdge,
        const MAP& pmap, const int period, const PARAM& params, std::vector< FP >& fpGuesses,
        std::set<SRDATA>& cache, std::set<SRDATA>& backCache )
{
    typedef typename MAP::lvec_type         vec_type;
    typedef EdgeRotationFailure<vec_type>   ERotFailure;
    const int& pmax = params.max_period;
    
    //Section separations locations
    std::list<ERotFailure> sepPtsList; //isValid = false means backward time
    std::set<SRDATA> edgeMapDisplacements; //Using the sorting ability
    
    //Run adaptive rotation algorithm with AdaptiveEdge data structure
    ntAdaptiveRotation(xEdge[0],vEdge[0],xEdge.back(),vEdge.back(),
                       pmap, period, params, edgeLength, cache, backCache,
                       edgeMapDisplacements, sepPtsList  );
                       
    //Sum through points, restarting at each section sep or singularity location (backward OR forward):
    double theta = 0.0, tmp = 0.0;
    sepPtsList.sort(); //Sort to work from left to right or bottom to top along edge
    
    //Find stop points for independent summation segments
    typename std::list<ERotFailure>::iterator sepListIT = sepPtsList.begin();
    for (; sepListIT!=sepPtsList.end(); sepListIT++) {
        ERotFailure& sepFail = (*sepListIT);
        sepFail.fixed = true;
        e.rotFailures.push_back( sepFail );
        //If this is a fixed point, add to guesses
        if ( (sepFail.type == ERotFailure::FIXED_POINT_SUSPECTED ) ||
                (sepFail.type == ERotFailure::DOUBLE_PERIOD_FIXED_POINT && std::abs(sepFail.period)<=pmax) ) {
            FP aGuess;
            aGuess.pos = sepFail.failurePos;
            aGuess.K = std::abs(sepFail.period);
            fpGuesses.push_back(aGuess);
            
        }
        //Insert blank mapDisplacements into data set for sep points
        std::vector<vec_type> dispValues;
        dispValues.push_back( vec_type(50.0) );
        //dispValues.push_back( vec_type(50.0) );
        //dispValues.push_back( vec_type(50.0) );
        edgeMapDisplacements.insert( SRDATA(sepFail.failurePos,dispValues) );
        //Inserting into the set should NOT overwrite nodes with discontinuties
    }
    
    //Make the list unique
    sepPtsList.unique(); //Checking if periods and positions match
    
    ///Sum through the map displacements along the edge, compute rotation and add to running total
    if(params.debug) {
        std::cout << " Summing rotations along Edge: [" << xEdge[0] << "," << xEdge.back() << "]\n";
        std::cout << "**********************************************************************\n";
    }
    typename std::set<SRDATA>::iterator dataIT = edgeMapDisplacements.begin();
    sepListIT = sepPtsList.begin();
    typename std::set<SRDATA>::iterator lastIT = dataIT;
    int seg = 0;
    for(++dataIT; dataIT!=edgeMapDisplacements.end(); ++lastIT,++dataIT) {
        //Evaluate rotation in the forward map displacement vectors
        vec_type delta0 = lastIT->returns[0];
        vec_type delta1 = dataIT->returns[0];
        //Compute the signed rotation from last step to next step
        tmp = xavier::signed_angle(delta0,delta1);
        if(params.debug) {
            std::cout << " Seg " << seg << " : x0 = " << lastIT->x0 << " x1 = " << dataIT->x0 << "\n";
            std::cout << "   delta0 = " << delta0 << "  delta1 = " << delta1 << "\n";
            std::cout << "   Current theta = " << theta << "\n";
            std::cout << "   Signed Angle = " << tmp << "\n";
            std::cout << "   Current Discontinuity = " << sepListIT->failurePos << " " << sepListIT->what() << "\n";
        }
        
        /*//Evaluate rotation in the map displacement tangent vectors
        vec_type eta0 = lastIT->returns[2];
        vec_type eta1 = dataIT->returns[2];
        //Compute the signed rotation from last step to next step
        tmp = xavier::signed_angle(eta0,eta1);*/
        
        //Check if point is a separation point
        bool currentIsSepPt = nvis::all(dataIT->x0 == sepListIT->failurePos);
        bool lastIsSepPt = nvis::all(lastIT->x0 == sepListIT->failurePos);
        if(params.debug) {
            std::cout << "   IsLastPtSep = " << lastIsSepPt
                      << "  IsCurrentPtSep = " << currentIsSepPt << "\n";
        }
        if ( !currentIsSepPt && !lastIsSepPt ) {
            //Add rotation to sum if not accross separation conditions
            theta += tmp;
        }
        //If type is stable manifold crossing, assume the rotation is ok
        else if (!lastIsSepPt && currentIsSepPt
                 && (sepListIT->type == ERotFailure::STABLE_MANIFOLD_CROSS) ) {
            //count the rotation stored
            theta += sepListIT->theta; //Stored as the rotation between nodes
        } else if (lastIsSepPt) {
            //If the trailing point is now on the sep Point, switch to next sep point
            //Increment to next separation point
            sepListIT++;
            if (sepListIT == sepPtsList.end()) {
                sepListIT = sepPtsList.begin();
            }
        }
        if(params.debug) {
            std::cout << "   Total theta = " << theta << "\n";
            std::cout << "-------------------------------------------------------------------\n";
        }
        seg++;
    }
    if(params.debug) {
        std::cout << " Total theta = " << theta << "\n";
        std::cout << "**********************************************************************\n";
    }
    
    //If everything worked, return the rotation value
    return theta;
}


/* -- THIS IS THE RECURSIVE-CALL VERSION OF THIS FUNCTION --
/// Function to recursively subdivide along an edge for use in
/// Non-Transverse version of adaptiveRotationAngle()
template<class MAP,class PARAM,class SRDATA,class EROTFAILURE>
double ntAdaptiveRotation(
        const vec_type& x0, const vec_type& v0, const vec_type& x1, const vec_type& v1,
        const MAP& pmap, const int period, const PARAM& params, const double edgeLength,
        std::set<SRDATA>& cache, std::set<SRDATA>& backCache,
        std::set<SRDATA>& dispSet,
        std::list<EROTFAILURE>& sepFailureList
        )
{
    typedef typename MAP::lvec_type vec_type;
    //typedef EdgeRotationFailure<vec_type> ERotFailure;
    typedef EROTFAILURE ERotFailure;

    //Ideally, we would like to continue to subdivide segments to max dx (=lmin)
    //for every section separation location, but this is time-prohibitive.

    //Evaluate the Angle displacement
    double theta = xavier::signed_angle(v0, v1);  //Computing theta here means nothing
    bool verbose = params.verbose;
    const int& ntEdgeDivisions = params.ntEdgeDivisions;
    const double lmin = edgeLength / ((double)ntEdgeDivisions-1.0); //Subdivide based on current edge
    //In theory, we don't need as many subdivisions knowing we have Transversality Violations


    //Compute the forward and backward mappings for the endpoints (they should exist)
    vec_type delta0 = mapDisplacementUsingCache(x0, pmap, period, params, cache);
    vec_type beta0  = mapDisplacementUsingCache(x0, pmap, -period, params, backCache);
    vec_type delta1 = mapDisplacementUsingCache(x1, pmap, period, params, cache);
    vec_type beta1  = mapDisplacementUsingCache(x1, pmap, -period, params, backCache);

    //Check if we have fixed point conditions
    ERotFailure fpTestFail;
    if ( isFixedPointSuspected(x0,delta0,beta0,v0,x1,delta1,beta1,v1,period,params,fpTestFail) ) {
      //throw fpTestFail;
      sepFailureList.push_back(fpTestFail);
    }

    //Check for successful tracking condition
    bool rotConditionOK = getEdgeRotationTrackingCondition(x0,v0,x1,v1,pmap,period,params,cache,backCache);
    if  (rotConditionOK ) {return theta;} //  Tracking Conditions are met
    //If not successful and at smallest allowed step, throw corresponding error
    else if (nvis::norm(x1-x0) <= lmin ) {
        //Check if we have fixed point conditions first
        ERotFailure fperr;
        if( isFixedPointSuspected(x0,delta0,beta0,v0,x1,delta1,beta1,v1,period,params,fperr) ) {
            if (verbose) {
                std::cerr << "Period (" << period << "): However, a fixed_point is suspected on the edge\n";
            }
            //throw fperr;
            sepFailureList.push_back(fpTestFail);
        }

        /// If not successful rotation nor a fixed point suspected,
        /// check if it's a section separation condition (have to test both forward and backward)
        //Check for singularity conditions
        if (detectSingularityIntersection(x0,x1,backCache,pmap,-period,params) ) {
          // If we are at lmin, store the backward singular point
          vec_type intersection = x0+x1;
          intersection /= 2.0;
          ERotFailure singErr(ERotFailure::BACKWARD_SINGULARITY,intersection,-period);
          sepFailureList.push_back( singErr );
          return 0.0; //Force a completion
        }
        else if (detectSingularityIntersection(x0,x1,cache,pmap,period,params)) {
          // If we are at lmin, store the forward singular point
          vec_type intersection = x0+x1;
          intersection /= 2.0;
          ERotFailure singErr(ERotFailure::FORWARD_SINGULARITY,intersection,period);
          sepFailureList.push_back( singErr );
          return 0.0; //Force a completion
        }
        else if (detectSectionSeparation(x0,x1,backCache,pmap,-period,params) ) {
          // If we are at lmin, store the backward separation point
          vec_type intersection = x0+x1;
          intersection /= 2.0;
          ERotFailure sepErr(ERotFailure::BACKWARD_SECTION_SEP,intersection,-period);
          sepFailureList.push_back( sepErr );
          return 0.0; //Force a completion
        }
        else if (detectSectionSeparation(x0,x1,cache,pmap,period,params) ) {
          // If we are at lmin, store the forward separation point
          vec_type intersection = x0+x1;
          intersection /= 2.0;
          ERotFailure sepErr(ERotFailure::FORWARD_SECTION_SEP,intersection,period);
          sepFailureList.push_back( sepErr );
          return 0.0; //Force a completion
        }
        //At this point, we don't know what went wrong, so should throw 'unknown'
        // HOWEVER, it is likely that these are still separation conditions (likely singularity
        // intersections) that are not detected by the heuristics because the subdivision isn't
        // refined enough.  Thus, lets pass these as singularity conditions
        else {
          vec_type intersection = x0+x1;
          intersection /= 2.0;
          bool isDeltaSmooth = runAngleTest(delta0,delta1,pmap,params);
          ERotFailure singErr(ERotFailure::FORWARD_SINGULARITY,intersection,period);
          if (isDeltaSmooth) {
            singErr.setType(ERotFailure::BACKWARD_SINGULARITY);
            singErr.period = -period;
          }
          sepFailureList.push_back( singErr );
          return 0.0; //Force a completion
        }

    }

    //Otherwise, we run the mapping, storing only separation errors, but throwing error if
#ifndef SECANT_METHOD
    vec_type x = 0.5 * (x0 + x1);
#else
    double u = secant<vec_type>(v0, v1);
    if (std::isnan(u) || std::isinf(u)) {
        std::cerr << "secant method returned NaN" << std::endl;
    }
    vec_type x = (1-u)*x0 + u*x1;
#endif

    //Displacement vectors
    vec_type eta,delta,beta;
    ///Forward and Backward Map Calls
    ERotFailure pf;
    //Forward map without section sep failures
    try {
      delta = mapDisplacementUsingCache(x,pmap,period,params,cache,pf);
    } catch(ERotFailure& ef) {
      sepFailureList.push_back(ef);
      return 0.0; //Force a completion of subdivision
    }
    //Backward map without section sep failures
    try {
      beta  = mapDisplacementUsingCache(x,pmap,-period,params,backCache,pf);
    } catch(ERotFailure& ef) {
      sepFailureList.push_back(ef);
      return 0.0; //Force a completion of subdivision
    }
    //Map Tangent vector:  Like a low-pass filter on the map, makes displacement field more like canonical saddle
    eta = (delta - beta) / 2.0;


    //Add displacement to cache
    std::vector<vec_type> dispValues;
    dispValues.push_back( delta );
    dispValues.push_back( beta );
    dispValues.push_back( eta );
    dispSet.insert( SRDATA(x,dispValues) );

    //Adaptively subdivide through recursive call:  returned value means nothing
    return ntAdaptiveRotation(x0,v0,x,eta,pmap,period,params,edgeLength,cache,backCache,dispSet,sepFailureList) +
           ntAdaptiveRotation(x,eta,x1,v1,pmap,period,params,edgeLength,cache,backCache,dispSet,sepFailureList);
}*/

} // end topology


#endif