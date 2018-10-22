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


//Section Transversality Detection :
//  Check if there is a section separation between two close initial points
//Author : Wayne Schlei (Purdue University)
// Note:  This may have to be tailored to a specific problem...

#ifndef __SECTION_TRANSVERSALITY_HPP__
#define __SECTION_TRANSVERSALITY_HPP__

#include <iostream>
#include <iterator>
#include <set>
#include <algorithm>
#include <math/fixed_vector.hpp>
#include <maps/definitions.hpp>
#include <maps/misc.hpp>
#include <topology/SortableReturnData.hpp>

using namespace xavier;

namespace topology {



/** Detecting Singularity crossing conditions: Returns 'true' if a singularity intersection
 *  exists between two subsequent points on a line on the map.
  *    Note:  It is assumed that the two points are fairly close together when operating this test.
  *    SRDATA is a template handle for SortableReturnData<>
  */
template <class MAP, class PARAM, class SRDATA>
bool detectSingularityIntersection(
    const vec_type& x0, const vec_type& x1,  //Positions of adjacent points
    std::set<SRDATA>& cache,                 //Cache storing mappings
    const MAP& pmap, const int& period, const PARAM& params)
{
    typedef typename MAP::lmat_type mat_type;
    typedef typename MAP::xstate_type xstate_type;
    typedef std::pair<vec_type,mat_type>  ProjPair;
    //Threshold parameters
    const double& MAX_REL_TIME = params.maxDeltaTimeFactor;
    const double& MAX_ABS_TIME = params.maxTimeChange;
    const vec_type& MAX_MAPDISP = params.maxDeltaMapDisplacement;
    static const int numSing = MAP::rhs_type::numSingularities;
    
    //Check if the nodes are singular within the cache
    if ( isNodeSingularInCache(x0,cache,period) ) {
        return true;
    }
    if ( isNodeSingularInCache(x1,cache,period) ) {
        return true;
    }
    
    //Gather data values from subsequent edge nodes
    typename std::set<SRDATA>::iterator cit0,cit1;
    cit0 = cache.find( SRDATA(x0) ); //Note:  Should be found
    cit1 = cache.find( SRDATA(x1) ); //Note:  Should be found
    vec_type p0 = (*cit0).returns[std::abs(period)-1]; //P^p(x0)
    vec_type p1 = (*cit1).returns[std::abs(period)-1]; //P^p(x1)
    const xstate_type* singPtr = pmap.rhs().singularities();
    //This has to fail the time test (otherwise will be on opposite sides of P1)
    //Check if propagation time is jumping too much
    //  - Time is last component of SRDATA elements
    //  - This is the BEST check for separation without singularities
    double dt0 = (*cit0).data[std::abs(period)-1][numSing];
    double dt1 = (*cit1).data[std::abs(period)-1][numSing];
    if (fabs(dt1-dt0) >= MAX_ABS_TIME) {
        return false;
    }
    if (fabs((dt1-dt0)/dt0) >= MAX_REL_TIME) {
        return false;    //Tangency not singularity
    }
    
    //Check if subsequent points possess a singularity intersection between them
    for(int k=0; k<numSing; k++) {
        double singularitySafeDistance = pmap.rhs().getSingularitySafeDistance(k);
        double minCloseApproachDist0 = (*cit0).data[std::abs(period)-1][k];
        double minCloseApproachDist1 = (*cit1).data[std::abs(period)-1][k];
        //Check for sign change in position coordinate from singularity
        ProjPair singPair = pmap.section().project( singPtr[k] );
        vec_type s = singPair.first;
        vec_type r0 = p0 - s;
        vec_type r1 = p1 - s;
        //If a point is within singularity limits
        if(minCloseApproachDist0<=singularitySafeDistance ||
                minCloseApproachDist1<=singularitySafeDistance) {
            int j=0; //Position Coord (rx changes sign=CR3BP HARD CODED!)
            if( r0[j]*r1[j]<0 ) { //sign change indicating not on same side of primary
                return true;
            }
            int i=1; //Velocity coord (xdot changes sign=CR3BP HARD CODED!)
            if( p0[i]*p1[i]<0 ) { //sign change on actual map coordinate
                return true;
            }
        }
        //The tricky hybrid of Singularity and Tangency
        int j=0; //Position Coord
        if( (r0[j]*r1[j]<0) &&  //Opposite sides of the primary
                ((fabs(r0[j])+fabs(r1[j]))<= MAX_MAPDISP[j])) { // AND in close proximity
            return true;
        }
    }
    
    //Otherwise, there is no signularity crossing
    return false;
}


/** Detecting if a node within a cache has a singularity for the specified period.
 *  ->Nodes possessing singularities will not have the appropriate data, and thus,
 *    produce an addressing error when referenced.
 *  Note: This is a simple work-around for these detection functions...
 */
template <class SRDATA>
bool isNodeSingularInCache( const vec_type& x, std::set<SRDATA>& cache, const int& period)
{
    typename std::set<SRDATA>::const_iterator cit;
    bool isSingular = false;
    cit = cache.find( SRDATA(x) );
    if (cit != cache.end() && cit->isThere(period) ) {
        isSingular = false;
    } else {
        isSingular = true;    //i.e., not found in cache
    }
    
    return isSingular;
}

/** Detecting Section Transversality violations:  Returns 'true' if transversality is violated.
  *    Note:  It is assumed that the two points are fairly close together when operating this test.
  *    SRDATA is a template handle for SortableReturnData<>
  */
template <class MAP, class PARAM, class SRDATA>
bool detectSectionSeparation(
    const vec_type& x0, const vec_type& x1,  //Positions of adjacent points in Section-coords
    std::set<SRDATA>& cache,     //Cache storing mappings
    const MAP& pmap, const int& period, const PARAM& params)
{
    //Threshold parameters
    //const double& MAX_TRANSSPEED = params.maxDeltaTransSpeed;
    const vec_type& MAX_MAPDISP = params.maxDeltaMapDisplacement;
    const double& MAX_REL_TIME = params.maxDeltaTimeFactor;
    const double& MAX_ABS_TIME = params.maxTimeChange;
    static const int numSing = MAP::rhs_type::numSingularities;
    //static const int stateDim = MAP::rhs_type::dimension;
    //const bool& debug = params.debug;
    const bool& debug = false;
    const int mapDim = (int)x0.size();
    
    //SECTION TRANSVERSALITY:  We also must check if the section is transverse to the flow
    // ->If transversality assumption is violated, we must indicate failure location and
    //   follow a slightly different procedure.
    if(debug) {
        std::cout << "-----------------------------------------------------\n";
        std::cout << " Detection Criteria :\n";
        std::cout << "  period = " << period << "\n";
        std::cout << "  x0 = " << x0 << " x1 = " << x1 << "\n";
    }
    //Check if either node is a singularity
    if ( isNodeSingularInCache(x0,cache,period) ) {
        return true;
    }
    if(debug) {
        std::cout << "  Node x0 is NOT singular in cache...\n";
    }
    if ( isNodeSingularInCache(x1,cache,period) ) {
        return true;
    }
    if(debug) {
        std::cout << "  Node x1 is NOT singular in cache...\n";
    }
    
    //Gather data values from subsequent edge nodes
    typename std::set<SRDATA>::iterator cit0,cit1;
    cit0 = cache.find( SRDATA(x0) ); //Note:  Should be found
    cit1 = cache.find( SRDATA(x1) ); //Note:  Should be found
    vec_type p0 = (*cit0).returns[std::abs(period)-1]; //P^p(x0)
    vec_type p1 = (*cit1).returns[std::abs(period)-1]; //P^p(x1)
    
    //Check if propagation time is jumping too much
    //  - Time is last component of SRDATA elements
    //  - This is the BEST check for separation without singularities
    double dt0 = (*cit0).data[std::abs(period)-1][numSing];
    double dt1 = (*cit1).data[std::abs(period)-1][numSing];
    if(debug) {
        std::cout << "  p0 = " << p0 << " p1 = " << p1 << "\n";
        std::cout << "  dt0 = " << dt0 << " dt1 = " << dt1 << "\n";
        std::cout << "  Relative Time Test: DeltaT/dt0 = " << fabs((dt1-dt0)/dt0) << " >= " << MAX_REL_TIME << " ? = "
                  << (fabs((dt1-dt0)/dt0) >= MAX_REL_TIME) << "\n";
    }
    if (fabs((dt1-dt0)/dt0) >= MAX_REL_TIME) {
        return true;
    }
    if(debug) {
        std::cout << "  Absolute Time Test: DeltaT = " << fabs(dt1-dt0) << " >= " << MAX_ABS_TIME << " ? = "
                  << (fabs(dt1-dt0) >= MAX_ABS_TIME) << "\n";
                  
    }
    if (fabs(dt1-dt0) >= MAX_ABS_TIME) {
        return true;
    }
    
    //Check if map displacment coordinate jumps too much
    vec_type mapDisp0 = p0 - x0; //P^p(x0)-x0
    vec_type mapDisp1 = p1 - x1; //P^p(x1)-x1
    //for(int i=0;i<(int)mapDim;i++) {//Should be this, but CR3BP HARD CODED!
    int i=0; //x-coord only!
    if(debug) {
        std::cout << "  Displacement Test: " << fabs(mapDisp1[i]-mapDisp0[i]) << " >= " << MAX_MAPDISP[i] << " ? = "
                  << (fabs(mapDisp1[i]-mapDisp0[i]) >= MAX_MAPDISP[i]) << "\n";
    }
    if (fabs(mapDisp1[i]-mapDisp0[i]) >= MAX_MAPDISP[i]) {
        return true;
    }
    //}
    
    
    //Also need a LOW-TO-HIGH Map Displacment Magnitude Test
    double mapDispMag0 = nvis::norm(mapDisp0);
    double mapDispMag1 = nvis::norm(mapDisp1);
    //If have a low displacement on one node,
    if( (mapDispMag0 < 0.5) || (mapDispMag1 < 0.5) ) {
        //see if map displacement magnitude jumps too much between nodes
        if(debug) {
            std::cout << "  LowToHigh Mag Test: " << fabs(mapDispMag1 - mapDispMag0) << " >= " << 0.5 << " ? = "
                      << ( fabs(mapDispMag1 - mapDispMag0) > 0.5 ) << "\n";
        }
        if ( fabs(mapDispMag1 - mapDispMag0) > 0.5 ) {
            return true;
        }
        //Note: high mag to high mag jumps are likely around singularities, but it
        //does not mean there is a singularity between them!  (Saddle at singularity)
    }
    
    //Map Coord Sign Change Tests
    for(int i=0; i<mapDim; i++) {
        //If sign flip in map coords AND too large of a change, then likely separation condition
        if(debug) {
            std::cout << "  MapCoord Sign Test: i=" << i << " : p0[i]*p1[i] = "
                      << p0[i]*p1[i] << " < 0.0 ? = " << (p0[i]*p1[i] < 0.0) << "\n";
            std::cout << "                      i=" << i << " : |p1[i]-p0[i]| = "
                      << fabs(p1[i]-p0[i]) << " >= 1.0 ? = " << (fabs(p1[i] - p0[i]) >= 1.0) << "\n";
                      
        }
        if( (p0[i]*p1[i] < 0.0) &&
                (fabs(p1[i] - p0[i]) >= 1.0) ) {
            return true;
        }
        //Note: Singularities may return positive in this test too!
    }
    
    //Compute the transverse velocity components of the mapping
    //double dy0 = pmap.section().transverseVelocity(p0);
    //double dy1 = pmap.section().transverseVelocity(p1);
    //Check that the transverse velocity component jumps too much
    //if (fabs(dy1-dy0) >= MAX_TRANSSPEED) return true;
    //Note: ydot and xdot go to +/- inf rapidly near primaries, causing this check to yield false positives!
    
    //Check if subsequent points possess a singularity intersection between them
    bool crossesSingularity = detectSingularityIntersection(x0,x1,cache,pmap,period,params);
    if(debug) {
        std::cout << "  Singularity Test : Crossing Singularity = " << crossesSingularity << "\n";
    }
    if (crossesSingularity) {
        return true;
    }
    
    //The two points passed the test, and their is likely not a Section Transversality Violation
    if(debug) {
        std::cout << "  Passed All Separation Tests\n";
        std::cout << "-----------------------------------------------------\n";
    }
    return false;
}

/** Detecting Section Transversality violations at a given point:
 *  Returns 'true' if transversality is violated.
 *    Note:  At a propagation point, the only post-mapping check is
 *           transverseVelocity(P^p(x0)) = 0.0;
 *    SRDATA is a template handle for SortableReturnData<>
 */
template <class MAP, class PARAM, class SRDATA>
bool detectSectionSeparation(
    const vec_type& x0,  //Positions of adjacent points in Section-coords
    std::set<SRDATA>& cache,     //Cache storing mappings
    const MAP& pmap, const int& period, const PARAM& params)
{
    //SECTION TRANSVERSALITY:  We must check if the section is transverse to the flow
    // ->If transversality assumption is violated, we must indicate failure location and
    //   follow a slightly different procedure.
    
    //std::cout << __FILE__ << ":detectSectionSep() x0 = " << x0 << " p= " << period << "\n";
    //std::cout << __FILE__ << ":isNodeSingularInCache() = " << isNodeSingularInCache(x0,cache,period) << "\n";
    
    //Check if node is a singularity based on empty cache data
    if ( isNodeSingularInCache(x0,cache,period) ) {
        return true;
    }
    
    //Gather data values from edge node
    typename std::set<SRDATA>::iterator cit0;
    cit0 = cache.find( SRDATA(x0) ); //Note:  Should be found
    vec_type p0 = (*cit0).returns[std::abs(period)-1]; //P^p(x0)
    
    //Transverse velocity component of mapping
    double dy0;
    try {
        dy0 = pmap.section().transverseVelocity(p0);
    } catch(...) {
        //Any error indicates that we are on the ZVC (at zero to numerical precision)
        dy0 = 0.0;
    }
    //std::cout << __FILE__ << ":transVel(x0) = " << dy0 << " <1.e-5 =" << (fabs(dy0) <= 1.e-5) << "\n";
    if(fabs(dy0) <= 1.e-5) {
        return true;
    }
    
    //Passed transversality test as long as mapping was ok:
    return false;
}

} //end topology

#endif // __SECTION_TRANSVERSALITY_HPP__