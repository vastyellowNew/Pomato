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


//////////////////////////////////////////////////////////////////////////////
// STHManifold.hpp
// Author:  Wayne Schlei & Xavier Tricoche (Purdue University)
// Date:    3/12/2015
// Purpose:  Utilize the ManCurve algorithm for computing manifolds on a
// Poincare section. This method is loosely based on "Computing 1D Global
// Manifolds By Continuation" by England et al. (SIAM Dyn Sys., 2005)
// but has been modified to utilize curve-refinement ideas.  Also, this method
// includes ways of handling non-transverse sections.
//
//
// Future Work:
// 1)This propagates manifolds on a 2D surface of section; for work with higher
// dimensional sections, be sure to conform to a template that specifies dimensions.
///////////////////////////////////////////////////////////////////////////////

#ifndef STH_MANIFOLD_HPP
#define STH_MANIFOLD_HPP

#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <topology/ManifoldClasses.hpp>
#include <topology/SortableReturnData.hpp>
#include <topology/AdaptiveEdge.hpp>
#include <topology/EdgeRotationFailure.hpp> //Discontinuities
#include <topology/SectionTransversality.hpp>
#include <topology/FixedPointTests.hpp>
#include <topology/EdgeRotationMapCalls.hpp>

using namespace xavier;

namespace topology {

/// Function to determine if downstream map points are too close together based on delta_min setting
template<class MAP,class PARAM>
bool isSegmentTooShort(
    const MAP& theMap, const PARAM& params,
    const typename MAP::lvec_type& p0,       // downstream point 0
    const typename MAP::lvec_type& p1,       // downstream point 1
    const ManifoldSettings& settings        // quality settings
)
{
    typedef typename MAP::lvec_type   VecType;
    const metric_type& theMetric = params.the_metric;
    
    //Evaluate the curve parameter delta (distance)
    double delta = theMetric.distance(p1,p0);
    
    //Downstream points are too close together => Done!
    if (fabs(delta) < settings.delta_min) {
        return true;
    }
    
    //Passed
    return false;
}

/// Function to determine if downstream map points are too far apart based on delta_max setting
template<class MAP,class PARAM>
bool isSegmentTooLong(
    const MAP& theMap, const PARAM& params,
    const typename MAP::lvec_type& p0,       // downstream point 0
    const typename MAP::lvec_type& p1,       // downstream point 1
    const ManifoldSettings& settings        // quality settings
)
{
    typedef typename MAP::lvec_type   VecType;
    const metric_type& theMetric = params.the_metric;
    
    //Evaluate the curve parameter delta (distance)
    double delta = theMetric.distance(p1,p0);
    
    //Downstream points are too close together => Done!
    if (fabs(delta) > settings.delta_max) {
        return true;
    }
    
    //Passed
    return false;
}

/// Determine if two segments are vastly different in magnitude (true = more than 3 orders)
template<class MAP,class PARAM>
bool segmentMagnitudeComparison(
    const MAP& theMap, const PARAM& params,
    const typename MAP::lvec_type& p0,       // downstream point 0
    const typename MAP::lvec_type& p1,       // downstream point 1
    const typename MAP::lvec_type& p2,       // downstream point 2
    const ManifoldSettings& settings        // quality settings
)
{
    typedef typename MAP::lvec_type   VecType;
    const metric_type& theMetric = params.the_metric;
    
    //Evaluate the curve parameter delta (distance)
    double delta0 = theMetric.distance(p1,p0);
    double delta1 = theMetric.distance(p2,p1);
    
    //Compare the relative magnitudes
    double lm = std::min(delta0,delta1);
    double rel = std::fabs(delta1 - delta0) / lm;
    
    //Downstream segments are more than 3 orders of magnitude different
    if (rel > 1000.0) {
        return true;
    }
    
    //Passed
    return false;
}

/// Checking if three successive points (prev,y,next) of downstream manifold (working segment)
/// meet ManifoldSettings (curve-refinement criteria)
template<class MAP,class PARAM,class SRDATA>
bool isManifoldCurveOK(
    const MAP& theMap, int period, const PARAM& params,
    const typename MAP::lvec_type& xp,    // upstream, prior to testing seg
    const typename MAP::lvec_type& x0,    // upstream point of current seg->pt0
    const typename MAP::lvec_type& x1,    // upstream point of current seg->pt1
    const typename MAP::lvec_type& prev,  // (downstream) point prior to testing seg
    const typename MAP::lvec_type& p0,    // (downstream) point of current seg->pt0
    const typename MAP::lvec_type& p1,    // (downstream) point of current seg->pt1
    const ManifoldSettings& settings,     // quality settings
    std::set<SRDATA>& cache,              // cache holding mapping data
    bool invert = false                   // invert segs so (prev,p0) is the current segment
)
{
    typedef typename MAP::lvec_type   VecType;
    typedef typename MAP::gvec_type   FullState;
    const metric_type& theMetric = params.the_metric;
    
    //Note:  Testing segment and comparison segment:
    // If (invert)
    //    TestSeg       = [prev, p0]
    //    ComparisonSeg = [  p0, p1]
    // Else
    //    TestSeg       = [  p0, p1]
    //    ComparisonSeg = [prev, p0]
    // End
    //
    //Note:  Comparison segment must be BOTH
    //     1) Transverse (i.e., no section separation between nodes)
    //     2) NOT Too Short (delta < delta_min) => Leads to erroneous angles
    
    //Evaluate the curve parameters of tuple
    double delta, alpha;
    alpha = theMetric.angle(prev,p0,p1);
    if (invert) {
        delta = theMetric.distance(p0,prev);
    } else {
        delta = theMetric.distance(p1,p0);
    }
    
    /*//Display refinement conditions
    std::cout << "    Delta = " << delta << "\n";
    std::cout << "    [DeltaMin = " << settings.delta_min << "] [DeltaMax = "
              << settings.delta_max << "]\n";
    std::cout << "    Alpha = " << alpha << " [Alpha_Max = " << settings.alpha_max << "]\n";
    std::cout << "    |Delta*Alpha| = " << fabs(alpha*delta) << " [(da)_Max = "
              << settings.delta_alpha_max << "]\n";*/
    
    
    //Section Transversality : Both segments have to be OK
    //--------------------------------------------------------------------------------------
    //if ( detectSectionSeparation(x0,x1,cache,theMap,period,params) )
    //    std::cout << "     Separation Detected On Segment!\n";
    if ( invert ) {
        //Test on 'current' segment
        if ( detectSectionSeparation(xp,x0,cache,theMap,period,params) ) {
            return false;
        }
        //Test on 'comparison' (or other) segment
        if ( detectSectionSeparation(x0,x1,cache,theMap,period,params) ) {
            return false;
        }
    } else {
        //Test on 'current' segment
        if ( detectSectionSeparation(x0,x1,cache,theMap,period,params) ) {
            return false;
        }
        //Test on 'comparison' (or other) segment
        if ( detectSectionSeparation(xp,x0,cache,theMap,period,params) ) {
            return false;
        }
    }
    
    
    //Check for violations of ZVC at 1/4,1/2,3/4 points of downstream test segment:
    VecType q, h, q3;
    if (invert) {
        double t = 0.25;
        q = (1.0 - t)*prev + t*p0;
        t = 0.5;
        h = (1.0 - t)*prev + t*p0;
        t = 0.75;
        q3 = (1.0 - t)*prev + t*p0;
    } else {
        double t = 0.25;
        q = (1.0 - t)*p0 + t*p1;
        t = 0.5;
        h = (1.0 - t)*p0 + t*p1;
        t = 0.75;
        q3 = (1.0 - t)*p0 + t*p1;
    }
    FullState y;
    //If error in reverse projection, violates bounds of section
    try {
        y = theMap.section().unproject( q );
    } catch(...) {
        return false;
    }
    try {
        y = theMap.section().unproject( h );
    } catch(...) {
        return false;
    }
    try {
        y = theMap.section().unproject( q3 );
    } catch(...) {
        return false;
    }
    
    
    
    //Curve-refinement criteria:
    //--------------------------------------------------------------------------------------
    //Downstream points are too close together => Done!
    if (fabs(delta) < settings.delta_min) {
        return true;
    }
    
    //Companion segment must not be too small  => Use Magnitude comparison, NOT generic shortness...
    if ( segmentMagnitudeComparison(theMap,params, prev, p0, p1, settings) ) {
        return false;
    }
    /*if (invert) {
      if ( isSegmentTooShort(theMap,params,p0,p1,settings) ) return false;
    } else {
      if ( isSegmentTooShort(theMap,params,prev,p0,settings) ) return false;
    }*/
    
    //Downstream point distance (A bit tricky, easier to shut off)
    //if (fabs(delta) > settings.delta_max ) return false; //Strict adherence has trouble handling asymptotes
    //Try delta_max ONLY when a segment's returns are both within the ORIGINAL domain
    if( (!invert && params.bounds.inside(p0)   && params.bounds.inside(p1) ) ||
            ( invert && params.bounds.inside(prev) && params.bounds.inside(p0) ) ) {
        if (fabs(delta) > settings.delta_max ) {
            return false;
        }
    }
    
    //Angle of points
    if (fabs(alpha) > settings.alpha_max) {
        return false;
    }
    
    //Arc length -> Control sampling distance
    if (fabs(alpha*delta) > settings.delta_alpha_max) {
        return false;
    }
    
    
    //Passed
    return true;
}

/// Checking if three successive points (prev,y,next) of downstream manifold (working segment)
/// meet ManifoldSettings (curve-refinement criteria)
/// -> Assumes NO section separation, which is useful for weakStrength manifolds
template<class MAP,class PARAM>
bool isManifoldCurveOK(
    const MAP& theMap, const PARAM& params,
    const typename MAP::lvec_type& prev,  // (downstream) point prior to testing seg
    const typename MAP::lvec_type& p0,    // (downstream) point of current seg->pt0
    const typename MAP::lvec_type& p1,    // (downstream) point of current seg->pt1
    const ManifoldSettings& settings,     // quality settings
    bool invert = false                   // invert segs so (prev,p0) is the current segment  (Don't need here!)
)
{
    typedef typename MAP::lvec_type   VecType;
    typedef typename MAP::gvec_type   FullState;
    const metric_type& theMetric = params.the_metric;
    
    //Note:  Testing segment and comparison segment:
    // If (invert)
    //    TestSeg       = [prev, p0]
    //    ComparisonSeg = [  p0, p1]
    // Else
    //    TestSeg       = [  p0, p1]
    //    ComparisonSeg = [prev, p0]
    // End
    //
    //Note:  Comparison segment must be BOTH
    //     1) Transverse (i.e., no section separation between nodes)
    //     2) NOT Too Short (delta < delta_min) => Leads to erroneous angles
    
    //Evaluate the curve parameters of tuple
    double delta, alpha;
    alpha = theMetric.angle(prev,p0,p1);
    if (invert) {
        delta = theMetric.distance(p0,prev);
    } else {
        delta = theMetric.distance(p1,p0);
    }
    
    /*//Display refinement conditions
    std::cout << "    Delta = " << delta << "\n";
    std::cout << "    [DeltaMin = " << settings.delta_min << "] [DeltaMax = "
              << settings.delta_max << "]\n";
    std::cout << "    Alpha = " << alpha << " [Alpha_Max = " << settings.alpha_max << "]\n";
    std::cout << "    |Delta*Alpha| = " << fabs(alpha*delta) << " [(da)_Max = "
              << settings.delta_alpha_max << "]\n";*/
    
    
    //Section Transversality : Both segments have to be OK => Assumed with weakStrength Manifold
    //--------------------------------------------------------------------------------------
    //if ( detectSectionSeparation(x0,x1,cache,theMap,period,params) )
    //    std::cout << "     Separation Detected On Segment!\n";
    /*if ( invert ) {
      //Test on 'current' segment
      if ( detectSectionSeparation(xp,x0,cache,theMap,period,params) )
        return false;
      //Test on 'comparison' (or other) segment
      if ( detectSectionSeparation(x0,x1,cache,theMap,period,params) )
        return false;
    } else {
      //Test on 'current' segment
      if ( detectSectionSeparation(x0,x1,cache,theMap,period,params) )
        return false;
      //Test on 'comparison' (or other) segment
      if ( detectSectionSeparation(xp,x0,cache,theMap,period,params) )
        return false;
    }*/
    
    //Check for violations of ZVC at 1/4,1/2,3/4 points of downstream test segment:
    VecType q, h, q3;
    if (invert) {
        double t = 0.25;
        q = (1.0 - t)*prev + t*p0;
        t = 0.5;
        h = (1.0 - t)*prev + t*p0;
        t = 0.75;
        q3 = (1.0 - t)*prev + t*p0;
    } else {
        double t = 0.25;
        q = (1.0 - t)*p0 + t*p1;
        t = 0.5;
        h = (1.0 - t)*p0 + t*p1;
        t = 0.75;
        q3 = (1.0 - t)*p0 + t*p1;
    }
    FullState y;
    //If error in reverse projection, violates bounds of section
    try {
        y = theMap.section().unproject( q );
    } catch(...) {
        return false;
    }
    try {
        y = theMap.section().unproject( h );
    } catch(...) {
        return false;
    }
    try {
        y = theMap.section().unproject( q3 );
    } catch(...) {
        return false;
    }
    
    //Curve-refinement criteria:
    //--------------------------------------------------------------------------------------
    //Downstream points are too close together => Done!
    if (fabs(delta) < settings.delta_min) {
        return true;
    }
    
    //Companion segment must not be too small => Use Magnitude comparison, NOT generic shortness...
    if ( segmentMagnitudeComparison(theMap,params, prev, p0, p1, settings) ) {
        return false;
    }
    /*if (invert) {
      if ( isSegmentTooShort(theMap,params,p0,p1,settings) ) return false;
    } else {
      if ( isSegmentTooShort(theMap,params,prev,p0,settings) ) return false;
    }*/
    
    //Downstream point distance (A bit tricky, easier to shut off)
    //if (fabs(delta) > settings.delta_max ) return false; //Strict adherence has trouble handling asymptotes
    //Try delta_max ONLY when a segment's returns are both within the ORIGINAL domain
    if( (!invert && params.bounds.inside(p0)   && params.bounds.inside(p1) ) ||
            ( invert && params.bounds.inside(prev) && params.bounds.inside(p0) ) ) {
        if (fabs(delta) > settings.delta_max ) {
            return false;
        }
    }
    
    //Angle of points
    if (fabs(alpha) > settings.alpha_max) {
        return false;
    }
    
    //Arc length -> Control sampling distance
    if (fabs(alpha*delta) > settings.delta_alpha_max) {
        return false;
    }
    
    
    //Passed
    return true;
}

/// Comparison struct for a std::pair<ivec2,DATA> set or list
template<typename DATA>
struct DataPairFirstCompare {
    typedef DataPairFirstCompare<DATA>   SelfType;
    typedef std::pair<nvis::ivec2,DATA>  PairType;
    bool operator()(const PairType& lhs, const PairType& rhs) const
    {
        //Use AdaptiveEdge::LineIDCompare
        LineIDCompare<nvis::ivec2> compare;
        return compare(lhs.first,rhs.first);
    }
};

/// Condition for checking on a DataPair if DataPair.ID > given ID
template<typename DATAPAIR>
struct DataPairBeyondID {
    typedef DataPairBeyondID<DATAPAIR>   SelfType;
    typedef nvis::ivec2              IdType;
    
    DataPairBeyondID(const IdType& input) : threshold(input) {}
    
    bool operator()(const DATAPAIR& other) const
    {
        //True when otherID > threshold
        LineIDCompare<IdType> compare; //Less than
        return compare(threshold,other.first); //True when threshold < otherID
    }
    
    IdType threshold;
};

/// Comparison function for an On-Edge compare for an EdgeRotationFailure object
template<typename MAPDISCONT>
bool mapDiscontOnEdgeCompare(const MAPDISCONT& lhs, const MAPDISCONT& rhs)
{
    //Use Adaptiveedge::LineIDCompare
    LineIDCompare<nvis::ivec2> compare;
    return compare(lhs.edgeID,rhs.edgeID);
}


//------------------------------------------------------------------------------------------------------
//Note: The next functions are for a crude version of the STHManifold generation algorithm that
//utilizes curve-refinement and heuristics to compute/advect manifolds.  It has no processing
//that allows for extraction of arcs and data.  Also, it is generally slow with one manifold per thread.
//------------------------------------------------------------------------------------------------------
/// Function to subdivide seed segment such that working segment meets specified ManifoldSettings (with discontinuities)
template<typename MAP,typename PARAM,typename SRDATA,typename MAPDISCONT>
void adaptiveManfioldSubdivision(
    const vec_type& x0, const vec_type& x1,
    const MAP& theMap, const int period, const PARAM& params,
    ManifoldProgress& progress, const ManifoldSettings& settings,
    std::set<SRDATA>& cache,
    std::vector<vec_type>& downstreamPoints,
    std::list<int>& downstreamSepPtIdx,
    std::list<MAPDISCONT>& sepFailureList
)
{
    typedef typename MAP::lvec_type            vec_type;
    typedef typename MAP::lmat_type            mat_type;
    typedef typename MAP::xstate_type          xstate_type;
    typedef typename MAP::lbox_type            BoundsType;
    typedef std::pair<vec_type,mat_type>       ProjPair;
    //typedef EdgeRotationFailure<vec_type>      MapDiscont;
    typedef MAPDISCONT                         MapDiscont;
    //Adaptive Edge data structures
    typedef AdaptiveEdge<vec_type,vec_type>    EdgeType;
    //typedef AdaptiveEdge<ManifoldPoint,vec_type>  EdgeType; //After upgrades
    typedef typename EdgeType::IdType          IdType;
    typedef std::pair<IdType,IdType>           Segment;
    typedef typename EdgeType::DataType        DataPairs;
    
    static const int numSing = MAP::rhs_type::numSingularities;
    
    //Ideally, we would like to continue to subdivide segments to max dx (=lmin)
    //for every section separation location, but this is time-prohibitive.
    
    //Evaluate the Angle displacement
    bool verbose = params.verbose;
    const double edgeLength = nvis::norm(x1-x0);
    const double dTauMin = settings.delta_tau_min;
    //const double dTauMax = edgeLength / 8.0; //Force at least 2 subdivisions
    bool fwd = (period > 0);
    //In theory, we don't need as many subdivisions knowing we have Discontinuities
    
    //Adaptive Edge Data structure
    EdgeType adaptiveEdge(x0,x1); //This is on the current seeding segment
    
    //By default, we have to add the midpoint of the seeding segment
    adaptiveEdge.refine(IdType(0,0));
    vec_type xm = (x0+x1) / 2.0;
    adaptiveEdge.setVertex(IdType(1,1), xm);
    
    //List of segments to try
    std::vector<Segment>  segsToCheck;
    segsToCheck.push_back( Segment(IdType(0,1),IdType(1,1)) );
    segsToCheck.push_back( Segment(IdType(1,1),IdType(1,0)) );
    
    
    //Compute the mappings for the endpoints (they should exist)
    MapDiscont eFail;
    vec_type p0, pm, p1;
    //Left point
    try {
        p0 = mapUsingCache(x0, theMap, period, params, cache, eFail);
    } catch(MapDiscont& ef) {
        //Singularity caught
        ef.edgeID = IdType(0,1);
        adaptiveEdge.insertDiscontinuityNode(IdType(0,1),ef);
        p0 = vec_type(50);
    }
    //Midpoint
    try {
        pm = mapUsingCache(xm, theMap, period, params, cache, eFail);
    } catch(MapDiscont& ef) {
        //Singularity caught
        ef.edgeID = IdType(1,1);
        adaptiveEdge.insertDiscontinuityNode(IdType(1,1),ef);;
        pm = vec_type(50);
    }
    //Right point
    try {
        p1 = mapUsingCache(x1, theMap, period, params, cache, eFail);
    } catch(MapDiscont& ef) {
        //Singularity caught
        ef.edgeID = IdType(1,0);
        adaptiveEdge.insertDiscontinuityNode(IdType(1,0),ef);;
        p1 = vec_type(50);
    }
    
    //Add the map information into AdaptiveEdge
    adaptiveEdge.setValue( IdType(0,1), p0 );
    adaptiveEdge.setValue( IdType(1,1), pm );
    adaptiveEdge.setValue( IdType(1,0), p1 );
    
    
    //Edge subdivision algorithm for seeding segment
    int emptyNodes = 1;
    int numSubdivisions = 0;
    while ( emptyNodes > 0 ) {
        std::vector<Segment> segs;
        typename std::vector<Segment>::iterator segIT, segITm1, segITp1;
        //Check each segment for subdivision    ----------------------------
        adaptiveEdge.getAllSegments(segs);
        segIT = segITm1 = segs.begin();
        for(; segIT!=segs.end(); ++segIT) {
            segITm1 = segIT;
            --segITm1;
            segITp1 = segIT;
            ++segITp1;
            //Work on previous and current segment
            bool firstSeg = (segIT==segs.begin());
            bool lastSeg = nvis::all(segIT->second == IdType(1,0));
            
            //std::cout << "Testing segment : " << segIT->first << " , " << segIT->second << "\n";
            
            //Checking if we have to subdivide based on three nodes
            vec_type xA = adaptiveEdge.getVertex(segIT->first);
            vec_type xB = adaptiveEdge.getVertex(segIT->second);
            vec_type xC = adaptiveEdge.getVertex(
                              (firstSeg)? segITp1->second : segITm1->first);
            //Get parameters
            double tauA = adaptiveEdge.getParameter(segIT->first);
            double tauB = adaptiveEdge.getParameter(segIT->second);
            //double tauC = adaptiveEdge.getParameter(
            //               (firstSeg)? segITp1->second : segITm1->first);
            bool subdivide = false;
            //If testing segment is too far apart, subdivide
            double dTau = tauB - tauA;
            //if(dTau > dTauMax) subdivide = true; //Don't add any extra subdivision
            //Check for a discontinuity point, subdivide if still above dTauMin
            if(!subdivide && (dTau > dTauMin)) {
                //Are Node 0 and Node 1 both discontinuities? stop
                if(adaptiveEdge.isDiscontinuity(segIT->first) &&
                        adaptiveEdge.isDiscontinuity(segIT->second)) {
                    //Stop
                    continue;
                }
                //Is Node 0 or Node 1 a discontinuity point? subdivide?
                else if(adaptiveEdge.isDiscontinuity(segIT->first) ||
                        adaptiveEdge.isDiscontinuity(segIT->second)) {
                    //subdivide the node
                    subdivide = true;
                }
            }
            //Check the downstream manifold for curve-refinement parameters
            if(!subdivide) {
                //Get map data (downstream)
                vec_type pA = adaptiveEdge.getValue(segIT->first);
                vec_type pB = adaptiveEdge.getValue(segIT->second);
                vec_type pC = adaptiveEdge.getValue(
                                  (firstSeg)? segITp1->second : segITm1->first);
                /*std::cout << " Data Points:\n";
                std::cout << "   Seg: x0    = " << xA << " , x1    = " << xB << "\n";
                std::cout << "   Seg: P(x0) = " << pA << " , P(x1) = " << pB << "\n";
                std::cout << "   Prev: xC = " << xC << " , P(xC) = " << pC << "\n";*/
                
                
                //Run the downstream curve-refinement checks
                bool curveOK = false;
                //If primary segment (A-B) has separation, then curve is not ok
                if( !detectSectionSeparation(xA,xB,cache,theMap,period,params) ) {
                    if(firstSeg) {
                        //On first seg, we have to use the next's segments results to test curve
                        curveOK = isManifoldCurveOK(
                                      theMap,period,params,xA,xB,xC,pA,pB,pC,settings,cache);
                        //otherwise, curve is not ok
                    } else {
                        //Test the other segment (C-A) for separation
                        bool otherSegHasSep =
                            detectSectionSeparation(xC,xA,cache,theMap,period,params);
                        if(otherSegHasSep) {
                            //if (lastSeg) curveOK = false;
                            if(!lastSeg) {
                                //Try to use the next point to test segment
                                vec_type xD = adaptiveEdge.getVertex(segITp1->second);
                                vec_type pD = adaptiveEdge.getValue(segITp1->second);
                                curveOK = isManifoldCurveOK(
                                              theMap,period,params,xA,xB,xD,pA,pB,pD,settings,cache);
                            }
                        } else {
                            //Standard segment test
                            curveOK = isManifoldCurveOK(
                                          theMap,period,params,xC,xA,xB,pC,pA,pB,settings,cache);
                        }
                    }
                }
                //std::cout << "   Curve-Refinement = " << curveOK << "\n";
                
                //Decide to subdivide segment based on manifold curve params
                if (curveOK) {
                    continue;   //Done with segment
                } else if ( dTau <= dTauMin ) {
                    //std::cout << "  Segment at dTauMin: " << dTau << " < dTauMin=" << dTauMin << "\n";
                    MapDiscont rotErr; //Empty discontinuity for checking errors
                    //First, check if nodes are discontinuities (any type),
                    //  if yes, no further errors to report
                    if ( adaptiveEdge.isDiscontinuity(segIT->first) ||
                            adaptiveEdge.isDiscontinuity(segIT->second) ) {
                        continue; //Done
                    }
                    //Check for fixed point conditions (of same period)
                    else if(
                        isFixedPointSuspected(xA,pA-xA,xB,pB-xB,period,params,rotErr)
                    ) {
                        if (verbose) {
                            std::cerr << "Period (" << period << "): Fixed point is suspected on the manifold segment!\n";
                        }
                        rotErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                        adaptiveEdge.insertDiscontinuity(rotErr);
                        continue; //Done
                    }
                    /// If at lmin without successful rotation conditions nor a fixed point suspected,
                    /// check if it's a section-separation condition (test both forward and backward):
                    bool discontinuityFound = false;
                    //Separation
                    if (detectSectionSeparation(xA,xB,cache,theMap,period,params) ) {
                        // If we are at lmin, store the forward separation point
                        vec_type intersection = xA+xB;
                        intersection /= 2.0;
                        MapDiscont sepErr((fwd)? MapDiscont::FORWARD_SECTION_SEP : MapDiscont::BACKWARD_SECTION_SEP,
                                          intersection,period);
                        sepErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                        if (detectSingularityIntersection(xA,xB,cache,theMap,period,params) ) {
                            sepErr.setType((fwd)? MapDiscont::FORWARD_SINGULARITY : MapDiscont::BACKWARD_SINGULARITY);
                            //Data output for debugging
                            typename std::set<SRDATA>::iterator cit0,cit1;
                            cit0 = cache.find( SRDATA(xA) ); //Note:  Should be found
                            cit1 = cache.find( SRDATA(xB) ); //Note:  Should be found
                            vec_type p0 = (*cit0).returns[std::abs(period)-1]; //P^p(x0)
                            vec_type p1 = (*cit1).returns[std::abs(period)-1]; //P^p(x1)
                            double dt0 = (*cit0).data[std::abs(period)-1][numSing];
                            double dt1 = (*cit1).data[std::abs(period)-1][numSing];
                            const xstate_type* singPtr = theMap.rhs().singularities();
                            std::cout << " Singularity detected :\n";
                            std::cout << "  period = " << period << "\n";
                            std::cout << "  x0 = " << xA << " x1 = " << xB << "\n";
                            std::cout << "  p0 = " << p0 << " p1 = " << p1 << "\n";
                            std::cout << "  dt0 = " << dt0 << " dt1 = " << dt1 << "\n";
                            for(int k=0; k<numSing; k++) {
                                std::cout << "  Singularity " << k << ":\n";
                                double singularitySafeDistance = theMap.rhs().getSingularitySafeDistance(k);
                                double minCloseApproachDist0 = (*cit0).data[std::abs(period)-1][k];
                                double minCloseApproachDist1 = (*cit1).data[std::abs(period)-1][k];
                                //Check for sign change in position coordinate from singularity
                                ProjPair singPair = theMap.section().project( singPtr[k] );
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
                    }
                    //At this point, we don't know what went wrong, so should throw 'unknown'
                    // HOWEVER, it is likely that these are still separation conditions (likely singularity
                    // intersections) that are not detected by the heuristics because the subdivision isn't
                    // refined enough.  Thus, lets pass these as singularity conditions
                    if (!discontinuityFound) {
                        vec_type intersection = xA+xB;
                        intersection /= 2.0;
                        MapDiscont singErr(MapDiscont::UNRESOLVED,intersection,period);
                        singErr.edgeID = adaptiveEdge.computeChildID(segIT->first,1);
                        adaptiveEdge.insertDiscontinuity( singErr );
                        std::cout << " Unknown Discontinuity Found: " << singErr.what() << "\n";
                        std::cout << " Data Points:\n";
                        std::cout << "   Seg: x0    = " << xA << " , x1    = " << xB << "\n";
                        std::cout << "   Seg: P(x0) = " << pA << " , P(x1) = " << pB << "\n";
                        std::cout << "   Prev: xC = " << xC << " , P(xC) = " << pC << "\n";
                        
                    }
                    continue; //Force a completion of this segment since at dTauMin
                    
                } else {
                    subdivide = true;    //Curve needs refinement
                }
            }
            
            //Run subdivision step
            if (subdivide) {
                numSubdivisions++;
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
        }
        //------------------------------------------------------------------
        //Propagate new nodes
        std::vector<IdType> seedIDs;
        adaptiveEdge.getEmptyLeaves(seedIDs);
        //Leaves that need data (and thus testing on round)
        emptyNodes = (int) seedIDs.size();
        typename std::vector<IdType>::iterator seedIT;
        for(seedIT=seedIDs.begin(); seedIT!=seedIDs.end(); ++seedIT) {
            //std::cout << " Propagation of New Nodes:\n";
            //std::cout << "   SeedID = " << (*seedIT) << "\n";
            //Get position node from structure
            vec_type x = adaptiveEdge.getVertex((*seedIT));
            //std::cout << "   Seed x = " << x << "\n";
            vec_type px;
            bool detectedSeparation = false;
            try {
                px = mapUsingCache(x, theMap, period, params, cache, eFail);
            } catch(MapDiscont& ef) {
                //Singularity caught
                ef.edgeID = (*seedIT);
                adaptiveEdge.insertDiscontinuityNode((*seedIT),ef);;
                px = vec_type(50);
                detectedSeparation = true;
                std::cout << "   Manifold Mapping Encountered Error = " << ef.what() << "\n";
            }
            //Check for a fixed point at this new point
            //if( !detectedSeparation && isFixedPointSuspected(x,px-x,period,params,eFail) ) {
            //    adaptiveEdge.insertDiscontinuityNode((*seedIT),eFail);
            //    detectedSeparation = true;
            //    std::cout << "   Fixed Point Detected = " << eFail.what() << "\n";
            //}
            //Run test for transversality condition (dydot(P^p(x)) != 0)
            if( !detectedSeparation ) {
                //Forward map
                if(detectSectionSeparation(x,cache,theMap,period,params)) {
                    MapDiscont nodeTransFail(
                        (fwd)? MapDiscont::FORWARD_SECTION_SEP_NODE : MapDiscont::BACKWARD_SECTION_SEP_NODE,
                        x,period);
                    nodeTransFail.edgeID = (*seedIT);
                    adaptiveEdge.insertDiscontinuityNode((*seedIT),nodeTransFail);
                    detectedSeparation = true;
                    std::cout << "   Section Separation Node = " << nodeTransFail.where() << "\n";
                }
            }
            
            //Store data into adaptive structure
            adaptiveEdge.setValue( (*seedIT), px );
        }
        
    } //End subdivision algorithm
    
    /// Output Data Ops: ------------------------------------------------------------------
    downstreamSepPtIdx.clear();
    adaptiveEdge.getDataInOrder(downstreamPoints);
    //Erase the first point since it's already stored in manifold from previous segment
    downstreamPoints.erase(downstreamPoints.begin());
    //If no subdivision, also erase next point.  (Usually on slow manifolds or weak saddles)
    //if(noSubdivision) downstreamPoints.erase(downstreamPoints.begin());
    adaptiveEdge.getDiscontinuityList(sepFailureList);//Upstream on working seg
    
    /// Sum through segments to update progress (i.e., arclength)
    std::vector<DataPairs> dataPairVec; //id->px pairs
    typename std::vector<DataPairs>::iterator dpIT;
    adaptiveEdge.getData(dataPairVec); //Comes in on-edge order
    dataPairVec.erase(dataPairVec.begin()); //To match downstreamPoints
    
    //Set that holds all downstream points and sep points in ID->Data pairs
    std::set<DataPairs,DataPairFirstCompare<vec_type> > wsDataSet;
    
    //Set analysis bounds for not counting large growth near singularities
    BoundsType bbox = params.bounds;
    //Config coords should never have singularity bounds
    //const double LARGE = std::numeric_limits<double>::max();
    //bbox.min()[0] = -LARGE; bbox.max()[0] = LARGE;  //Expanding x bounds?
    bool hasSingularities = (numSing > 0);
    
    //Pass through downstream data
    for(dpIT=dataPairVec.begin(); dpIT!=dataPairVec.end(); ++dpIT) {
        //Fill out a working segment data set <x,px>
        wsDataSet.insert((*dpIT));
        //Look for places near singularity that are out of analysis region
        if (hasSingularities && !bbox.inside( dpIT->second ) ) {
            //These points are going to xdot= +/-\inf
            //Don't visualize these points
            downstreamSepPtIdx.push_back( dpIT - dataPairVec.begin() );
            //Also, don't use in arclength computation
        }
    }
    
    
    //Should pass in STOPPER but it has inconsistent assumptions
    
    //Approaching (Homoclinic connection) test
    
    //Sort the sepFailureList
    sepFailureList.sort( mapDiscontOnEdgeCompare<MapDiscont> );
    
    //Add in separation failures to wsDataSet & find downstream IDs
    typename std::list<MapDiscont>::iterator sepListIT = sepFailureList.begin();
    for(; sepListIT!=sepFailureList.end(); sepListIT++) {
        //Find the downstreamPoints index just after a separation condition
        dpIT = std::find_if(dataPairVec.begin(),dataPairVec.end(),
                            DataPairBeyondID<DataPairs>(sepListIT->edgeID) );
        if (dpIT != dataPairVec.end()) {
            downstreamSepPtIdx.push_back( dpIT-dataPairVec.begin() );
        }
        //Add to set
        wsDataSet.insert( DataPairs(sepListIT->edgeID, sepListIT->failurePos) );
    }
    
    //Computing Arclength
    // 1) Include parts that are Transverse regions
    // 2) Exclude pieces that leave the analysis region
    double totalNewArcLength = 0.0;
    int numSepPts = (int) sepFailureList.size();
    sepListIT = sepFailureList.begin();
    typename std::set<DataPairs,DataPairFirstCompare<vec_type> >::iterator wsIT,wsLastIT;
    wsIT = wsLastIT = wsDataSet.begin();
    //For each working segment portion
    for(wsIT++; wsIT!=wsDataSet.end(); ++wsIT,++wsLastIT) {
        //Compute arclength downstream
        vec_type p1 = wsIT->second;
        vec_type p0 = wsLastIT->second;
        double aLength = params.the_metric.distance(p0,p1);
        
        //First check for leaving bounds
        bool inBounds = true;
        //Don't add arc lengths outside bounds (as they tend toward +/- inf)
        if ( !bbox.inside(p0) || !bbox.inside(p1) ) {
            inBounds=false;
        }
        
        
        //Check for separation points
        if(numSepPts > 0 ) {
            //Add to total if not a separation condition (compare ID's)
            bool currentIsSepPt = nvis::all( wsIT->first == sepListIT->edgeID );
            bool lastIsSepPt = nvis::all( wsLastIT->first == sepListIT->edgeID );
            if ( !currentIsSepPt && !lastIsSepPt ) {
                //Normal segment - if not near singularity
                if(inBounds) {
                    totalNewArcLength += aLength;
                }
            } else if (lastIsSepPt) {
                //If the trailing point is now on the sep point, switch to next sep point
                sepListIT++;
                if (sepListIT == sepFailureList.end()) {
                    sepListIT = sepFailureList.begin();
                }
            }
        } else {
            //Otherwise, just always add up arclengths a
            if(inBounds) {
                totalNewArcLength += aLength;
            }
        }
    }
    
    //Update Manifold Progress
    progress.length += totalNewArcLength;
    progress.seedingLength += edgeLength;
    
    //Sort the downstreamPtIdx list
    downstreamSepPtIdx.sort(); //sorting ints
    downstreamSepPtIdx.unique(); //no duplicates
}



/// The ManCurve algorithm for computing 1D manifolds as a curve-refinement process (one-manifold per thread)
template<typename MAP, typename STOP, typename PARAM>
void ManCurve(Separatrix& theManifold, ManifoldProgress& progress,
              const MAP& theMap, STOP& stopper, const PARAM& params,
              const fixpoint& saddle, const Perturbation& mType,
              const ManifoldSettings& settings)
{
    typedef typename MAP::section_type::lvec_type            VecType;
    typedef typename MAP::rhs_type                           RHStype;
    static const int s = RHStype::numSingularities;
    typedef nvis::fixed_vector<double,s+1>                   ExtendedMapDataVec;
    typedef SortableReturnData<vec_type,ExtendedMapDataVec>  SortableData;
    typedef EdgeRotationFailure<VecType>                     MapDiscont;
    
    //Reserve some space for manifold points
    theManifold.manifold.reserve(100);
    
    const metric_type& theMetric = params.the_metric;
    unsigned int i0, i1;
    double length;
    int period = saddle.K;
    MapDiscont eFail;
    //Propagation direction (stable or unstable)
    bool fwd = (mType == UNSTABLE_PLUS || mType == UNSTABLE_MINUS) ? true : false;
    int thePeriod = (fwd)? period : -period; //For Map calls
    //Plus or minus perturbation
    bool p_or_m = (mType == STABLE_PLUS || mType == UNSTABLE_PLUS) ? true : false;
    //Data set holding all mappings for easy reference
    std::set<SortableData> cache;
    
    //Find the initial step => Use function relating max eigenvalue and eps
    double eps = settings.eps;
    double epsMax = settings.sdelta_max;
    double lambda = fabs(saddle.eval[1]); //Unstable eigval(i.e., max, indicates saddle strength)
    //Piecewise function - Hard codeded for CR3BP
    if (!settings.manualStep) {
        if (lambda > 150.0) {
            double a = 1e-9 - 1e-7;
            double c = 1e-7;
            eps = a*(1.0-exp( (150.0-lambda)/200 )) + c;
        } else {
            double a = 1e-3 - 1e-7;
            double c = 1e-7;
            eps = a*(1.0-exp( 2.0*(150.0-lambda)/(1.0-lambda) )) + c;
        }
        epsMax = 20*eps;
    }
    if (params.verbose) {
        std::cerr << " ManCurve: lambdaMax = " << lambda << " eps = " << eps << "\n";
    }
    
    //if (progress.length == 0) {
    if(params.verbose) {
        theManifold.manifold.clear();
    }
    theManifold.manifold.push_back(saddle.pos);
    if (params.verbose) {
        std::cerr << " ManCurve:  Computing phi1 (intial guess phi1=phi0+eps*v0)\n";
    }
    // move along eigenvector until p-step aligns within alpha_max with the eigenvector
    VecType evec = (fwd ? saddle.evec[1] : saddle.evec[0]);
    if (params.verbose)
        std::cerr << "        Perturbation: type = " << (fwd ? "unstable" : "stable")
                  << ", period = " << thePeriod << "\n";
    if (params.verbose)
        std::cerr << "                    :  evec = " << evec << ",  step = " << eps
                  << " with (" << (p_or_m ? "+":"-") << ")\n";
    VecType p = saddle.pos;
    bool already_aligned = false;
    int tempCount = 0;
    while (true) {
        p += eps * (p_or_m ? 1 : -1) * evec;
        if (params.verbose)
            std::cerr << "  Step " << tempCount
                      << " (eps = " << (tempCount+1)*eps << "):  p = " << p;
        if ((tempCount+1)*eps > epsMax) {
            break;    //Reached max eps
        }
        VecType map_of_p;
        try {
            map_of_p = mapUsingCache(p,theMap,thePeriod,params,cache,eFail);
        } catch(MapDiscont& ef) {
            //Likely this won't happen unless fixed point is close to singularity
            std::cout << " STHMandifold : Encountered mapping error in initial perturbation\n";
            //Move to next step
            tempCount++;
            continue;
        }
        if (params.verbose) {
            std::cerr << " map(p) = " << map_of_p << "\n";
        }
        VecType dir = theMetric.displacement(p, map_of_p);
        // std::cerr << "     at [u1-u0] = " << dir << " norm = " << nvis::norm(dir)
        //     << " ( where delta_min = " << settings.delta_min << ") \n";
        // std::cerr << "     (test : u1-u0 = " << map_of_p-p
        //     << " and theMetric.displacement(u0,u1) = " << dir << ")\n";
        dir /= nvis::norm(dir);
        VecType p0u0 = theMetric.displacement(saddle.pos,p);
        //std::cerr << "     (test : u0-p0 = " << p-saddle.pos
        //     << " and theMetric.displacement(p0,u0) = " << p0u0 << ")\n";
        p0u0 /= nvis::norm(p0u0);
        double cosalpha = nvis::inner( p0u0 , dir );
        //std::cerr << "     at cosAngle = "
        // << cosalpha << " so alpha = " << acos(cosalpha)
        // << " (alpha_max = " << settings.alpha_max << ") \n";
        if (!already_aligned && cosalpha > cos(settings.alpha_max)) {
            already_aligned = true;
            break; //break if aligned
        } else if (already_aligned) {
            break;
        }
        tempCount++;
    }
    theManifold.manifold.push_back(p);
    length = theMetric.distance(theManifold.manifold[0], theManifold.manifold[1]);
    if (length > 2.0*epsMax) { //Max distance allowed (50x eps was in original code)
        theManifold.manifold[1] = epsMax * (p_or_m ? 1 : -1)*evec;
        progress.length = epsMax;
        std::cerr << " ManCurve:  Initial step larger than allowed first step (2*sdelta_max).\n";
        return;
    }
    i0 = 0;
    i1 = 1;
    
    //Total manifold members
    bool breakPointsExist = false;
    
    //Process the manifold
    int segmentCount = 0;
    //while (!stopper(theManifold.manifold.back()) ) {
    while (
        progress.length < settings.max_arc_length &&          //Total downstream arclength
        progress.seedingLength < settings.maxSeedArclength    //Seeding arclength
    ) {
        if ( settings.enableMaxSeg &&
                segmentCount >= settings.maxNumSeg ) {
            break;    //Max Number of segments
        }
        
        //Gather the seeding segment
        const VecType& p0 = theManifold.manifold[i0];
        const VecType& p1 = theManifold.manifold[i1];
        int numPts = (int) theManifold.manifold.size();
        double seedLength = progress.seedingLength;
        double curLength = progress.length;
        int currentNumBreakPts = (int) theManifold.breakIDs.size();
        
        //Call function to subdivide the seeding segment and generate new points
        manifold_type          downstreamPoints;
        std::list<int>         downstreamSepPtIdx;
        std::list<MapDiscont>  sepFailureList;
        adaptiveManfioldSubdivision<MAP,PARAM,SortableData,MapDiscont>
        (
            p0,p1,theMap,thePeriod,params,
            progress,settings,cache,
            downstreamPoints,downstreamSepPtIdx,sepFailureList
        );
        
        //Observe output
        segmentCount++;
        progress.segment += 1;
        std::cout << " In ManCurve: Step " << segmentCount << " computed :\n";
        std::cout << "   i0 = " << i0 << " i1 = " << i1 << "\n";
        std::cout << "   phi0 = " << p0 << " phi1 = " << p1 << "\n";
        std::cout << "   current SeedLength = " << seedLength << "\n";
        std::cout << "   current ArcLength = " << curLength << "\n";
        std::cout << "   current Num Sep Pts = " << currentNumBreakPts << "\n";
        std::cout << "   Computed " << (int) downstreamPoints.size() << " downstream points\n";
        std::cout << "   Found " << (int) downstreamSepPtIdx.size() << " downstream breaks\n";
        std::cout << "   Found " << (int) sepFailureList.size() << " New Sep Conditions\n";
        std::cout << "   new SeedLength = " << progress.seedingLength << "\n";
        std::cout << "   new ArcLength = " << progress.length << "\n";
        
        //Move temp data into manifold data structure
        typename manifold_type::iterator it;
        for(it=downstreamPoints.begin(); it!=downstreamPoints.end(); it++) {
            theManifold.manifold.push_back( (*it) );
        }
        std::list<int>::iterator idxIT;
        breakPointsExist = ((currentNumBreakPts+(int)downstreamSepPtIdx.size())>0);
        for(idxIT=downstreamSepPtIdx.begin(); idxIT!=downstreamSepPtIdx.end(); ++idxIT) {
            theManifold.breakIDs.insert( (*idxIT) + numPts );
        }
        typename std::list<MapDiscont>::iterator lit;
        for(lit=sepFailureList.begin(); lit!=sepFailureList.end(); lit++) {
            theManifold.discontList.push_back( *lit );
        }
        theManifold.lastSeed = i1;
        std::cout << "   Total Points = " << theManifold.manifold.size()
                  << " SepPts = " << theManifold.discontList.size()
                  << " BreakPts = " << theManifold.breakIDs.size() << "\n";
                  
                  
        /// Breaking at other conditions
        //If last point is nearing the homoclinic connection to another point in the chain, exit
        if(period > 1) { //Only run if more than p=1
            if(stopper.homoclinicTest(theManifold.manifold.back())) {
                break;
            }
        }
        
        //Increment to next seeding segment
        i0++;
        i1++;
        //Keep incrementing until past separation points
        if (breakPointsExist) {
            bool i1FoundInBreakIDs = true;
            i0--;
            i1--; //Reset for loop
            while( i1FoundInBreakIDs ) {
                //Increment at start of loop
                i0++;
                i1++;
                //BreakIDs are first index beyond a sep
                std::set<int>::iterator bIT;
                //bIT = theManifold.breakIDs.find( i0 );
                //i0Found = (bIT!=theManifold.breakIDs.end())? true : false;
                bIT = theManifold.breakIDs.find( i1 );
                i1FoundInBreakIDs = (bIT!=theManifold.breakIDs.end())? true : false;
            }
        }
        
    }
    
}


} //end topology


#endif //STH_MANIFOLD_HPP
