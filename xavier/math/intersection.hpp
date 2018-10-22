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


#ifndef INTERSECTION_HPP
#define INTERSECTION_HPP

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>


namespace xavier {

/** \brief Class that computes intersection between 2D line segments
 * A Class for computing the intersections between two 2D Segments
 * Wayne Schlei - Purdue University
 *
 * The templates "Segment" and "Segment2" can be different, BUT
 * both types require that the end points can be accessed with
 * the operator [] where seg[0] and seg[1] are the first and
 * second points respectively.
 * It is also assumed that the "VectorType" of each Segment object
 * matches.  This class shouldn't compile unless they do match.
 *
 * Note: Up to user to adjust input segmetns to match space
 * metric (displacement via periodic space).
 */
template <typename Segment, typename Segment2 = Segment>
class LineSegmentIntersection2D {
public :
    typedef typename Segment::value_type          VecType;
    typedef nvis::bounding_box<VecType>           BoundsType;
    
    /// Constructor which requires the two segments
    LineSegmentIntersection2D(const Segment& s0, const Segment2& s1) :
        eps(1e-15),
        seg0(s0), tau(-1.0), t0(-1.0), t1(-1.0),
        seg1(s1), nu(-1.0), solType(NONE)
    {
        if (isIntersectionPossible()) {
            compute();
            //Otherwise, there is no intersection.
        } else {
            solType = NONE;
        }
    }
    
    /// Compute a bbox from end points (i=0 is seg0, i=1 is seg1)
    BoundsType getSegmentBBox(const int i) const
    {
        //Make sure i is 0 or 1
        int ii = i;
        if (ii>1) {
            ii = 1;
        }
        if (ii<0) {
            ii = 0;
        }
        
        
        
        //Switch based on which point
        VecType x0(0), x1(0);
        VecType xMin(0), xMax(0);
        if (ii==0) {
            //Set up end points
            x0 = seg0[0];
            x1 = seg0[1];
            xMin = x0;
            xMax = x1;
        } else {
            //Set up end points
            x0 = seg1[0];
            x1 = seg1[1];
            xMin = x0;
            xMax = x1;
        }
        //Find the min/max values
        for (int j=0; j<(int)x0.size(); j++) {
            if(x1[j] < xMin[j]) {
                xMin[j] = x1[j];
            }
            if(x0[j] > xMax[j]) {
                xMax[j] = x0[j];
            }
        }
        return BoundsType(xMin,xMax);
    }
    
    /// Function to see if we have to run the intersection protocol (useful to expedite processing)
    bool isIntersectionPossible() const
    {
        BoundsType bbox0 = getSegmentBBox(0);
        BoundsType bbox1 = getSegmentBBox(1);
        
        //Evaluate if centers are in other bboxes
        if (bbox0.inside(bbox1.center()) || bbox1.inside(bbox0.center()) ) {
            return true;
        }
        //Evaluate if any end points are within the bounding boxes
        if (bbox0.inside(seg1[0]) || bbox0.inside(seg1[1]) ||
                bbox1.inside(seg0[0]) || bbox1.inside(seg0[1]) ) {
            return true;
        }
        //If either segment is segment is vertical or horizontal,
        // have to test the ranges and modified centers for overlap.
        // Note: There may be a more elegant way to do this...
        VecType range0 = bbox0.size();
        VecType range1 = bbox1.size();
        VecType c0 = bbox0.center();
        VecType c1 = bbox1.center();
        //If vertical lines:   HARD-CODED AS 2D!!!!
        if ((std::fabs(range0[0]) <= eps) && (std::fabs(range1[0]) > eps)) { // Seg0 is vertical
            VecType testPt( c0[0], c1[1] );
            if ( bbox1.inside( testPt ) ) {
                return true;
            }
        }
        if ((std::fabs(range0[0]) > eps) && (std::fabs(range1[0]) <= eps)) { // Seg1 is vertical
            VecType testPt( c1[0], c0[1] );
            if ( bbox0.inside( testPt ) ) {
                return true;
            }
        }
        //If horizontal lines:    HARD-CODED AS 2D!!!!
        if ((std::fabs(range0[1]) <= eps) && (std::fabs(range1[1]) > eps)) { // Seg0 is horizontal
            VecType testPt( c1[0], c0[1] );
            if ( bbox1.inside( testPt ) ) {
                return true;
            }
        }
        if ((std::fabs(range0[1]) > eps) && (std::fabs(range1[1]) <= eps)) { // Seg1 is horizontal
            VecType testPt( c0[0], c1[1] );
            if ( bbox0.inside( testPt ) ) {
                return true;
            }
        }
        //Otherwise, these don't overlap and there are no intersections:
        return false;
    }
    
    /// Set the tolerance
    void setTolerance(const double e)
    {
        eps = e;
    }
    /// Get the tolerance
    double getTolerance()
    {
        return eps;
    }
    
    /// Enumerator indicating solution types
    enum SolutionType {
        NONE,                //No intersection exists
        PARALLEL,            //Segments are parallel (no intersection)
        COLLINEAR_DISJOINT,  //Segments are collinear but don't overlap
        COLLINEAR,           //Segments are collinear and DO overlap
        FOUND                //Intersection found
    };
    
    /// Get the line segment interesection solution type
    SolutionType getSolutionType() const
    {
        return solType;
    }
    
    /** Get the linear parameter indicating the crossing
     *  Note:  There is no check on whether or not a solution is generated.
     */
    double getIntersectionParameter(const int i=0) const
    {
        if (i>0) {
            return nu;
        } else {
            return tau;
        }
    }
    
    /// Get the intersection based on Segment::VecType
    VecType getIntersection() const
    {
        VecType sol(0);
        VecType r = seg0[1] - seg0[0];
        if (solType == FOUND) {
            sol = seg0[0] + tau*r;
        } else if (solType == COLLINEAR) {
            if (t0 >= 0.0 && t0 <= 1.0) {
                sol = seg0[0] + t0*r;
            } else {
                sol = seg0[0] + t1*r;
            }
        }
        return sol;
    }
    
    /// Compute the interesection
    void compute();
    
    /// Retrieve the angle of separation between the two segments
    double getAngle() const;
    
    /// Return a string indicating the solution type
    std::string getSolutionTypeString() const;
    
    /// Compute the distance between parallel lines
    double getParallelDistance() const;
    
    /// Compute distance from point to segment(id) with the boolean representing if the point is on the inside of the segment
    std::pair<bool,double> getPointToSegmentDistance(const VecType& pt, const int segID = 0) const;
    
    /// Compute the min distance between these two segments (bool indicates overlap)
    std::pair<bool,double> getMinDistance() const;
    
private :
    /// Tolerance on zero evaluation
    double eps;
    
    /// The first segment
    const Segment& seg0;
    /// The first segment's linear parameter indicating the intersection
    double tau, t0, t1;
    /// The second segment
    const Segment2& seg1;
    /// The second segment's linear parameter indicating the intersection
    double nu;
    /// The solution type information
    SolutionType solType;
    
};

/// Computing an intersection between line segments
template<typename Segment, typename Segment2>
void LineSegmentIntersection2D<Segment,Segment2>::
compute()
{
    //Construct vectors for testing
    VecType p = seg0[0];
    VecType r = seg0[1] - seg0[0];
    VecType q = seg1[0];
    VecType s = seg1[1] - seg1[0];
    VecType qmp = q - p;
    
    // Eval assuming we have 2D vectors : FIX in future versions!!!
    double r_x_s = r[0]*s[1] - r[1]*s[0];
    double qmp_x_r = qmp[0]*r[1] - qmp[1]*r[0];
    double qmp_x_s = qmp[0]*s[1] - qmp[1]*s[0];
    double qmp_d_r = nvis::inner( qmp, r);
    double s_d_r = nvis::inner( s, r);
    double r_d_r = nvis::inner( r, r);
    
    // Solution types :
    if (fabs(r_x_s) < eps && fabs(qmp_x_s) < eps ) {
        // Collinear
        solType = COLLINEAR_DISJOINT;
        t0 = qmp_d_r / r_d_r;
        t1 = t0 + s_d_r / r_d_r;
        //Check if it is within the interval to create an overlap
        if ( (t0>=0.0 && t0 <= 1.0) || (t1>=0.0 && t1 <= 1.0)) {
            solType = COLLINEAR;
            //Have to output something:
            if (t0>=0.0 && t0 <= 1.0) {
                tau = t0;
                nu = 0.0;
            } else {
                tau = 0.0;
                nu = t1;
            }
        }
    } else if (fabs(r_x_s) < eps) {
        // Segments are parallel, and thus, have no intersection.
        solType = PARALLEL;
    } else {
        //Evaluate tau and nu now that r_x_s is non-zero
        tau = qmp_x_s / r_x_s;
        nu =  qmp_x_r / r_x_s;
        // Solution exists if both values are on linear range.
        if ( (tau>=0.0 && tau <=1.0) && (nu>=0.0 && nu<=1.0) ) {
            solType = FOUND;
        } else {
            //No intersection exists on given segments.
            solType = NONE;
        }
    }//Parallel option
    
}

/// Retrieve the angle of separation between the two segments
template<typename Segment, typename Segment2>
double LineSegmentIntersection2D<Segment,Segment2>::
getAngle() const
{
    VecType r = seg0[1] - seg0[0];
    VecType s = seg1[1] - seg1[0];
    double rMag = nvis::norm( r );
    double sMag = nvis::norm( s );
    double r_d_s = nvis::inner( r, s);
    return acos( r_d_s / (rMag*sMag) );
}

/// Get a string with solution type
template<typename Segment, typename Segment2>
std::string LineSegmentIntersection2D<Segment,Segment2>::
getSolutionTypeString() const
{
    std::string output;
    switch(solType) {
        case NONE :
            output = "NONE";
            break;
        case PARALLEL :
            output = "PARALLEL";
            break;
        case COLLINEAR_DISJOINT :
            output = "COLLINEAR_DISJOINT";
            break;
        case COLLINEAR :
            output = "COLLINEAR";
            break;
        case FOUND :
            output = "FOUND";
            break;
    }
    return output;
}

/// Compute the distance (perpendicular separation) between two parallel segments
template<typename Segment, typename Segment2>
double LineSegmentIntersection2D<Segment,Segment2>::
getParallelDistance() const
{
    //Warn user that you are calling this incorrectly
    if(solType != PARALLEL) {
        std::cerr << "LineSegmentIntersection2D:  Inappropriate call to getParallelDistance().  [Not PARALLEL!]\n";
        return 1000.0;
    }
    
    //Test x to see if lines are too vertical
    double dx0 = seg0[1][0] - seg0[0][0];
    double dx1 = seg1[1][0] - seg1[0][0];
    //Distance parameters
    double a=0.0,b=0.0,c0=0.0,c1=0.0,d=1000.0;
    if (dx0 <= 1.e-12) {
        //Compute the distance as if the lines were vertical
        double m = dx0 / (seg0[1][1] - seg0[0][1]);
        c0 = seg0[0][0] - m*seg0[0][1];
        c1 = seg1[0][0] - m*seg1[0][1];
        a = -1.0;
        b = m;
        d = abs(c1 - c0) / sqrt(a*a + b*b);
    } else {
        //Compute the distance
        double m = (seg0[1][1] - seg0[0][1]) / dx0; //slope
        c0 = seg0[0][1] - m*seg0[0][0];
        c1 = seg1[0][1] - m*seg1[0][0]; //Same slope
        a = m;
        b = -1.0;
        d = abs(c1 - c0) / sqrt(a*a + b*b);
    }
    return d;
}

/// Get the distance between a point and one of the two segments.
template<typename Segment, typename Segment2>
std::pair<bool,double> LineSegmentIntersection2D<Segment,Segment2>::
getPointToSegmentDistance(
    const typename LineSegmentIntersection2D<Segment,Segment2>::VecType& pt,
    const int segmentID
) const
{
    std::pair<bool,double> output(false,1000.0);
    VecType v, w, u, pb;
    if (segmentID == 0) {
        v = seg0[1] - seg0[0];
        w = pt - seg0[0];
        u = pt - seg0[1];
        pb = seg0[0];
    } else {
        v = seg1[1] - seg1[0];
        w = pt - seg1[0];
        u = pt - seg1[1];
        pb = seg1[0];
    }
    
    //Evaluate dot products
    double wdotv = nvis::inner(w,v);
    double vdotv = nvis::inner(v,v);
    if (wdotv <= 0.0 ) {
        //outside so use vector distance
        output.first = false;
        if (wdotv == 0.0) {
            output.first = true;
        }
        output.second = nvis::norm(w); //before
    } else if (vdotv <= wdotv) {
        //outside so use vector distance
        output.first = false;
        if (vdotv == wdotv) {
            output.first = true;
        }
        output.second = nvis::norm(u); //after
    } else {
        pb +=  (wdotv / vdotv) * v;
        output.first = true;
        output.second = nvis::norm( pb - pt );
    }
    //Returning output pair
    return output;
}

/// Get the minimum distance between two segments
template<typename Segment, typename Segment2>
std::pair<bool,double> LineSegmentIntersection2D<Segment,Segment2>::
getMinDistance() const
{

    //If collinear
    if(solType == COLLINEAR) {
        return std::pair<bool,double>(true,0.0);
    }
    
    //Otherwise, we run the point to segment tests
    double minDist = 1000.0;
    std::vector<std::pair<bool,double> > ptTests;
    ptTests.push_back( getPointToSegmentDistance(seg0[0],1) );
    ptTests.push_back( getPointToSegmentDistance(seg0[1],1) );
    ptTests.push_back( getPointToSegmentDistance(seg1[0],0) );
    ptTests.push_back( getPointToSegmentDistance(seg1[1],0) );
    
    //Overlap check
    bool overlap = false;
    for(int i=0; i<4; i++) {
        //Minimum distance
        if (ptTests[i].second <= minDist) {
            minDist = ptTests[i].second;
        }
        //Overlap
        if (ptTests[i].first) {
            overlap = true;
        }
    }
    //Deliver output (overlap, minDist)
    return std::pair<bool,double>(overlap,minDist);
}


} //End xavier

#endif
