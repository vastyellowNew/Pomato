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


#ifndef MAP_MANIFOLD_SEGMENT_HPP
#define MAP_MANIFOLD_SEGMENT_HPP

#include <cstdio>
#include <cassert>
#include <fstream>
#include <vector>
#include <set>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/intersection.hpp>
#include <data/edge.hpp>
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <orbital/progenitor.hpp>
#include <topology/EdgeRotationFailure.hpp>

using namespace xavier;

namespace topology {

/** Container for Invariant manifold segments
 *
 */
template <typename MAPSTATE, typename FULLSTATE, typename DATA>
class MapManifoldSegment : public nvis::fixed_vector<MAPSTATE,2> {
public:
    typedef MAPSTATE                                     VecType;
    typedef FULLSTATE                                    StateType;
    typedef nvis::fixed_vector<MAPSTATE,2>               BaseType;
    typedef nvis::bounding_box<VecType>                  BoundsType;
    typedef MapManifoldSegment<MAPSTATE,FULLSTATE,DATA>  SelfType;
    typedef EdgeRotationFailure<MAPSTATE>                MapDiscont;
    typedef orbital::ProgenitorState<StateType,VecType>  ProgState;
    typedef LineSegmentIntersection2D<SelfType>          IntersectionFinder;
    
    ///Constructor for first segment call
    MapManifoldSegment(const MAPSTATE& x0, const MAPSTATE& x1, const DATA& d1) :
        BaseType(x0,x1),
        previous(-1),
        next(-1),
        segID(0),
        parent(-1),
        tau0(0.0),
        tau1(1.0),
        data0(d1),  //Note:  This is copied to both the first steps and should be found from orbit.
        data1(d1),
        ps0(), ps1(), //Assume blank, we don't know them.
        transverseSection(true)
    {}
    
    ///Constructor for generic new segment
    MapManifoldSegment( const int idx,        //Segment Index (in vector container)
                        const MAPSTATE& x0, const MAPSTATE& x1, //End points
                        const double& t0, const double& t1,     //Linear parameters (tau) at end nodes
                        const DATA& d0, const DATA& d1,         //Data from CATtracker at end nodes
                        const SelfType& source, const SelfType& streamPrev) :
        BaseType(x0,x1),
        previous(streamPrev.segID),
        next(-1),
        segID(idx),
        parent(source.segID),
        tau0(t0),
        tau1(t1),
        data0(d0),
        data1(d1),
        ps0(), ps1(), //Don't know ProgenitorStates on generic case.
        transverseSection(true)
    {}
    
    /// Copy constructor
    MapManifoldSegment(const SelfType& other) :
        BaseType(other[0],other[1]),
        previous( other.previous ),
        next( other.next ),
        segID( other.segID ),
        parent( other.parent ),
        children( other.children ),
        tau0( other.tau0 ),
        tau1( other.tau1 ),
        data0( other.data0 ),
        data1( other.data1 ),
        ps0( other.ps0 ),
        ps1( other.ps1 ),
        transverseSection( other.transverseSection ),
        discontinuityList( other.discontinuityList )
    {}
    
    /// Read constructor
    MapManifoldSegment(FILE* f) :
        BaseType(MAPSTATE(0),MAPSTATE(0)),
        previous(-1),
        next(-1),
        segID(-1),
        parent(-1),
        tau0(0.0),
        tau1(1.0),
        data0(),
        data1(),
        ps0(),
        ps1(),
        transverseSection(true)
    {
        read(f);
    }
    
    
    //Connection in Advection Process ---------------------------------------------
    /// Int for position Previous segment on the Poincare section
    int previous;
    /// Int for position of the Next segment
    int next;
    /// Get reference to next segment
    SelfType& getNext(std::vector<SelfType>& segs)
    {
        assert( (next>=0 && next<(int)segs.size()) );
        return segs[next];
    }
    /// Get reference to previous segment
    SelfType& getPrevious(std::vector<SelfType>& segs)
    {
        assert( (previous>=0 && previous<(int)segs.size()) );
        return segs[previous];
    }
    
    //Interaction Information: ----------------------------------------------------
    /// Current index counter for easy access
    int segID;
    /// Parent (source) int with -1 representing the starting object
    int parent;
    /// Does this segment have a parent segment
    bool hasParent() const
    {
        return (parent<0)? false : true;
    }
    /// Container for children (downstream segments)
    std::vector<int> children;
    /// Adding a child
    void addChild(int kidID)
    {
        children.push_back(kidID);
    }
    /// Does this segment have child segments
    bool hasKids() const
    {
        return ((int)children.size()>0)? true : false;
    }
    /// Does an indicated parent tau track to this segment
    bool isOnSegment(const double& parentTau) const
    {
        return (parentTau >= tau0 && parentTau <= tau1);
    }
    
    /// Container for cousins (sub-p iterates) ??? - likely in cache??
    
    /// Segment length
    template <class PARAM>
    double length(const PARAM& params) const
    {
        MAPSTATE tmp = params.the_metric.displacement((*this)[1],(*this)[0]);
        return nvis::norm(tmp);
    }
    
    /// Linear parameters from parent segment that generates the end points
    double tau0, tau1;
    /** Additional Data about path from parent segment to node at given period
     *  which is required for various purposes (filtering and visualization)
     *  (e.g., CloseApproach distances and time from CATtracker)
     */
    DATA data0, data1;
    /// Progenitor State Data for each point of this segment
    ProgState ps0, ps1;
    /// Is the section transverse for this segment
    bool transverseSection;
    /// Map Discontinuities from downstream propagation if non-transverse segment
    std::list< MapDiscont > discontinuityList;
    
    ///Insert a child manifold segment
    void insertChild(SelfType* kid)
    {
        children.push_back(kid);
    }
    ///Remove all children
    void removeChildren()
    {
        children.clear();
    }
    ///Set Parent Pointer
    void setParent(SelfType* mom)
    {
        parent = mom->segID;
    }
    
    
    /// Write to file under a MapManifold object
    void write(FILE* f) const
    {
        //SegID transverseSection parent->segID numData numKids numSeps
        int numData = (int) data0.size(); //Assume it's a vector?
        int numKids = (int) children.size();
        int numSeps = (int) discontinuityList.size();
        fprintf(f,"%d %d %d %d %d %d\n", segID, (transverseSection)? 1 : 0, parent, numData, numKids, numSeps);
        //Write point 0 : tau0, data0
        fprintf(f,"%.15f %.15f %.15f ", (*this)[0][0], (*this)[0][1], tau0);
        for(int k=0; k<numData; k++) {
            fprintf(f,"%.15f ",data0[k]);
        }
        fprintf(f,"\n");
        //Write point 0 : Progenitor State 0
        fprintf(f,"%d %.15f %.15f %.15f %.15f ",
                (ps0.found)? 1 : 0, ps0.alpha, ps0.dv[0], ps0.dv[1], ps0.tof
               );
        for(int k=0; k<6; k++) {
            fprintf(f,"%.15f ", ps0.orbitState[k]);
        }
        fprintf(f,"\n");
        //Write point 1 : tau1, data1
        fprintf(f,"%.15f %.15f %.15f ", (*this)[1][0], (*this)[1][1], tau1);
        for(int k=0; k<numData; k++) {
            fprintf(f,"%.15f ",data1[k]);
        }
        fprintf(f,"\n");
        //Write point 1 : Progenitor State 1
        fprintf(f,"%d %.15f %.15f %.15f %.15f ",
                (ps1.found)? 1 : 0, ps1.alpha, ps1.dv[0], ps1.dv[1], ps1.tof
               );
        for(int k=0; k<6; k++) {
            fprintf(f,"%.15f ", ps1.orbitState[k]);
        }
        fprintf(f,"\n");
        
        //Previous,Next are easily reconstructed in container object
        
        //Discontinuities
        typename std::list<MapDiscont>::const_iterator lit = discontinuityList.begin();
        for(; lit!=discontinuityList.end(); ++lit) {
            unsigned int t = lit->getTypeNumber(lit->type);
            int p = lit->period;
            VecType pos = lit->failurePos;
            fprintf(f,"%d %d %.15f %.15f\n",t,p,pos[0],pos[1]);
        }
        
        //Children - segment ID's that spawn from this segment as working segment
        for(int k=0; k<numKids; k++) {
            fprintf(f,"%d ",children[k]);
        }
        if(numKids>0) {
            fprintf(f,"\n");
        }
    }
    
    /** Read a MapManifoldSegment from ManifoldData file:
     *  Note: have to set next and previous ids externally based on container
     */
    void read(FILE* f)
    {
        int tempBool = 0;
        int numData = 0, numSeps = 0, numKids = 0;
        //Base data
        fscanf(f, "%d %d %d %d %d %d", &segID, &tempBool, &parent, &numData, &numKids, &numSeps);
        transverseSection = (tempBool == 1);
        //Point 0 - HardCoded!
        DATA temp(numData,50.0);
        fscanf(f, "%lf %lf %lf %lf %lf %lf",
               &((*this)[0][0]), &((*this)[0][1]), &tau0,
               &(temp[0]),&(temp[1]),&(temp[2]));
        //for(int k=0;k<numData;k++) fscanf(f,"%lf",&(temp[k]));
        data0 = temp;
        //ProgenitorState 0 - HARD CODED!
        fscanf(f, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &tempBool, &(ps0.alpha), &(ps0.dv[0]), &(ps0.dv[1]), &(ps0.tof),
               &(ps0.orbitState[0]), &(ps0.orbitState[1]), &(ps0.orbitState[2]),
               &(ps0.orbitState[3]), &(ps0.orbitState[4]), &(ps0.orbitState[5])
              );
        ps0.found = (tempBool==1);
        //Point 1 - HardCoded!
        fscanf(f, "%lf %lf %lf %lf %lf %lf",
               &((*this)[1][0]), &((*this)[1][1]), &tau1,
               &(temp[0]),&(temp[1]),&(temp[2]));
        //for(int k=0;k<numData;k++) fscanf(f,"%lf",&(temp[k]));
        data1 = temp;
        //ProgenitorState 1 - HARD CODED!
        fscanf(f, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &tempBool, &(ps1.alpha), &(ps1.dv[0]), &(ps1.dv[1]), &(ps1.tof),
               &(ps1.orbitState[0]), &(ps1.orbitState[1]), &(ps1.orbitState[2]),
               &(ps1.orbitState[3]), &(ps1.orbitState[4]), &(ps1.orbitState[5])
              );
        ps1.found = (tempBool==1);
        
        //Section Discontinuities List
        for(int k=0; k<numSeps; k++) {
            int p=0, t=0;
            VecType x;
            fscanf(f,"%d %d %lf %lf",&t,&p,&(x[0]),&(x[1]));
            MapDiscont thisDiscont(MapDiscont::UNKNOWN,x,p);
            thisDiscont.setType((unsigned int) t);
            discontinuityList.push_back( thisDiscont );
        }
        
        //Children segIDs
        for(int i=0; i<numKids; i++) {
            int kid = 0;
            fscanf(f,"%d",&kid);
            children.push_back(kid);
        }
        
    }
    
    /// Interpolate data based on linear parameter on [0,1]
    double getData(const double& u, const int dataID) const
    {
        if (!(u>=0.0 && u<=1.0)) {
            throw std::out_of_range("MapManifoldSegment: Parameter out of range [0,1] in getData()");
        }
        if (!(dataID < (int)data0.size())) {
            throw std::invalid_argument("MapManifoldSegment: Data index is invalid in getData()");
        }
        return (1.0-u)*data0[dataID] + u*data1[dataID];
    }
    
    
    /// Get a bounding box for this segment
    //template <class PARAM>  //technically, this needs a template
    BoundsType getSegmentBounds() const
    {
        //TODO:  adjust to match space metric (displacement, periodic, etc.)
        //First check if points are even close to intersecting by testing bounding boxes
        VecType x0 = (*this)[0];
        VecType x1 = (*this)[1];
        VecType xMin(x0), xMax(x1);
        //Find the min/max values
        for (int i=0; i<(int)x0.size(); i++) {
            if(x1[i] < xMin[i]) {
                xMin[i] = x1[i];
            }
            if(x0[i] > xMax[i]) {
                xMax[i] = x0[i];
            }
        }
        
        //Make a bounding box with the two end points of this segment
        return BoundsType(xMin,xMax);
        
    }
    
    
    /// Get the linear parameter for a given point (assert checks if point is within bbox of segment)
    //template <class PARAM>  //technically, this needs a template
    double getLinearParam(const MAPSTATE& x) const
    {
        //Make sure inside bounds
        BoundsType bbox = getSegmentBounds();
        //Scale up by 50% to avoid borderline cases that are caused by single-precision picking
        bbox.scale(1.5);
        //assert( bbox.inside(x) ); //Temporary exit for debugging
        //if ( !bbox.inside(x) ) throw std::out_of_range("MapManifoldSegment: State is outside bbox in getLinearParam()");
        
        
        //Run linear parameter computation (TODO: Really need a template for metric.displacement here)
        MAPSTATE p0pt = x - (*this)[0];
        MAPSTATE p0p1 = (*this)[1] - (*this)[0];
        double u = nvis::norm(p0pt) / nvis::norm(p0p1);
        //Check if inside bbox
        if(!bbox.inside(x) ) {
            //Project the point to the line if this occurs
            u = nvis::inner( p0pt, p0p1 ) / nvis::norm(p0p1);
        }
        //Map to ends if conditions met:
        if ( u > 1.0 ) {
            u = 1.0;
        }
        if ( nvis::inner( p0pt, p0p1 ) < 0.0 ) {
            u = 0.0;
        }
        return u;
    }
    
    /// Get a linear parameter on [0,1] for this segment given an upstream linear parameter that should be on this segment
    double getLinearParam(const double& parentTau) const
    {
        //Test to make sure this is on this segment
        assert( isOnSegment(parentTau) );
        //Compute the linear parameter based on weighting between end nodes of segment
        double u = (parentTau - tau0) / (tau1 - tau0);
        if ( u > 1.0 ) {
            u = 1.0;
        }
        if ( u < 0.0 ) {
            u = 0.0;
        }
        return u;
    }
    
    /// Get Point given a linear parameter for a given segment
    MAPSTATE getPoint(const double& u) const
    {
        if (!(u>=0.0 && u<=1.0)) {
            throw std::out_of_range("MapManifoldSegment: Parameter out of range [0,1] in getPoint()");
        }
        return (1.0-u)*(*this)[0] + u*(*this)[1];
    }
    
    /// Get Point (double) given an approximate input state (may be slightly off line or in float  )
    MAPSTATE getPoint(const MAPSTATE& xApprox) const
    {
        double u = getLinearParam(xApprox);
        return (1.0-u)*(*this)[0] + u*(*this)[1];
    }
    
    
    /// Does this segment intersect another?
    bool intersect(const SelfType& other) const
    {
        //Create an intersection object to run the test:
        IntersectionFinder lsx((*this),other); //Computes in constructor
        //If segments intersect or overlap:
        if(lsx.getSolutionType() == IntersectionFinder::FOUND ||
                lsx.getSolutionType() == IntersectionFinder::COLLINEAR ) {
            return true;
        }
        
        //Otherwise No intersection
        return false;
    }
    
    /// Test between two (complementary) segments to determine if they are along a saddle-loop (e.g., KAM manifold)
    bool isSaddleLoop(const SelfType& other, const double eps = 1.e-6) const
    {
        //Create an intersection object to run the test:
        IntersectionFinder lsx((*this),other); //Computes in constructor
        
        //The first test is to see if the segments are pointing towards one another
        double theta = lsx.getAngle(); //Angle between vectors
        double alpha = M_PI - theta;
        //Condition indicates that they are not pointed towards each other:
        if ( fabs(alpha) >= 3.0*M_PI/180.0 ) {
            return false;
        }
        
        //If segments COLLINEAR AND overlap:
        if(lsx.getSolutionType() == IntersectionFinder::COLLINEAR ) {
            return true;
        }
        
        //Let's also check if their is an intersection with a very small angle
        //Note: we don't want this to be too big because it will trip a lot of erroneous stops.
        if(lsx.getSolutionType() == IntersectionFinder::FOUND ) {
            return true;
            //Note: this is highly an unlikely case since manifolds are the separatrix of the flow.
        }
        
        //Also, lines could be parallel but really close together.
        if(lsx.getSolutionType() == IntersectionFinder::PARALLEL ) {
            double distance = lsx.getParallelDistance();
            if (distance <= eps) {
                return true;
            }
        }
        
        //And lastly, also possible for no intersection but with overlap due to numerical roundoff
        if(lsx.getSolutionType() == IntersectionFinder::NONE ) {
            std::pair<bool,double> test = lsx.getMinDistance();
            //Check for overlap AND small map separation
            if(test.first && (test.second <= eps) ) {
                return true;
            }
        }
        
        //Otherwise, unable to determine if KAM manifold
        return false;
    }
    
    /// Get intersection conditions
    bool getIntersection(const SelfType& other, double& u, double& uOther) const
    {
        //Create an intersection object to run the computation:
        IntersectionFinder lsx((*this),other);
        
        //Always gather the linear parameters: (returns -1.0 if no intersection)
        //Evaluate the linear parameter for THIS segment
        u = lsx.getIntersectionParameter(0);
        //Evaluate the linear parameter for the OTHER segment
        uOther = lsx.getIntersectionParameter(1);
        
        //If segments overlap:
        if(lsx.getSolutionType() == IntersectionFinder::FOUND ||
                lsx.getSolutionType() == IntersectionFinder::COLLINEAR ) {
            return true;
        } else {
            //Otherwise No intersection
            return false;
        }
    }
    
    /// Is this segment vertical based on the angle tolerance
    bool isNearlyVertical(const double eps = 3.0*M_PI/180.0) const
    {
        VecType s = (*this)[1] - (*this)[0];
        VecType v(0);
        v[1] = 1.0;
        double theta = acos( nvis::inner(s,v) / nvis::norm(s) );
        double alpha = M_PI - theta;
        if (theta < eps || alpha < eps) {
            return true;
        }
        return false;
    }
    
};

} //end topology

#endif  // MAP_MANIFOLD_SEGMENT_HPP
