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


/* Base class of trajectory-type objects
 * Author: Wayne Schlei
 */
#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
//CR3BP objects
//#include <State.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/intersection.hpp>

using namespace xavier;

namespace pmateDesign {

template<typename STATE>
class Trajectory {
public :
    typedef STATE                 value_type;
    typedef Trajectory<STATE>     SelfType;
    
    /// Construction of a blank Trajectory object
    Trajectory() :
        initialized(false), valid(true), ic(0.0), epoch(0.0), deltaT(0.0)
    {}
    
    /// Construction of Trajectory object given input states
    Trajectory(std::vector<STATE>& in, std::vector<double>& t) :
        initialized(false), valid(true), ic(0.0), epoch(0.0), deltaT(0.0)
    {
        if ((int) in.size() == 0 ) {
            return;
        }
        initialized = true;
        ic = in[0];
        epoch = t[0];
        states = in;
        times = t;
        deltaT = t.back() - t[0];
    }
    
    Trajectory(const std::vector<STATE>& in, const std::vector<double>& t) :
        initialized(false), valid(true), ic(0.0), epoch(0.0), deltaT(0.0)
    {
        if ((int) in.size() == 0 ) {
            return;
        }
        initialized = true;
        ic = in[0];
        epoch = t[0];
        setStates(in,t);
        deltaT = t.back() - t[0];
    }
    
    /// Copy Constructor
    Trajectory(const SelfType& other) :
        initialized(other.initialized), valid(other.valid),
        ic(other.ic), epoch(other.epoch), deltaT(other.deltaT)
    {
        if (other.getNumStates() == 0) {
            return;
        }
        setStates( other.states, other.times );
    }
    
    /// Construction of a Trajectory object given an initial state, time of flight, and map engine???
    //template <class PMAP>
    //Trajectory(PMAP& theMap, const STATE& in, const double t0, const double tf);
    
    //Destructor
    //~Trajectory();
    
    /// Get the initial condition
    STATE getIC() const
    {
        return ic;
    }
    
    /// Setting the epoch time: sets all times accordingly
    void setEpoch(const double t0);
    /// Get the epoch time
    double getEpoch() const
    {
        return epoch;
    }
    /// Get the time at an index
    virtual double getTime(const int i) const;
    /// Get the time of flight for the current trajectory (result in whatever units are input)
    double getTimeOfFlight() const
    {
        return deltaT;
    }
    
    /// Set the initial condition (creates a need to repropagate)
    void setIC(const STATE& x0)
    {
        ic = x0;
    }
    /// Set states assuming t is already at the desired epochs
    void setStates(const std::vector<STATE>& in,const std::vector<double>& t,const bool reverseTime = false);
    /// Set states with t0 as the desired epoch
    void setStates(const double t0, const std::vector<STATE>& in, const std::vector<double>& t, const bool reverseTime=false);
    
    /// Get the number of states
    int getNumStates() const
    {
        return states.size();
    }
    /// Get a reference to the state vector
    std::vector<STATE>& getStates() const
    {
        return states;
    }
    /// Get a particular state at an index i
    void getState(const int& i, STATE& out) const;
    /// Get a state at a time using linear interpolation
    void getState(const double& t, STATE& out) const;
    ///Returns the index before t
    int getStateAndIndex(const double& t, STATE& out) const;
    
    /// Clear the current data
    void clear();
    /// Check if the states/times are initialized
    bool isInit() const
    {
        return initialized;
    };
    
    /// The states of the trajectory
    std::vector<STATE> states;
    /// The times of the states (Note: all time is stored in forward direction by convention)
    std::vector<double> times;
    
    /// A global epoch time (Useful for ephemeris and design)
    //static double globalEpoch;
    
    //---------------------------------------------------------------------------------------
    // Visualization functions
    //---------------------------------------------------------------------------------------
    /// Write this trajectory to an HxLineSet file for use with Avizo
    bool writeToHxLineSet(const char* filename) const;
    ///Get the minimum values
    //nvis::vec3 getMinValues();
    ///Get the maximum values
    //nvis::vec3 getMaxValues();
    ///Get 2D (Position) Bounding box
    //nvis::bbox2 getBoundingBox2D();
    ///Get 3D (Position) bounding box
    //nvis::bbox3 getBoundingBox();
    
    //---------------------------------------------------------------------------------------
    // Extraction functions:  (Comparing Trajectory objects)
    //---------------------------------------------------------------------------------------
    /// Find the intersection (in position space) between this trajectory and another
    void findIntersection(const SelfType& other,
                          std::vector<STATE>& xStates, std::vector<STATE>& otherXStates,
                          std::vector<double>& t, std::vector<double>& ot) const;
                          
    /// Find the Nearest (min distance in position) pair of nodes between this Trajectory and another
    void findNearestNode(const SelfType& other, STATE& tState, double& t,
                         STATE& oState, double& ot) const;
                         
    /// Find the Nearest node pair (min position distance) with parallel processing
    //void findNearestNodeParallel(const SelfType& other, STATE& tState, double& t, STATE& oState, double& ot) const;
    
    /** Find the Nearest (min distance in position) pair of points between
     *  this Trajectory and another using interpolation in time to find the closest points.
     */
    //void findNearestPoint(const SelfType& other, STATE& tState, double& t, STATE& oState, double& ot) const;
    
    /// Find the Nearest node pair (min position distance) with parallel processing
    //void findNearestNodeParallel(const SelfType& other, STATE& tState, double& t, STATE& oState, double& ot) const;
    //---------------------------------------------------------------------------------------
    
    /// Is this a valid trajectory (flag that is set by user)
    bool isValid() const
    {
        return valid;
    }
    /// Set the validity flag (user option)
    void setValidity(const bool v)
    {
        valid = v;
    }
    
private :
    //Orbit Info
    /// Does this object contain data or is the IC set.
    bool initialized;
    /// Validity (for ease of use)
    bool valid;
    /// The initial condition [typically for use with failed propagation]
    STATE ic;
    /// Epoch time (usually 0)
    double epoch;
    /// Time of flight from epoch
    double deltaT;
};

template<typename STATE>
void Trajectory<STATE>::clear()
{
    initialized = false;
    states.clear();
    times.clear();
}

/// Get Time at an index
template<typename STATE>
double Trajectory<STATE>::
getTime(const int i) const
{
    if (i<0 || i>=(int)times.size()) {
        std::cerr << "Access to wrong element of trajectory data\n";
        throw std::out_of_range("Trajectory::getTime() - Index out of range");
        return 0.0;
    }
    return times[i];
}

//This version assume epoch is already applied to the t states
template<typename STATE>
void Trajectory<STATE>::
setStates(const std::vector<STATE>& in,const std::vector<double>& t, const bool reverseTime)
{
    clear();
    initialized = true;
    if (reverseTime) {
        //Have to flip the direction
        int numStates = in.size();
        for (int i = numStates-1; i>=0; i--) {
            states.push_back( in[i] );
            times.push_back( t[i] );
        }
        ic = states[0];
    } else {
        ic = in[0];
        states = in;
        times = t;
    }
    epoch = t[0];
    deltaT = t.back() - t[0];
}

//This version will set all times to be relative to the provided epoch
template<typename STATE>
void Trajectory<STATE>::
setStates(const double t0, const std::vector<STATE>& in, const std::vector<double>& t, const bool reverseTime)
{
    setStates(in,t,reverseTime);
    setEpoch(t0);
}

//Sets all the times to follow new t0
template<typename STATE>
void Trajectory<STATE>::setEpoch(const double t0)
{
    double currentEpoch = epoch;
    epoch = t0;
    std::vector<double>::iterator it;
//Reset all the times
    for(it=times.begin(); it!=times.end(); it++) {
        (*it) += epoch - currentEpoch;
    }
}

///Get state based on index
template<typename STATE>
void Trajectory<STATE>::getState(const int& i, STATE& out) const
{
    if (i<0 || i>=(int)states.size()) {
        std::cerr << "Access to wrong element of trajectory data\n";
        throw std::out_of_range("Trajectory::getState() - Index out of range");
        return;
    }
    out = states[i];
}

///Get a state at a time based on linear interpolation in time
template<typename STATE>
void Trajectory<STATE>::getState(const double& t, STATE& out) const
{
    // Check to see if indicated time is within time range
    double eps = 1e-3;
    int index = 0;
    if ( (t < epoch - 1e-3) || (t > epoch + deltaT + 1e-3) ) {
        //Throw error if time is way off
        std::cerr << "Trajectory:  Failed access to states, time (" << t << ") is out of range\n";
        throw std::out_of_range("Trajectory::getState() - Time is out of range");
        return;
    } else if ( t <= epoch ) {
        out = states[0];
    } else if ( t >= times.back() ) {
        out = states.back();
    } else {
        //Loop through times to find index
        while(times[index] <= t) {
            index++;
        }
        index--;
        //Use the index for linear interpolation
        double u = (t - times[index]) / (times[index+1] - times[index]);
        out = (1.0-u)*states[index] + u*states[index+1];
    }
}


///Get state based on linear interpolation and return the data index just prior to 't'
template<typename STATE>
int Trajectory<STATE>::getStateAndIndex(const double& t, STATE& out) const
{
    //Get the trajectory state through interpolation
    getState(t,out);
    //Loop through times to get the index
    int index = 0;
    while(times[index] <= t) {
        index++;
    }
    //when found, backtrack index
    return index--;
}


/// Find the nearest-neighbor node (including time) with respect to another Trajectory object
template <typename STATE>
void Trajectory<STATE>::
findNearestNode(
    const typename Trajectory<STATE>::SelfType& other, //Second trajectory with which to find the closest point
    STATE& tState, double& t,                //Closest point and time on this trajectory
    STATE& oState, double& ot                //Closest point and time on the second trajectory
) const
{
    double minDist = 1000;
    int minIdx = 0, minOtherIdx = 0;
    for (int i=0; i<(int)times.size(); i++) {
        //Check the distance between other state and base state
        nvis::vec2 trajPos( states[i][0], states[i][1]);  // HARD CODED to Planar CR3BP!
        //TODO: template the extraction of position coordinates?
        for(int j=0; j<(int)other.times.size(); j++) {
            nvis::vec2 otherPos( other.states[j][0], other.states[j][1] ); //Position - planar
            double dist = nvis::norm( trajPos - otherPos );
            if(dist < minDist) {
                minDist = dist;
                minIdx = i;
                minOtherIdx = j;
            }
        }
    }
    //Store
    tState = states[minIdx];
    t = times[minIdx]; //+ globalEpoch;
    other.getState(minOtherIdx,oState);
    ot = other.getTime(minOtherIdx);
    
}

/// Find the intersections (in PLANAR position space) between this trajectory and another
template <typename STATE>
void Trajectory<STATE>::
findIntersection(
    const typename Trajectory<STATE>::SelfType& other,
    std::vector<STATE>& xStates,
    std::vector<STATE>& otherXStates,
    std::vector<double>& t,
    std::vector<double>& ot
) const
{
    //NOTE: HARD-CODED to 2D position vector!
    typedef  nvis::vec2                                          PosVec;
    typedef  nvis::fixed_vector< PosVec, 2 >                     PositionSegment;
    typedef  xavier::LineSegmentIntersection2D<PositionSegment>  IntersectionFinder;
    
    //Search for intesection on a Single Thread:
    int numOther = (int) other.states.size();
    int n = (int) states.size();
    //For every line segment on the other trajectory
    PosVec p0, p1, q0, q1;
    for(int j=0; j<(numOther-1); j++) {
        for(int k=0; k<2; k++) {
            q0[k] = other.states[j][k];
            q1[k] = other.states[j+1][k];
        }
        PositionSegment qSeg(q0,q1);
        //For every line segment on this
        for(int i=0; i<(n-1); i++) {
            for(int k=0; k<2; k++) {
                p0[k] = states[i][k];
                p1[k] = states[i+1][k];
            }
            PositionSegment pSeg(p0,p1);
            
            //Build the Intersection computation object
            IntersectionFinder lsx(pSeg,qSeg); //Computes in constructor
            
            //If found, store:
            if (lsx.getSolutionType() == IntersectionFinder::FOUND ||
                    lsx.getSolutionType() == IntersectionFinder::COLLINEAR ) {
                double tau = lsx.getIntersectionParameter(0); //On pSeg
                double nu = lsx.getIntersectionParameter(1); //On qSeg
                //Segment time bounds
                double dtp = times[i+1] - times[i];
                double dtq = other.times[j+1] - other.times[j];
                //Approximate states with linear interpolation
                STATE pState, qState;
                getState( times[i] + tau*dtp, pState );
                xStates.push_back( pState );
                t.push_back( times[i] + tau*dtp );
                other.getState( other.times[j] + nu*dtq, qState );
                otherXStates.push_back( qState );
                ot.push_back( other.times[j] + nu*dtq );
                
            } //Otherwise, no intersection. Don't store.
            
        }
    }
    
}



//Write to HxLineSet Object for Avizo
template <typename STATE>
bool Trajectory<STATE>::
writeToHxLineSet(const char* filename) const
{
    //Check sizes
    int numPts = (int) states.size();
    int numTimes = (int) times.size();
    if (numPts != numTimes) {
        return false;    //Shouldn't happen
    }
    
    //Open file
    std::ofstream theFile;
    theFile.open( filename );
    
    //Write Header:
    theFile << "# AmiraMesh ASCII 1.0\n";
    theFile << "# Generated by Trajectory.cpp\n\n";
    
    char buffer[200];
    int n=0;
    
    //Definitions
    n = sprintf(buffer, "define Lines %d\n",numPts+1);
    theFile << buffer;
    n = sprintf(buffer, "define Vertices %d\n\n",numPts);
    theFile << buffer;
    
    theFile << "Parameters {\n ";
    theFile << "    ContentType \"HxLineSet\"\n";
    theFile << "}\n\n";
    theFile << "Vertices { float[3] Coordinates } = @1\n";
    theFile << "Vertices { float Data1 } = @2\n";
    theFile << "Vertices { float Data2 } = @3\n";
    theFile << "Vertices { float Data3 } = @4\n";
    theFile << "Vertices { float Data4 } = @5\n";
    theFile << "Lines { int LineIdx } = @6\n\n";
    
    //NOTE: This is tempored to STATE being a poincare_map<>::state_type or xstate_type
    //      and assumes you have a [x,y,z,xdot,ydot,zdot] elements as the first
    //      6 elements of your vector.
    
    //x y z per point
    theFile << "@1\n";
    for (int i=0; i<numPts; i++) {
        theFile << states[i][0] << " " << states[i][1] << " " << states[i][2] << "\n";
    }
    
    //xd;yd;zd;t in separate lines
    theFile << "\n@2\n";
    for (int i=0; i<numPts; i++) {
        //theFile << states[i].xd << " ";
        theFile << states[i][3] << " ";
    }
    theFile << "\n";//End of data Line
    theFile << "\n@3\n";
    for (int i=0; i<numPts; i++) {
        //theFile << states[i].yd << " ";
        theFile << states[i][4] << " ";
    }
    theFile << "\n";//End of data Line
    theFile << "\n@4\n";
    for (int i=0; i<numPts; i++) {
        //theFile << states[i].zd << " ";
        theFile << states[i][5] << " ";
    }
    theFile << "\n";//End of data Line
    theFile << "\n@5\n";
    for (int i=0; i<numPts; i++) {
        theFile << times[i] << " ";
    }
    theFile << "\n";//End of data Line
    
    //LineIndex
    theFile << "\n@6\n";
    for (int i=0; i<numPts; i++) {
        theFile << i << " ";
    }
    theFile << "-1\n"; //close line
    
    //Close file
    theFile.close();
    
    return true;
}



} //End pmateDesign


#endif //TRAJECTORY_HPP
