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


/// Edge Rotation Exceptions:
/// Errors that occur during map displacement vector rotation tracking
// Author: Wayne Schlei

#ifndef __EDGE_ROTATION_FAILURE_HPP__
#define __EDGE_ROTATION_FAILURE_HPP__

#include <exception>
#include <string>


namespace xavier {

template <class VEC>
class EdgeRotationFailure : public std::exception {
    typedef EdgeRotationFailure<VEC> self_type;
public :
    /// Failure Types common in tracking "Map-Tangent Vector" along an edge
    enum FailureType {
        UNKNOWN = 0,                      //Unknown failure
        BACKWARD_MAP,                     //Poincare Map failed in backward time
        FORWARD_MAP,                      //Poincare Map failed in forward time
        BACKWARD_SINGULARITY,             //Intersecting an EOM singularity in backward time
        FORWARD_SINGULARITY,              //Intersecting an EOM singularity in forward time
        BACKWARD_SECTION_SEP,             //Trajectory separated from defined section in backward time
        FORWARD_SECTION_SEP,              //Trajectory separated from defined section in forward time
        BACKWARD_SECTION_SEP_NODE,        //Separated from defined section in backward time at Node
        FORWARD_SECTION_SEP_NODE,         //Separated from defined section in forward time at Node
        DOUBLE_PERIOD_FIXED_POINT,        //Double-period fixed point is interfering (store guess if applicable)
        FIXED_POINT_SUSPECTED,            //Fixed point of given period is on the edge
        STABLE_MANIFOLD_CROSS,            //Special unresolved case that shows up in non-separation cases
        UNRESOLVED                        //Unable to resolve rotation
    };
    
    ///Constructor
    EdgeRotationFailure() :
        type(UNKNOWN), failurePos( VEC(0) ), edgeID(-1), period(0), failIterate(0), theta(0.0), fixed(false), useMetaCell(false),
        mesg("Unknown Edge Rotation Failure")
    {}
    ///Constructor given FailureType, position, period
    EdgeRotationFailure(const FailureType& t, const VEC& pos, const int p) :
        type(t), failurePos(pos), edgeID(-1), period(p), failIterate(p), theta(0.0), fixed(false), useMetaCell(false)
    {
        setMessage(t);
    }
    ///Constructor given a name, position, period
    EdgeRotationFailure(const char* m, const VEC& pos, const int p) :
        type(UNKNOWN), failurePos(pos), edgeID(-1), period(p), failIterate(p), theta(0.0), fixed(false), useMetaCell(false),
        mesg(m)
    {}
    ///Destructor
    ~EdgeRotationFailure() throw() {}
    
    
    ///Set the message based on the failure type
    void setMessage(const FailureType& t)
    {
        switch(t) {
            case BACKWARD_MAP :
                mesg = std::string("Backward Map Failed");
                break;
            case FORWARD_MAP :
                mesg = std::string("Forward Map Failed");
                break;
            case BACKWARD_SINGULARITY :
                mesg = std::string("Backward Map Hit Singularity");
                break;
            case FORWARD_SINGULARITY :
                mesg = std::string("Forward Map Hit Singularity");
                break;
            case BACKWARD_SECTION_SEP :
                mesg = std::string("Separation from section in backward time");
                break;
            case FORWARD_SECTION_SEP :
                mesg = std::string("Separation from section in forward time");
                break;
            case BACKWARD_SECTION_SEP_NODE :
                mesg = std::string("NODE separated from section in backward time");
                break;
            case FORWARD_SECTION_SEP_NODE :
                mesg = std::string("NODE separated from section in forward time");
                break;
            case DOUBLE_PERIOD_FIXED_POINT :
                mesg = std::string("Interference from double period fixed point");
                useMetaCell = true;
                break;
            case FIXED_POINT_SUSPECTED :
                mesg = std::string("Fixed Point Suspected");
                break;
            case STABLE_MANIFOLD_CROSS :
                mesg = std::string("Suspected crossing of stable manifold");
                break;
            case UNRESOLVED :
                mesg = std::string("Unable to resolve separation condition");
                break;
            default :
                //Note:  This indicates an unknown reason for why the rotation tracking conditions are not met
                mesg = std::string("Unknown Edge Rotation Failure");
                break;
        }
    }
    
    ///Set a custom failure message
    void setMessage(const char* msg)
    {
        mesg = std::string( msg );
    }
    
    ///Set the error type
    void setType(const FailureType& t)
    {
        type = t;
        setMessage(type);
    }
    /// Set type based on unsigned int
    void setType(const unsigned int& i)
    {
        setType( getType(i) );
    }
    
    //Members:
    /// The type of edge rotation failure
    FailureType type;
    /// The location along an edge that the failure occurs
    VEC failurePos;
    /// Edge ID - used in conjunction with AdaptiveEdge
    nvis::ivec2 edgeID;
    /// The period associated with the edge failure
    int period;
    /// The achieved amount of iterates a failure occurs on
    int failIterate;
    /// The rotation angle value at error
    double theta;
    /// Flag (usually 'false' and not used)
    bool fixed;
    /// Suggestion of whether or not to use a meta cell to compensate for the failure
    bool useMetaCell;
    /// The error message
    std::string mesg;
    
    ///Return the section position of this failure
    VEC where() const
    {
        return failurePos;
    }
    
    ///Reset to an Unknown,unfixed error
    void reset()
    {
        setType( UNKNOWN );
        fixed = false;
    }
    
    
    ///Overwritten what() statement
    virtual const char* what() const throw()
    {
        return mesg.c_str();
    }
    
    /// Operator for less than (used in sorting)
    bool operator<(const self_type& rhs) const
    {
        //First by period
        if( this->period < rhs.period ) {
            return true;
        }
        //Next by position
        nvis::lexicographical_order compare;
        return compare(this->failurePos,rhs.failurePos);
        //Finally by type
        //if( this->type < rhs.type ) return true;
        return false;
    }
    
    /// Operator for "equals to" (used in sorting/uniquify)
    bool operator==(const self_type& rhs) const
    {
        if( this->period != rhs.period ) {
            return false;
        }
        //if( this->type != rhs.type ) return false;
        //Position
        return nvis::all(this->failurePos == rhs.failurePos);
    }
    
    /// Get an unsigned int representing the failure type
    unsigned int getTypeNumber(const FailureType& t) const
    {
        unsigned int n;
        switch(t) {
            case BACKWARD_MAP :
                n=1;
                break;
            case FORWARD_MAP :
                n=2;
                break;
            case BACKWARD_SINGULARITY :
                n=3;
                break;
            case FORWARD_SINGULARITY :
                n=4;
                break;
            case BACKWARD_SECTION_SEP :
                n=5;
                break;
            case FORWARD_SECTION_SEP :
                n=6;
                break;
            case BACKWARD_SECTION_SEP_NODE :
                n=7;
                break;
            case FORWARD_SECTION_SEP_NODE :
                n=8;
                break;
            case DOUBLE_PERIOD_FIXED_POINT :
                n=9;
                break;
            case FIXED_POINT_SUSPECTED :
                n=10;
                break;
            case STABLE_MANIFOLD_CROSS :
                n=11;
                break;
            case UNRESOLVED :
                n=12;
                break;
            default :
                n=0;
                break;
        }
        return n;
    }
    /// Get type given number
    FailureType getType(const unsigned int& i) const
    {
        FailureType t;
        switch(i) {
            case 1 :
                t = BACKWARD_MAP;
                break;
            case 2 :
                t = FORWARD_MAP;
                break;
            case 3 :
                t = BACKWARD_SINGULARITY;
                break;
            case 4 :
                t = FORWARD_SINGULARITY;
                break;
            case 5 :
                t = BACKWARD_SECTION_SEP;
                break;
            case 6 :
                t = FORWARD_SECTION_SEP;
                break;
            case 7 :
                t = BACKWARD_SECTION_SEP_NODE;
                break;
            case 8 :
                t = FORWARD_SECTION_SEP_NODE;
                break;
            case 9 :
                t = DOUBLE_PERIOD_FIXED_POINT;
                break;
            case 10 :
                t = FIXED_POINT_SUSPECTED;
                break;
            case 11 :
                t = STABLE_MANIFOLD_CROSS;
                break;
            case 12 :
                t = UNRESOLVED;
                break;
            default :
                t = UNKNOWN;
                break;
        }
        return t;
    }
    
};

} //end xavier


#endif
