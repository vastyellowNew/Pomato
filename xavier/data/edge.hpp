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


#ifndef __EDGE_HPP__
#define __EDGE_HPP__

#include <exception>
#include <functional>
#include <set>
#include <vector>
#include <map>
#include <list>
#include <math/fixed_vector.hpp>
#include <topology/EdgeRotationFailure.hpp>

namespace xavier {
/** An Edge Class - Stores all periods that need to be computed for angle rotation
   Author: Xavier Tricoche & Wayne Schlei - Purdue University
   Date: 7/16/2013
   Keeps a list of possible periods.
   Inherits from fixed_vector, so access to endpoints is simply edge[0], edge[1]

   The type (T) is usually some form of vector like nvis::ivec2
 */
template<typename T, typename Compare = std::less<T> >
class edge : public nvis::fixed_vector<T, 2> {
    void swap(T& t0, T& t1)
    {
        T tmp = t1;
        t1 = t0;
        t0 = tmp;
    }
    
public:
    typedef nvis::fixed_vector<T, 2>  base_type;
    typedef edge<T, Compare>          self_type;
    typedef T                         index_type;
    
    edge() : nvis::fixed_vector<T, 2>() {}
    edge(const T& t0, const T& t1) : nvis::fixed_vector<T, 2>(t0, t1)
    {
        Compare Lt;
        if (Lt(t1,t0)) {
            swap((*this)[0], (*this)[1]);
        }
    }
    ///Create from "edge" type (Copy constructor just to be sure)
    edge(const self_type& e) :
        nvis::fixed_vector<T,2>(e[0],e[1])
    {
        periods = e.periods;
        rotationAngleMap = e.rotationAngleMap;
        mapError = e.mapError;
        //rotFailurePos = e.roFailurePos;
        rotFailures = e.rotFailures;
        transverseSection = e.transverseSection;
        xMinDeltaMap = e.xMinDeltaMap;
        minDeltaMap = e.minDeltaMap;
    }
    bool operator<(const self_type& e) const
    {
        Compare Lt;
        if (Lt((*this)[0], e[0])) {
            return true;
        } else if (Lt(e[0], (*this)[0])) {
            return false;
        }
        return Lt((*this)[1], e[1]);
    }
    
    /// Set of periods for angle rotation computation
    std::set<int> periods;
    /// Map between period and rotation Value
    std::map<int,double> rotationAngleMap;
    /// Error flags : flag per period.  Set to true if could NOT evaluate rotation correctly
    std::map<int,bool> mapError;
    /// Map Failure locations on this edge (period->position)
    //std::list< std::pair<int,nvis::vec2> > rotFailurePos;
    /// Edge Rotation Failure data on this edge (locations, periods, types, fixed or not)
    std::list< EdgeRotationFailure<nvis::vec2> > rotFailures;
    /// Map indicating if transverse section assumption is violated at a given period on this edge
    std::map<int,bool> transverseSection;
    /// Location of minimum delta along edge per period
    std::map<int,nvis::vec2> xMinDeltaMap;
    /// Minimum delta value per period
    std::map<int,nvis::vec2> minDeltaMap;
};


/** Comparison structure : compares a nodeID as a subID given preset maximum depth
    - Useful for AdaptiveGrid edge information.
    - Author:  Wayne Schlei
*/
struct OnEdgeCompare {
    ///Maximum depth for comparison parameters
    static int maxDepth;
    
    ///Comparison for an nvis::ivec3 vector (=Vec)
    template<typename Vec>
    bool operator() (const Vec& lhs, const Vec& rhs) const
    {
        //Convert to sub_id first at max depth
        int lhsD = lhs[2];
        int rhsD = rhs[2];
        int lhsShift = (int)pow((float)2,maxDepth-lhsD);
        int rhsShift = (int)pow((float)2,maxDepth-rhsD);
        Vec lhsMaxD(lhsShift*lhs[0],lhsShift*lhs[1],maxDepth);
        Vec rhsMaxD(rhsShift*rhs[0],rhsShift*rhs[1],maxDepth);
        //Run comparison with lexicographical_order
        nvis::lexicographical_order baseCompare;
        bool result = baseCompare(lhsMaxD,rhsMaxD);
        return result;
    }
};

/** An AdaptiveGridTerminalEdge class representing a "bottom-level"
 *  edge within the adaptive grid that cannot be subdivided
 *  Future:  try to replace pow() calls with a bit shift
 *  Author: Wayne Schlei (Purdue University)
 */
template<typename T>
class TerminalEdge : public edge<T,OnEdgeCompare> {
public :
    typedef edge<T,OnEdgeCompare> edge_type;
    typedef TerminalEdge<T>     self_type;
    
    
    ///Constructors
    TerminalEdge() : edge_type() {}
    TerminalEdge(const T& t0, const T& t1) : edge_type(t0,t1) {}
    TerminalEdge(const edge_type& e) : edge_type(e) {}
    
    /// Return the unit difference vector (determines vertical or horizontal)
    T deltaID()
    {
        int depth0 = (*this)[0][2];
        int depth1 = (*this)[1][2];
        T delta;
        if (depth0 == depth1) {
            //Same depth
            delta = (*this)[1] - (*this)[0];
        } else if( depth0 < depth1 ) {
            //Convert first element to depth1 and take difference
            int shift = (int)pow((float)2,depth1-depth0);//this should probably be a bit shift
            T id0( shift*(*this)[0][0], shift*(*this)[0][1],depth1);
            delta = (*this)[1] - id0;
        } else {
            //Convert second element to depth0 and take difference
            int shift = (int)pow((float)2,depth0-depth1); //this should probably be a bit shift...
            T id1( shift*(*this)[1][0], shift*(*this)[1][1],depth0);
            delta = id1 - (*this)[0];
        }
        return delta;
    }
    
    /// Compute the midpoint index for this Edge in the adaptive grid
    T midpointID()
    {
        //Compute the new depth
        int newDepth = 1+std::max((*this)[0][2],(*this)[1][2]);
        //(*this)[0] will be the left or bottom of edge
        int shift = (int) pow((float)2,newDepth - (*this)[0][2]); //This should probably be a bit shift...
        T midpoint = T(shift*(*this)[0][0],shift*(*this)[0][1],newDepth) + deltaID();
        return midpoint;
    }
    
};


/** A subclass of Edge that represents edges composed of other edges
 *  employed in adaptive sampling technique.
 *
 *  T - essentially an nvis::ivec3 for the 2D map case
 *
 *  Edge_Rotation will work through a set of base-level TerminalEdges (or on terminating leaves)
 *  The CompositeEdges exist in a separate set and are constructed from the information
 *  in the base-level set.
 *
 *  Future:  try to replace pow() calls with a bit shift
 *  Author: Wayne Schlei (Purdue University)
 */
template<typename T >
class CompositeEdge : public TerminalEdge<T> {
public :
    typedef TerminalEdge<T>     edge_type;
    typedef CompositeEdge<T>     self_type;
    typedef typename std::set<T, OnEdgeCompare> NodeSetType;
    
    
    ///Constructors
    CompositeEdge() : edge_type() {}
    CompositeEdge(const T& t0, const T& t1) : edge_type(t0,t1)
    {
        //Insert the edge nodes to the nodeSet
        this->nodeSet.insert( t0 );
        this->nodeSet.insert( t1 );
    }
    CompositeEdge(const edge_type& e) : edge_type(e)
    {
        //Insert the edge nodes to the nodeset
        this->nodeSet.insert( e[0]);
        this->nodeSet.insert( e[1]);
    }
    
    ///Node Set - Nodes that make this edge (sorted by their order on edge)
    NodeSetType  nodeSet;
    
};


//Datastructure namespace - Unused in map topology code
namespace datastructure {

struct Edge {
    Edge() : _i(0), _j(0) {}
    Edge( int i, int j )
    {
        if ( i > j ) {
            _i = j;
            _j = i;
        } else {
            _i = i;
            _j = j;
        }
    }
    Edge( const Edge& e ) : _i(e._i), _j(e._j) {}
    
    Edge& operator=( const Edge& e )
    {
        _i = e._i;
        _j = e._j;
        return *this;
    }
    
    int operator==( const Edge& e ) const
    {
        return ( _i == e._i && _j == e._j );
    }
    
    int operator<( const Edge& e ) const
    {
        return ( _i < e._i ||
                 ( _i == e._i && _j < e._j) );
    }
    
    int _i, _j;
};

struct Segment {
    Segment() : _e0(), _e1() {}
    Segment( const Edge& e0, const Edge& e1 )
    {
        if ( e1 < e0 ) {
            _e0 = e1;
            _e1 = e0;
        } else {
            _e0 = e0;
            _e1 = e1;
        }
    }
    Segment( const Segment& s ) : _e0(s._e0), _e1(s._e1) {}
    
    const Edge& val1() const
    {
        return _e0;
    }
    const Edge& val2() const
    {
        return _e1;
    }
    
    int operator<( const Segment& s ) const
    {
        return ( _e0 < s._e0 ||
                 ( _e0 == s._e0 && _e1 < s._e1 ) );
    }
    
    Edge _e0, _e1;
};

};

};

#endif
