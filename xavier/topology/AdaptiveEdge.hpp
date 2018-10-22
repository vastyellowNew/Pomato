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


/** AdaptiveEdge - Data structure for adaptively computing the rotation of the map displacement
    vector along a specified edge that also accommodates the possible discontinuities.
    Author: Wayne Schlei & Xavier Tricoche (Purdue University)

    Notes:
    1) This works a lot like AdaptiveEdge but with one less dimension
    2) When this is called/setup/run, it is assumed that there are NO invalid
       points in the field.  In other words, the map can be run everywhere on this edge
       apart from singularity conditions.
    2) Have to supply checking/refinement function externally
*/

#ifndef __ADAPTIVE_EDGE_HPP__
#define __ADAPTIVE_EDGE_HPP__

#include <list>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include <iostream>
#include <memory>
#if defined(HX_HAS_STD) || defined(C_0X)
#include <boost/shared_ptr.hpp>
#endif
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
//PMap API
#include <data/edge.hpp>
#include <topology/EdgeRotationFailure.hpp>

namespace topology {



/// Line parameter (0,1) comparison function for AdaptiveEdge Ids
template<typename IDTYPE>
struct LineIDCompare {
    ///Lessthan Comparison of two IdTypes along an edge
    bool operator() (const IDTYPE& lhs, const IDTYPE& rhs) const
    {
        //Find the maximum depth value
        unsigned int maxDepth = std::max(lhs[1],rhs[1]);
        unsigned int lhs0 = toDepth(lhs,maxDepth);
        unsigned int rhs0 = toDepth(rhs,maxDepth);
        return (lhs0 < rhs0);
    }
    
    unsigned int toDepth(const IDTYPE& id, const unsigned int& depth) const
    {
        unsigned int d = id[1];
        int i = id[0];
        if(d==depth) {
            return i;
        }
        while(d<depth) {
            i *= 2;
            d++;
        }
        return i;
    }
};

///AdaptiveEdge Base class : Make a derived class that defines refinement check function
template<typename T, typename VEC>
class AdaptiveEdge {
public :
    typedef AdaptiveEdge<T,VEC>                     SelfType;
    //NonTransverse Edge_Rot:  Map Displacement vector
    //STH Manifold:            Downstream (p) mapping data
    typedef T                                       ValueType;
    typedef VEC                                     PosType;
    //Id's are the tau values of a segment (positions from 0 to 1)
    typedef nvis::ivec2                             IdType;
    typedef std::pair<IdType,IdType>                SegType;
    typedef nvis::lexicographical_order             LexOrder;
    typedef LineIDCompare<IdType>                   LeafCompare;
    typedef std::pair<IdType,ValueType>             DataType;
    typedef EdgeRotationFailure<VEC>                MapDiscont;
    typedef std::map<IdType,MapDiscont,LexOrder>    SepMap;  //(id,MapDiscont)
    typedef std::list<MapDiscont>                   SepList;
    typedef std::map<IdType,ValueType,LexOrder>     DataMap; //(id,ValueType)
    typedef std::map<IdType,PosType,LexOrder>       VertexMap; //(id,PosValue)
    
    ///Constructor with edge end points
    AdaptiveEdge(const VEC& x0, const VEC& x1) :
        _vertexMap(), _dataMap(), _sepDataMap(), _additionalSepData(), _raster(2)
    {
        _endPoints[0] = x0;
        _endPoints[1] = x1;
        _edgeLength = nvis::norm(x1-x0);
        _leaves.clear();
        for(int j=0; j<2; j++) {
            IdType id(j,0);
            _raster[j].reset(new Node(id));
            _leaves.insert(id);
            _vertexMap.insert( std::pair<IdType,PosType>(id,_endPoints[j]) );
            //std::cout << " Added VertexMap: id=" << id << " -> Vertex = " << _vertexMap[id] << "\n";
            //No data to add yet
        }
        
    }
    
    
    ///Reset data arrays (assumes that end points are already set AND you will insert end point values)
    void reset()
    {
        _raster.clear();
        //Clear current information
        _leaves.clear();
        _vertexMap.clear();
        _dataMap.clear();
        _sepDataMap.clear();
        _additionalSepData.clear();
        //Initialize
        _raster.resize(2);
        for(int j=0; j<2; j++) {
            IdType id(j,0);
            _raster[j].reset(new Node(id));
            _leaves.insert(id);
            _vertexMap.insert(std::pair<IdType,PosType>(id,_endPoints[j]));
            // std::cout << " Added VertexMap: id=" << id << " -> Vertex = " << _vertexMap[id] << "\n";
            //No data to add yet
        }
        
    }
    
    ///Reset object and data arrays with new end points
    void reset(const VEC& x0, const VEC& x1)
    {
        _endPoints[0] = x0;
        _endPoints[1] = x1;
        _edgeLength = nvis::norm(x1-x0);
        
        //Data clear & initialize
        reset();
    }
    
    //Access (What are the node values) -----------------------------------------
    ///Return a line parametric value (on [0,1]) for a given node [Assumes each entry is a Midpoint!]
    double getParameter(const IdType& id) const; //assumes nodes are midpoints
    
    ///Return a vector of all pos->value data stored on the edge (returned in on-edge order)
    void getData(std::vector<DataType>& data) const;
    ///Return a vector of all data values
    void getData(std::vector<ValueType>& data) const;
    ///Return a set of all data values
    void getData(std::set<ValueType>& dataSet) const;
    ///Return ordered data output based On-Edge order
    void getDataInOrder(std::vector<ValueType>& orderedData) const;
    ///Flag to see if a point has data
    bool hasData(const IdType& id) const;
    ///Return the value stored at a given index
    const ValueType& getValue(const IdType& id) const;
    ///Insert a id->Value pair into the corresponding node
    void setValue(const IdType& id, const ValueType& v);
    ///Get data map ( id->data list within cell)
    const DataMap& getDataMap() const
    {
        return _dataMap;
    }
    
    ///Get the all the terminating leaves of this grid
    void  getLeaves(std::vector<IdType>& leaves) const;
    ///Get the number of leaves
    int getNumLeaves() const
    {
        return (int) _leaves.size();
    }
    ///Get the leaf indexes of all the terminating leaves
    void getNodes(std::vector<IdType>& nodeIDs) const;
    ///Get the positions of all the Nodes in the terminating leaves
    void getNodes(std::vector<PosType>& verts) const;
    ///Get empty leaves (without data yet)
    void getEmptyLeaves(std::vector<IdType>& emptyLeaves) const;
    ///Get number of empty leaves (nodes with no data)
    int getNumEmptyLeaves() const;
    ///Get new segments for testing
    void getNewSegments(std::vector<SegType>& newSegs) const;
    ///Get all segments for testing
    void getAllSegments(std::vector<SegType>& segs) const;
    ///Get a vertex given an id
    PosType getVertex(const IdType& id) const
    {
        return nodePoint(id);
    }
    ///Set a vertex value given an id
    void setVertex(const IdType& id, const PosType& x);
    ///Get the vertex map container
    const VertexMap& getVertexMap() const
    {
        return _vertexMap;
    }
    /** Perform refinement on a particular node
     *    New Vertex,Data,Sep must be entered after refine() operation
     */
    void refine(const IdType& id);
    ///Find the id of a Node's specified child
    IdType getChildID(const IdType& parentID, const int& child) const;
    ///Compute child id (Not really there, but what it would be | child = 0 or 1)
    IdType computeChildID(const IdType& parentID, const int& child) const
    {
        return sub_id(parentID) + IdType(child,0);
    }
    
    ///Return the leaf id(i,d) of the terminating leaf with the given index (as long as it's ok)
    IdType getLeaf(const int& idx) const;
    
    ///Spacing at a given level (for external access
    double spacing(const unsigned int& d) const
    {
        return resolution(d)*_edgeLength;
    }
    
    
    //Node ID operations ---------------------------------------------------------------------
    ///Return the id of a node one depth down (d+1 or one node-level up) on the tree structure
    IdType sub_id(const IdType& id) const;
    ///Return the id of a node at depth=depth with node id(i,depth)
    IdType raise_id(const IdType& id, unsigned int depth) const;
    ///Compute top-level index for a node id(i,d)
    unsigned int id_to_index(const IdType& id) const;
    
    //Separation/Singularity point commands --------------------------------------------------
    /// Is this node a separation point
    const bool isDiscontinuity(const IdType& node_id) const;
    /// Is this node a separation point
    bool isDiscontinuity(const IdType& node_id);
    /// Is this node a singularity point
    const bool isSingularityPoint(const IdType& node_id) const;
    /// Is this node a singularity point
    bool isSingularityPoint(const IdType& node_id);
    ///Get separation map (id->MapDiscont)
    const SepMap& getDiscontinuityMap() const
    {
        return _sepDataMap;
    }
    ///Get list of discontinuities
    void getDiscontinuityList(std::list<MapDiscont>& theDisconts) const;
    /// Get a map discontinuity (returns an 'Unknown' if not found)
    MapDiscont& getDiscontinuity(const IdType& id) const;
    /// Set whether the node triggers separation/singularity conditions
    void insertDiscontinuityNode(const IdType& id, const MapDiscont& v);
    /// Add a separation/singularity condition location (usually between nodes)
    void insertDiscontinuity(const MapDiscont& v);
    
private:
    ///Leaves (terminal nodes) of the edge (ordered from left->right or bottom->top)
    std::set<IdType,LeafCompare> _leaves;
    ///Map container for ids (key) to vertex values (value)
    VertexMap _vertexMap;
    ///Map container for ids (key) to data values (value)
    DataMap _dataMap;
    ///Map container for ids (key) to map discontinuities (value)
    SepMap  _sepDataMap; //Inserted only if detected ON NODE
    ///List containing additional separation data (usually between nodes)
    SepList _additionalSepData;
    /// End point positions
    PosType  _endPoints[2];
    /// Length of the edge
    double   _edgeLength;
    ///Get a vertex (position vector) of a node
    PosType nodePoint(const IdType& id) const;
    ///Get resolution at a given depth
    double resolution(unsigned int depth) const;
    
public:

    ///Tree-node structure
    struct Node {
        ///Kids
#if defined(HX_HAS_STD) || defined(C_0X)
        boost::shared_ptr<Node>   _child[2];
#else
        std::shared_ptr<Node>   _child[2];
#endif
        ///(i,d)
        IdType _id;
        ///Constructor:
#if defined(HX_HAS_STD) || defined(C_0X)
        Node() : _id(0,-1)               {  } //Pointers are auto-set to NULL in boost
        Node(const IdType& id) : _id(id) {  }
        
        //Access
        const IdType&  get_id() const
        {
            return _id;
        }
        
        const boost::shared_ptr<Node>    get_child(int i) const
        {
            return _child[i];
        }
        boost::shared_ptr<Node>          get_child(int i)
        {
            return _child[i];
        }
#else
        Node() : _id(0,-1)
        {
            _child[0] = NULL;
        }
        Node(const IdType& id) : _id(id)
        {
            _child[0] = NULL;
        }
        //Access
        const IdType&  get_id() const
        {
            return _id;
        }
        
        const std::shared_ptr<Node>    get_child(int i) const
        {
            return _child[i];
        }
        std::shared_ptr<Node>          get_child(int i)
        {
            return _child[i];
        }
#endif
        ///Get a child's id
        IdType get_id(int i) const
        {
            switch (i) {
                case  0:
                    return _id;
                case  1:
                    return _id + IdType(1,0);
                default:
                    return IdType(2*_id[0]+1, _id[1]+1);
            }
        }
        
        bool is_leaf() const
        {
            return (_child[0]==NULL);
        }
        ///Create children nodes when refinement needed (Left Node and midpoint of segment)
        void split()
        {
            for (int i=0 ; i<2 ; ++i) {
                _child[i].reset(new Node());
            }
            IdType new_id(2*_id[0], _id[1]+1);
            _child[0]->_id = new_id;
            _child[1]->_id = new_id + IdType(1,0);
        }
    };
    ///Refine a node
#if defined(HX_HAS_STD) || defined(C_0X)
    void refine(boost::shared_ptr<Node> n);
#else
    void refine(std::shared_ptr<Node> n);
#endif
    
private :
#if defined(HX_HAS_STD) || defined(C_0X)
    std::vector< boost::shared_ptr<Node> > _raster;
#else
    std::vector< std::shared_ptr<Node> > _raster;
#endif
};


///Return resolution at a given depth
template<typename T,typename VEC>
inline double AdaptiveEdge<T,VEC>::
resolution(unsigned int depth) const
{
    int f = 1;
    f = f << depth;
    return 1./(double)f;
}


///Get a vertex (position vector) of a node
template<typename T, typename VEC>
inline typename AdaptiveEdge<T, VEC>::PosType
AdaptiveEdge<T, VEC>::
nodePoint(const typename AdaptiveEdge<T,VEC>::IdType& id) const
{
    typename VertexMap::const_iterator it = _vertexMap.find(id);
    if ( it == _vertexMap.end() ) {
        throw std::runtime_error("No vertex data for node");
    }
    return it->second;
    //Old format:  Not consistent evaluation!
    /*double tau = getParameter(id);
    //Position from linear parameterization
    return (1.0-tau)*_endPoints[0] + tau*_endPoints[1];*/
}

///Return a line parametric value (on [0,1]) for a given node [Assumes each entry is a Midpoint!]
template<typename T,typename VEC>
inline double AdaptiveEdge<T, VEC>::
getParameter(const typename AdaptiveEdge<T,VEC>::IdType& id) const
{
    double dtau = resolution(id[1]);
    return dtau*id[0]; //Assuming this is a midpoint!
}

///Return a vector of all posID->value data stored on the edge (returned in on-edge order)
template<typename T,typename VEC>
inline void AdaptiveEdge<T, VEC>::
getData(std::vector<typename AdaptiveEdge<T, VEC>::DataType>& data) const
{
    data.clear();
    typename std::set<IdType,LeafCompare>::const_iterator leafIT;
    typename DataMap::const_iterator dataIT;
    for(leafIT=_leaves.begin(); leafIT!=_leaves.end(); ++leafIT) {
        dataIT = _dataMap.find( (*leafIT) );
        if (dataIT == _dataMap.end()) {
            //No data ->Probably not the best idea to call this without full data
            data.push_back( DataType((*leafIT),ValueType(50.0)) );
        } else {
            data.push_back( DataType((*leafIT),dataIT->second) );
        }
    }
}

///Return a vector of all data values
template<typename T,typename VEC>
inline void AdaptiveEdge<T, VEC>::
getData(std::vector<typename AdaptiveEdge<T, VEC>::ValueType>& data) const
{
    data.clear();
    typename DataMap::const_iterator vDataIT;
    for(vDataIT=_dataMap.begin(); vDataIT!=_dataMap.end(); ++vDataIT) {
        data.push_back( vDataIT->second );
    }
}

///Set return assuming ValueType has some inherent sorting
template<typename T,typename VEC>
inline void AdaptiveEdge<T, VEC>::
getData(std::set<typename AdaptiveEdge<T, VEC>::ValueType>& dataSet) const
{
    dataSet.clear();
    typename DataMap::const_iterator vDataIT;
    for(vDataIT=_dataMap.begin(); vDataIT!=_dataMap.end(); ++vDataIT) {
        dataSet.insert( vDataIT->second );
    }
}

///Return ordered data output based On-Edge compare (Call with full data)
template<typename T,typename VEC>
inline void AdaptiveEdge<T, VEC>::
getDataInOrder(std::vector<typename AdaptiveEdge<T, VEC>::ValueType>& orderedData) const
{
    orderedData.clear();
    //Work through leaves which are in order along edge
    typename std::set<IdType,LeafCompare>::const_iterator leafIT;
    typename DataMap::const_iterator dataIT;
    for(leafIT=_leaves.begin(); leafIT!=_leaves.end(); ++leafIT) {
        dataIT = _dataMap.find( (*leafIT) );
        if (dataIT == _dataMap.end()) {
            //No data ->Probably not the best idea to call this without full data
            orderedData.push_back( ValueType(50.0) );
        } else {
            orderedData.push_back( dataIT->second );
        }
    }
    
}

///Flag to see if a point has data
template<typename T, typename VEC>
bool AdaptiveEdge<T, VEC>::
hasData(const IdType& id) const
{
    typename DataMap::const_iterator it = _dataMap.find(id);
    return (it != _dataMap.end());
}
///Return the value stored at a given index
template<typename T, typename VEC>
const typename AdaptiveEdge<T, VEC>::ValueType&
AdaptiveEdge<T, VEC>::
getValue(const typename AdaptiveEdge<T, VEC>::IdType& id) const
{
    typename  DataMap::const_iterator it = _dataMap.find(id);
    if (it == _dataMap.end()) {
        throw std::runtime_error("Provided (i,d) has no data in AdaptiveEdge");
    }
    return it->second;
}

///Insert a id->Value pair into the corresponding node
template<typename T,typename VEC>
void AdaptiveEdge<T, VEC>::
setValue(const typename AdaptiveEdge<T, VEC>::IdType& id,
         const typename AdaptiveEdge<T, VEC>::ValueType& v)
{
    std::pair<typename DataMap::iterator,bool> ret;
    ret = _dataMap.insert(std::pair<IdType,ValueType>(id,v));
    //If Data value already exists, overwrite
    if(ret.second == false) {
        _dataMap[id] = v;
    }
    
}

///Set a vertex value given an id
template<typename T,typename VEC>
void AdaptiveEdge<T, VEC>::
setVertex(const typename AdaptiveEdge<T, VEC>::IdType& id,
          const typename AdaptiveEdge<T, VEC>::PosType& x)
{
    std::pair<typename VertexMap::iterator,bool> ret;
    ret = _vertexMap.insert( std::pair<IdType,PosType>(id,x) );
//if vertex value already exists, overwrite
    if(ret.second == false) {
        _vertexMap[id] = x;
    }
//std::cout << " Added VertexMap: id=" << id << " -> Vertex = " << _vertexMap[id] << "\n";

}

///Get the all the terminating leaves of this grid
template<typename T,typename VEC>
inline void AdaptiveEdge<T, VEC>::
getLeaves(std::vector<typename AdaptiveEdge<T, VEC>::IdType>& leaves) const
{
    leaves.clear();
    std::copy(_leaves.begin(), _leaves.end(), std::back_inserter(leaves));
}


///Get the positions of all the Nodes in the terminating leaves
template<typename T,typename VEC>
inline void AdaptiveEdge<T, VEC>::
getNodes(std::vector<typename AdaptiveEdge<T,VEC>::PosType>& verts) const
{
    verts.clear();
    typename std::set<IdType,LeafCompare>::iterator leafIT;
    for( leafIT=_leaves.begin(); leafIT!=_leaves.end(); ++leafIT) {
        verts.push_back( nodePoint( (*leafIT) ) );
    }
}

/// Is this node a separation point
template<typename T,typename VEC>
const bool AdaptiveEdge<T, VEC>::
isDiscontinuity(const typename AdaptiveEdge<T, VEC>::IdType& id) const
{
    typename SepMap::const_iterator sepIT = _sepDataMap.find(id);
    if (sepIT != _sepDataMap.end()) {
        return true;
    } else {
        return false;
    }
}

/// Is this node a separation point
template<typename T,typename VEC>
bool AdaptiveEdge<T, VEC>::
isDiscontinuity(const typename AdaptiveEdge<T, VEC>::IdType& id)
{
    typename SepMap::iterator sepIT = _sepDataMap.find(id);
    if (sepIT != _sepDataMap.end()) {
        return true;
    } else {
        return false;
    }
}

/// Is this node a singularity point
template<typename T,typename VEC>
const bool AdaptiveEdge<T, VEC>::
isSingularityPoint(const typename AdaptiveEdge<T, VEC>::IdType& id) const
{
    typename SepMap::const_iterator sepIT = _sepDataMap.find(id);
    if (sepIT != _sepDataMap.end()) {
        if((sepIT->second.type == MapDiscont::BACKWARD_SINGULARITY) ||
                (sepIT->second.type == MapDiscont::FORWARD_SINGULARITY) ) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}
/// Is this node a singularity point
template<typename T,typename VEC>
bool AdaptiveEdge<T, VEC>::
isSingularityPoint(const typename AdaptiveEdge<T, VEC>::IdType& id)
{
    typename SepMap::const_iterator sepIT = _sepDataMap.find(id);
    if (sepIT != _sepDataMap.end()) {
        if((sepIT->second.type == MapDiscont::BACKWARD_SINGULARITY) ||
                (sepIT->second.type == MapDiscont::FORWARD_SINGULARITY) ) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

/// Get a map discontinuity (returns an 'Unknown' with p=0 if not found)
template<typename T,typename VEC>
inline typename AdaptiveEdge<T,VEC>::MapDiscont&
AdaptiveEdge<T,VEC>::
getDiscontinuity(const typename AdaptiveEdge<T,VEC>::IdType& id) const
{
    typename SepMap::const_iterator sepIT = _sepDataMap.find(id);
    if (sepIT != _sepDataMap.cend()) {
        return (*sepIT);
    } else {
        return MapDiscont(MapDiscont::UNKNOWN,nodePoint(id),0);
    }
}

/// Set whether the node triggers separation/singularity conditions
template<typename T,typename VEC>
inline void AdaptiveEdge<T,VEC>::
insertDiscontinuityNode(const typename AdaptiveEdge<T,VEC>::IdType& id,
                        const typename AdaptiveEdge<T,VEC>::MapDiscont& v)
{
    std::pair<typename SepMap::iterator,bool> ret;
    ret = _sepDataMap.insert( std::pair<IdType,MapDiscont>(id,v) );
    //if sep value already exists, overwrite
    if(ret.second == false) {
        _sepDataMap[id] = v;
    }
    //std::cout << " Inserted Sep Data [" << id << "]->[" << v.what() << "]\n";
}

template<typename T,typename VEC>
inline void AdaptiveEdge<T,VEC>::
insertDiscontinuity(const typename AdaptiveEdge<T,VEC>::MapDiscont& v)
{
    //std::cout << " Adding Sep Data ->[" << v.what() << "] at " << v.failurePos << "\n";
    _additionalSepData.push_back(v);
    //std::cout << "  Back of _additionalSepData = " << _additionalSepData.back().what()
    //<< " at " << _additionalSepData.back().failurePos << "\n";
}

/// Get a list of discontinuities
template<typename T,typename VEC>
inline void AdaptiveEdge<T,VEC>::
getDiscontinuityList(std::list<typename AdaptiveEdge<T,VEC>::MapDiscont>& theDisconts) const
{
    theDisconts.clear();
    //Loop through the nodes
    typename SepMap::const_iterator sepIT;
    for(sepIT=_sepDataMap.begin(); sepIT!=_sepDataMap.end(); ++sepIT) {
        theDisconts.push_back( sepIT->second );
    }
    //Add the additional supplement of discontinuities (midpoints of remaining segments)
    typename SepList::const_iterator sIT;
    for(sIT=_additionalSepData.begin(); sIT!=_additionalSepData.end(); ++sIT) {
        theDisconts.push_back( (*sIT) );
    }
}


/// Get the leaf corresponding to a set index (like a vector look-up)
template<typename T, typename VEC>
inline typename AdaptiveEdge<T, VEC>::IdType
AdaptiveEdge<T,VEC>::
getLeaf(const int& idx) const
{
    std::set<IdType,LeafCompare>::iterator lit;
    if (idx >= (int)_leaves.size() || idx<0) {
        return IdType(0,0);
    }
    lit = _leaves.begin();
    std::advance(lit,idx);
    return (*lit);
}

///Return the id of a cell one depth down (or one node-level up) on the tree structure
template<typename T, typename VEC>
inline typename AdaptiveEdge<T, VEC>::IdType
AdaptiveEdge<T, VEC>::
sub_id(const typename AdaptiveEdge<T, VEC>::IdType& id) const
{
    return IdType(id[0]*2, id[1]+1);
}


///Return the id of a cell at depth=depth containing cell id(i,d)
template<typename T, typename VEC>
inline typename AdaptiveEdge<T, VEC>::IdType
AdaptiveEdge<T, VEC>::
raise_id(const typename AdaptiveEdge<T, VEC>::IdType& id, unsigned int depth) const
{
    const int& d = id[1];
    if (!d) {
        return id;
    }
    assert((int)depth < d);
    int f = 1;
    f = f << (d-depth);
    return IdType(id[0]/f, depth);
}

///Get the index of the starting grid (MxN) given an id (i,depth)
template<typename T, typename VEC>
inline unsigned int AdaptiveEdge<T, VEC>::
id_to_index(const typename AdaptiveEdge<T, VEC>::IdType& id) const
{
    IdType top_id = raise_id(id, 0);
    return top_id[0];
}



///Refine a node
template<typename T, typename VEC>
inline void AdaptiveEdge<T, VEC>::
#if defined(HX_HAS_STD) || defined(C_0X)
refine(boost::shared_ptr<typename AdaptiveEdge<T,VEC>::Node> n)
#else
refine(std::shared_ptr<typename AdaptiveEdge<T,VEC>::Node> n)
#endif
{
    n->split();
    const IdType& baseID = n->get_id();
    //Pull out the base ID from _leaves container
    _leaves.erase(baseID);
    // create (empty) entries for new nodes
    for (int i=0; i<2; i++) {
        const IdType& id = n->get_child(i)->get_id();
        _leaves.insert(id);
    }
    
    //Copy existing data to next level (Only Left point)
    const IdType& newID = n->get_child(0)->get_id();
    //Values (only know parent node value)
    typename DataMap::iterator dataIT = _dataMap.find(baseID);
    if (dataIT != _dataMap.end()) {
        _dataMap.insert( std::pair<IdType,ValueType>(newID, dataIT->second) );
    }
    
    //Vertex of parent node
    typename VertexMap::iterator vertIT = _vertexMap.find(baseID);
    if (vertIT != _vertexMap.end()) {
        _vertexMap.insert( std::pair<IdType,PosType>(newID, vertIT->second) );
    }
    
    //Discontinuities (only if on parent node)
    typename SepMap::iterator sepIT = _sepDataMap.find(baseID);
    if (sepIT != _sepDataMap.end()) {
        _sepDataMap.insert( std::pair<IdType,MapDiscont>(newID, sepIT->second) );
    }
    
    //Note:  You must add data with setVertex(),setValue(),insertDiscontinuityNode()
    //       commands AFTER running the refine() command for the new points!
}

///Perform refinement on a particular node
template<typename T, typename VEC>
inline void AdaptiveEdge<T, VEC>::
refine(const typename AdaptiveEdge<T,VEC>::IdType& idIN)
{
    // (i,d) -> (i/2,d-1) -> ... -> (i/2^d, 0)
    IdType theID = raise_id(idIN,0);
    const int& depth = idIN[1];
    
    // Get the node for the depth=0 cell (should be left/bottom point)
#if defined(HX_HAS_STD) || defined(C_0X)
    boost::shared_ptr<Node> n = _raster[id_to_index(idIN)];
#else
    std::shared_ptr<Node> n = _raster[id_to_index(idIN)];
#endif
    
    //Travel down levels to find node to refine (through structure in _raster)
    for (int d=1; d<=depth; ++d) {
        IdType baseID = sub_id(theID); //Find ID at next level
        IdType actualID = idIN;
        if(d!=depth) {
            actualID = raise_id(idIN,d);
        }
        IdType diff = actualID - baseID;
        if (nvis::all(diff == IdType(0,0))) {
            n = n->get_child(0);
        } else {
            n = n->get_child(1);
        }
        //Get the new node id
        theID = n->get_id();
    }
    //Run Node refine()
    refine(n);
    
}

///Find the id of a Node's specified child
template<typename T, typename VEC>
inline typename AdaptiveEdge<T,VEC>::IdType AdaptiveEdge<T,VEC>::
getChildID(const typename AdaptiveEdge<T,VEC>::IdType& parentID, const int& child) const
{
    // (i,d) -> (i/2,d-1) -> ... -> (i/2^d, 0)
    if(std::abs(child)>1) {
        throw std::runtime_error("Child of node should be 0 or 1");
    }
    IdType childID(parentID[0]*2+child, parentID[1]+1);
    return childID;
}

///Get empty leaves (without data yet)
template<typename T, typename VEC>
inline void AdaptiveEdge<T,VEC>::
getEmptyLeaves(std::vector<typename AdaptiveEdge<T,VEC>::IdType>& emptyLeaves) const
{
    emptyLeaves.clear();
    typename std::set<IdType,LeafCompare>::const_iterator leafIT;
    typename DataMap::const_iterator dataMapIT;
    for(leafIT=_leaves.begin(); leafIT!=_leaves.end(); ++leafIT) {
        dataMapIT = _dataMap.find( (*leafIT) );
        if (dataMapIT == _dataMap.end()) {
            emptyLeaves.push_back( (*leafIT) );
        }
    }
}

///Get number of empty leaves (nodes with no data)
template<typename T, typename VEC>
int AdaptiveEdge<T,VEC>::getNumEmptyLeaves() const
{
    std::vector<IdType> emptyLeaves;
    getEmptyLeaves(emptyLeaves);
    return (int) emptyLeaves.size();
}

///Get new segments for testing after a call to refine()
template<typename T, typename VEC>
inline void AdaptiveEdge<T,VEC>::
getNewSegments(std::vector<typename AdaptiveEdge<T,VEC>::SegType>& newSegs) const
{
    newSegs.clear();
    getAllSegments(newSegs);
    typename std::vector<SegType>::iterator segIT = newSegs.begin();
    int numSegs = (int) newSegs.size();
    for(int i=0; i<numSegs; i++) {
        IdType id0 = (*segIT).first;
        IdType id1 = (*segIT).second;
        //std::cout << "  Segment : " << id0 << " -> " << id1 << " \n";
        //std::cout << "    Node0 - hasData=" << hasData(id0) << " isDiscontinuity=" << isDiscontinuity(id0) << "\n";
        //std::cout << "    Node1 - hasData=" << hasData(id1) << " isDiscontinuity=" << isDiscontinuity(id1) << "\n";
        //Check for separation nodes or data
        if ( (!hasData(id0) && !isDiscontinuity(id0) ) ||  //Node 0 is empty and is not a discontinuity
                (!hasData(id1) && !isDiscontinuity(id1) ) ) { //Node 1 is empty and is not a discontinuity
            //We need to test this segment, so move to the next
            ++segIT;
        } else {
            //If both nodes have data or are separation points,
            // we don't need to check this segment
            segIT = newSegs.erase(segIT);
            /*if(segIT!=newSegs.end()) {
              std::cout << "    Erasing segment! Move iterator to " << (*segIT).first << " -> " << (*segIT).second << "\n";
            } else {
              std::cout << "    Erasing segment! End of loop (i=" << i << " of " << numSegs << ")\n";
            }*/
        }
    }
    
}

///Get all segments for testing
template<typename T, typename VEC>
inline void AdaptiveEdge<T,VEC>::
getAllSegments(std::vector<typename AdaptiveEdge<T,VEC>::SegType>& segs) const
{
    segs.clear();
    //Loop through all leaves (i,i-1) to find segments
    typename std::set<IdType,LeafCompare>::const_iterator leafIT, lastLeafIT;
    leafIT = _leaves.begin();
    lastLeafIT = leafIT;
    //Generate all segments
    for(++leafIT; leafIT!=_leaves.end(); ++leafIT,++lastLeafIT) {
        segs.push_back( SegType((*lastLeafIT),(*leafIT)) );
    }
}


} // end topology


#endif //End __ADAPTIVE_EDGE_HPP__
