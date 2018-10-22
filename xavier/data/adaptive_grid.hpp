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


/** AdaptiveGrid - Data structure for adaptively modifying grid via a quad-tree structure
 based on creating a convex hull with winding number.
 Author:  Xavier Tricoche & Wayne Schlei (Purdue University)

 Notes:
  1) The criteria to refine a cell should be evaluated before a call to AdaptiveGrid::refine()
  2) Evaluate based on the convex hull of the winding number(s) concept in map topology analysis
     -Max change in winding number over cell
     -Min cell size (maximum depth parameter)
     -Forbidden cells (with node in forbidden region) will be refined to max depth
  3) (QuadTree ONLY) This version ONLY works in 2D section space!!!  Will need upgrade to general ND case...
*/

#ifndef __ADAPTIVE_GRID_HPP__
#define __ADAPTIVE_GRID_HPP__

#include <list>
#include <vector>
#include <map>
#include <set>
#include <exception>
#include <stdexcept>
#include <memory>
#if defined(HX_HAS_STD) || defined(C_0X)
#include <boost/shared_ptr.hpp>
#endif
#include <assert.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif



namespace xavier {

/** Default ValueTraits mechanism to indicate unknown or invalid values
 */
struct DefaultValueTraits {
public:
    DefaultValueTraits() {}
    
    //Static members signifying "invalid" values and "unknown" values
    static const double invalid;
    static const double unknown;
    //Initialize independently in a staticInit.cpp file
    //const double xavier::DefaultValueTraits::invalid = 1000.;
    //const double xavier::DefaultValueTraits::unknown = 0.0;
};


/** Position type paired with depth*/
struct posWithDepth : public nvis::vec2 {
public:
    size_t depth;
    posWithDepth(const nvis::vec2& p, const size_t& d) :
        nvis::vec2(p), depth(d) {}
    posWithDepth(const nvis::vec2& p) :
        nvis::vec2(p), depth(0) {}
};

/** Comparison structure : compares a nodeID as a subID given preset maximum depth
    - Author:  Wayne Schlei
*/
struct ID_Compare {
    ///Maximum depth for comparison parameters
    static int maxDepth;
    
    ///Comparison for an nvis::ivec3 vector (=Vec)
    template<typename Vec>
    bool operator() (const Vec& lhs, const Vec& rhs) const
    {
        //Convert to sub_id first at max depth
        int lhsD = lhs[2];
        int rhsD = rhs[2];
        int lhsShift = (int) pow((float)2.0,maxDepth-lhsD);
        int rhsShift = (int) pow((float)2.0,maxDepth-rhsD);
        Vec lhsMaxD(lhsShift*lhs[0],lhsShift*lhs[1],maxDepth);
        Vec rhsMaxD(rhsShift*rhs[0],rhsShift*rhs[1],maxDepth);
        //Run comparison with lexicographical_order
        nvis::lexicographical_order baseCompare;
        bool result = baseCompare(lhsMaxD,rhsMaxD);
        return result;
    }
};

/** Class responsible for adaptively refining the poincare map sampling mesh and housing
 * relevant container information. (i.e., quad tree info)
 *  T - winding number type (e.g., vec3 for planar cr3bp)
 *  ValueTraits - Mechanism to describe a value as unknown or invalid
 *  VertexInfo - Data objects stored per vertex (doubles up storage?)
 */
template<typename T, typename ValueTraits = DefaultValueTraits>
class AdaptiveGrid {
public:
    typedef AdaptiveGrid<T, ValueTraits>                                    self_type;
    typedef T                                                               value_type;
    typedef ValueTraits                                                     traits_type;
    typedef nvis::vec2                                                      pos_type;
    typedef std::pair<posWithDepth, value_type>                             data_type;
    typedef nvis::bbox2                                                     bounds_type;
    typedef nvis::ivec3                                                     id_type;
    typedef std::list<data_type>                                            cell_data_type;
    typedef std::map<id_type, value_type, nvis::lexicographical_order>      vertex_container_type;
    typedef std::map<id_type, bool, nvis::lexicographical_order>            validity_container_type;
    typedef std::map<id_type, cell_data_type, nvis::lexicographical_order>  cell_container_type;
    
    ///Constructor:
    AdaptiveGrid(const bounds_type& bounds, size_t sizeX, size_t sizeY)
        : _bounds(bounds), _size(sizeX, sizeY), _step(bounds.size() / _size),
          _v_data(), _c_data(), _c_valid(), _raster(sizeX*sizeY)
    {
        _leaves.clear();
        // initialize cells' data with empty lists
        for (int i=0 ; i<(((int)sizeX)*((int)sizeY)) ; ++i) {
            int x = i % (int) sizeX;
            int y = i / (int) sizeX;
            id_type id(x,y,0);
            _raster[i].reset( new Node(id) );
            _leaves.insert(id);
            _c_data[id_type(x,y,0)] = cell_data_type();
            _c_valid[id_type(x,y,0)] = true;
        }
        // initialize vertexes's data with "unknowns" values
        const value_type unknown = ValueTraits::unknown;
        for (int i=0 ; i<(((int)sizeX+1)*((int)sizeY+1)) ; ++i) {
            int x = i % ((int)sizeX+1);
            int y = i / ((int)sizeX+1);
            _v_data[id_type(x,y,0)] = unknown;
        }
    }
    
    /// Reset the grid to new bounds/sizes
    void reset(const bounds_type& bounds, size_t sizeX, size_t sizeY)
    {
        _bounds = bounds;
        _size = pos_type(sizeX,sizeY);
        _raster.clear();
        _raster.resize(sizeX*sizeY);
        //Initialize
        _step = _bounds.size() / _size;
        _leaves.clear();
        //Clear current information
        _v_data.clear();
        _c_data.clear();
        _c_valid.clear();
        // initialize cells' data with empty lists
        for (int i=0 ; i<(int)sizeX*(int)sizeY ; ++i) {
            int x = i % (int) sizeX;
            int y = i / (int) sizeX;
            id_type id(x,y,0);
            _raster[i].reset( new Node(id) );
            _leaves.insert(id);
            _c_data[id_type(x,y,0)] = cell_data_type();
            _c_valid[id_type(x,y,0)] = true;
        }
        // initialize vertexes's' data with "unknown" values
        const value_type unknown = ValueTraits::unknown;
        for (int i=0 ; i<((int)sizeX+1)*((int)sizeY+1) ; ++i) {
            int x = i % ((int)sizeX+1);
            int y = i / ((int)sizeX+1);
            _v_data[id_type(x,y,0)] = unknown;
        }
    }
    
    //Cell access (what is in a cell or cell container) ------------------------------------
    
    ///Insert a position->paramValue pair into the corresponding cell
    id_type insertCellValue(const pos_type& x, const value_type& v, const int& depth);
    
    ///Return the list of cell data (pair: pos->parameter)
    const cell_data_type& getCellData(const id_type& cell_id) const;
    ///Return a vector of all pos->value data stored in the grid
    void getData(std::vector<data_type>& data);
    ///Get cell map ( id->data list within cell)
    const cell_container_type& getCellDataMap() const
    {
        return _c_data;
    }
    ///Get cell map (id->validity of cell)
    const validity_container_type& getCellValidityMap() const
    {
        return _c_valid;
    }
    ///Get cells that need data inside
    void  getEmptyCells(std::vector<id_type>& empty) const;
    ///Get the all the terminating leaves of this grid
    void  getLeaves(std::vector<id_type>& leaves) const;
    ///Get the number of leaves
    int getNumLeaves() const
    {
        return (int) _leaves.size();
    }
    ///Get the leaf indexes of all the corners in the terminating leaves
    void getGridNodes(std::vector<id_type>& cornerIDs) const;
    ///Get the positions of all the corners in the terminating leaves
    void getGridNodes(std::vector<pos_type>& corners) const;
    
    ///Return a cell's bounding box given id(i,j,d)
    bounds_type cellBounds(const id_type& id) const;
    
    ///Return the cell id(i,j,d) that surrounds a given coordinate (returns (-1,-1,-1) if outside grid bounds)
    id_type getCell(const pos_type& x) const;
    ///Return the leaf id(i,j,d) of the terminating leaf with the given index (as long as it's ok)
    id_type getLeaf(const int& idx) const;
    
    ///Spacing at a given level (for external access
    pos_type spacing(const unsigned int& d)
    {
        return resolution(d);
    }
    pos_type spacing(const unsigned int& d) const
    {
        return resolution(d);
    }
    
    
    //Node ID operations ---------------------------------------------------------------------
    ///Return the id of a cell one depth down (d+1 or one node-level up) on the tree structure
    id_type sub_id(const id_type& id) const;
    ///Return the id of a cell at depth=depth containing cell id(i,j,d)
    id_type raise_id(const id_type& id, unsigned int depth) const;
    
    ///Compute top-level grid index for a cell id(i,j,d)
    unsigned int id_to_index(const id_type& id) const;
    
    
    
    //Vertex access -------------------------------------------------------------------------
    ///Get the vertex given an id(i,j,d)
    pos_type  getVertex(const id_type& vert_id) const;
    ///Set the parameter value at a vertex given an id(i,j,d)
    void setVertexValue(const id_type& vert_id, const value_type& v);
    
    ///Get the parameter value at a vertex given an id(i,j,d)
    const value_type& getVertexValue(const id_type& vert_id) const;
    ///Get all the grid vertexes
    const vertex_container_type&  getvertexes() const
    {
        return _v_data;
    }
    ///Get vertexes that are undefined
    void  getUndefinedVertexes(std::vector<id_type>& undef) const;
    
    ///Command to performs refinement.
    void refine(const id_type& cell_id);
    
    //Cell validity commands ----------------------------------------------------------------
    /// Is this node valid (i.e., no "invalid" value)
    const bool valid(const id_type& node_id) const;
    /// Is this node valid (i.e., no "invalid" value)
    bool valid(const id_type& node_id);
    /// Set whether the cell is valid or not
    void setCellValidity(const id_type& cell_id, const bool v);
    /// Update cell validity values after data is inserted
    void updateCellValidity();
    /// Is this cell valid (i.e., specified "valid" flag per cell is ok)
    const bool isCellValid(const id_type& cell_id) const;
    /// Are all the corners of a cell invalid
    const bool allCornersInvalid(const id_type& cell_id) const;
    /// Are all the corners of a cell valid
    const bool allCornersValid(const id_type& cell_id) const;
    
    ///Remove the leaves with an "invalid" value at a corner
    void removeInvalidLeaves();
    
    
private:
    bounds_type                                     _bounds;
    pos_type                                        _size;
    pos_type                                        _step;
    ///Map container for vertex ids (key) to value at vertex (value)
    vertex_container_type                           _v_data;
    ///Map container for cell ids (key) to list of data in cell (value)
    cell_container_type                             _c_data; //Map cellID->List of data in cell
    ///Map container for cell ids (key) to cell validity (value) which can be set externally
    validity_container_type                         _c_valid;
    std::set<id_type, nvis::lexicographical_order>  _leaves;
    
    ///Get a vertex (position vector) of a node given index i (0-LL,1-LR,2-UR,3-UL,4-C)
    pos_type   node_point(const id_type&, int i) const;
    ///Return resolution at a given depth
    pos_type   resolution(unsigned int depth) const;
    
public:
    ///Tree-node structure
    struct Node {
        ///Kids
#if defined(HX_HAS_STD) || defined(C_0X)
        boost::shared_ptr<Node> _child[4];
#else
        std::shared_ptr<Node>  _child[4];
#endif
        //Node*   _child[4];
        ///(i,j,d)
        id_type _id;
#if defined(HX_HAS_STD) || defined(C_0X)
        ///Constructor:
        Node() : _id(0, 0,-1)             {  } //Pointers are auto-set to NULL in boost
        Node(const id_type& id) : _id(id) {  }
        //Access
        const id_type& get_id() const
        {
            return _id;
        }
        const boost::shared_ptr<Node>   get_child(int i) const
        {
            return _child[i];
        }
        boost::shared_ptr<Node>         get_child(int i)
        {
            return _child[i];
        }
#else
        ///Constructor:
        Node() : _id(0, 0,-1)
        {
            _child[0] = NULL;
        }
        Node(const id_type& id) : _id(id)
        {
            _child[0] = NULL;
        }
        //Access
        const id_type& get_id() const
        {
            return _id;
        }
        const std::shared_ptr<Node>   get_child(int i) const
        {
            return _child[i];
        }
        std::shared_ptr<Node>         get_child(int i)
        {
            return _child[i];
        }
#endif
        
        //const Node*    get_child(int i) const { return _child[i]; }
        //Node*          get_child(int i)       { return _child[i]; }
        
        id_type get_id(int i) const
        {
            switch (i) {
                case  0:
                    return _id;
                case  1:
                    return _id + id_type(1,0,0);
                case  2:
                    return _id + id_type(1,1,0);
                case  3:
                    return _id + id_type(0,1,0);
                default:
                    return id_type(2*_id[0]+1, 2*_id[1]+1, _id[2]+1);
            }
        }
        
        bool is_leaf() const
        {
            return (_child[0]==NULL);
        }
        ///Create children nodes when refinement needed
        void split()
        {
            //for (int i=0 ; i<4 ; ++i) _child[i] = new Node();
            for (int i=0 ; i<4 ; ++i) {
                _child[i].reset(new Node());
            }
            id_type new_id(2*_id[0], 2*_id[1], _id[2]+1);
            _child[0]->_id = new_id;
            _child[1]->_id = new_id + id_type(1,0,0);
            _child[2]->_id = new_id + id_type(1,1,0);
            _child[3]->_id = new_id + id_type(0,1,0);
        }
    };
#if defined(HX_HAS_STD) || defined(C_0X)
    void refine(boost::shared_ptr<Node> cell);
#else
    void refine(std::shared_ptr<Node> cell);
#endif
    
private:
#if defined(HX_HAS_STD) || defined(C_0X)
    /// Smart Pointers to each node (auto deletes objects when references are gone)
    std::vector< boost::shared_ptr<Node> >   _raster;
#else
    /// Smart Pointers to each node (auto deletes objects when references are gone)
    std::vector< std::shared_ptr<Node> >   _raster;
#endif
};

///Return resolution at a given depth
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::pos_type
AdaptiveGrid<T, ValueTraits>::
resolution(unsigned int depth) const
{
    int f = 1;
    f = f << depth;
    return 1./(double)f * _step;
}

///Get the vertex given an id(i,j,d)
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::pos_type
AdaptiveGrid<T, ValueTraits>::
getVertex(const typename AdaptiveGrid<T, ValueTraits>::id_type& id) const
{
    pos_type delta = resolution(id[2]);
    return _bounds.min() + delta*nvis::vec2(nvis::subv<0,2,int,3>(id));
}

///Return a bounding box of a cell indicated by id(i,j,d)
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::bounds_type
AdaptiveGrid<T, ValueTraits>::
cellBounds(const typename AdaptiveGrid<T, ValueTraits>::id_type& id) const
{
    pos_type delta = resolution(id[2]);
    pos_type _min = _bounds.min() + delta*nvis::vec2(nvis::subv<0,2,int,3>(id));
    return bounds_type(_min, _min + delta);
}

///Return the cell that contains a given coordinate (returns (-1,...,-1) if outside grid bounds)
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::id_type
AdaptiveGrid<T,ValueTraits>::
getCell(const typename AdaptiveGrid<T,ValueTraits>::pos_type& x) const
{
    //Return (-1,-1,-1) if outside grid (other version throws a runtime_error)
    if ( !_bounds.inside(x) ) {
        return id_type(-1);
    }
    
    pos_type tmp = (x-_bounds.min())/_step;
    int idx = floor(tmp[1])*_size[0] + floor(tmp[0]);
    
    //Search down leaves
#if defined(HX_HAS_STD) || defined(C_0X)
    boost::shared_ptr<Node> n = _raster[idx];
#else
    std::shared_ptr<Node> n = _raster[idx];
#endif
    while (!n->is_leaf()) {
        const id_type& id = n->get_id();
        pos_type c = node_point(id, -1);
        if (nvis::all(x <= c)) {
            n = n->_child[0];
        } else if (nvis::all(x > c)) {
            n = n->_child[2];
        } else if (x[0] < c[0]) {
            n = n->_child[3];
        } else {
            n = n->_child[1];
        }
    }
    //At terminal leaf
    return n->get_id();
}

/// Get the leaf corresponding to a set index (like a vector look-up)
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::id_type
AdaptiveGrid<T,ValueTraits>::
getLeaf(const int& idx) const
{
    std::set<id_type, nvis::lexicographical_order>::iterator lit;
    if (idx >= (int)_leaves.size() || idx<0) {
        return id_type(0,0,0);
    }
    lit = _leaves.begin();
    std::advance(lit,idx);
    return (*lit);
}

///Return the id of a cell one depth down (or one node-level up) on the tree structure
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::id_type
AdaptiveGrid<T, ValueTraits>::
sub_id(const typename AdaptiveGrid<T, ValueTraits>::id_type& id) const
{
    return id_type(id[0]*2, id[1]*2, id[2]+1);
}


///Return the id of a cell at depth=depth containing cell id(i,j,d)
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::id_type
AdaptiveGrid<T, ValueTraits>::
raise_id(const typename AdaptiveGrid<T, ValueTraits>::id_type& id, unsigned int depth) const
{
    const int& d = id[2];
    if (!d) {
        return id;
    }
    assert((int)depth < d);
    int f = 1;
    f = f << (d-depth);
    return id_type(id[0]/f, id[1]/f, depth);
}

///Get the index of the starting grid (MxN) given an id (i,j,depth)
template<typename T, typename ValueTraits>
inline unsigned int AdaptiveGrid<T, ValueTraits>::
id_to_index(const typename AdaptiveGrid<T, ValueTraits>::id_type& id) const
{
    id_type top_id = raise_id(id, 0);
    return top_id[0] + _size[0]*top_id[1];
}

///Get a vertex (position vector) of a node given local_index i (0-LL,1-LR,2-UR,3-UL,4-C)
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::pos_type
AdaptiveGrid<T, ValueTraits>::
node_point(const typename AdaptiveGrid<T, ValueTraits>::id_type& id, int i) const
{
    bounds_type b = cellBounds(id);
    switch (i) {
        case 0:
            return b.min();
        case 1:
            return pos_type(b.max()[0], b.min()[1]);
        case 2:
            return b.max();
        case 3:
            return pos_type(b.min()[0], b.max()[1]);
        default:
            return b.center();
    }
}

///Return the list of cell data (pair: pos->parameter)
template<typename T, typename ValueTraits>
const std::list<typename AdaptiveGrid<T, ValueTraits>::data_type>&
AdaptiveGrid<T, ValueTraits>::
getCellData(const typename AdaptiveGrid<T, ValueTraits>::id_type& id) const
{
    typename cell_container_type::const_iterator it = _c_data.find(id);
    if (it == _c_data.end()) {
        throw std::runtime_error("unknown cell coordinates");
    }
    return it->second;
}

///Return a vector of all pos->value data stored in the grid
template<typename T, typename ValueTraits>
void AdaptiveGrid<T, ValueTraits>::
getData(std::vector<typename AdaptiveGrid<T, ValueTraits>::data_type>& data)
{
    data.clear();
    //Grab all values at vertexes
    typename vertex_container_type::const_iterator vDataIT;
    for(vDataIT=_v_data.begin(); vDataIT!=_v_data.end(); ++vDataIT) {
        data.push_back(
            data_type( posWithDepth(getVertex(vDataIT->first),vDataIT->first[2]),
                       vDataIT->second ) );
    }
    
    //For each leaf, grab cell data
    std::set<id_type, nvis::lexicographical_order>::iterator it;
    for(it=_leaves.begin(); it!=_leaves.end(); ++it) {
        //get cell data (list)
        cell_data_type cDataList = getCellData( *it );
        //Add to input data vector
        typename cell_data_type::iterator cdIter;
        for(cdIter = cDataList.begin(); cdIter != cDataList.end(); ++cdIter) {
            data.push_back( *cdIter );
        }
    }
}

///Get the parameter value at a vertex given an id(i,j,d)
template<typename T, typename ValueTraits>
const typename AdaptiveGrid<T, ValueTraits>::value_type&
AdaptiveGrid<T, ValueTraits>::
getVertexValue(const typename AdaptiveGrid<T, ValueTraits>::id_type& id) const
{
    typename vertex_container_type::const_iterator it = _v_data.find(id);
    if (it == _v_data.end()) {
        throw std::runtime_error("unknown vertex coordinates");
    }
    return it->second;
}

///Set the parameter value at a given vertex id(i,j,d)
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
setVertexValue(const typename AdaptiveGrid<T, ValueTraits>::id_type& id,
               const typename AdaptiveGrid<T, ValueTraits>::value_type& v)
{
    _v_data[id] = v;
    
}

///Insert a position->paramValue pair into the corresponding cell
template<typename T, typename ValueTraits>
inline typename AdaptiveGrid<T, ValueTraits>::id_type
AdaptiveGrid<T, ValueTraits>::
insertCellValue(const typename AdaptiveGrid<T, ValueTraits>::pos_type& x,
                const typename AdaptiveGrid<T, ValueTraits>::value_type& v,
                const int& depth)
{
    if (!_bounds.inside(x)) {
        throw std::runtime_error("invalid position");
    }
    
    pos_type tmp = (x-_bounds.min())/_step;
    int idx = floor(tmp[1])*_size[0] + floor(tmp[0]);
    
#if defined(HX_HAS_STD) || defined(C_0X)
    boost::shared_ptr<Node> n = _raster[idx];
#else
    std::shared_ptr<Node> n = _raster[idx];
#endif
    while (!n->is_leaf()) {
        const id_type& id = n->get_id();
        pos_type c = node_point(id, -1);
        if (nvis::all(x <= c)) {
            n = n->_child[0];
        } else if (nvis::all(x > c)) {
            n = n->_child[2];
        } else if (x[0] < c[0]) {
            n = n->_child[3];
        } else {
            n = n->_child[1];
        }
    }
    const id_type& id = n->get_id();
    _c_data[id].push_back(data_type(posWithDepth(x,depth), v));
    return id;
}

///Get all the ids of leaves in the tree  (All valid leaves after removeInvalidLeaves() is called)
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
getLeaves(std::vector<typename AdaptiveGrid<T, ValueTraits>::id_type>& leaves) const
{
    leaves.clear();
    std::copy(_leaves.begin(), _leaves.end(), std::back_inserter(leaves));
}

///Get all the undefined (i.e., ValueTraits::unknown) vertexes as indexes
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
getUndefinedVertexes(std::vector<typename AdaptiveGrid<T, ValueTraits>::id_type>& undef) const
{
    undef.clear();
    typename vertex_container_type::const_iterator it;
    const value_type unknown = ValueTraits::unknown;
    for (it=_v_data.begin() ; it!=_v_data.end() ; ++it) {
        if (it->second == unknown) {
            undef.push_back(it->first);
        }
    }
}

///Get cells that need info within
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
getEmptyCells(std::vector<typename AdaptiveGrid<T, ValueTraits>::id_type>& empty) const
{
    empty.clear();
    typename cell_container_type::const_iterator it;
    for (it=_c_data.begin(); it!= _c_data.end() ; ++it) {
        if (it->second.empty()) {
            empty.push_back(it->first);
        }
    }
}

template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
#if defined(HX_HAS_STD) || defined(C_0X)
refine(boost::shared_ptr<typename AdaptiveGrid<T, ValueTraits>::Node> n)
{
#else
refine(std::shared_ptr<typename AdaptiveGrid<T, ValueTraits>::Node> n)
{
#endif
    static const int shift [][2] = { {1,0}, {1,1}, {2,1}, {1,2}, {0,1} };
    
    n->split();
    
    // create (empty) entries for newly created cells
    typedef typename cell_container_type::iterator cc_iterator_type;
    cc_iterator_type sub_it[4];
    for (int i=0 ; i<4 ; ++i) {
        const id_type& id = n->get_child(i)->get_id();
        _leaves.insert(id);
        typename cell_container_type::value_type v(id, cell_data_type());
        std::pair<cc_iterator_type, bool> r = _c_data.insert(v);
        sub_it[i] = r.first;
        //Assume the cell is valid to start
        _c_valid.insert(std::pair<id_type,bool>(id, true));
    }
    _leaves.erase(n->get_id());
    
    // assign contents of old cell to subcells
    const id_type& cell_id = n->get_id();
    cell_data_type& data = _c_data[cell_id];
    const pos_type center = getVertex(n->get_id(-1));
    for (typename cell_data_type::const_iterator it=data.begin() ; it!=data.end() ; ++it) {
        const pos_type& x = it->first;
        if      (nvis::all(x <= center)) {
            sub_it[0]->second.push_back(*it);
        } else if (nvis::all(x > center)) {
            sub_it[2]->second.push_back(*it);
        } else if (x[0] <= center[0]) {
            sub_it[3]->second.push_back(*it);
        } else {
            sub_it[1]->second.push_back(*it);
        }
    }
    data.clear(); // empty parent cell
    
    // create (unknown) entries for newly created vertexes
    const id_type base_id = sub_id(cell_id);
    typedef typename vertex_container_type::iterator vc_iterator_type;
    const value_type unknown = ValueTraits::unknown;
    for (int i=0 ; i<5 ; ++i) {
        id_type id = base_id + id_type(shift[i][0], shift[i][1], 0);
        typename vertex_container_type::value_type v(id, unknown);
        // make sure not to overwrite existing vertexes, if any
        //std::pair<vc_iterator_type, bool> r = _v_data.insert(v);
        _v_data.insert(v);
    }
    // copy existing entries to next resolution level
    _v_data[base_id                 ] = _v_data[cell_id                 ];
    _v_data[base_id + id_type(2,0,0)] = _v_data[cell_id + id_type(1,0,0)];
    _v_data[base_id + id_type(2,2,0)] = _v_data[cell_id + id_type(1,1,0)];
    _v_data[base_id + id_type(0,2,0)] = _v_data[cell_id + id_type(0,1,0)];
}

template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
refine(const typename AdaptiveGrid<T, ValueTraits>::id_type& cell_id)
{
    // (i,j,d) -> (i/2, j/2, d-1) -> ... -> (i/2^d, j/2^d, 0)
    assert(_c_data.find(cell_id) != _c_data.end());
    id_type id = raise_id(cell_id, 0);
    const int& depth = cell_id[2];
    
    // Get the node for the depth=0 cell
#if defined(HX_HAS_STD) || defined(C_0X)
    boost::shared_ptr<Node> n = _raster[id_to_index(id)];
#else
    std::shared_ptr<Node> n = _raster[id_to_index(id)];
#endif
    
    //Travel down levels to find node to refine (through structure in _raster)
    for (int d=1 ; d<=depth; ++d) {
        id_type base_id   = sub_id(id); //Find ID at next level (lower-left)
        id_type actual_id = cell_id;
        if(d!=depth) {
            actual_id = raise_id(cell_id, d);
        }
        id_type diff      = actual_id-base_id;
        if      (nvis::all(diff == id_type(0,0,0))) {
            n = n->get_child(0);
        } else if (nvis::all(diff == id_type(1,0,0))) {
            n = n->get_child(1);
        } else if (nvis::all(diff == id_type(1,1,0))) {
            n = n->get_child(2);
        } else {
            n = n->get_child(3);
        }
        //Get the new node id
        id = n->_id;
    }
    //Run Node refine (split)
    refine(n);
}

/// Is this node valid (i.e., no "invalid" value)
template<typename T, typename ValueTraits>
inline const bool AdaptiveGrid<T, ValueTraits>::
valid(const typename AdaptiveGrid<T, ValueTraits>::id_type& node_id) const
{
    T cornerValue = getVertexValue(node_id);
    if ( cornerValue == ValueTraits::invalid ) {
        return false;
    }
    return true;
}
/// Is this node valid (i.e., no "invalid" value)
template<typename T, typename ValueTraits>
inline bool AdaptiveGrid<T, ValueTraits>::
valid(const typename AdaptiveGrid<T, ValueTraits>::id_type& node_id)
{
    T cornerValue = getVertexValue(node_id);
    if ( cornerValue == ValueTraits::invalid ) {
        return false;
    }
    return true;
}

// Is this cell valid
template<typename T, typename ValueTraits>
inline const bool AdaptiveGrid<T, ValueTraits>::
isCellValid(const typename AdaptiveGrid<T, ValueTraits>::id_type& cell_id) const
{

    typename validity_container_type::const_iterator it = _c_valid.find(cell_id);
    if (it == _c_valid.end()) {
        throw std::runtime_error("unknown cell coordinates");
    }
    return it->second;
}

// Set the cell validity
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
setCellValidity(const typename AdaptiveGrid<T, ValueTraits>::id_type& cell_id, const bool v)
{
    _c_valid[cell_id] = v;
}

// Update cell validity after information is added to the grid
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
updateCellValidity()
{
    //For all leaves:
    std::set<id_type, nvis::lexicographical_order>::iterator lit;
    for(lit=_leaves.begin(); lit!=_leaves.end(); ++lit) {
        bool v = allCornersValid( (*lit) );
        _c_valid[(*lit)] = v;
    }
}



template<typename T, typename ValueTraits>
inline const bool AdaptiveGrid<T, ValueTraits>::
allCornersValid(const typename AdaptiveGrid<T, ValueTraits>::id_type& cell_id) const
{

    T cornerValue = getVertexValue(cell_id);
    if ( cornerValue == ValueTraits::invalid ) {
        return false;
    }
    cornerValue  = getVertexValue(cell_id + id_type(1,0,0));
    if ( cornerValue == ValueTraits::invalid ) {
        return false;
    }
    cornerValue  = getVertexValue(cell_id + id_type(1,1,0));
    if ( cornerValue == ValueTraits::invalid ) {
        return false;
    }
    cornerValue  = getVertexValue(cell_id + id_type(0,1,0));
    if ( cornerValue == ValueTraits::invalid ) {
        return false;
    }
    
    return true;
}

// Are all corners of a cell invalid
template<typename T, typename ValueTraits>
inline const bool AdaptiveGrid<T, ValueTraits>::
allCornersInvalid(const typename AdaptiveGrid<T, ValueTraits>::id_type& cell_id) const
{

    T cornerValue = getVertexValue(cell_id);
    if ( cornerValue != ValueTraits::invalid ) {
        return false;
    }
    cornerValue  = getVertexValue(cell_id + id_type(1,0,0));
    if ( cornerValue != ValueTraits::invalid ) {
        return false;
    }
    cornerValue  = getVertexValue(cell_id + id_type(1,1,0));
    if ( cornerValue != ValueTraits::invalid ) {
        return false;
    }
    cornerValue  = getVertexValue(cell_id + id_type(0,1,0));
    if ( cornerValue != ValueTraits::invalid ) {
        return false;
    }
    
    return true;
}


///Remove the leaves with an "invalid" value at a corner ->pull from set
template<typename T, typename ValueTraits>
inline void AdaptiveGrid<T, ValueTraits>::
removeInvalidLeaves()
{
    //Loop through the leaves set to find invalid values
    std::vector<id_type> invalidLeaves;
    std::vector<id_type>::iterator invIT;
    std::set<id_type, nvis::lexicographical_order>::iterator setIT;
    
    //Need Special Integration Failure case
    
    //Remove any cell with at least 1 invalid corner or marked as invalid
    for(setIT=_leaves.begin(); setIT!=_leaves.end(); ++setIT) {
        if ( !allCornersValid( *setIT ) || (!isCellValid( *setIT )) ) {
            invalidLeaves.push_back( *setIT );
        }
    }
    //Loop through invalid values to erase from leaves set
    for(invIT=invalidLeaves.begin(); invIT!=invalidLeaves.end(); ++invIT) {
        _leaves.erase( *invIT );
    }
}

///Get the leaf indexes of all the corners in the terminating leaves
template<typename T, typename ValueTraits>
void AdaptiveGrid<T,ValueTraits>::
getGridNodes(std::vector<typename AdaptiveGrid<T,ValueTraits>::id_type>& cornerIDs) const
{
    cornerIDs.clear();
    std::vector<id_type> leaves;
    typename std::vector<id_type>::iterator it;
    getLeaves(leaves);
    //Get an isolated corner id using a set
    std::set<id_type,ID_Compare> cornerSet;
    for(it=leaves.begin(); it!=leaves.end(); ++it) {
        //Try to insert all 4 corners into the node set
        cornerSet.insert( *it );
        cornerSet.insert( *it + id_type(1,0,0) );
        cornerSet.insert( *it + id_type(1,1,0) );
        cornerSet.insert( *it + id_type(0,1,0) );
    }
    //Copy the set to the output vector
    std::copy(cornerSet.begin(),cornerSet.end(),std::back_inserter(cornerIDs));
}

///Get the positions of all the corners in the terminating leaves
template<typename T, typename ValueTraits>
void AdaptiveGrid<T,ValueTraits>::
getGridNodes(std::vector<typename AdaptiveGrid<T,ValueTraits>::pos_type>& corners) const
{
    corners.clear();
    std::vector<id_type> cornerIDs;
    getGridNodes(cornerIDs);
    for(int i=0; i<(int)cornerIDs.size(); i++) {
        corners.push_back( getVertex(cornerIDs[i]) );
    }
}

} // namespace xavier

#endif
