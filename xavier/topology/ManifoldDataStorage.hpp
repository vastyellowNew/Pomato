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


//Manifold Data Stroage structure
// Author: Wayne Schlei
// Date: 2/7/2013
#ifndef __MANIFOLD_DATA_STORAGE_HPP__
#define __MANIFOLD_DATA_STORAGE_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <maps/fixpoints.hpp>
#include <topology/invariant_manifold.hpp>

using namespace xavier;

namespace topology {

typedef std::vector<fp_chain>          chain_type;
typedef std::vector<Separatrix>        sep_type;

///Data structure for keeping data available without recomputing (Global Variable Hack)
struct ManifoldDataStorage {
    ManifoldDataStorage() {}
    
    ~ManifoldDataStorage()
    {
        clear();
    }
    void clear()
    {
        chains.clear();
        separatrices.clear();
        broken_manifolds.clear();
    }
    
    bool isEmpty()
    {
        if ((int) chains.size() == 0) {
            return true;
        }
        return false;
    }
    
    bool isComplete()
    {
        if ((int) separatrices.size() == 0) {
            return false;
        } else {
            return true;
        }
    }
    chain_type            chains;
    sep_type              separatrices;
    std::vector< std::vector<nvis::vec2> > broken_manifolds;
};

} //orbital

#endif //__MANIFOLD_DATA_STORAGE_HPP__
