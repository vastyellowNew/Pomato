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


// Sortable Return Data : Useful for temporary storage in EdgeRotation map evaluations
// Author:  Wayne Schlei (Purdue University)

#ifndef __SORTABLE_RETURN_DATA_HPP
#define __SORTABLE_RETURN_DATA_HPP

#include <vector>
#include <iostream>
#include <math/fixed_vector.hpp>
#include <maps/definitions.hpp>
#include <maps/misc.hpp>
#include <maps/index.hpp>

namespace topology {

/// Structure for storing pts and their returns in a set
template<class VEC,class MAPDATA>
struct SortableReturnData {
    typedef SortableReturnData<VEC,MAPDATA>  SelfType;
    typedef VEC                              VecType;
    typedef MAPDATA                          MapDataType;
    
    SortableReturnData() : x0(0.,0.), returns(), data(), isValid(false)
    {}
    SortableReturnData(const VEC& ic) : x0(ic), returns(), data(), isValid(true)
    {}
    SortableReturnData(const VEC& ic, const std::vector<VEC>& iterates) :
        x0(ic),
        returns(iterates),
        data(), //Note: won't be the right size
        isValid(true)
    {}
    SortableReturnData(const VEC& ic, const std::vector<VEC>& iterates,
                       const std::vector<MAPDATA>& mapData) :
        x0(ic),
        returns(iterates),
        data(mapData),
        isValid(true)
    {}
    /// For creating points from a dataGrid Node (dataField->orbit_data->steps)
    template<class DATASET, class IVEC>
    SortableReturnData(const VEC& ic, const DATASET& dataField, const IVEC& index, const int& p)
    {
        //Note:  raster_data stores first point in orbit_data, but adaptive_grid does NOT
        //x0 = dataField(index).steps[0]; //raster_data
        x0 = ic; //adaptive_grid
        isValid = false;
        if (std::abs(p)<(int) dataField(index).steps.size()) {
            isValid = true;
            //for (int i=1;i<(p+1);i++) returns.push_back(dataField(index).steps[i]);//raster_data
            for (int i=0; i<std::abs(p); i++) {
                returns.push_back(dataField(index).steps[i]);    //adpative_grid
            }
        }
        
    }
    /// Check if p exists
    inline bool isThere(const int& p) const
    {
        return (isValid && std::abs(p)<=(int) returns.size());
    }
    /// Get the pth return
    VEC getReturn(const int& p) const
    {
        if (isThere(p)) {
            return returns[std::abs(p)-1];
        } else {
            //This was invalid! ->Force this case to never be called!
            return VEC(0,0);
        }
    }
    
    /// Get the pth return data values
    MAPDATA getDataAtReturn(const int& p) const
    {
        if (isThere(p)) {
            return data[std::abs(p)-1];
        } else {
            //This was invalid! ->Force this case to never be called!
            return MAPDATA(0);
        }
    }
    
    ///Overloaded < operator for working with std::set<>
    bool operator<(const SelfType& other) const
    {
        nvis::lexicographical_order compare;
        return compare(this->x0,other.x0);
    }
    
    /*///Print data
    template<class VEC,class MAPDATA>
    std::ostream& operator<<(std::&ostream os,
                           const SortableReturnData<VEC,MAPDATA>& rdata) {
      os << "ic: " << rdata.x0 << '\n';
      os << "returns:\n";
      for (int i=0 ; i<(int)rdata.returns.size() ; ++i) {
          os << i << ": " << rdata.returns[i] << "\t->\t" << rdata.data[i] << '\n';
      }
      return os;
    }*/
    
    
    //Members:
    ///IC
    VEC x0;
    ///Returns
    std::vector< VEC > returns;
    ///Map Data - any other necessary data per return
    std::vector< MAPDATA > data;
    ///isValid
    bool isValid;
    
    
};






} //end xavier

#endif