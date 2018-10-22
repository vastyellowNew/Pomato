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


#ifndef __DISPLAY_TOPOLOGY_HPP__
#define __DISPLAY_TOPOLOGY_HPP__

#include <limits>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <vector>

// nvis
#include <math/fixed_vector.hpp>
// xavier
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/planar_section.hpp>
#include <cr3bp/multipleAngleTracker.hpp>
#include <maps/definitions.hpp>
#include <pmate/ManifoldData.hpp>


/// A class for assigning a render priority for ManifoldData
class ManifoldDataRenderPriority
{
public:
    typedef xavier::fixpoint                                          FPType;
    typedef std::vector<xavier::fixpoint>                             FPChainType;
    typedef pmate::FixedPointData<xavier::fixpoint>                   FPData;
    typedef orbital::cr3bp                                            RHStype; //Hard-coded (template?)
    static const int numSing = RHStype::numSingularities;
    typedef nvis::fixed_vector<double,numSing+1>                      ExtendedMapDataVec;
    typedef nvis::fixed_vector<double,2>                              VecType;
    typedef topology::SortableReturnData<VecType,ExtendedMapDataVec>  SortableData;
    typedef xavier::EdgeRotationFailure<VecType>                      MapDiscont;
    typedef pmate::ManifoldData<SortableData,FPType>                  ManifoldDataType;
    typedef ManifoldDataType::ManifoldType                            ManifoldType;
    typedef ManifoldType::ManifoldSeg                                 ManifoldSeg;

    /// Constructor
    ManifoldDataRenderPriority(ManifoldDataType& mData) :
      data(mData) {}

    /// Comparison for priority_queue
    bool operator()(const int& mID1, const int& mID2) const {
      //First, test the orbit ID (long, highly unstable orbits have highest priority)
      if(data.mapManifolds[mID1].fpdOrbitIdx != data.mapManifolds[mID2].fpdOrbitIdx) {
        return (data.mapManifolds[mID1].fpdOrbitIdx < data.mapManifolds[mID2].fpdOrbitIdx);
      }
      //Next, test the fixed point index
      if(data.mapManifolds[mID1].fpdPointIdx != data.mapManifolds[mID2].fpdPointIdx) {
        return (data.mapManifolds[mID1].fpdPointIdx < data.mapManifolds[mID2].fpdPointIdx);
      }
      //Next, test the type (stable < unstable)
      if(data.mapManifolds[mID1].isForward() != data.mapManifolds[mID2].isForward()) {
        return (data.mapManifolds[mID1].isForward() < data.mapManifolds[mID2].isForward());
      }
      //And finally, manifold ID
      return (mID1 < mID2);
    }

    /// Reference to ManifoldData
    ManifoldDataType &data;
};

#endif
