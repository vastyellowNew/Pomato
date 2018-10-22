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


/*
 *  ManifoldData.h - a Header for stable/unstable manifold data class
 *
 *  Author: Wayne Schlei
 *          Purdue University
 *
 *  Date:  03/2/2015
 *
 */

#ifndef MANIFOLD_DATA_H
#define MANIFOLD_DATA_H

//Avizo & GUI Headers
#include <hxcore/HxData.h>
#include <hxcore/HxPortInfo.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxcore/HxPortSeparator.h>
#include <hxcore/HxConnection.h>
#include <hxcore/HxPortFilename.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortButtonList.h>
#include <hxcore/HxPortTabBar.h>
#include <hxcore/HxPortRadioBox.h>
#include <hxcore/HxPortToggleList.h>
#include <hxcore/HxPortText.h>
#include <hxtime/HxPortTime.h>

//Containers
#include <QMap>
#include <vector>
#include <list>

//API
#include <math/fixed_vector.hpp>
#include <maps/definitions.hpp>
#include <maps/map_analysis_param.hpp>
#include <maps/fixpoints.hpp>
#include <cr3bp/cr3bp.hpp>
#include <topology/ManifoldDataStorage.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/SortableReturnData.hpp>
#include <pmate/ManifoldData.hpp>
#include <poincare/poincareAPI.h>



using namespace std;

/** Manifold data storage class
*/
class POINCARE_API ManifoldData : public HxData
{
    HX_HEADER(ManifoldData);

  public :
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
    ManifoldData();


    //Ports
    HxConnection fpxDataConnection;
    HxPortInfo portInfo;
    HxPortInfo portSeps;
    /// FixedPointData filename 
    HxPortFilename portFPDataFile;

    /// Compute: modifies ports and will update time values
    virtual void compute();


    /// Save certain port states
    //virtual void savePorts(FILE* fp);

    /// Tcl Interpreter for functions
    //virtual int parse(Tcl_Interp* t, int argc, char **argv);

    //The Old way...
    /*///Manifold Data container
    topology::ManifoldDataStorage manifoldData;

    /// Get number of separtrices
    int getNumManifolds() const
    { return (int) manifoldData.separatrices.size(); }

    /// Get number of fixed points
    int getNumOrbits() const
    { return (int) manifoldData.chains.size(); }*/
    
    /// Manifold Data Object - (note: blank read, missing pointers to be assigned later).
    ManifoldDataType theManifoldData;  
    
    
    /// Number of Manifolds within data
    int getNumManifolds() const
    { return (int) theManifoldData.mapManifolds.size(); }

  protected :
    /// Destructor
    ~ManifoldData();

};

/// A class for assigning a render priority for ManifoldData
class POINCARE_API ManifoldDataRenderPriority
{
  public :
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

#endif // ManifoldData_H
