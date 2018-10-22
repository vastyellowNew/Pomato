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


/**  AdaptiveSampler
  *    Author : Wayne Schlei
  *    Purpose:  Class object that adaptively samples a Poincare section
  *    using the adaptive grid
  */

#ifndef ADAPTIVE_SAMPLER_HPP
#define ADAPTIVE_SAMPLER_HPP


//Standard Libs
#include <vector>
#include <ostream>
#include <list>
#include <map>
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/array.hpp>

//API - chris
#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

//API - xavier
#include <data/grid.hpp>
#include <data/adaptive_grid.hpp>
#include <data/cell.hpp>
#include <data/edge.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/definitions.hpp>
#include <maps/adaptiveGridMapInterface.hpp>
#include <maps/map_analysis.hpp>

//Point cloud file output
#include <cr3bp/psiWrite.hpp>


#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;

namespace pmate {
/// Class for adaptively sampling a Poincare map with respect to smoothness in winding number
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
class AdaptiveSampler {
public :
    typedef typename PMAP::value_type            value_type;
    typedef typename PMAP::rhs_type              rhs_type;
    typedef typename PMAP::section_type          section_type;
    typedef AdaptiveGrid<WVEC,WTRAITS>           AdaptiveGridType;
    typedef AdaptiveGridNodeData<DATATYPE,AdaptiveGridType,PARAM> AGNodeData;
    static const int N = rhs_type::dimension; //EOM dimension
    static const int M = section_type::local_dimension; //Section dimension
    typedef grid<value_type,M>                   PlanarGridType;
    typedef typename PlanarGridType::ivec_type   DimType;
    typedef typename PlanarGridType::bounds_type BoundsType;
    typedef typename AdaptiveGridType::data_type DataPairType;
    typedef typename AdaptiveGridType::id_type   IdType;
    
    ///Constructor
    AdaptiveSampler(const BoundsType& bbox, const DimType& res, PMAP& mapIn, PARAM& params) :
        gridBounds(bbox),
        gridRes(res),
        needsUpdate(true),
        theMap(&mapIn),
        theMapParams(&params),
        thePlanarGrid(gridRes,gridBounds),
        theGrid(bbox,res[0]-1,res[1]-1),
        theData(theGrid),
        cellChecker(theGrid,params.winding_convexity_tols,params.winding_cell_maxDist)
    {}
    
    
    /// Compute the sampling to given maximum depth
    void sample();
    
    /// Write the grid data to file
    void write(const char* filename) const;
    
    /// Write internal grid data as puncture plot representation with winding numbers (for Avizo)
    void writeToPSI(const char* psiFile) const;
    
    /// Read grid data from file
    void read(const char* filename);
    
    /// Remove the invalid leaves
    void removeInvalidLeaves();
    
    /// Reset the data - Changes the map params
    void reset();
    /// Reset from change in Map Parameters
    void resetFromParams();
    
    /// Does the sampling data need updating
    bool isUpdateRequired() const
    {
        return needsUpdate;
    }
    
    //Functions to set Computation parameters from outside
    /// Set the grid bounding box
    void setGridBounds(const BoundsType& bbox)
    {
        gridBounds = bbox;
        reset();
    }
    /// Set the grid initial resolution
    void setGridRes(const DimType& res)
    {
        gridRes = res;
        reset();
    }
    /// Set the maximum depth level (Default is 3)
    void setMaxDepth(const int d)
    {
        maxDepth = d;
        reset();
    }
    /// Set the map parameters object
    void setMapParams(PARAM& params)
    {
        theMapParams = &params;
        resetFromParams();
    }
    
    /// Set the winding number tolerances (for cell convexity checks)
    void setWindingTols(const WVEC& tols)
    {
        for(int i=0; i<(int)tols.size(); i++) {
            theMapParams->winding_convexity_tols[i] = tols[i];
        }
        reset();
    }
    
    /// Set the winding number maximum distances (for cell convexity checks)
    void setWindingMaxDist(const WVEC& mDist)
    {
        for(int i=0; i<(int)mDist.size(); i++) {
            theMapParams->winding_cell_maxDist[i] = mDist[i];
        }
        reset();
    }
    
    /// Convexity Functor
    CONVEXITY cellChecker;
    
    /// Planar grid (to access grid points easily with thePlanarGrid(i,j) )
    PlanarGridType thePlanarGrid;
    /// The adaptive grid object
    AdaptiveGridType theGrid;
    /// The data corresponding to the adaptive grid
    AGNodeData theData;
    
private :
    PMAP* theMap;
    PARAM* theMapParams;
    BoundsType gridBounds;
    DimType    gridRes;
    int maxDepth;
    bool needsUpdate;
    
};

/// Compute the sampling to given maximum depth
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
sample()
{
    //Perform adaptive sampling with respect to winding number
    //  (from maps/map_analysis.hpp)
    adaptiveSampleRaster<PMAP,TRACKER,AdaptiveGridType,AGNodeData,WVEC,CONVEXITY>
    (theGrid,theData,cellChecker,(*theMap),(*theMapParams),true);
    
    //Ready for additional processing
    needsUpdate = false;
}

/// Write the data to file
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
write(const char* filename) const
{
    //Make sure the data is there
    if (needsUpdate) {
        throw("AdaptiveSampler write() Error:  requesting data that is unavailable");
        return;
    }
    
    //Write to file using AdaptiveGridNodeData<>
    DimType tempGridRes = gridRes;
    for(int j=0; j<(int)tempGridRes.size(); j++) {
        tempGridRes[j] -= 1;
    }
    theData.writeDataToFile( filename, gridBounds, tempGridRes,
                             theMap->rhs().desired_hamiltonian(),(*theMapParams), maxDepth);
}

/// Read from file
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
read(const char* filename)
{
    //Note:  the read can change parameters & C value
    double ham = theMap->rhs().desired_hamiltonian();
    // Call the "read" function and set DataSet & Adaptive Grid
    DimType baseDim;
    adaptiveRasterFromFile<PMAP,AdaptiveGridType,AGNodeData,WVEC,CONVEXITY>
    (filename, theGrid, theData, gridBounds, baseDim, maxDepth,
     cellChecker, (*theMap), ham, (*theMapParams), true);
     
    gridRes = baseDim + DimType(1); //Add 1 to all entry to accomodate shift
    // Set the new pmap components
    rhs_type* theNewRHS = new rhs_type(theMap->rhs());
    section_type* theNewSection = new section_type(theMap->section());
    theNewRHS->setJacobi( ham ); //**CR3BP**//
    theNewSection->set_rhs( *theNewRHS );
    theMap->setNewComponents( *theNewRHS, *theNewSection );
    thePlanarGrid.setNewBounds(gridRes,gridBounds);
    delete theNewRHS;
    delete theNewSection;
    
    //Ready to do other stuff
    needsUpdate = false;
}

/// Reset from new Map Parameters
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
resetFromParams()
{
    //The map parameters are now holding
    gridBounds = theMapParams->bounds;
    gridRes = theMapParams->resolution;
    thePlanarGrid.setNewBounds(gridRes,gridBounds);
    theGrid.reset(gridBounds,gridRes[0]-1,gridRes[1]-1);
    theData.nodeDataMap.clear(); //clear the nodes
    //Max Depth
    maxDepth = theMapParams->max_depth;
    OnEdgeCompare::maxDepth = maxDepth;
    ID_Compare::maxDepth = maxDepth;
    //Convexity checker
    cellChecker.tols = theMapParams->winding_convexity_tols;
    cellChecker.maxDists = theMapParams->winding_cell_maxDist;
    
    //We will need to recompute the grid if this is called
    needsUpdate = true;
    
}

/// Reset the data - also changes Map Parameters
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
reset()
{
    thePlanarGrid.setNewBounds(gridRes,gridBounds);
    theGrid.reset(gridBounds,gridRes[0]-1,gridRes[1]-1);
    theData.nodeDataMap.clear(); //clear the nodes
    theMapParams->bounds = gridBounds;
    theMapParams->resolution = gridRes;
    //Map Params
    OnEdgeCompare::maxDepth = maxDepth;
    ID_Compare::maxDepth = maxDepth;
    theMapParams->max_depth = maxDepth;
    //Convexity checker
    cellChecker.tols = theMapParams->winding_convexity_tols;
    cellChecker.maxDists = theMapParams->winding_cell_maxDist;
    
    //Assumed that these cellChecker elements are set outside of this:
    //testCellChecker->setMu(sys->mup);
    //double lstar = sys->lstar;
    //cellChecker->setBodyRadii(sys->radius1 / lstar, sys->radius2 / lstar);
    
    //We will need to recompute the grid if this is called
    needsUpdate = true;
    
}

/// Remove the invalid leaves
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
removeInvalidLeaves()
{
    //Make sure the data is there
    if (needsUpdate) {
        throw("AdaptiveSampler removeInvalidLeaves() Error:  requesting data that is unavailable");
        return;
    }
    
    //Remove invalids if data is already computed
    theGrid.removeInvalidLeaves();
}

/// Write internal grid data as puncture plot representation with winding numbers (for Avizo)
template <typename PMAP, typename DATATYPE,
          typename TRACKER, typename WVEC, typename WTRAITS,
          typename CONVEXITY, typename PARAM>
void AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeToPSI(const char* psiFile) const
{
    //Make sure the data is there
    if (needsUpdate) {
        throw("AdaptiveSampler writeToPSI() Error:  requesting data that is unavailable");
        return;
    }
    
    //Point data cloud (for Avizo input)
    std::vector<PointPSI> points;
    std::vector<int> ids;
    std::vector< std::vector<double> > totalData;
    std::vector< std::string > psi_labels, psi_symbols;
    psi_labels.push_back( std::string("Winding Number x-xd") );
    psi_labels.push_back( std::string("Winding Number x-yd") );
    psi_labels.push_back( std::string("Winding Number xd-yd") );
    psi_labels.push_back( std::string("Computation Depth") );
    //psi_labels.push_back( std::string("Poloidal Period Approx") );
    //psi_labels.push_back( std::string("Crossing Number") );
    //psi_labels.push_back( std::string("Time") );
    psi_symbols.push_back( std::string("w"));
    psi_symbols.push_back( std::string("u"));
    psi_symbols.push_back( std::string("v"));
    psi_symbols.push_back( std::string("d"));
    //psi_symbols.push_back( std::string("p"));
    //psi_symbols.push_back( std::string("c"));
    //psi_symbols.push_back( std::string("t"));
    
    //Organize the output data into the appropriate vectors
    std::vector<double> wNums, uNums, vNums, dNums;
    std::vector<DataPairType>  allCellDataPairs;
    theGrid.getData( allCellDataPairs );
    int numPoints = (int) allCellDataPairs.size();
    for(int i=0; i<numPoints; i++) {
        nvis::vec2 pt = allCellDataPairs[i].first;
        points.push_back( PointPSI(pt[0],0.0,pt[1]) ); //unscaled
        ids.push_back(i);
        wNums.push_back( allCellDataPairs[i].second[0] );
        uNums.push_back( allCellDataPairs[i].second[1] );
        vNums.push_back( allCellDataPairs[i].second[2] );
        dNums.push_back( (double) allCellDataPairs[i].first.depth );
    }
    
    //Write information to PSI file
    std::vector< std::vector<double> > dataVecs;
    dataVecs.push_back(wNums);
    dataVecs.push_back(uNums);
    dataVecs.push_back(vNums);
    dataVecs.push_back(dNums);
    //dataVecs.push_back(pNums);
    //dataVecs.push_back(cNums);
    //dataVecs.push_back(tVec);
    bool ok = PSI_WriteToFile(points,ids,dataVecs,psi_labels,psi_symbols,psiFile);
    if(!ok) {
        throw("AdaptiveSampler::Error in writing to PSI file!");
    }
    
}

}

#endif