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


/** Adaptive Grid Interface objects for Topology analysis calls
 *  Author: Wayne Schlei (Purdue University)
 */
#ifndef __ADAPTIVE_GRID_MAP_INTERFACE_HPP
#define __ADAPTIVE_GRID_MAP_INTERFACE_HPP

#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math/fixed_vector.hpp>
#include <data/adaptive_grid.hpp>


namespace xavier {

/** Data structure holding a windingNumber<->OrbitReturns pair
 */
template<typename WINDING_VEC, typename POS, typename IDTYPE>
struct OrbitWindingPair {
public:
    OrbitWindingPair() : seed_id(IDTYPE(0))
    {
        orbit.first = WINDING_VEC(DefaultValueTraits::invalid);
        orbit.second.clear();
    }
    OrbitWindingPair(const IDTYPE& id) : seed_id(id)
    {
        orbit.first = WINDING_VEC(DefaultValueTraits::invalid);
        orbit.second.clear();
    }
    IDTYPE seed_id;
    std::pair< WINDING_VEC, std::vector<POS> > orbit;
    bool isICValid;
};

/** Data structure that holds the node->orbit_data interface for adaptive grid
 *  - This is the DATASET type in adaptiveSampleRaster().  [DATATYPE = orbit_data]
 *  - Note:  Data is stored at id it first appears as in the tree, so (4,4,4)
 *    is actually not in the nodeDataMap, but (1,1,2) is it's first appearance
 *    index.
 */
template<typename DATATYPE, typename GRID, typename PARAMS>
class AdaptiveGridNodeData {
public :
    typedef typename GRID::pos_type           vec_type;
    typedef typename GRID::id_type            id_type;
    typedef typename GRID::value_type         value_type;
    typedef typename GRID::bounds_type        bounds_type;
    typedef DATATYPE                          data_type;
    typedef std::map<id_type, DATATYPE, nvis::lexicographical_order>        DataMap;
    typedef typename DataMap::iterator        DataMapIterator;
    typedef typename DataMap::const_iterator  DataMapConstIterator;
    
    //Invalid point runtime error
    struct InvalidPointError : public std::runtime_error {
        InvalidPointError() : std::runtime_error("Invalid Point") {}
    };
    
    ///Constructor
    AdaptiveGridNodeData(const GRID& g) : theGrid(&g) {}
    
    /// Add data to map
    void addData(const id_type& id, const DATATYPE& data)
    {
        nodeDataMap.insert( std::pair<id_type,DATATYPE>(id,data) );
    }
    
    /// Check to see if this is in the grid
    bool inGrid(const id_type& id) const
    {
        DATATYPE v;
        try {
            v = (*this)(id); //Call Operator()
        } catch(...) {
            return false;
        }
        return true;
    }
    
    ///Is this point valid
    const bool valid(const id_type& id) const
    {
        return theGrid->valid(id);
    }
    bool valid(const id_type& id)
    {
        return theGrid->valid(id);
    }
    
    ///Cell spacing
    vec_type getSpacing(int depth=0)
    {
        return theGrid->spacing(depth);
    }
    vec_type getSpacing(int depth=0) const
    {
        return theGrid->spacing(depth);
    }
    
    ///Get vertex of a grid node
    vec_type getVertex(const id_type& id)
    {
        return theGrid->getVertex(id);
    }
    ///Get vertex of a grid node
    vec_type getVertex(const id_type& id) const
    {
        return theGrid->getVertex(id);
    }
    
    ///Get the value (winding number) at a vertex
    const value_type& getVertexValue(const id_type& id) const
    {
        return theGrid->getVertexValue(id);
    }
    
    ///Operator returns the data structure at a given point
    DATATYPE& operator()(const id_type& id)
    {
        typename DataMap::iterator mit;
        //See if it's in the data map
        mit = nodeDataMap.find(id);
        if (mit != nodeDataMap.end()) {
            return mit->second; //It's in there
        } else if (mit == nodeDataMap.end() && id[2]!=0) {
            //Counter a call to some high depth duplicate
            //e.g., (2,2,4) may not exist but (1,1,2) does
            //Recursively track up the tree to find the first-appearance index
            for(int d=id[2]-1; d>=0; d--) {
                id_type temp_id = theGrid->raise_id(id,d);
                mit = nodeDataMap.find(temp_id);
                if ( mit != nodeDataMap.end()) {
                    break;    //We found the point
                }
                if ( (temp_id[0] % 2 != 0) || (temp_id[1] % 2 != 0) ||(d==0) ) {
                    //We have nowhere left to go since we are at the origin point
                    throw InvalidPointError();
                }
            }
        } else if (id[2]==0) {
            throw InvalidPointError(); //It's not here
        } else {
            throw ("Unknown Error in Data lookup");
        }
        //We found it if it makes it here
        return mit->second;
    }
    ///Const operator returns the data structure at a given point
    const DATATYPE& operator()(const id_type& id) const
    {
        //See if this is in the map
        typename DataMap::const_iterator mit = nodeDataMap.find(id);
        if (mit != nodeDataMap.end()) {
            return mit->second; //It's in there
        } else if (mit == nodeDataMap.end() && id[2]!=0) {
            //Counter a call to some high depth duplicate
            //e.g., (2,2,4) may not exist but (1,1,2) does
            //Recursively track up the tree to find the first-appearance index
            for(int d=id[2]-1; d>=0; d--) {
                id_type temp_id = theGrid->raise_id(id,d);
                mit = nodeDataMap.find(temp_id);
                if ( mit != nodeDataMap.end()) {
                    break;    //We found the point
                }
                if ( (temp_id[0] % 2 != 0) || (temp_id[1] % 2 != 0) ||(d==0) ) {
                    //We have nowhere left to go since we are at the origin point
                    throw InvalidPointError();
                }
            }
        } else if (id[2]==0) {
            throw InvalidPointError(); //It's not here
        } else {
            throw ("Unknown Error in Data lookup");
        }
        //We found it if it makes it here
        return mit->second;
        
    }
    
    /// Write Node Data to file with other parameters for use of reconstructing a grid
    void writeDataToFile(const char* filename, const bounds_type& bbox, const nvis::ivec2& dim0,
                         const double& ham, const PARAMS& mapParams, const int& maxDepth) const;
                         
    /// Read Node Data from file -> Sets data map and can be used to refill AdaptiveGrid
    void readData(const char* filename, GRID& newGrid, bounds_type& bbox, nvis::ivec2& dim0, double& ham, PARAMS& mapParams);
    
    /// Get the terminating leaves in the adaptive grid
    void getLeaves(std::vector<id_type>& leaves) const
    {
        theGrid->getLeaves(leaves);
    }
    
    ///Map container between ID and data
    std::map<id_type, DATATYPE, nvis::lexicographical_order> nodeDataMap;
    
private :
    /// Adaptive grid pointer
    const GRID* theGrid;
    mutable bool _verbose;
};

} //end xavier

/// Write Node Data to file with other parameters for use of reconstructing a grid
template<typename DATATYPE, typename GRID, typename PARAMS>
void xavier::AdaptiveGridNodeData<DATATYPE,GRID,PARAMS>::
writeDataToFile(const char* filename,
                const typename GRID::bounds_type& bbox,//Grid Bounds on section
                const nvis::ivec2& dim0,               //Depth=0 dimension
                const double& ham,                     //Hamiltonian Value (constant)
                const PARAMS& mapParams,               //Map parameters
                const int& maxDepth                    //Maximum depth
               ) const
{
    //Open the file
    FILE* f = fopen(filename, "w");
    
    if (!f) {
        std::cerr << "Problem occurred in opening write file!\n";
        throw("Problem opening file for writing...");
        return;
    }
    
    //Write a header
    fprintf(f, "# Adaptive Grid Node Data\n");
    
    //Write information regarding grid
    int spaceDim = 2; //Fixed for now
    const DATATYPE& firstData = nodeDataMap.begin()->second;
    int wDim = (int) firstData.wn.size(); //Winding number size
    int numDataNodes = (int) nodeDataMap.size(); //Number of data nodes
    fprintf(f,"%d %d %d %d %d %d\n",spaceDim,dim0[0],dim0[1],maxDepth,numDataNodes,wDim); //Dims
    fprintf(f,"%.15f %f %f %f %f\n", ham,bbox.min()[0],bbox.min()[1],bbox.max()[0],bbox.max()[1]);
    //Tolerances
    for(int j=0; j<wDim; j++) {
        fprintf(f,"%f ", mapParams.winding_convexity_tols[j]);
    }
    fprintf(f,"\n");
    //Cell Max Distances
    for(int j=0; j<wDim; j++) {
        fprintf(f,"%f ", mapParams.winding_cell_maxDist[j]);
    }
    fprintf(f,"\n");
    
    
    //Iterate through all data
    typename std::map<id_type, DATATYPE, nvis::lexicographical_order>::const_iterator mapit;
    for (mapit=nodeDataMap.begin(); mapit!=nodeDataMap.end(); mapit++) {
        const id_type id = mapit->first;
        const DATATYPE& data = mapit->second;
        int numSteps = (int) data.steps.size();
        fprintf(f,"%d %d %d %d\n", id[0],id[1],id[2],numSteps);
        //HARDCODED for CR3BP maps!!!!
        //Print the points to file
        for(int i=0; i<numSteps; i++) {
            fprintf(f,"%.15f %.15f \n",data.steps[i][0],data.steps[i][1]);
        }
        //Print the winding numbers and isICValid
        int isICValid = (data.isICValid)? 1 : 0;
        fprintf(f,"%.15f %.15f %.15f %d\n", data.wn[0], data.wn[1], data.wn[2], isICValid);
        
    }
    
    //Close file
    fclose(f);
    
}


/// Read Node Data from file, also set some parameters of the grid
template<typename DATATYPE, typename GRID, typename PARAMS>
void xavier::AdaptiveGridNodeData<DATATYPE,GRID,PARAMS>::
readData(    const char* filename,
             GRID& newGrid,                //The adaptive grid (which is reset)
             typename GRID::bounds_type& bbox,    //Grid Bounds on section
             nvis::ivec2& dim0,            //Depth=0 dimension
             double& ham,                  //Hamiltonian value loaded from file
             PARAMS& mapParams             //Map parameters - will be modified on read()
        )
{
    //Open the File
    FILE* f = fopen(filename, "r"); //Read only
    
    if (!f) {
        std::cerr << "Error opening file!\n";
        throw("Cannot open file for read...");
        return;
    }
    
    //Buffer for skipping
    char buf[80];
    //Read header
    fgets(buf, 80, f); //Skip first line
    
    //Read the dimensions and bounding box
    int spaceDim, numDataNodes, wDim;
    int maxDepth,numIters;
    fscanf(f,"%d %d %d %d %d %d", &spaceDim, &(dim0[0]), &(dim0[1]), &maxDepth, &numDataNodes, &wDim);
    nvis::vec2 minBB, maxBB;
    fscanf(f,"%lf %lf %lf %lf %lf", &ham, &(minBB[0]),&(minBB[1]),&(maxBB[0]),&(maxBB[1]));
    bbox = typename GRID::bounds_type(minBB,maxBB);
    
    //Load the winding number tolerances
    for(int j=0; j<wDim; j++) {
        fscanf(f,"%lf ", &(mapParams.winding_convexity_tols[j]) );
    }
    //Load the winding number Max distances
    for(int j=0; j<wDim; j++) {
        fscanf(f,"%lf ", &(mapParams.winding_cell_maxDist[j]) );
    }
    
    //Check that file is ok
    if (spaceDim != 2) {
        std::cerr << "Error:  Hard-coded to 2D section space...\n";
        throw("Cannot support section space dimension");
        return;
    }
    
    
    //Reset the grid to the new parameters read from file
    newGrid.reset(bbox,dim0[0],dim0[1]);
    theGrid = &newGrid; //Assign const GRID * theGrid;
    
    //Assign data to nodeDataMap as we read
    nodeDataMap.clear(); //Start fresh
    numIters = 0;
    //HARD-CODED to CR3BP maps!!!!
    for (int i=0; i<numDataNodes; i++) {
        //Get the ID & number of steps
        int numSteps = 0;
        id_type id;
        fscanf(f,"%d %d %d %d", &(id[0]),&(id[1]),&(id[2]),&numSteps);
        if(numSteps > numIters) {
            numIters = numSteps;
        }
        DATATYPE dataValue;
        //Read data.steps in
        for(int j=0; j<numSteps; j++) {
            nvis::vec2 pos;
            fscanf(f,"%lf %lf",&(pos[0]),&(pos[1]));
            dataValue.steps.push_back(pos);
        }
        //Read winding number and data completeness
        double* w = new double[3];
        int isICValid = 1; //True by default
        fscanf(f,"%lf %lf %lf %d",w,w+1,w+2,&isICValid);
        for(int j=0; j<3; j++) {
            dataValue.wn.push_back( w[j] );
        }
        dataValue.isICValid = (isICValid==1)? true : false;
        
        //Add to the map container
        addData(id, dataValue);
    }
    
    mapParams.max_depth = maxDepth;
    mapParams.nb_iterations = numIters;
    
    //Close the file
    fclose(f);
}


#endif
