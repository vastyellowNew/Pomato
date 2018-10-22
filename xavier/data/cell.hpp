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


/* Cell class for a 2D cell in topology analysis.
 * Author: Wayne Schlei
 * Date: 7/18/2013
 * Purpose:  For building a list of edges, storing periods, and computing the Poincare Index
 */

#ifndef __CELL_HPP__
#define __CELL_HPP__

#include <iostream>
#include <vector>
#include <list>
#include <boost/static_assert.hpp>
#include <stdexcept>
//Chris
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
//Map Topology
#include <data/grid.hpp>
#include <data/adaptive_grid.hpp>
#include <data/edge.hpp>
#include <maps/metric.hpp>
#include <maps/index.hpp>
#include <maps/adaptiveGridMapInterface.hpp>
#include <maps/map_analysis.hpp>
//OpenMP
#if _OPENMP
#include <omp.h>
#endif
//Rounding access
#ifndef _WIN32
#include <fenv.h>
#pragma STDC FENV_ACCESS ON
#endif

namespace xavier {

/** A Cell class for 2D cells with N edges
      (usually N=4, but there could be special cases for future work)

    Poincare Index Call Sequence:
    1) create cells
    2) setCellPeriods()
    3) buildEdges()
    4) Call xavier::Edge_Rotation() outside this
    5) poincareIndex()
    6) detection(p) -> tells you if there is a fp in this cell for p
    7) getFixedPointGuess() -> will only find a guess if Poincare Index exists and is non-zero

    In the AdaptiveGrid case, I want DATASET to be the AdaptiveGridNodeData....
 */
template<class DATASET, typename IVEC = nvis::ivec3>
class Cell {
public :
    typedef Cell<DATASET, IVEC>                    self_type;
    //Vertex positions - on section 2D only
    typedef nvis::vec2                             vec_type;
    //Cell index (ivec2 for regular grid or ivec3 for adaptive_grid)
    typedef IVEC                                   ivec_type;
    //Edges
    typedef xavier::TerminalEdge<ivec_type>        edge_type;
    typedef xavier::CompositeEdge<ivec_type>       comp_edge_type;
    //DataSet access - AdaptiveGridNodeData<>
    typedef typename DATASET::DataMapConstIterator DataMapConstIterator;
    typedef typename DATASET::data_type            data_type;
    /// Edge types
    enum EdgeKind {
        BASE_EDGE = 0,
        COMPOSITE_EDGE
    };
    //Constructors:
    Cell() : dataField(0), numEdges(4), id(ivec_type(0,0))
    {
        //Default - just create a unit square
        nodeIDs.push_back(id); //node 0
        nodeIDs.push_back(id+ivec_type(1,0)); //node 1 - Need to be able to add IVEC for corners
        nodeIDs.push_back(id+ivec_type(1,1)); //node 2
        nodeIDs.push_back(id+ivec_type(0,1)); //node 3
        //Edge kinds
        std::vector<EdgeKind> ek(4,self_type::BASE_EDGE);
        edgeKinds = ek;
        verbose = false;
    }
    Cell(const int& n) : dataField(0), numEdges(n), id(ivec_type(0,0)), verbose(false) {}
    
    ///Construct from a DATASET and cell_id
    Cell(const DATASET& dataGrid, const ivec_type& cell_id) :
        dataField(&dataGrid),
        numEdges(4),
        id(cell_id)
    {
        //Fill nodes assuming a 4-node cell
        nodeIDs.push_back(cell_id); //node 0
        nodeIDs.push_back(cell_id+ivec_type(1,0)); //node 1
        nodeIDs.push_back(cell_id+ivec_type(1,1)); //node 2
        nodeIDs.push_back(cell_id+ivec_type(0,1)); //node 3
        //Edge kinds
        std::vector<EdgeKind> ek(4,self_type::BASE_EDGE);
        edgeKinds = ek;
        //... Future Work ...
        // -> Build a 3-node (triangle) cell
        // -> Build a 5-node (pentagon) cell
        // -> Build a n-node cell (for approximating curved edges like on ZVC, AHHH NO TIME!!!!)
        // -> Work in higher N-dimensions
        verbose = false;
    }
    
    
    //Member variables:
    /// DataField Pointer
    const DATASET* dataField;
    /// Sub-sampling grid data (for storing initial guess sub-grid)
    std::vector< std::vector<vec_type> > subGridDeltaSamples, subGridBetaSamples;
    /// Number of edges
    int numEdges;
    /// Cell ID
    ivec_type id;
    
    /// Cell Edge type
    std::vector<EdgeKind> edgeKinds;
    /// Cell node IDs (should be numEdges nodeIDs)
    std::vector<ivec_type> nodeIDs;
    /// Cell periods
    std::set<int> periods;
    /// Cell Poincare Index per period
    std::map<int,long int> pIndex;
    /// Map indicating if the section is transverse throughout cell at given period
    std::map<int,bool> isSectionTransverse;
    /// Map for minimum delta locations per period (and edge)
    std::map<int,std::vector<vec_type> > minDeltaPtsMap;
    /// Map for the delta values at minimum delta locations per period (and edge)
    std::map<int,std::vector<vec_type> > minDeltaValueMap;
    /// Turn on the output
    bool verbose;
    
    //Member functions:
    ///Subdivide cell based on fraction (u,v) and return actual point (parameters u and v are from 0 to 1)
    vec_type getSubPoint(const vec_type& uvVec)
    {
        vec_type spacing = dataField->getSpacing();
        return vec_type(uvVec[0]*spacing[0],uvVec[1]*spacing[1]) + dataField->getVertex(nodeIDs[0]);
    }
    ///Set the list of possible periods
    void setCellPeriods(std::vector<int>& periodVector)
    {
        std::vector<int>::iterator it;
        for (it=periodVector.begin(); it!=periodVector.end(); ++it) {
            //Sort to set
            periods.insert( *it );
        }
    }
    
    /// Set the set of periods from a set
    inline void setCellPeriods(const std::set<int>& periodSet)
    {
        periods = periodSet;
    }
    
    
    /// Insert an edge_type edge into the baseEdgeSet
    void insertEdgeToBaseEdgeSet(edge_type& newEdge, std::set<edge_type>& baseEdgeSet)
    {
        //See if edge already exists in list (all members of a set are const)
        typename std::set<edge_type>::iterator edgeIt;
        edgeIt = baseEdgeSet.find( newEdge );
        //If edge already exists, just insert this cell's periods
        if (edgeIt != baseEdgeSet.end()) {
            //Make a copy of the edge
            edge_type edgeCopy(*edgeIt);
            //Remove the edge from the set
            baseEdgeSet.erase(edgeIt);
            //Insert the additional periods to the edge
            std::set<int>::iterator pit;
            for (pit = periods.begin(); pit!=periods.end(); ++pit) {
                edgeCopy.periods.insert(*pit); //Checks if it is already there
            }
            //Insert the edge back into the master list
            baseEdgeSet.insert(edgeCopy);
        } else {
            //Edge wasn't found so we create it and add it to the list
            //Add periods
            std::set<int>::iterator pit;
            for (pit = periods.begin(); pit!=periods.end(); ++pit) {
                newEdge.periods.insert(*pit);
            }
            //Add to list
            baseEdgeSet.insert(newEdge);
        }
    }
    
    ///Build a set of edges after this cell is constructed
    void buildEdges(std::set<edge_type>& baseEdgeSet, std::set<comp_edge_type>& compositeEdgeSet)
    {
        //Don't add anything if cell is invalid
        if (!isCellValid()) {
            return;
        }
        //Loop through edge CCW
        for (int i=0; i<4; i++) {
            //Edge is connecting node i->node j
            int j= (i+1>3)? 0 : i+1;
            edge_type newEdge(nodeIDs[i],nodeIDs[j]);
            //Get the midpoint index, look up!
            ivec_type midID = newEdge.midpointID();
            //Is this a base-level edge or a CompositeEdge
            DataMapConstIterator mit = dataField->nodeDataMap.find(midID);
            if (mit == dataField->nodeDataMap.end() ) {
                // TerminalEdge - Base-level, no midpoint is found in nodeDataMap
                //  - Business as usual, add edge (or periods) to baseEdgeSet
                insertEdgeToBaseEdgeSet(newEdge, baseEdgeSet);
            } else {
                edgeKinds[i] = COMPOSITE_EDGE; //Flag for index
                // Composite Edge - made of subdivision edges
                comp_edge_type newCompEdge(newEdge);
                // Find the terminal segments and add to set of subnodes
                int numMidpoints = 1;
                //Recursive until found all sub edges (no midpoints left)
                while (numMidpoints > 0) { //Note: a parent->child setup might work better here
                    std::vector<ivec_type> nodesToAdd;
                    typename std::vector<ivec_type>::iterator nodeIT;
                    typename comp_edge_type::NodeSetType::iterator nodeSetIT, nodeSetIT2;
                    nodeSetIT = newCompEdge.nodeSet.begin();
                    nodeSetIT2 = nodeSetIT;
                    ++nodeSetIT2; //next node set member
                    for(; nodeSetIT2 != newCompEdge.nodeSet.end(); ++nodeSetIT,++nodeSetIT2 ) {
                        edge_type tempEdge((*nodeSetIT), (*nodeSetIT2));
                        midID = tempEdge.midpointID();
                        mit = dataField->nodeDataMap.find(midID);
                        //Insert as a subnode if found in data
                        if (mit != dataField->nodeDataMap.end()) {
                            nodesToAdd.push_back(midID);
                        }
                    }
                    numMidpoints = (int) nodesToAdd.size();
                    //Insert the new nodes
                    for(nodeIT=nodesToAdd.begin(); nodeIT!=nodesToAdd.end(); ++nodeIT) {
                        newCompEdge.nodeSet.insert( *nodeIT );
                    }
                }
                // Assign periods to composite edge
                std::set<int>::iterator pit;
                for (pit = periods.begin(); pit!=periods.end(); ++pit) {
                    newCompEdge.periods.insert(*pit);
                }
                // Look to see if edge already exists - I don't think they do at this point
                typename std::set<comp_edge_type>::iterator cEdgeIT;
                cEdgeIT = compositeEdgeSet.find(newCompEdge);
                if (cEdgeIT != compositeEdgeSet.end() ) {
                    //Composite Edges should be unique (not shared by cells)
                    compositeEdgeSet.erase(cEdgeIT);
                }
                // Insert the sub edges (TerminalEdges) or new data to the baseEdgeSet
                typename comp_edge_type::NodeSetType::iterator nodeSetIT, nodeSetIT2;
                nodeSetIT = newCompEdge.nodeSet.begin();
                nodeSetIT2 = nodeSetIT;
                ++nodeSetIT2; //next node in set
                for(; nodeSetIT2 != newCompEdge.nodeSet.end(); ++nodeSetIT,++nodeSetIT2) {
                    edge_type subEdge((*nodeSetIT),(*nodeSetIT2));
                    //Insert subEdge into base set (or just corresponding periods)
                    insertEdgeToBaseEdgeSet(subEdge, baseEdgeSet);
                }
                // Add completed CompositeEdge to structure
                compositeEdgeSet.insert(newCompEdge);
            }
        }//End for each edge loop
    }
    
#ifdef _WIN32
    ///Custom rint() function if using Visual Studio
    long int custom_rint(double x)
    {
        //middle value point test
        if (ceil(x+0.5) == floor(x+0.5)) {
            long int a = (long int)ceil(x);
            if (a%2 == 0) {
                return a;
            } else {
                return (long int) floor(x);
            }
        }
        
        else {
            return (long int) floor(x+0.5);
        }
    }
#endif
    
    /// Compute index for each period from edges (note: call topology/Edge_Rotation() before this)
    void poincareIndex(std::set<edge_type>& baseEdgeSet,std::set<comp_edge_type>& compositeEdgeSet)
    {
        //clear the current value
        pIndex.clear();
        //Don't
        if (!isCellValid()) {
            if(verbose) {
                std::cout << "This cell is INVALID! Exiting Poincare Index Computation!\n";
            }
            return;
        }
        
        std::set<int>::iterator periodIt;
        typename std::set<edge_type>::iterator eit;
        typename std::set<comp_edge_type>::iterator cEdgeIT;
        //For each period
        for(periodIt = periods.begin(); periodIt!=periods.end(); ++periodIt) {
            bool allEdgesOK = true;
            int p = *periodIt;
            double theta = 0.0;
            bool sectionTransversality = true;
            if(verbose) {
                std::cout << "P.I. for period = " << p << ":\n";
            }
            //Vector indicating minimum delta positions
            std::vector<vec_type> minDeltaPts, minDeltaValues;
            //Loop through edges CCW
            for (int i=0; i<4; i++) {
                //Find the ith edge
                int j= (i+1>3) ? 0 : i+1;
                edge_type thisEdge(nodeIDs[i],nodeIDs[j]);
                if (edgeKinds[i] == BASE_EDGE) { //Standard, terminal edge
                    eit = baseEdgeSet.find( thisEdge );
                    //Double-check that edge is in the list (it should be!!)
                    if (eit == baseEdgeSet.end()) {
                        theta = 0.0;
                        allEdgesOK = false;
                        if (verbose) {
                            std::cout << "   ERROR:Edge (" << i << "," << j << ") is not in the master edge list!!!\n";
                        }
                    }
                    //Compute total theta if no errors in edge
                    bool errorInMapRun = eit->mapError.find(p)->second;
                    if (errorInMapRun) {
                        //An Edge failed in computing rotation for this period
                        theta = 0.0;
                        allEdgesOK = false;
                        if (verbose) std::cout << "   ERROR:Edge (" << i << "," << j
                                                   << ") encountered a mapping error!  Setting theta = 0.0\n";
                    } else {
                        //Rotation is counted CCW, but edges are always stored
                        // left->right and bottom->top
                        std::map<int,double>::const_iterator mit = eit->rotationAngleMap.find(p);
                        if (mit == eit->rotationAngleMap.end()) {
                            theta = 0.0;
                            allEdgesOK = false;
                            if (verbose) std::cout << "    ERROR:Edge (" << i << "," << j
                                                       << ") does not have a have a rotation angle for period = " << p << "\n";
                        } else {
                            theta += ((i>1)?-1:1)*(mit->second);
                            if (verbose) {
                                std::cout << "   Edge (" << i << "," << j << ") computes rotation component as:\n";
                                std::cout << "       dTheta = " <<  ((i>1)?-1:1)*(eit->rotationAngleMap.find(p)->second)
                                          << " and now, theta = " << theta << "\n";
                            }
                            //Store the min delta points
                            minDeltaPts.push_back( eit->xMinDeltaMap.find(p)->second );
                            minDeltaValues.push_back( eit->minDeltaMap.find(p)->second );
                        }
                    }
                    //Lookup the transversality condition
                    if (!(eit->transverseSection.find(p)->second)) {
                        sectionTransversality = false;
                    }
                } else {
                    //This is a composite edge
                    cEdgeIT = compositeEdgeSet.find( comp_edge_type(thisEdge) );
                    //Double-check that this is an edge (it should be)
                    if (cEdgeIT == compositeEdgeSet.end() ) {
                        theta = 0.0;
                        allEdgesOK = false;
                        if (verbose) std::cout << "   ERROR: CompositeEdge (" << i << "," << j
                                                   << ") is not in the CompositeEdgeSet!!!\n";
                    }
                    //Compute total theta if no errors in edge
                    bool errorInMapRun = cEdgeIT->mapError.find(p)->second;
                    if (errorInMapRun) {
                        //An Edge failed in computing rotation for this period
                        theta = 0.0;
                        allEdgesOK = false;
                        if (verbose) std::cout << "   ERROR: CompositeEdge (" << i << "," << j
                                                   << ") has a failed subEdge due to a mapping error!  Setting theta = 0.0\n";
                    } else {
                        //Rotation is counted CCW, but edges are always stored
                        // left->right and bottom->top
                        std::map<int,double>::const_iterator mit = cEdgeIT->rotationAngleMap.find(p);
                        if (mit == cEdgeIT->rotationAngleMap.end()) {
                            theta = 0.0;
                            allEdgesOK = false;
                            if (verbose) std::cout << "    ERROR: CompositeEdge (" << i << "," << j
                                                       << ") does not have a have a rotation angle for period = " << p << "\n";
                        } else {
                            theta += ((i>1)?-1:1)*(mit->second);
                            if (verbose) {
                                std::cout << "   CompositeEdge (" << i << "," << j << ") computes rotation component as:\n";
                                std::cout << "       NumSeg = " << cEdgeIT->nodeSet.size() << "\n";
                                std::cout << "       dTheta = " <<  ((i>1)?-1:1)*(mit->second)
                                          << " and now, theta = " << theta << "\n";
                            }
                            
                            //Store the min delta points
                            minDeltaPts.push_back( cEdgeIT->xMinDeltaMap.find(p)->second );
                            minDeltaValues.push_back( cEdgeIT->minDeltaMap.find(p)->second );
                        }
                    }
                    //Lookup the transversality condition
                    if (!(cEdgeIT->transverseSection.find(p)->second)) {
                        sectionTransversality = false;
                    }
                }//End edge kind check
            }//End loop through edges CCW
            
            //Add the minimum delta points to storage
            minDeltaPtsMap.insert( std::pair<int,std::vector<vec_type> >(p,minDeltaPts) );
            minDeltaValueMap.insert( std::pair<int,std::vector<vec_type> >(p,minDeltaValues) );
            
            //Compute the index
            long int idx = 0;
#ifndef _WIN32
            fesetround(FE_TONEAREST);
            idx = lrint(0.5*theta/M_PI);
#else
            idx = custom_rint( 0.5*theta / M_PI );
#endif
            if (allEdgesOK) {
                if (verbose) {
                    std::cout << "Inserting period-index pair of <" << p << "," << idx << ">\n";
                }
                pIndex.insert( std::pair<int,long int>(p, idx) );
            }
            //Add the Transversality condition for the section
            isSectionTransverse.insert( std::pair<int,bool>(p,sectionTransversality) );
            
        }//End period loop
    }
    
    /// Is this edge valid
    bool isEdgeValid(const int& node0, const int& node1)
    {
        //No data so invalid cell/edge
        if (!dataField) {
            return false;
        }
        return (dataField->valid(nodeIDs[node0]) &&
                dataField->valid(nodeIDs[node1]));
    }
    
    /// Is this cell valid
    bool isCellValid()
    {
        //No data so invalid cell/edge
        if (!dataField) {
            return false;
        }
        //Check dataset for validity
        for(int i=0; i<4; i++) {
            bool valid = dataField->valid(nodeIDs[i]);
            if (!valid) {
                return false;
            }
        }
        return true;
    }
    
    /// Is the Poincare index non-zero
    inline bool detection(const int& p)
    {
        int isThere = (int) pIndex.count(p);
        if (isThere<1) {
            if (verbose) {
                std::cout << " Period " << p << " is not found for this cell.  detection(p) = false\n";
            }
            return false;
        }
        
        std::map<int,long int>::iterator it;
        it = pIndex.find(p);
        if (it==pIndex.end()) { //It didn't find it
            if(verbose) {
                std::cout << "Cell (" << id << ") has NO INDEX for period = " << p << "\n";
            }
            return false;
        } else if (it->second != 0 ) {
            if (verbose) {
                std::cout << "Cell (" << id << ") has DETECTED a fixed point with P.I. = " << it->second << " for period = " << p << "\n";
            }
            return true;
        } else {
            if (verbose) {
                std::cout << "Cell (" << id << ") has a ZERO INDEX for period = " << p << "\n";
            }
            return false;
        }
    }
    
    
    /** Computing a guess with a brute force approach - should store iterates for reuse with smaller period.
    *   - This MUST BE employed in a REVERSE loop on period.  Start from the largest and go to smallest for a
    *     cell.Traverse two leaves of a quadtree to find the best guess for a fixed point at a given period
    *  Note: This is not really used... */
    /*template <class MAP, class PARAMS>
    bool getFixedPointGuess(const MAP& theMap, const PARAMS& mapParams, const int& p, vec_type& guess)
    {
      //Check if p has a fixed point here (i.e., index is nonzero)
      bool ok = detection(p);
      if (!ok) return ok;
    
      //Run the map at the midpoint
      double delta = 0.5;
      vec_type parametricPoint(delta,delta);
      vec_type xMidpoint = getSubPoint(parametricPoint);
      vec_type gMidpoint(0.0);
      try {
        gMidpoint = mapDisplacement(xMidpoint,theMap,p,mapParams);
      } catch(...) {}
    
      //Storage for checking mins
      vec_type xMin(xMidpoint), gMin(gMidpoint);
      double gMinNorm = nvis::norm(gMin);
      if (gMinNorm < 1.e-6) {
          //This should be good enough so store it
          guess = xMin;
          return ok;
      }
      //Quadtree nodes - Unscaled, so multiply by delta
      std::vector<vec_type> qtNodes(4,vec_type(0,0));
      qtNodes[0] = vec_type(-1,-1);
      qtNodes[1] = vec_type(1,-1);
      qtNodes[2] = vec_type(1,1);
      qtNodes[3] = vec_type(-1,1);
      //Quadtree first level - always run
      delta *= 0.5;
      int minLeaf = -1;
      for (int leaf=0;leaf<4; leaf++) {
          vec_type gLeaf = mapDisplacement( getSubPoint(parametricPoint + delta*qtNodes[leaf]),
                    theMap,p,mapParams);
          if (nvis::norm(gLeaf) < nvis::norm(gMin)) {
        gMin = gLeaf;
        minLeaf = leaf;
          }
      }
      //If midpoint is still lowest, export midpoint as guess
      if (minLeaf == -1) {
          guess = xMin;
          return ok;
      }
      //Otherwise we reset and try second level
      parametricPoint += delta*qtNodes[minLeaf];
      xMin = getSubPoint(parametricPoint);
      //Start second level
      delta *= 0.5; minLeaf = -1;//Reset to new "center"
      for (int leaf=0;leaf<4; leaf++) {
          vec_type x = getSubPoint(parametricPoint + delta*qtNodes[leaf]);
          vec_type gLeaf = mapDisplacement(x,theMap,p,mapParams);
          if (nvis::norm(gLeaf) < nvis::norm(gMin)) {
        gMin = gLeaf;
        xMin = x;
        minLeaf = leaf;
          }
      }
      //Output best of second level
      guess = xMin;
      return ok;
    }*/
    
    /** Computing a guess employing a Model-fitting approach from sample sites within the
      *   analysis cell.  The initial guess is a linear-form saddle, but use Lev-Mar function to compute
      *   the best fit model.
    
      NOTE:  This is definitely a spot where GPU application could help!!!
      FUTURE:  Employ an option to use the LinearFit vs QuadraticFit
    */
    template <class MAP, class PARAMS>
    bool getFixedPointGuess(const MAP& theMap, const PARAMS& mapParams, const int& p,
                            const nvis::ivec2& gridRes, vec_type& guess)
    {
        //Check if p has a fixed point here (i.e., index is nonzero)
        bool ok = detection(p);
        if (!ok) {
            return ok;
        }
        
        //Check if grid data exists
        bool subSampleAvailable = false;
        if ( ( (int) subGridDeltaSamples.size() > 0 ) &&
                ((int) subGridDeltaSamples[0].size() > p) ) {
            //Sub grid data exists and covers the requested period
            subSampleAvailable = true;
        } else {
            //start fresh
            subGridDeltaSamples.clear();
            subGridBetaSamples.clear();
        }
        
        // If non-transverse section, try a guess strategy that utilizes stable-manifold behavior
        if ( !isSectionTransverse[p] ) {
            //Call non-transverse section sampling function
            ok = ntGuessFinder(theMap,mapParams,p,gridRes,guess);
            return ok;
        }
        
        //Create a gridRes[0]xgridRes[1] sampling grid inside edges (4x4 means a 6x6 cell)
        int numPts = gridRes[0]*gridRes[1];
        vec_type x2 = dataField->getVertex(nodeIDs[2]);
        vec_type x0 = dataField->getVertex(nodeIDs[0]);
        vec_type spacing( x2 - x0 );
        spacing[0] = spacing[0]/((double) gridRes[0]+1.);
        spacing[1] = spacing[1]/((double) gridRes[1]+1.);
        vec_type xMin(x0);
        vec_type dispMin(1000.); //Forward displacement minimum
        std::vector<vec_type> xVals, deltaVals, betaVals;
        bool integOK = true;
        double minNorm = 1000.;
        //Propagate through grid and collect best guess (note runs on single thread, but cells are on parallel threads)
        for (int n=0; n<numPts; n++) {
            int i = n % gridRes[0];
            int j = n / gridRes[0];
            vec_type delta(spacing[0]*((double)i+1),spacing[1]*((double)j+1));//Start inside edge
            vec_type x = x0 + delta;
            vec_type deltaVec(0.0), betaVec(0.0);
            
            // If we have subsample data, use it!  Otherwise, compute and store
            if (subSampleAvailable) {
                deltaVec = mapParams.the_metric.displacement(x,subGridDeltaSamples[n][p-1]); // P^p(x) - x
                betaVec = mapParams.the_metric.displacement(x,subGridBetaSamples[n][p-1]); // P^-p(x) - x
            } else {
                //Compute the map - Note this is on a single thread
                MAP* amap = theMap.clone();
                std::vector<nvis::vec2> tmp, btmp;
                try {
                    amap->map(x, tmp, p);
                    deltaVec = mapParams.the_metric.displacement(x, tmp[p-1]);
                    amap->map(x, btmp, -p); //Backward map
                    betaVec = mapParams.the_metric.displacement(x, btmp[p-1]);
                } catch(...) {
                    std::cout << "   Unable to compute the subsampled point " << x << "\n";
                    tmp.clear();
                    btmp.clear();
                    for (int i=0; i<p; i++) {
                        tmp.push_back(nvis::vec2(1000,1000));
                        btmp.push_back(nvis::vec2(1000,1000));
                    }
                    deltaVec = betaVec = nvis::vec2(1000.,1000.);
                    integOK = false;
                }
                if (integOK) {
                    //Store info for model fitting
                    xVals.push_back(x);
                    deltaVals.push_back(deltaVec);
                    betaVals.push_back(betaVec);
                }
                //Store point
                subGridDeltaSamples.push_back( tmp ); //Store nth point's P^p(x)
                subGridBetaSamples.push_back( btmp );
            }
            
            //Check if minimum
            double mapDisp = nvis::norm(deltaVec);
            if (mapDisp < minNorm) {
                xMin = x;
                dispMin = deltaVec;
                minNorm = mapDisp;
            }
        }
        
        //Set the guess as the guess with the smallest error
        guess = xMin;
        
        //If cell edges indicate that section transversality assumption is violated, then we must stop here
        //if ( !isSectionTransverse[p] ) return ok; //Brute force guesses are poor!
        
        //Run the Saddle Model-fitting if the index indicates this is a saddle
        std::map<int,long int>::iterator it;
        it = pIndex.find(p);
        if (it->second < 0) {
            //Saddle guess generator using Levenberg-Marquardt
            vec_type center = (x2+x0)/2.0, deltaVec;
            VectorXd coeffs; //The coefficient values
            
            //LINEAR Fit first!
            vec_type mfGuess = saddleGuessFromModelFit<vec_type> (center,xVals,deltaVals,betaVals,coeffs);
            int nLinParams = (int) coeffs.size();
            //Map the result to be sure it's better than subgrid
            MAP* amap = theMap.clone();
            std::vector<vec_type> tmp;
            try {
                amap->map(mfGuess, tmp, p);
                deltaVec = mapParams.the_metric.displacement(mfGuess, tmp[p-1]);
                if (nvis::norm(deltaVec) < minNorm) {
                    guess = mfGuess;
                }
            } catch (...) {}
            
            //QUADRATIC Fit - Uses Linear Solution as guess
            VectorXd qdCoeffs;
            qdCoeffs.setZero( nLinParams + 2*2*2 ); //Add Tensor (2x2x2)
            qdCoeffs.head(nLinParams) = coeffs; //Use Linear result
            vec_type qmfGuess = saddleGuessFromModelFit<vec_type> (mfGuess,xVals,deltaVals,betaVals,qdCoeffs,true);
            //Map the result to see if it is superior to the subgrid points
            tmp.clear();
            try {
                amap->map(qmfGuess, tmp, p);
                deltaVec = mapParams.the_metric.displacement(qmfGuess, tmp[p-1]);
                if (nvis::norm(deltaVec) < minNorm) {
                    guess = qmfGuess;
                }
            } catch (...) {}
            
            //LINEAR FIT ON NORMALIZED MAP TANGENT FIELD
            coeffs.setZero( nLinParams );
            vec_type nmfGuess = saddleGuessFromMapTangentNorm<vec_type>(center,xVals,deltaVals,betaVals,coeffs);
            //Map the result to see if it is superior to the subgrid points
            tmp.clear();
            try {
                amap->map(nmfGuess, tmp, p);
                deltaVec = mapParams.the_metric.displacement(nmfGuess, tmp[p-1]);
                if (nvis::norm(deltaVec) < minNorm) {
                    guess = nmfGuess;
                }
            } catch (...) {}
        }
        
        return ok;
    }
    
    ///Non-transverse Section sampling to find a guess (utilizes stable manifold behavior to find guess)
    template <class MAP, class PARAMS>
    bool ntGuessFinder(const MAP& theMap, const PARAMS& mapParams, const int& p,
                       const nvis::ivec2& gridRes, vec_type& guess)
    {
        //Num sample points per edge
        int numSamples = gridRes[0];
        
        //Grab the 4 edge minimum delta points
        std::vector<vec_type> minDeltaPts = minDeltaPtsMap.find(p)->second;
        std::vector<vec_type> minDeltaValues = minDeltaValueMap.find(p)->second;
        vec_type xMin;
        double deltaMin = 1000.0;
        double dmu = 1.0/((double) numSamples + 1);
        //For each sampling line (0:1,2,3; 1:2,3; 2:3)
        // ->Recall, Cells are processed in parallel
        for(int i=0; i<4; i++) {
            //Check the min point values
            if (nvis::norm(minDeltaValues[i]) < deltaMin) {
                xMin = minDeltaPts[i];
                deltaMin = nvis::norm(minDeltaValues[i]);
            }
            for(int j=i+1; j<4; j++) {
                //Form the line
                //-> Assumes that one line that connects minimum delta points is the stable manifold
                vec_type& ptI = minDeltaPts[i];
                vec_type& ptJ = minDeltaPts[j];
                //Check the sampling points delta for best guess
                for (int k=1; k<numSamples+1; k++) {
                    double mu = dmu*k;
                    //Weighted average
                    vec_type ptK = ptI*(1.0-mu) + ptJ*mu;
                    //Run the map
                    vec_type delta;
                    try {
                        delta = mapDisplacement(ptK,theMap,std::abs(p),mapParams);
                        if (nvis::norm(delta) < deltaMin) {
                            xMin = ptK;
                            deltaMin = nvis::norm(delta);
                        }
                    } catch(...) {
                        //We just have to skip this point if there's a mapping error
                    }
                }
            }
        }
        
        //Set the guess as the guess with the smallest error
        guess = xMin;
        return true;
    }
};


} //end xavier


#endif //__CELL_HPP__
