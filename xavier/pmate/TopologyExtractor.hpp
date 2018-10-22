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


/**  TopologyExtractor
  *    Author : Wayne Schlei
  *    Purpose:  Class object that autonomously extract topology
  *    from a Poincare section
  *    Utilizes:
  *     AdaptiveSampler (if not input) : Runs the sampling
  *
  *    Inputs:
  *
  *      Note: This class only samples the intial grid; however, you
  *      can specify ahead of time to read AdaptiveGridNodeData by using:
  *         TopologyExtractor *topoExtractor = new TopologyExtractor(...);
  *         topoExtractor->theSampler.read(agndFile);
  *      before the 'compute()' function call.
  *    Outputs:
  */

#ifndef TOPO_ANALYSIS_HPP
#define TOPO_ANALYSIS_HPP

#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip> //C++11
#include <memory> //C++11
#include <ctime>
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>
//#include <boost/shared_ptr.hpp> //boost alternative???

//API - chris
#include <math/rational.hpp>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>
//API - xavier
#include <math/rational.hpp>
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/cell.hpp>
#include <data/raster_data.hpp>
#include <data/adaptive_grid.hpp>
#include <maps/definitions.hpp>
#include <maps/DP45wrapper.hpp>
#include <maps/adaptiveGridMapInterface.hpp>
#include <maps/map_analysis.hpp>
#include <maps/index.hpp>
#include <maps/poincare_map.hpp>
//API - orbital
#include <orbital/corrections.hpp>
#include <orbital/monodromy.hpp>
//API - topology
#include <topology/ManifoldDataStorage.hpp>
#include <topology/Edge_Rotation.hpp>
#include <topology/EdgeRotationFailure.hpp>
#include <topology/invariant_manifold.hpp>
//API - pmate
#include <pmate/AdaptiveSampler.hpp>
#include <pmate/FixedPointData.hpp>
#include <pmate/FixedPointDataFilters.hpp>

//Point cloud file output (for Avizo)
#include <cr3bp/psiWrite.hpp>


#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;

namespace pmate {

/// Class for adaptively sampling a Poincare map with respect to smoothness in winding number
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
class TopologyExtractor {
    typedef AdaptiveSampler<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>  AdaptiveSamplerType;
    typedef orbit_data                               OrbitData;
    typedef std::vector<xavier::fixpoint>            FixedPointChain;
    typedef FixedPointData<xavier::fixpoint>         FPData;
    // Note: 2D Map analysis objects
#ifdef _WIN32
    typedef AdaptiveSamplerType::AdaptiveGridType    AdaptiveGridType;
    typedef AdaptiveSamplerType::BoundsType          BoundsType;
    typedef AdaptiveSamplerType::DimType             ResType;
    typedef AdaptiveSamplerType::AGNodeData          DataSet;
    typedef AdaptiveGridType::id_type                IdType;
    typedef DataSet::DataMapIterator                 DataMapIterator;
    typedef PMAP::rhs_type                           rhs_type;
    typedef PMAP::section_type                       section_type;
    typedef PMAP::lvec_type                          lvec_type;
    typedef xavier::EdgeRotationFailure<lvec_type>   EdgeRotFailure;
    typedef xavier::TerminalEdge<IdType>             EdgeType;
    typedef xavier::CompositeEdge<IdType>            CompEdgeType;
    typedef std::set<EdgeType>                       EdgeSet;
    typedef std::set<CompEdgeType>                   CompositeEdgeSet;
    typedef std::set<EdgeType>::const_iterator       EdgeSetIter;
    typedef std::set<CompEdgeType>::const_iterator   CompositeEdgeSetIter;
    typedef AdaptiveGridType::data_type              DataPairType;
    typedef xavier::Cell<DataSet,IdType>             MapCell;
#else
    typedef typename AdaptiveSamplerType::AdaptiveGridType  AdaptiveGridType;
    typedef typename AdaptiveSamplerType::BoundsType        BoundsType;
    typedef typename AdaptiveSamplerType::DimType           ResType;
    typedef typename AdaptiveSamplerType::AGNodeData        DataSet;
    typedef typename AdaptiveGridType::id_type              IdType;
    typedef typename DataSet::DataMapIterator               DataMapIterator;
    typedef typename PMAP::rhs_type                         rhs_type;
    typedef typename PMAP::section_type                     section_type;
    typedef typename PMAP::lvec_type                        lvec_type;
    typedef typename xavier::EdgeRotationFailure<lvec_type> EdgeRotFailure;
    typedef typename xavier::TerminalEdge<IdType>           EdgeType;
    typedef typename xavier::CompositeEdge<IdType>          CompEdgeType;
    typedef typename std::set<EdgeType>                     EdgeSet;
    typedef typename std::set<CompEdgeType>                 CompositeEdgeSet;
    typedef typename std::set<EdgeType>::const_iterator     EdgeSetIter;
    typedef typename std::set<CompEdgeType>::const_iterator CompositeEdgeSetIter;
    typedef typename AdaptiveGridType::data_type            data_pair_type;
    typedef typename xavier::Cell<DataSet,IdType>           MapCell;
#endif
    
    
public :
    ///Constructor
    TopologyExtractor(PMAP& mapIn, PARAM& params) :
        theMap(&mapIn),
        theMapParams(&params),
        theSampler(params.bounds, params.resolution, mapIn, params),
        enableWatchCell(false),
        theWatchCell(0),
        runName("none"),
        agndFileStr("none"),
        mpFileStr("none"),
        //erotFileStr("none"),
        mapDiscontFileStr("none"),
        fpGuessFileStr("none"),
        fpRefineFileStr("none"),
        fpFileStr("none"),
        nbthreads(1)
    {}
    
    
    /// Execute the autonomous Poincare map topology extraction process
    void compute()
    {
        //Phase A: Domain, Initialize based on settings
        //Note: Mostly completed in constructor.
        
        //Phase B: Sampling & Characterization based on winding number
        if(theSampler.isUpdateRequired()) {
            theSampler.sample();
        }
        // Clear out all invalid cells from the adaptive grid
        theSampler.removeInvalidLeaves();
        //Create cells and edges and assign information
        cellPeriods();
        generateEdges();
        
        //Phase C: Detection via Poincare Index
        edgeRotation();
        cellPoincareIndex();
        //Phase D: Guess generation for fixed points (in cellPoincareIndex())
        writeFixedPointGuesses();
        
        //Phase E: Fixed Point Refinement
        fixedPointRefinement();
        if ((int) fp_chains.size() > 0) {
            FPData* fps = getFixedPointData();
            fixedPointDataFilter((*fps)); //Check for doubles and mirror theorem
            fps->write(fpFileStr.c_str()); //Write results to file
            writeFailedFixedPoints();
        } else {
            std::cerr << "PMATE: No Fixed Points detected! Exiting...\n";
        }
        //Phase F:  Manifold Extraction (separate class)
        // -> Uses FixedPointData, theMap, theMapParams as input.
    }
    
    
    // Functions for individual phases of extraction process: -------------
    /// Compute per-cell Period information
    void cellPeriods();
    /// Generate Edge list from cells
    void generateEdges();
    /// Compute the per-edge map displacement vector rotation
    void edgeRotation();
    /// Compute the Poincare Index for cell from edges and formulate fixed point guesses
    void cellPoincareIndex();
    /// Fixed point refinement from guesses
    void fixedPointRefinement();
    /// Filtering the fixed points results (step 2)
    void fixedPointDataFilter(FPData& fpd);
    //---------------------------------------------------------------------
    
    /// The Adaptive Sampling Object (stores thePlanarGrid,theAdaptiveGrid,theAdaptiveGridData)
    AdaptiveSamplerType theSampler;
    
    
    /// Reset the grid information
    void resetGrid()
    {
        theSampler.reset();
    }
    
    
    //Functions to set Computation parameters from outside
    /// Set the grid bounding box
    void setGridBounds(const BoundsType& bbox)
    {
        theSampler.setGridBounds(bbox);
    }
    /// Set the grid initial resolution
    void setGridRes(const ResType& res)
    {
        theSampler.setGridRes(res);
    }
    /// Set the maximum depth level (Default is 3)
    void setMaxDepth(const int d)
    {
        theSampler.setMaxDepth(d);
    }
    
    /// Set the winding number tolerances (for cell convexity checks)
    void setWindingTols(const WVEC& tols)
    {
        for(int i=0; i<(int)tols.size(); i++) {
            theMapParams->winding_convexity_tols[i] = tols[i];
        }
        resetGrid();
    }
    
    /// Set the winding number maximum distances (for cell convexity checks)
    void setWindingMaxDist(const WVEC& mDist)
    {
        for(int i=0; i<(int)mDist.size(); i++) {
            theMapParams->winding_cell_maxDist[i] = mDist[i];
        }
        resetGrid();
    }
    
    /// Set the map parameters object (with grid-based info)
    void setMapParams(PARAM& params)
    {
        theSampler.setMapParams(params);
    }
    
    /// Enabling a specific cell monitor (useful for debugging)
    bool enableWatchCell;
    /// Watch Cell:  Outputs information for a provided watch cell if enabled
    IdType theWatchCell;
    
    
    // Writing members and functions -----------------------------------------------------------------
    /// String indicating the name of the run (If not "none", this will be included in every filename)
    std::string runName;
    /// String indicating Adaptive Grid Node Data file
    std::string agndFileStr; //.agnd
    /// String indicating map parameters file name
    std::string mpFileStr; //.param
    /// Name of Edge Rotation results file
    //std::string erotFileStr; //.erot
    /// Name of Map Discontinuities file
    std::string mapDiscontFileStr; //.psi (.discont?)
    /// Name of Fixed Point Guesses file
    std::string fpGuessFileStr; //.fpguess
    /// Name of Fixed Point Faiures file
    std::string fpFailFileStr; //.fpfail
    /// Name of Fixed Point Refinement results file
    std::string fpRefineFileStr; //.fprefine
    /// Name of Fixed Point Data file
    std::string fpFileStr; //.fpx
    /// Compose all the names off of runName member
    void composeFileNames();
    /// Compose all names off of input
    void composeFileNames(const char* baseName)
    {
        runName = baseName;
        composeFileNames();
    }
    
    /// Write adaptive grid node data to file
    void writeSamplingData();
    /// Write Map Analysis Parameters to file
    void writeMapParams();
    /// Write Edge Rotation Results to file
    //void writeEdgeRotation() const;
    /// Write Map Discontinuities to file
    void writeMapDiscont();
    /// Write Fixed Point Guesses to file
    void writeFixedPointGuesses() const;
    /// Write Fixed Point Refinement Info to file
    //void writeFixedPointRefine() const; //Handled within output strings
    /// Write the Failed Fixed Point Guesses
    void writeFailedFixedPoints() const;
    /// Get Fixed Point Data object of results
    FPData* getFixedPointData() const;
    /// Write Fixed Point Data object of results (assuming they are computed)
    void writeFixedPointResults(const char* filename) const;
    //--------------------------------------------------------------------------------------------------
    /// Number of computational threads used in parallel (openMP)
    int nbthreads;
    
    
private :
    /// Pointer to the map engine
    PMAP* theMap;
    /// Pointer to the map parameters object
    PARAM* theMapParams;
    //Data objects indicating elements of topology analysis:
    //   Analysis
    EdgeSet                             edges;
    CompositeEdgeSet                    compositeEdges;
    std::vector<MapCell>                cellsVector;
    //   Fixed Points
    std::vector<FixedPointChain>       fp_chains;
    FixedPointChain                    fixed_points;
    FixedPointChain                    fpGuesses, failed_guesses;
    
};


/// Compute per-cell Period information
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
cellPeriods()
{
    std::cout << "PMATE:  Processing Cells and Periods ...\n";
    std::vector<IdType> leaves;
    theSampler.theData.getLeaves(leaves);
    int ncells = (int) leaves.size();
    int pmax = theMapParams->max_period;
    //Per-thread cell cache
    std::vector<MapCell>* cached_cells = new std::vector<MapCell>[nbthreads];
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<ncells ; ++n) {
            //int i = n % cellres[0];
            //int j = n / cellres[0];
            IdType cell_id = leaves[n];
            int i = cell_id[0];
            int j = cell_id[1];
            int d = cell_id[2];
            
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Check if cell is being watched...
            bool cellWatcher = false;
            if ( enableWatchCell && nvis::all(cell_id == theWatchCell) ) {
                cellWatcher = true;
            }
            
            //Cell Processing: ------------------------------------------------------------
            //Compute the set of viable periods for a cell
            std::vector<int> periods;
            std::set<int> periodSet;
            //Call an algorithm using the best rational approximation based on pmax
            xavier::period_range<DataSet,IdType>(periods, theSampler.theData, cell_id, pmax, cellWatcher);
            //Convert periods to a set
            for (int k=0 ; k<(int)periods.size() ; ++k) {
                if (periods[k]<1 || periods[k]>pmax) {
                    std::cerr << "PMATE Warning: Wrong period (" << periods[k] << ") in cell (" << cell_id << ")\n";
                    //theMsg->printf("Wrong Period ( %d ) in cell (%d,%d,%d)",periods[k], i,j,d);
                }
                //Create a set for the periods (unique) in this cell
                periodSet.insert( periods[k] );
            }
            //Also run best period approximation at each vertex
            int bestPApprox = 0;
            std::vector<int> bestPvec;
            IdType cornerID[4] = {cell_id, cell_id+IdType(1,0,0), cell_id+IdType(1,1,0), cell_id+IdType(0,1,0)};
            for (int ii=0; ii<4; ii++) {
                //Call best_period() as long as data is available
                int numSteps = (int) theSampler.theData(cornerID[ii]).steps.size();
                int testStepsNum = std::min(pmax,numSteps-1);
                //Run only if valid
                if (testStepsNum > 0)
                    bestPApprox = xavier::best_period(theSampler.theData(cornerID[ii]).steps,
                                                      testStepsNum, theMapParams->the_metric);
                                                      
                if(cellWatcher) {
                    std::cout << "PMATE:  Best Period = " << bestPApprox << " for watch cell (" << cell_id << "), corner " << ii << "\n";
                    //theMsg->printf(" Best period = %d in cell (%d,%d,%d), corner %d ", bestPApprox,i,j,d,ii);
                }
                if (bestPApprox >= 1 && bestPApprox <= pmax) { //in valid range
                    //Add to set (if not already there)
                    periodSet.insert( bestPApprox );
                    periods.push_back( bestPApprox ); //Add to cell input
                }
                bestPvec.push_back(bestPApprox);
            }
            
            //Create Cell and Add to cache
            MapCell thisCell(theSampler.theData,cell_id);
            if (cellWatcher) {
                thisCell.verbose = true;
            }
            thisCell.setCellPeriods( periods );
            cached_cells[thread_id].push_back(thisCell);
            //Display ----------------------------------------------------------------------------------------------
            //Watching a cell - Avizo Output of periods
            /*if (cellWatcher) {
                QString cellPeriodsString("  Cell (");
                QString iText, jText;
                iText.setNum(i); jText.setNum(j);
                cellPeriodsString += iText + "," + jText + "): IsCellValid() = ";
                cellPeriodsString += ( (thisCell.isCellValid()) ? "true ->" : "false ->");
                cellPeriodsString += " Possible Periods [";
                QString periodString;
                std::set<int>::iterator it;
                std::set<int>::iterator ePSIT = periodSet.end();
                --ePSIT;
                for ( it = periodSet.begin(); it!=ePSIT; ++it) {
                        periodString.setNum(*it);
                        cellPeriodsString += periodString + ",";
                }
                periodString.setNum(*ePSIT);
                cellPeriodsString += periodString + "] with BestP = (";
                for(int ii=0;ii<3;ii++) {
                  periodString.setNum(bestPvec[ii]);
                  cellPeriodsString += periodString + ",";
                }
                periodString.setNum(bestPvec[3]);
                cellPeriodsString += periodString + ")";
                theMsg->printf(cellPeriodsString);
            }*/
            
            std::cerr << "\rPMATE:  Computed period of " << n << " / " << ncells
                      << " (" << 100*n/ncells << "%)          \r" << std::flush;
                      
            //------------------------------------------------------------------------------------------------------
            
        } //End loop
    }//End parallel Statement
    std::cerr << "PMATE:  Computed cell periods for " << ncells << "\n";
    //std::cerr << "\nperiod range computation took " << _timer.elapsed() << '\n';
    
    for (int i=0 ; i<nbthreads ; ++i) {
        std::copy(cached_cells[i].begin(), cached_cells[i].end(), std::back_inserter(cellsVector));
    }
    for (int i=0 ; i<nbthreads ; ++i) {
        cached_cells[i].clear();
    }
    delete[] cached_cells;
    
}

/// Generate Edge list from cells
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
generateEdges()
{
    edges.clear();
    compositeEdges.clear();
    //For each cell (valid)
    for(int cid=0; cid<(int)cellsVector.size(); cid++) {
        //Create edges without duplicates
        cellsVector[cid].buildEdges(edges, compositeEdges);
    }
}

/// Compute the per-edge map displacement vector rotation
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
edgeRotation()
{
    //Update the map engine with lower (more precise) integration tolerance
    //theMap->setPrecision(theMapParams->refinementIntegTol);
    // NOT worth the added computation time!
    
    //Use Edge_Rotation.hpp code: per-edge computation that cuts down on pMap calls
    std::cerr << "PMATE:  Evaluating rotation of map displacement along "
              << (int)edges.size() << " base edges forming cells\n";
    topology::Edge_Rotation<DataSet,PMAP,PARAM,EdgeSet,xavier::fixpoint>
    (theSampler.theData,(*theMap),(*theMapParams),edges,fpGuesses);
    
    //Compute the composite Edges after the terminal edges are computed (No pMap calls)
    topology::compositeEdgeRotation<CompositeEdgeSet,EdgeSet>( compositeEdges, edges);
    
    //Printing Edge_Rotation results:
    //???
}

/// Compute the Poincare Index for cell from edges and formulate fixed point guesses
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
cellPoincareIndex()
{
    //Call the Poincare index computation
    std::cout << "PMATE: Computing the Poincare Index for " << cellsVector.size() << " cells\n";
    //(Note: single thread because I don't want to make copies, just data lookup)
    int ncells = (int) cellsVector.size();
    for(int cid=0; cid<ncells; cid++) {
        cellsVector[cid].poincareIndex( edges , compositeEdges );
    }
    
    // ---------------------------------------------------------------------------------------
    //Gather the fixpoint guesses for each cell in parallel
    std::cout << "PMATE: Gathering guesses for fixed points for cells with non-zero Poincare Index\n";
    int numERguess = (int) fpGuesses.size();
    std::cout << "PMATE: There are " << fpGuesses.size() << " fixed point guesses are from Edge_Rotation\n";
    std::cout << "       (i.e., found on cell edges during Poincare Index evaluation)\n";
    
    //Note:  We may want to up the precision here for a more-precise estimate
    //theMap->setPrecision( theMapParams->refinementIntegTol );
    
    //Per-thread fixpoint cache
    FixedPointChain* cached_fps = new FixedPointChain[nbthreads];
    //Run cell-internal sampling grid to find the guess
    // -> For use in Least displacement, LevMar model fitting, modified sampling in NT cells
    nvis::ivec2 bruteForceGrid = theMapParams->subcellResolution;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<ncells ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Loop through all periods
            std::set<int>::iterator pit;
            for( pit=cellsVector[n].periods.begin(); pit != cellsVector[n].periods.end(); ++pit ) {
                //Gather the fixed point guess using subsampling (and LevMar ModelFit)
                nvis::vec2 theGuess(0,0);
                bool found = cellsVector[n].getFixedPointGuess((*theMap),(*theMapParams),*pit,bruteForceGrid,theGuess);
                if (found) {
                    xavier::fixpoint aFP;
                    aFP.K = *pit;
                    aFP.pos = theGuess;
                    cached_fps[thread_id].push_back( aFP );
                }
            }
            if (thread_id == 0) {
                float prog = (float) n/(float)ncells;
                std::cerr << "\rPMATE: Computed Guesses for " << n << "/" << ncells
                          << " cells = " << prog*100 << "\% complete\r" << std::flush;
            }
            
            
        }//End parallel for
    }//End parallel statment
    
    //Add to guess vector (note, some guesses might be there already from Edge_Rotation()
    for(int i=0; i<nbthreads; i++) {
        for(int j=0; j<(int)cached_fps[i].size(); j++) {
            fpGuesses.push_back( cached_fps[i][j] );
        }
    }
    for(int i=0; i<nbthreads; i++) {
        cached_fps[i].clear();
    }
    delete[] cached_fps;
    
    //Prompt user
    std::cout << "PMATE: There are " << fpGuesses.size()-numERguess << " fixed-point seeds from Poincare Index evaulation.\n";
}

/// Store fixed point guesses
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeFixedPointGuesses() const
{
    if (fpGuessFileStr == "none") {
        std::cerr << "Bad filename for FixedPointGuessFile... Skipping guess storage.\n";
    } else {
        bool didItWork = writeFixpointVecBasic(fpGuessFileStr.c_str(),fpGuesses);
        if (!didItWork) {
            std::cerr << "FixedPointGuessFile write error... Skipping guess storage.\n";
        }
    }
    
    //We may also want to output these to a PSI file to load as an HxCluster for Avizo
}

/// Write the failed fixed point guesses to file
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeFailedFixedPoints() const
{
    if (fpFailFileStr == "none") {
        std::cerr << "Bad filename for FailedFixedPointFile... Skipping guess storage.\n";
    } else {
        bool didItWork = writeFixpointVecBasic(fpFailFileStr.c_str(),failed_guesses);
        if (!didItWork) {
            std::cerr << "FixedPointGuessFile write error... Skipping guess storage.\n";
        }
    }
    
    //We may also want to output these to a PSI file to load as an HxCluster for Avizo
}


/// Fixed point refinement from guesses
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
fixedPointRefinement()
{
    //Look through points and see if we can filter out obvious duplicates at higher periods
    int inputNumGuesses = (int) fpGuesses.size();
    //Convert to list<fixpoint> so we can do sorting
    std::list<xavier::fixpoint> guessList;
    std::list<xavier::fixpoint>::iterator gIT;
    for(int i=0; i<inputNumGuesses; i++) {
        guessList.push_back( fpGuesses[i] );
    }
    //Sort the list in order of higher periods to lower periods
    guessList.sort(xavier::fpPeriodCompare<xavier::fixpoint>);
    //Run "unique" command to remove duplicates at higher period (predicate in fixpoints.hpp)
    guessList.unique(xavier::is_highperiod_duplicate(1.e-3));
    int numGuesses = (int) guessList.size();
    int droppedNumGuesses = inputNumGuesses - numGuesses;
    if (droppedNumGuesses > 0) {
        std::cout << "PMATE: Pre-screening filtered out " << droppedNumGuesses << " guesses at higher periods.\n";
    }
    
    //Convert back to a vector for omp for
    FixedPointChain guessVec;
    for(gIT=guessList.begin(); gIT!=guessList.end(); ++gIT) {
        guessVec.push_back( *gIT );
    }
    
    //Update the map engine with the refinement integration tolerance
    theMap->setPrecision(theMapParams->refinementIntegTol);
    PMAP* amap = theMap->clone(); //Create new map copy in this scope
    
    //Call Multi-Tiered Targeting scheme within each thread to solve for a fixed point
    //int numGuesses = (int) guesses.size();
    std::vector<FixedPointChain>*  cached_fps = new std::vector<FixedPointChain>[nbthreads]; //Note: these are stored as chains of fixed points (vector<vector<fixpoint>>)
    FixedPointChain* cached_failed = new FixedPointChain[nbthreads];
    int counter = 0, foundCount=0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0; n<numGuesses; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            //Computation inputs
            nvis::vec2 x0 = guessVec[n].pos;
            int p = guessVec[n].K;
            xavier::fixpoint fp_result;
            FixedPointChain resultChain;
            //Parameters in solveFixedPoints():
            nvis::bbox2 bounds; //Stops points from going too far away from guess
            nvis::vec2 theGrid_spacing = theSampler.theGrid.spacing(theMapParams->max_depth);
            bounds.min() = x0 - theGrid_spacing;
            bounds.max() = x0 + theGrid_spacing; //Allow this to enter adjacent cells
            int maxIters = theMapParams->refinementMaxIters;
            int numPtsPerPeriod = theMapParams->refinementPointsPerPeriod;
            //Run multiple shooting with Periodicity, Jacobi constant, and Section constraints
            bool found = orbital::solveFixedPoints<PMAP>(
                             *amap, theMapParams->the_metric, bounds,
                             x0, amap->rhs().desired_hamiltonian(), 5, p,
                             fp_result, resultChain, theMapParams->refinementConvTol, theMapParams->verbose,
                             maxIters, numPtsPerPeriod, theMapParams->linearMonodromy);
                             
            //Store the result if converged or failed
            if (found) {
                cached_fps[thread_id].push_back( resultChain );
                std::cout << "PMATE:  Corrections found the " << (resultChain[0].saddle?"SADDLE ":"CENTER ") << resultChain[0].pos << " for p=" << p << "\n";
                #pragma omp atomic
                foundCount++;
            } else {
                cached_failed[thread_id].push_back( guessVec[n] );
            }
            
            #pragma omp atomic
            counter++;
            
            if (theMapParams->verbose) {
                std::ostringstream os;
                os << "\rCompleted Corrections process for " << counter << " / " << numGuesses << " ("
                   << 100.*(float)counter/(float)numGuesses << "%), "
                   << foundCount << " fixed points found.     \r" << std::flush;
                std::cout << os.str();
            }
            
        }//End parallel for
    }//End parallel statement
    
    //Reset Data
    fixed_points.clear();
    fp_chains.clear();
    failed_guesses.clear();
    //Data Process per thread
    for (int i=0; i<nbthreads; i++) {
        //Assemble the found fixed points
        for (int j=0; j<(int)cached_fps[i].size(); j++) {
            fp_chains.push_back( cached_fps[i][j] );
            //Store all fps to external cache for graphics
            for (int chainID =0; chainID < (int) cached_fps[i][j].size(); ++chainID) {
                fixed_points.push_back( cached_fps[i][j][chainID] );
            }
        }
        //Assemble the failed fixed points
        for (int k=0; k<(int)cached_failed[i].size(); k++) {
            failed_guesses.push_back( cached_failed[i][k] );
        }
    }
    //Clear thread storage
    for (int i=0; i<nbthreads; i++) {
        cached_fps[i].clear();
        cached_failed[i].clear();
    }
    delete[] cached_fps;
    delete[] cached_failed;
    delete amap;
    std::cout << "PMATE : Corrections process found " <<  (int)fp_chains.size() << " orbits.\n";
    
    //Print pre-filter fixed points???
}


/// Filtering the fixed points results - Remove duplicates (same and high-p) and apply symmetry property
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
fixedPointDataFilter(FPData& fpd)
{
    //Call function from FixedPointDataFilters.hpp
    filterFixedPointData<PMAP,PARAM,xavier::fixpoint>(fpd, (*theMap), (*theMapParams));
}

/// Get Fixed Point Data object of results
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
FixedPointData<xavier::fixpoint>* TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
getFixedPointData() const
{
    FixedPointData<xavier::fixpoint>* data = new FixedPointData<xavier::fixpoint>
    (theMap->rhs().desired_hamiltonian(),fp_chains);
    return data;
}

/// Write the fixed points to file
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeFixedPointResults(const char* filename) const
{
    FPData* fps = getFixedPointData();
    fixedPointDataFilter((*fps)); //Check for doubles and mirror theorem
    fps->write(filename); //Write results to file
}


/// Compose all the names off of runName member
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
composeFileNames()
{
    if (runName != "none") {
        //Rewrite the names of the files
        agndFileStr = runName + ".agnd";
        mpFileStr = runName + ".param";
        //erotFileStr = runName + ".erot";
        mapDiscontFileStr = runName + ".psi";
        fpGuessFileStr = runName + ".fpguess";
        fpFailFileStr = runName + ".fpfail";
        fpRefineFileStr = runName + ".fprefine";
        fpFileStr = runName + ".fpx";
    } //else, do nothing.
}

/// Write adaptive grid node data to file - Recommended: Set agndFileStr to desired name.
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeSamplingData()
{
    //Check if there's a filename
    if (agndFileStr == "none") {
        //If not, make one utilizing the current time
        agndFileStr = "samplingData.Run";
        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer,80,"%d-%m-%Y_%I:%M:%S",timeinfo);
        std::string tstr(buffer);
        agndFileStr = tstr + std::string(".agnd");
    }
    theSampler.write(agndFileStr.c_str());
}

/// Write Map Analysis Parameters to file - Recommended: Set mpFileStr to desired name.
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeMapParams()
{
    //Check if there's a filename
    if (mpFileStr == "none") {
        //If not, make one utilizing the current time
        mpFileStr = "mapParams.Run";
        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer,80,"%d-%m-%Y_%I:%M:%S",timeinfo);
        std::string tstr(buffer);
        mpFileStr = tstr + std::string(".param");
    }
    theMapParams->write(mpFileStr.c_str());
}

/// Write Map Discontinuities to file
template <class PMAP, class DATATYPE,
          class TRACKER, typename WVEC, typename WTRAITS,
          class CONVEXITY, class PARAM>
void TopologyExtractor<PMAP,DATATYPE,TRACKER,WVEC,WTRAITS,CONVEXITY,PARAM>::
writeMapDiscont()
{
    //Loop through all edges to collect the failure points
    EdgeSetIter eIT;
    std::vector<EdgeRotFailure> mapDisconts;
    for(eIT=edges.begin(); eIT!=edges.end(); eIT++) {
        typename std::list< EdgeRotFailure >::const_iterator failIT;
        for(failIT=eIT->rotFailures.begin(); failIT!=eIT->rotFailures.end(); ++failIT) {
            mapDisconts.push_back( *failIT );
        }
    }
    
    
    int numPts = (int) mapDisconts.size();
    
    //Build a Point Cloud for output to Avizo:
    std::vector<PointPSI> points;
    std::vector<int> ids;
    std::vector< std::vector<double> > totalData;
    std::vector< std::string > psi_labels, psi_symbols;
    psi_labels.push_back( std::string("DiscontinuityType") );
    psi_labels.push_back( std::string("Period") );
    psi_symbols.push_back( std::string("d"));
    psi_symbols.push_back( std::string("p"));
    
    //Organize the output data into the appropriate vectors
    std::vector<double> pNums, dNums;
    for(int i=0; i<numPts; i++) {
        nvis::vec2 pt = mapDisconts[i].failurePos;
        float p = mapDisconts[i].period;
        float typeVal = (float) mapDisconts[i].getTypeNumber(mapDisconts[i].type);
        points.push_back( PointPSI(pt[0],0.0,pt[1]) ); //unscaled
        ids.push_back(i);
        dNums.push_back( typeVal );
        pNums.push_back( p );
    }
    
    //Write information to PSI file
    std::vector< std::vector<double> > dataVecs;
    dataVecs.push_back( dNums );
    dataVecs.push_back( pNums );
    bool ok = PSI_WriteToFile(points,ids,dataVecs,
                              psi_labels,psi_symbols,mapDiscontFileStr.c_str());
    if(!ok) {
        throw("TopologyExtractor::Error in writing map Disconts to file!");
    }
    
    
}




}

#endif
