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


#ifndef __CELL_CONVEXITY_CHECK_HPP__
#define __CELL_CONVEXITY_CHECK_HPP__

#include <iostream>
#include <iomanip>
#include <limits>
#include <assert.h>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>

#include <data/adaptive_grid.hpp>

namespace xavier {

/** Functor for checking cell convexity with respect to winding numbers.  This tells
 *  the AdaptiveGrid class (in DATAGRID) whether or not to refine a given cell
 *  to the tolerances specified (per winding number). The "maxDist" tolerance
 *  indicates the maximum allowable gap in values at the vertices (or deltaW = Wmax-Wmin).
 *  The "tolFactor" tolerance applies to values within the cell fitting
 *  Wmin-tolFactor <= W_internal <= Wmax + tolFactor.  This function also determines
 *  if a cell has all valid values at the corners.  "false" is returned if it fails
 *  to meet the criteria mentioned.
 *
 *  Create a subclass of this for different problems
 */
template<typename DATAGRID, typename IDTYPE, typename WINDING_VEC>
class CellConvexityFunctor {
public :
    typedef CellConvexityFunctor<DATAGRID,IDTYPE,WINDING_VEC> self_type;
    typedef IDTYPE                                            id_type;
    ///Constructor
    CellConvexityFunctor(DATAGRID& grid, const WINDING_VEC& tolFactor, const WINDING_VEC& maxDist) :
        theGrid(&grid), tols(tolFactor), maxDists(maxDist), indexMarker(false) {}
    ///Copy constructor
    CellConvexityFunctor(const self_type& c) :
        theGrid((*c.theGrid)), tols(c.tols), maxDists(c.maxDists), indexMarker(c.indexMarker) {}
        
    /// Create a clone
    self_type* clone() const
    {
        return new self_type(*this);
    }
    
    /// Operator() runs the convexity check for a given cell (Overload for given problem)
    bool operator()(const IDTYPE& cell_id)
    {
        //Default -> you may have to overload this in your implementation
        typedef typename DATAGRID::cell_data_type     cell_data_type;
        typedef typename DATAGRID::traits_type        ValueTraits;
        
        const cell_data_type& data = theGrid->getCellData(cell_id);
        
        ///Check cell validity : All invalid cells are non-convex
        if ( !( theGrid->isCellValid(cell_id) ) ) {
            return false;
        }
        
        std::vector<WINDING_VEC> corner_values(4);
        corner_values[0] = theGrid->getVertexValue(cell_id);
        corner_values[1] = theGrid->getVertexValue(cell_id + IDTYPE(1,0,0));
        corner_values[2] = theGrid->getVertexValue(cell_id + IDTYPE(1,1,0));
        corner_values[3] = theGrid->getVertexValue(cell_id + IDTYPE(0,1,0));
        
        ///Will check the winding number(or vector) to see if the change is less than given tolerance
        int numWNValues = corner_values[0].size();
        for (int i=0; i<numWNValues; i++) {
            std::vector<double> wValues(4);
            for(int j=0; j<4; j++) {
                wValues[j] = corner_values[j][i];
            }
            double min = *std::min_element(wValues.begin(), wValues.end());
            double max = *std::max_element(wValues.begin(), wValues.end());
            
            //Distance check - encompasses too much variance
            if ((max-min) > maxDists[i]) {
                return false;
            }
            
            //Linearity check - internal points follow convex hull
            double eps = tols[i]*maxDists[i];
            min -= eps;
            max += eps;
            typedef typename cell_data_type::const_iterator iterator_type;
            for (iterator_type it=data.begin() ; it!=data.end() ; ++it) {
                //If the value is valid (sometimes invalids pop up within the valid grid cells)
                if (it->second != ValueTraits::invalid) {
                    double value = it->second[i];
                    if (value<min || value>max) {
                        return false;
                    }
                }
                //Otherwise throw an error
                else {
                    std::cerr << "Convexity check trying to reference an invalid value!\n";
                    //print output
                    std::cerr << " Cell ID : " << cell_id << "\n";
                    std::cerr << " Calling point : " << it->first << "\n";
                    std::cerr << " Winding Value : " << it->second << "\n";
                    //Throw error
                    //throw std::runtime_error("Invalid value referenced in AdaptiveGrid");
                }
            }
        }
        return true;
    }// end operator()
    
    /// Maximum depth validity check -> runs at max depth (end of refinement) to mark remaining invalid cells
    bool validityAfterRefinement(const IDTYPE& cell_id)
    {
        ///Return true for a valid cell, false for invalid
        return theGrid->isCellValid(cell_id);
        /// You can also custom define a validity setting here.  For example, in the CR3BP,
        /// cells that are entirely within a primary can be set as invalid here with the
        /// call:  theGrid->setCellValidity(cell_id,false);
    }
    
    /// Function to test if all input cells satisfy convexity constraints
    bool allConvexCells(const std::vector<IDTYPE>& cellIDs)
    {
        typedef typename std::vector<IDTYPE>::const_iterator IDConstIterator;
        IDConstIterator cit;
        //If any cell in list is NOT convex
        for(cit=cellIDs.begin(); cit!=cellIDs.end(); ++cit) if( !((*this)(*cit)) ) {
                return false;
            }
        //Otherwise they are all satisfy convexity constraints
        return true;
    }
    
    ///Grid Structure pointer (AdaptiveGrid<>)
    DATAGRID* theGrid;
    ///Tolerance factors on maxDistance parameters
    WINDING_VEC tols;
    ///Maximum allowable distance of winding vector between cell corners
    WINDING_VEC maxDists;
    ///An index marker that can be used as a flag per operator() call
    bool indexMarker;
};


/** Comparison operator for priority_queue processing during subdivision
 *  Purpose:  Essentially this sorts cells for processing such that cells
 *  at the top of the queue are the most in need of refinement.
 */
template <class CONVEXITY_FUNCTOR>
class CompareCellConvexity {
    typedef typename CONVEXITY_FUNCTOR::id_type   id_type;
    
public:
    /// Constructor
    CompareCellConvexity(CONVEXITY_FUNCTOR& cFunctor) :
        cfunc(&cFunctor)
    { }
    
    ///Set the functor object
    //void set(const CONVEXITY_FUNCTOR& cFunctor)
    //{ cfunc = &cFunctor; }
    
    /** Sorting operator for priority_queue
    *  Priority is (in order from high to low) invalid cells first, partially
    *  invalid cells, non-convex with low depth (like 0), non-convex with
    *  higher depth, convex cells with low depth (like 0), then
    *  convex cells at high depth.
    *
    */
    bool operator()(const id_type& cell1, const id_type& cell2) const
    {
        //Returns if cell1 has lower priority than cell2:
        
        //Completely invalid
        if ((cfunc->theGrid->isCellValid(cell1)) &&
                !(cfunc->theGrid->isCellValid(cell2))) {
            return true;
        } else if (!(cfunc->theGrid->isCellValid(cell1)) &&
                   (cfunc->theGrid->isCellValid(cell2))) {
            return false;
        }
        
        //Partially valid (some corners are outside bounds)
        else if ((cfunc->theGrid->allCornersValid(cell1)) &&
                 !(cfunc->theGrid->allCornersValid(cell2))) {
            return true;
        } else if (!(cfunc->theGrid->allCornersValid(cell1)) &&
                   (cfunc->theGrid->allCornersValid(cell2))) {
            return false;
        }
        
        //Convex winding number (yes or no)
        else if ( ((*cfunc)(cell1)) && !((*cfunc)(cell2)) ) {
            return true;
        } else if ( !((*cfunc)(cell1)) && ((*cfunc)(cell2)) ) {
            return false;
        }
        
        //Depth (or size of cell) - low depth is higher priority
        else {
            int dIdx = (int) cell1.size() -1;
            return ( cell1[dIdx] > cell2[dIdx] );
        }
    }
    
private:
    CONVEXITY_FUNCTOR* cfunc;
};


} //end xavier

#endif
