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


/** Functor for checking whether or not to subdivide a cell.
 *  - Figures out cell convexity wrt winding numbers but
 *    ignores the x-ydot number.
 *  - Also determines if a cell resides completely within a primary
 *    and marks it as invalid within the adaptive grid.
 *
 *  Author:  Wayne Schlei (Purdue University)
*/

#ifndef __CR3BP_CELL_CHECKER_HPP__
#define __CR3BP_CELL_CHECKER_HPP__

#include <vector>
#include <exception>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <map>
#include <iomanip>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <data/adaptive_grid.hpp>
#include <maps/cellConvexityChecker.hpp>

namespace xavier {

///Check if a cell has convexity with winding number (or is invalid)
//->return false means this cell must be subdivided
// This is outside of AdaptiveGrid because this is unique to map topology analysis
template<typename DATAGRID, typename IDTYPE, typename WINDING_VEC>
class CR3BP_Convexity : public CellConvexityFunctor<DATAGRID,IDTYPE,WINDING_VEC> {
    typedef CR3BP_Convexity<DATAGRID,IDTYPE,WINDING_VEC> self_type;
    typedef CellConvexityFunctor<DATAGRID,IDTYPE,WINDING_VEC> base_type;
public :
    CR3BP_Convexity(DATAGRID& grid, const WINDING_VEC& tolFactor, const WINDING_VEC& maxDist) :
        base_type(grid,tolFactor,maxDist),
        r1(0.0), r2(0.0)
    {}
    
    CR3BP_Convexity(const self_type& c) :
        base_type((*c.theGrid),c.tols,c.maxDists),
        r1(c.r1), r2(c.r2), mup(c.mup)
    {}
    
    /// Create a clone
    self_type* clone() const
    {
        return new self_type(*this);
    }
    
    /** Modulo function for w_{x,ydot}
     * With x-ydot number, we have a continuous spectrum where +inf maps to -inf.
     * If we limit the range to +/- wmax, we should also generate a moduluo
     * function that wraps around this limit.
     */
    double wxydModulo(const double& wxyd)
    {
        return wxyd + 250.0;
    }
    /** Range function for w_{x,ydot}
     *  Used to map value to an indicated allowed range
     */
    double wxydRanged(const double& wxyd, double wmax = 1000.0)
    {
        double wValue = wxyd;
        if( wValue > wmax ) {
            wValue = wmax;
        }
        if( wValue < -wmax ) {
            wValue = -wmax;
        }
        return wValue;
    }
    
    /// Returns false if subdivision is needed
    bool operator()(const IDTYPE& cell_id)
    {
        typedef typename DATAGRID::cell_data_type        cell_data_type;
        typedef typename DATAGRID::pos_type              pos_type;
        typedef typename DATAGRID::bounds_type           bounds_type;
        typedef typename DATAGRID::traits_type           ValueTraits;
        
        const cell_data_type& data = this->theGrid->getCellData(cell_id);
        
        ///For CR3BP, we want to ignore cells that are entirely within a primary.
        ///Run check and set the validity flag in grid before
        pos_type lowLeft = this->theGrid->getVertex(cell_id);
        pos_type upRight = this->theGrid->getVertex(cell_id + IDTYPE(1,1,0));
        //bounds_type bbox = this->theGrid->cellBounds(cell_id); //seems buggy
        this->indexMarker = false; //By default, cells are valid
        
        //Output for debugging
        //std::cerr << "Cell " << cell_id << ":\n";
        //std::cerr << "   Bounds : " << bbox << "\n";
        //std::cerr << "   mup = " << mup << " r1 = " << r1 << " r2 = " << r2 << "\n";
        //std::cerr << "   Upper right: " << upRight << " |<=? (-mup+r1,0) " << nvis::vec2(-mup+r1,0.0) << "\n";
        //std::cerr << "   Lower left : " << lowLeft << " |>=? (-mup-r1,0) " << nvis::vec2(-mup-r1,0.0) << "\n";
        
        //Inside P1
        if ( (upRight[0] <= (-mup+r1)) && (lowLeft[0] >= (-mup-r1)) ) {
            //Indicate we need to mark this cell as invalid (outside this functor)
            this->indexMarker = true;
            return true; //We don't want to subdivide invalid cells
        }
        //Inside P2
        if ( (upRight[0] <= (1.0-mup+r2)) && (lowLeft[0] >= (1.0-mup-r2)) ) {
            this->indexMarker = true; //Mark invalid
            return true; //Don't subdivide invalid cells
        }
        //Fails all corners check -> mark as invalid so we don't subdivide next time through
        if ( this->theGrid->allCornersInvalid(cell_id) ) {
            this->indexMarker = true; //Mark invalid
            return true; //Don't subdivide invalid cells
        }
        
        ///Check cell validity : All invalid cells are non-convex
        if ( !( this->theGrid->isCellValid(cell_id) ) ) {
            return false;
        }
        
        std::vector<WINDING_VEC> corner_values(4);
        corner_values[0] = this->theGrid->getVertexValue(cell_id);
        corner_values[1] = this->theGrid->getVertexValue(cell_id + IDTYPE(1,0,0));
        corner_values[2] = this->theGrid->getVertexValue(cell_id + IDTYPE(1,1,0));
        corner_values[3] = this->theGrid->getVertexValue(cell_id + IDTYPE(0,1,0));
        
        ///Will check the winding number(or vector) to see if the change is less than given tolerance
        int numWNValues = (int) corner_values[0].size();
        for (int i=0; i<numWNValues; i++) {
            std::vector<double> wValues(4);
            for(int j=0; j<4; j++) {
                wValues[j] = corner_values[j][i];
            }
            
            //Special Check for x-ydot number in CR3BP
            if (i==1) {
                //First, map values to prescribed range (-wmax,wmax)
                for(int j=0; j<4; j++) {
                    wValues[j] = wxydRanged(wValues[j]);
                }
                
                //Compute the max distance [Another way could use modulo function]
                double maxDistance = 0.0;
                for(int j=0; j<4; j++) {
                    for(int jp=j+1; jp<4; jp++) {
                        //Test the three possiblities for largest value
                        double dist = std::fabs(wValues[j]-wValues[jp]);
                        if (dist > maxDistance) {
                            maxDistance = dist;
                        }
                    }
                }
                //Distance check - encompasses too much variance
                if ( maxDistance > this->maxDists[i]) {
                    return false;
                }
                
                //Linearity check - internal points follow convex hull
                double eps = (this->tols[i])*(this->maxDists[i]);
                
                //Min and Max at corners
                double min = *std::min_element(wValues.begin(), wValues.end());
                double max = *std::max_element(wValues.begin(), wValues.end());
                min -= eps;
                max += eps;
                typedef typename cell_data_type::const_iterator iterator_type;
                for (iterator_type it=data.begin() ; it!=data.end() ; ++it) {
                    //If the value is valid
                    if (it->second != ValueTraits::invalid) {
                        double value = it->second[i];
                        //Normalize to range
                        value = wxydRanged(value);
                        if (value<min || value>max) {
                            return false;
                        }
                    }
                }
                
            } else {
                //Min and Max at corners
                double min = *std::min_element(wValues.begin(), wValues.end());
                double max = *std::max_element(wValues.begin(), wValues.end());
                
                //Distance check - encompasses too much variance
                if ((max-min) > this->maxDists[i]) {
                    return false;
                }
                
                //Linearity check - internal points follow convex hull
                double eps = (this->tols[i])*(this->maxDists[i]);
                min -= eps;
                max += eps;
                typedef typename cell_data_type::const_iterator iterator_type;
                for (iterator_type it=data.begin() ; it!=data.end() ; ++it) {
                    //If the value is valid
                    if (it->second != ValueTraits::invalid) {
                        double value = it->second[i];
                        if (value<min || value>max) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    } //End operator()
    
    /// Maximum depth validity check -> runs at max depth (end of refinement) to mark remaining invalid cells
    bool validityAfterRefinement(const IDTYPE& cell_id)
    {
        typedef typename DATAGRID::cell_data_type         cell_data_type;
        typedef typename DATAGRID::pos_type               pos_type;
        typedef typename DATAGRID::bounds_type            bounds_type;
        typedef typename DATAGRID::traits_type            ValueTraits;
        
        //const cell_data_type& data = this->theGrid->getCellData(cell_id);
        
        // Check if the cell is valid in the first place (And do nothing)
        if ( !(this->theGrid->isCellValid(cell_id)) ) {
            return false;
        }
        
        ///For CR3BP, we want to ignore cells that are entirely within a primary,
        /// or they have an edge inside P1 for sure (can also do P2 if you desire)
        pos_type lowLeft = this->theGrid->getVertex(cell_id);
        pos_type upRight = this->theGrid->getVertex(cell_id + IDTYPE(1,1,0));
        
        
        // Checks on cells inside primaries - mark invalid if found
        if ( (upRight[0] <= (-mup+r1)) && (upRight[0] >= (-mup-r1)) ) {
            //Indicate this cell as invalid
            this->theGrid->setCellValidity(cell_id, false);
            return false; //Cell crosses P1
        } else if ( (lowLeft[0] <=(-mup+r1)) && (lowLeft[0] >=(-mup-r1)) ) {
            //Indicate this cell as invalid
            this->theGrid->setCellValidity(cell_id, false);
            return false; //Cell crosses P1
        } else if ( (upRight[0] <= (1.0-mup+r2)) && (lowLeft[0] >= (1.0-mup-r2)) ) {
            //Indicate this cell as invalid
            this->theGrid->setCellValidity(cell_id, false);
            return false; //Whole cell inside P2
        }
        
        //Cell is valid
        return true;
    }
    
    /// Set the radii for P1,P2 in nondim values
    void setBodyRadii(const double rP1, const double rP2)
    {
        r1 = rP1;
        r2 = rP2;
    };
    /// Set system gravity parameter
    void setMu(const double mu)
    {
        mup = mu;
    }
    
    /// Body radii in nondim units
    double r1, r2;
    /// System mu
    double mup;
};

} //End xavier

#endif
