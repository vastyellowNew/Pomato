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
 *  FixedPointDataFilters.hpp -
 *    Filtering function for FixedPointData
 *     1) Symmetry
 *     2) Duplicate Removal
 *
 *  Author: Wayne Schlei
 *          Purdue University
 *
 *  Date:  02/01/2016
 *
 */

#ifndef PMATE_FIXED_POINT_DATA_FILTERS_HPP
#define PMATE_FIXED_POINT_DATA_FILTERS_HPP

#include <cstdio>
#include <vector>
#include <map>
#include <list>
#include <maps/fixpoints.hpp>
#include <pmate/FixedPointData.hpp>

namespace pmate {


/// Binary predicate implemented to check if two chains are the same orbit
struct IsSameFixedPointChain {
    typedef std::vector<xavier::fixpoint>  ChainType;
    IsSameFixedPointChain() : eps(9e-6) {}
    IsSameFixedPointChain(const double& e) : eps(e) {}
    
    ///Compare if two chains are the same
    bool operator() (const ChainType& first, const ChainType& second)
    {
        //If not the same period (AFTER High-Period reduction), it's not the same orbit.
        if (first[0].K != second[0].K) {
            return false;
        }
        //Now, check if any points overlap
        bool isSame = false;
        const nvis::vec2& x1 = first[0].pos;
        for(int k=0; k<(int)second.size(); k++) {
            const nvis::vec2& x2 = second[k].pos;
            nvis::vec2 dx = x2 - x1; //On-section position difference
            //Norm - Total error (which is not working in all cases)
            //if ( nvis::norm(dx) <= eps) isSame = true;
            //All components - each component is same to tolerance
            nvis::vec2 adx(fabs(dx[0]),fabs(dx[1]));
            if ( nvis::all(adx <= nvis::vec2(eps)) ) {
                isSame = true;
            }
        }
        return isSame;
    }
    
    ///Testing tolerance to see if two map-space points are the same (default=1e-6)
    double eps;
};


/// Binary predicate implemented to check if a fixed point is repeated within the same chain
struct HasLowPeriodFPDetector {
    typedef std::vector<xavier::fixpoint> ChainType;
    HasLowPeriodFPDetector() : p(-1), eps(1e-6) {}
    HasLowPeriodFPDetector(const double& e) : p(-1), eps(e) {}
    
    ///Check if fixed point chain has a lower-period version
    bool operator() (const ChainType& theChain)
    {
        const nvis::vec2& x0 = theChain[0].pos;
        int K = theChain[0].K; //Current period
        if (K==1) {
            return false;
        }
        
        for(int k=1; k<K; k++) {
            const nvis::vec2& xk = theChain[k].pos;
            if (nvis::norm(xk-x0) < eps) {
                p = k; //The actual period
                return true;
            }
        }
        return false;
    }
    ///Actual period of testing orbit (changes if low period is detected)
    int p;
    ///Testing tolerance to see if two map-space points are the same (default=1e-6)
    double eps;
};

/// Binary predicate for checking if fixed point has symmetry
struct FixedPointSymmetryDetector {
    typedef std::vector<xavier::fixpoint> ChainType;
    FixedPointSymmetryDetector() : eps(1e-6) {}
    FixedPointSymmetryDetector(const double& e) : eps(e) {}
    
    ///Check if fixed point chain is symmetric ([x,0] or symmetric pair)
    bool operator() (const ChainType& theChain)
    {
        int p = theChain[0].K;
        nvis::vec2 x0mirror = theChain[0].pos;
        x0mirror[1] *= -1.0; //Mirror Theorem for CR3BP
        for(int k=0; k<p; k++) {
            const nvis::vec2& xk = theChain[k].pos;
            //Perpendicular crossing test
            if (std::fabs(xk[1]) < eps) {
                return true;
            }
            
            //Symmetric pair test
            if( (k>0) && (nvis::norm(xk-x0mirror) < eps) ) {
                return true;
            }
            
        }
        return false;
    }
    
    ///Testing tolerance to see if two map-space points are the same (default=1e-6)
    double eps;
};


///Comparing fixed point chains (For sort: 1) period)
bool fpChainPeriodCompare(const std::vector<xavier::fixpoint>& lhs,
                          const std::vector<xavier::fixpoint>& rhs)
{
    return (lhs[0].K < rhs[0].K);
}

///Comparing fixed point chains (For sort: 1) period 2) type 3) |stability_index| 4) position)
bool fpChainAdvancedCompare(const std::vector<xavier::fixpoint>& lhs,
                            const std::vector<xavier::fixpoint>& rhs)
{
    //First by period
    if (lhs[0].K != rhs[0].K) {
        return lhs[0].K < rhs[0].K;
    } else {
        //Next by type (center then saddle)
        if (lhs[0].saddle != rhs[0].saddle ) {
            return ( !(lhs[0].saddle) );
        } else {
            //If saddles, sort by magnitude of stability index
            if (lhs[0].saddle && (fabs(lhs[0].si) != fabs(rhs[0].si)) ) {
                return (fabs(lhs[0].si) < fabs(rhs[0].si) );
            } else {
                //Otherwise, sort by position
                nvis::lexicographical_order compare;
                return compare(lhs[0].pos,rhs[0].pos);
            }
        }
    }
}


/// Custom 'uniquify' command to pull out repeated fixed point chains in a list
void uniquifyFPList( std::list< std::vector<xavier::fixpoint> >& chainList, const double& eps)
{
    typedef std::vector<xavier::fixpoint>  ChainType;
    typedef std::list<ChainType>::iterator ChainListIterator;
    ChainListIterator lit = chainList.begin();
    IsSameFixedPointChain sameChecker(eps);
    //Search entire list for copies
    while( lit != chainList.end() ) {
        //Examine each remaining chain for a duplicate
        ChainListIterator lit2 = lit;
        lit2++;
        for(; lit2 != chainList.end(); lit2++) {
            //Evaluate the comparitive struct
            if ( sameChecker((*lit),(*lit2)) ) {
                //If the same, erase the chain
                lit2 = chainList.erase(lit2);
                lit2--; //Move back to the last viable element
            }
        }
        //Increment chain
        lit++;
    }
}

/// Filter fixed points in FixedPointData<FP>
template <typename PMAP,typename PARAMS,typename FP>
void filterFixedPointData(
    FixedPointData<FP>& fpData,
    const PMAP& theMap,
    const PARAMS& theMapParams,
    const double eps = 5.e-6,
    const bool highPeriodCutoff = false,
    const bool enableBoxFilter = false
)
{
    typedef std::vector<FP>      ChainType;
    
    int numOrbits = fpData.getNumOrbits();
    std::vector< ChainType > theChains;
    fpData.getData( theChains );
    double jc1 = fpData.getJacobiConstant();
    
    //The 3 filtering steps:
    // 1) Symmetry - create the symmetric pairings if they don't exist
    // 2) Reduction - reduce high-p to the low-p equivalents
    // 3) Removal - remove duplicates
    // 4) Cutoff Period - drop to only specified period in theMapParams
    if( theMap.section().isSymmetric() ) { //Only do symmetry portion if it applies to section
        //Run Mirror theorem to create the missing Asymmetric pairs
        //Functor for finding symmetric orbits from fixed points
        FixedPointSymmetryDetector symDetector(eps);
        //Find the asym orbits
        std::vector<int> asymIDs;
        std::vector< ChainType > mirrors;
        int orbitIDTally = (int) theChains.size();
        for(int k=0; k<numOrbits; k++) {
            bool isSymmetric = symDetector( theChains[k] );
            if (!isSymmetric) {
                orbitIDTally++;
                asymIDs.push_back(k);
                ChainType counterPart = theChains[k];
                typename ChainType::iterator it;
                for(it=counterPart.begin(); it!=counterPart.end(); ++it) {
                    //Hard coded to CR3BP!   Future Work: create "mirror()" function for states in template.
                    it->pos[1] *= -1.0;
                    it->evec[0][1] *= -1.0;
                    it->evec[1][1] *= -1.0;
                    it->orbitID = orbitIDTally;
                }
                //Note: we need to flip the order and timing of points
                ChainType counterPartOrdered = counterPart;
                int period = counterPart[0].K;
                double tp = counterPart[0].timePeriod;
                if(period > 0) {
                    for(int j=period-1,newj=1; j>0; j--) {
                        //Put in reverse order and change the times
                        counterPartOrdered[newj] = counterPart[j];
                        counterPartOrdered[newj].t = tp - counterPart[j].t;
                        newj++;
                    }
                }
                //Store the ordered complimentary mirror orbit
                mirrors.push_back( counterPartOrdered );
            }
        }
        //Insert into full data
        typename std::vector<ChainType>::iterator mit;
        for(mit=mirrors.begin(); mit!=mirrors.end(); mit++) {
            theChains.push_back( (*mit) );
        }
    }
    
    //Sorting out duplicates
    // First, reduce high-order chains to the correct low-order p's
    // (e.g., p=9 -> p=3, if applicable)
    HasLowPeriodFPDetector lowPeriodChecker(eps);
    //Update total
    numOrbits = (int) theChains.size();
    for(int k=0; k<numOrbits; k++) {
        bool lowPeriodThere = lowPeriodChecker( theChains[k] );
        //If detected, reduce to low period version
        if(lowPeriodThere) {
            int p = lowPeriodChecker.p;
            int originalP = theChains[k][0].K;
            int scaleFactor = originalP/p;  //Should be integer multiple
            theChains[k].erase(theChains[k].begin()+p,theChains[k].end());
            //Reset the period,time, and timePeriod at each fp
            typename ChainType::iterator it;
            double t0 = theChains[k][0].t;
            for(it=theChains[k].begin(); it!=theChains[k].end(); ++it) {
                it->K = p;
                it->timePeriod /= (double) scaleFactor;
                it->t -= t0;
            }
        }
    }
    
    
    // Next, remove duplicates
    //Create a list of the fixed point chains
    std::list<ChainType> chainList;
    numOrbits = (int) theChains.size();
    for(int i=0; i<numOrbits; i++) {
        chainList.push_back(theChains[i]);
    }
    //Sort the list using custom advanced compare [period/type/si/pos]
    chainList.sort(fpChainAdvancedCompare);
    //Uniquify using comparitive structure
    chainList.unique( IsSameFixedPointChain(eps) );
    //Uniquify using custom function
    uniquifyFPList( chainList, eps );
    
    //Copy to new data while checking for period cutoff
    typename std::list<ChainType>::iterator cIT;
    theChains.clear();
    for(cIT=chainList.begin(); cIT!=chainList.end(); cIT++) {
        bool addChain = true;
        //High-period cutoff
        if(highPeriodCutoff) {
            if( (*cIT)[0].K > theMapParams.max_period ) {
                addChain = false;
            }
        }
        //Box filter (on whole chain)
        if(addChain && enableBoxFilter) {
            for(int i=0; i<(int)(*cIT).size(); i++) {
                if ( !(theMapParams.the_metric.bounds().inside( (*cIT)[i].pos )) ) {
                    addChain = false;
                }
            }
        }
        //Time period cutoff (Hard coded! - Prevent from running manifolds for Horseshoe orbits)
        if(addChain) {
            if( (*cIT)[0].timePeriod > 2000.0 ) {
                addChain = false;
            }
        }
        
        if(addChain) {
            //Add enabled
            theChains.push_back( (*cIT) );
        }
        
    }
    
    
    //Reset the orbitID data values to indicate the chain index :
    numOrbits = (int) theChains.size();
    for(int k=0; k<numOrbits; k++) {
        for(int i=0; i<(int)theChains[k].size(); i++) {
            theChains[k][i].orbitID = k;
        }
    }
    
    // Move data into FixedPointData object
    fpData.clearData();
    fpData.insertData(theChains);
    
    // Done with this info
    theChains.clear();
}


} //end pmate

#endif