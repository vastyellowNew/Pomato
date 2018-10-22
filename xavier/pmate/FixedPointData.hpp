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
 *  FixedPointData.hpp - a Header for a Fixed point data class
 *
 *  Author: Wayne Schlei
 *          Purdue University
 *
 *  Date:  02/01/2016
 *

 */

#ifndef PMATE_FIXED_POINT_DATA_HPP
#define PMATE_FIXED_POINT_DATA_HPP

#include <cstdio>
#include <cassert>
#include <vector>
#include <map>
#include <list>
#include <maps/fixpoints.hpp>

namespace pmate {

/** Fixed point data storage class
 *  Avizo Note: This is replicated in the 'poincare' package data class FixedPointData.
*/
template <typename FP>
class FixedPointData {
public :
    typedef FP                                       FixedPointType;
    typedef typename FP::VecType                     VecType;
    typedef typename FP::StateType                   StateType;
    typedef typename std::vector<FP>                 FixedPointChain;
    typedef typename std::vector<std::vector<FP> >   FPChainsContainer;
    /// Constructor (blank)
    FixedPointData();
    /// Constructor with input fixed point chains
    FixedPointData(const double& jc, const FPChainsContainer& chainsIN);
    /// Constructor with filename
    FixedPointData(const char* filename);
    
    
    /// Get the Jacobi Constant
    double getJacobiConstant() const
    {
        return jacobi;
    }
    
    /// Set the Jacobi Constant
    void setJacobiConstant(const double& jc);
    
    /// Clear the current data (fp_chains)
    void clearData();
    
    /// Insert a fixpoint chain into the data
    void insertOrbit(const FixedPointChain& orbit);
    
    /// Insert a whole set of data
    void insertData(const FPChainsContainer& data);
    
    /// Get a copy of the fixed point data
    void getData(FPChainsContainer& data) const;
    
    /// Returns the number of fixed points
    int getNumFixedPoints() const;
    
    /// Get the fixed point information given an index
    const FP& getFixedPoint(const int orbitID, const int fpID) const;
    
    /// Get a fixed point given a single index (works through all orbits)
    const FP& getFixedPoint(const int id) const;
    
    /// Returns the number of orbits
    int getNumOrbits() const
    {
        return (int) fp_chains.size();
    }
    
    /// Get the maximum (fp.orbitID) from data
    int getMaxOrbitID() const;
    
    /// Get the minimum (fp.orbitID) from data
    int getMinOrbitID() const;
    
    /// Get the orbit id classifier (fp.orbitID) from data
    int getOrbitID(const int orbitID, const int fpID) const;
    
    /// Get an orbit (chain) from set (if id>fp_chains.size(), returns 1st orbit)
    const FixedPointChain& getOrbit(const int orbitID) const;
    
    /// Get the period of an orbit indicated by the orbit index
    int getOrbitPeriod(const int orbitID) const;
    
    /// Get the orbit index from point index (not fp.orbitID)
    int getOrbitIndex(const int ptIdx) const;
    
    /// Is the selected orbit a saddle
    bool isSaddle(const int orbitID) const;
    
    /// Get an orbit (chain) from point index
    void getOrbitFromPointIndex(const int ptIdx, FixedPointChain& fp) const;
    
    ///Comparing fixed points by period
    bool fpPeriodCompare(const FP& lhs, const FP& rhs) const
    {
        return (lhs.K < rhs.K);
    }
    
    ///Comparing fixed point chains by period
    bool fpPeriodCompare(const FixedPointChain& lhs,
                         const FixedPointChain& rhs) const
    {
        return (lhs[0].K < rhs[0].K);
    }
    
    
    /// Reading data from file
    bool read(const char* filename);
    /// Writing data to file
    bool write(const char* filename) const;
    
private :
    ///Fixed Point data container
    std::vector< FixedPointChain > fp_chains;
    ///Jacobi constant of fixed points
    double jacobi;
    /// Hash table to link pointIndex to (orbitIndex,crossing)
    std::map<int,std::pair<int,int> > indexMap;
    /// Build the index maps AFTER new data is entered
    void buildIndexMap();
    
};

/// Constructor (blank)
template <typename FP>
FixedPointData<FP>::FixedPointData() :
    jacobi(2.96)
{}

/// Constructor with input fixed point chains
template <typename FP>
FixedPointData<FP>::
FixedPointData(const double& jc, const FPChainsContainer& chainsIN)
{
    fp_chains = chainsIN;
    jacobi = jc;
}

/// Constructor with filename
template <typename FP>
FixedPointData<FP>::
FixedPointData(const char* filename)
{
    read(filename);
}



/// Write this object to indicated file
template <typename FP>
bool FixedPointData<FP>::write(const char* filename) const
{
    FILE* f = fopen(filename, "w"); // open the file
    
    if (!f) {
        std::cerr << "FixedPointData:: Unable to write to file = " << filename << "!\n";
        return false; // indicate error
    }
    
    // Write header:
    fprintf(f, "# Fixed Point Data\n");
    
    // Write number of fixed point chains
    int numFPs = getNumOrbits();
    int spaceDim = 2; //Fixed for now
    fprintf(f, "%d %d %.15f\n", spaceDim, numFPs, getJacobiConstant());
    
    // Write fixed point data values:
    // period fixpoint0 fixpoint1 ... fixpointP
    for (int k=0; k<numFPs; k++) {
        FixedPointChain chain;
        chain = getOrbit(k);
        //9 Commons per fixed point
        int period = chain[0].K;
        int orbitID = chain[0].orbitID;
        double timePeriod = chain[0].timePeriod;
        double si = chain[0].si, otherSI = chain[0].otherSI;
        //Note, this is hard coded
        double evals[2];
        evals[0] = chain[0].eval[0];
        evals[1] = chain[0].eval[1];
        int saddle = (chain[0].saddle) ? 1 : 0;
        int isolated = (chain[0].isolated) ? 1 : 0;
        //Write commons to first line
        fprintf(f,"%d %d %d %d %.15f %.15f %.15f %.15f %.15f\n",
                period,saddle,orbitID,isolated,timePeriod,evals[0],evals[1],si,otherSI);
                
        //Point information - One line per indicated point in chain
        for (int i=0; i<period; i++) {
            nvis::vec2& pos = chain[i].pos; //hard coded
            //fixpoint print: t (x xdot) (evec0) (evec1)
            fprintf(f,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
                    chain[i].t, pos[0], pos[1],
                    chain[i].evec[0][0], chain[i].evec[0][1],
                    chain[i].evec[1][0], chain[i].evec[1][1]);
        }
    }
    
    // Close the file.
    fclose(f);
    
    return true; // indicate success
}

/// Read and popluate this object
template <typename FP>
bool FixedPointData<FP>::read(const char* filename)
{
    FILE* f = fopen(filename, "r"); // open the file
    
    if (!f) {
        std::cerr << "FixedPointData:: Unable to read from file = " << filename << "!\n";
        return false; // indicate error
    }
    
    // Skip header (first line). We could do some checking here:
    char buf[80];
    fgets(buf, 80, f);
    
    // Read size of fixed point array
    int spaceDim, numChains;
    fscanf(f, "%d %d %lf", &spaceDim, &numChains, &jacobi);
    
    // Do some consistency checking.
    if (spaceDim != 2) {
        std::cerr << "FixedPointData Read Error: Can't support section space other than dimension=2 at present.\n"
                  << "Error in file = " << filename << ".\n";
        fclose(f);
        return false;
    }
    
    
    // Now we have to read data values for each chain
    fp_chains.clear();
    for (int c=0; c<numChains; c++) {
        FixedPointChain chain;
        int p = 0, saddle=0, orbitID=0, isolated=0;
        double timePeriod,si,otherSI, eval[2];
        //Scan for 9 value line
        fscanf(f,"%d %d %d %d %lf %lf %lf %lf %lf",
               &p,&saddle,&orbitID,&isolated,&timePeriod,&(eval[0]),&(eval[1]),
               &si,&otherSI);
        //Scan for all fixpoint info
        for (int i=0; i<p; i++) {
            FP fp;
            //Commons
            fp.K = p;
            fp.saddle = (saddle==1)? true : false;
            fp.isolated = (isolated==1)? true : false;
            fp.eval[0] = eval[0];
            fp.eval[1] = eval[1];
            fp.orbitID = orbitID;
            fp.timePeriod = timePeriod;
            fp.si = si;
            fp.otherSI = otherSI;
            //Read rest from file
            double t = 0.0;
            nvis::vec2 pos, evec[2];
            fscanf(f,"%lf %lf %lf %lf %lf %lf %lf",
                   &t, &(pos[0]), &(pos[1]),
                   &(evec[0][0]), &(evec[0][1]),
                   &(evec[1][0]), &(evec[1][1]));
            fp.pos = pos;
            fp.t = t;
            fp.evec[0] = evec[0];
            fp.evec[1] = evec[1];
            chain.push_back(fp);
        }
        //Store the chain in data object
        fp_chains.push_back( chain );
    }
    // We are done with reading, close the file.
    fclose(f);
    
    return true; // indicate success
}


/// Set the Jacobi Constant
template <typename FP>
void FixedPointData<FP>::setJacobiConstant(const double& jc)
{
    jacobi = jc;
}

/// Build the index maps AFTER data is entered
template <typename FP>
void FixedPointData<FP>::buildIndexMap()
{
    //Update/build index maps
    indexMap.clear();
    typename std::vector< FixedPointChain >::iterator it;
    int chainID=0, idx = 0;
    for(it=fp_chains.begin(); it!=fp_chains.end(); ++it) {
        //Build the indexMap
        int p = (*it)[0].K;
        for(int i=0; i<p; i++) {
            std::pair<int,int> indexPair(chainID,i);
            indexMap.insert( std::pair<int,std::pair<int,int> >( idx, indexPair ) );
            idx++; //Total point index
        }
        chainID++; //increment to next fp_chain
    }
}


/// Clear the current data (fp_chains)
template <typename FP>
void FixedPointData<FP>::clearData()
{
    fp_chains.clear();
    indexMap.clear();
}

/// Insert a fixpoint chain into the data
template <typename FP>
void FixedPointData<FP>::insertOrbit(const FixedPointChain& orbit)
{
    fp_chains.push_back( orbit );
    buildIndexMap();
}

/// Insert a whole set of data
template <typename FP>
void FixedPointData<FP>::insertData(const FPChainsContainer& data)
{
    for (int k=0; k<(int) data.size(); k++) {
        fp_chains.push_back( data[k] );
    }
    buildIndexMap();
}

///Get the data as a copy
template <typename FP>
void FixedPointData<FP>::getData(FPChainsContainer& data) const
{
    data = fp_chains;
}

template <typename FP>
int FixedPointData<FP>::getNumFixedPoints() const
{
    typename std::vector< FixedPointChain >::const_iterator it;
    int numFP = 0;
    for(it=fp_chains.begin(); it!=fp_chains.end(); ++it) {
        numFP += (int) it->size();
    }
    return numFP;
}

/// Get the fixed point information given an index (Reference to actual storage location)
template <typename FP>
const FP& FixedPointData<FP>::getFixedPoint(const int orbitID, const int fpID) const
{
    assert( orbitID>=0 && orbitID < (int) fp_chains.size() );
    assert( fpID>=0 && fpID<(int) fp_chains[orbitID].size());
    return fp_chains[orbitID][fpID];
}

/// Get the period of an orbit indicated by the orbit index
template <typename FP>
int FixedPointData<FP>::getOrbitPeriod(const int orbitID) const
{
    std::vector<FP> orbit;
    orbit = getOrbit(orbitID);
    return (int)orbit.size();
}


/// Get a fixed point given a single index (works through all orbits
template <typename FP>
const FP& FixedPointData<FP>::getFixedPoint(const int id) const
{
    std::map<int,std::pair<int,int> >::const_iterator mit = indexMap.find( id );
    int theID = 0, cID = 0;
    if (mit != indexMap.end()) {
        cID = mit->second.first;
        theID = mit->second.second;
        //theMsg->printf("%s: indexMap[%d] = [%d,%d]", __FILE__,id,cID,theID);
    } else {
        std::cerr << __FILE__ << ": Failed to lookup point id " << id << " in Orbit List! Returning 0\n";
    }
    
    return fp_chains[cID][theID];
}

/// Is the selected orbit a saddle
template <typename FP>
bool FixedPointData<FP>::isSaddle(const int orbitID) const
{
    return fp_chains[orbitID][0].saddle;
}

/// Get the orbit index from point index (not fp.orbitID)
template <typename FP>
int FixedPointData<FP>::getOrbitIndex(const int ptIdx) const
{
    std::map<int,std::pair<int,int> >::const_iterator cit = indexMap.find( ptIdx );
    int cID = 0;
    if (cit != indexMap.end() ) {
        cID = cit->second.first;
        //theMsg->printf("%s: indexMap[%d] = [%d,%d]", __FILE__,ptIdx,cID,cit->second.second);
    } else {
        std::cerr << __FILE__ << ": Failed to lookup point id " << ptIdx << " in Orbit List! Returning 0\n";
    }
    
    return cID;
}

/// Get an orbit (chain) from set (if id>fp_chains.size(), returns 1st orbit)
template <typename FP>
const typename std::vector<FP>& FixedPointData<FP>::
getOrbit(const int orbitIdx) const
{
    int i = orbitIdx;
    if (orbitIdx<0 || orbitIdx >= (int) fp_chains.size()) {
        std::cerr << __FILE__ << ": Requested orbit index" << orbitIdx
                  << " which is larger than available number of fixed points ("
                  << (int) fp_chains.size() << ")\n";
        i=0;
    }
    return fp_chains[i];
}

/// Get an orbit (chain) from point index
template <typename FP>
void FixedPointData<FP>::getOrbitFromPointIndex(const int ptIdx, FixedPointChain& fp) const
{
    int cID = getOrbitIndex(ptIdx);
    fp = fp_chains[cID];
}

/// Get the maximum orbitID from data
template <typename FP>
int FixedPointData<FP>::getMaxOrbitID() const
{
    typename std::vector< FixedPointChain >::const_iterator it;
    int maxOID = -1e6;
    for(it=fp_chains.begin(); it!=fp_chains.end(); ++it) {
        int theID = (*it)[0].orbitID;
        if( theID > maxOID ) {
            maxOID = theID;
        }
    }
    return maxOID;
}

/// Get the minimum orbitID from data
template<typename FP>
int FixedPointData<FP>::getMinOrbitID() const
{
    typename FPChainsContainer::const_iterator it;
    int minOID = 1e6;
    for(it=fp_chains.begin(); it!=fp_chains.end(); ++it) {
        int theID = (*it)[0].orbitID;
        if( theID < minOID ) {
            minOID = theID;
        }
    }
    return minOID;
}

/// Get the orbit id classifier (fp.orbitID) from data
template<typename FP>
int FixedPointData<FP>::getOrbitID(const int orbitID, const int fpID) const
{
    assert( orbitID>=0 && orbitID < (int) fp_chains.size() );
    assert( fpID>=0 && fpID<(int) fp_chains[orbitID].size());
    //The data value (fp.orbitID) is a different identifier than the fp_chain position (orbitID)
    return fp_chains[orbitID][fpID].orbitID;
}


} // end pmate

#endif
