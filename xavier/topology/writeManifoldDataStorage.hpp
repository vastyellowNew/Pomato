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


///////////////////////////////////////////////////////////////////////
//
// Write routine ManifoldDataStorage
//
///////////////////////////////////////////////////////////////////////
#ifndef __WRITE_MANIFOLD_DATA_STORAGE_HPP
#define __WRITE_MANIFOLD_DATA_STORAGE_HPP

#include <cstdio>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <math/fixed_vector.hpp>
#include <maps/fixpoints.hpp>
#include <topology/ManifoldClasses.hpp>
#include <topology/invariant_manifold.hpp>
#include <topology/ManifoldDataStorage.hpp>

using namespace xavier;

namespace topology {

bool writeManifoldDataStorage(ManifoldDataStorage* data, const char* filename)
{
#if defined(_WIN32) || defined(C_0X)
    typedef Separatrix::MapDiscont MapDiscont;
#else
    typedef typename Separatrix::MapDiscont MapDiscont;
#endif
    FILE* f = fopen(filename, "w"); // open the file
    
    if (!f) {
        std::cerr << " Bad filename ...\n";
        return false; // indicate error
    }
    
    // Write header:
    fprintf(f, "# Manifold Data\n");
    
    ///-------------------------------------------------------------------------------
    fprintf(f,"# FP Chains\n");
    // Write chain information
    int numFPs = (int) data->chains.size();
    int spaceDim = 2; //Fixed for now
    fprintf(f, "%d %d\n", spaceDim, numFPs);
    
    // Write fixed point data values:
    // period fixpoint0 fixpoint1 ... fixpointP
    std::vector< fp_chain >::iterator it;
    for (it=data->chains.begin(); it!=data->chains.end(); it++) {
        //9 Commons per fixed point
        int period = (*it)[0].K;
        int orbitID = (*it)[0].orbitID;
        double timePeriod = (*it)[0].timePeriod;
        double si = (*it)[0].si, otherSI = (*it)[0].otherSI;
        double evals[2];
        evals[0] = (*it)[0].eval[0];
        evals[1] = (*it)[0].eval[1];
        int saddle = ((*it)[0].saddle) ? 1 : 0;
        int isolated = ((*it)[0].isolated) ? 1 : 0;
        //Write commons to first line
        fprintf(f,"%d %d %d %d %.15f %.15f %.15f %.15f %.15f\n",
                period,saddle,orbitID,isolated,timePeriod,evals[0],evals[1],si,otherSI);
                
        //Point information - One line per indicated point in chain
        for (int i=0; i<period; i++) {
            nvis::vec2 pos = (*it)[i].pos;
            //fixpoint print: t (x xdot) (evec0) (evec1)
            fprintf(f,"%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
                    (*it)[i].t, pos[0], pos[1],
                    (*it)[i].evec[0][0], (*it)[i].evec[0][1],
                    (*it)[i].evec[1][0], (*it)[i].evec[1][1]);
        }
    }
    
    ///--------------------------------------------------------------------------------
    //Separatrices
    //fprintf(f, "# Separatrix Data\n");
    int numSepx = (int) data->separatrices.size();
    fprintf(f,"%d\n", numSepx);
    std::vector<Separatrix>::iterator mit;
    for(mit=data->separatrices.begin(); mit!=data->separatrices.end(); ++mit) {
        //Common values
        int start0 = mit->start.first;
        int start1 = mit->start.second;
        int end0 = mit->end.first;
        int end1 = mit->end.second;
        int lastSeed = mit->lastSeed;
        int forward = (mit->forward)? 1 : -1; //Unstable = 1, stable = -1
        int numManPts = (int) mit->manifold.size();
        int numBreaks = (int) mit->breakIDs.size();
        int numSeps = (int) mit->discontList.size();
        double length = mit->length;
        fprintf(f,"%d %d %d %d %d %d %d %d %d %.15f\n",
                start0,start1,end0,end1,lastSeed,forward,numManPts,numBreaks,numSeps,length);
                
        //Manifold Points - Maybe modify later to hold more info?
        for(int j=0; j<numManPts; j++) {
            nvis::vec2& pos = mit->manifold[j];
            fprintf(f,"%.15f %.15f\n",pos[0],pos[1]);
        }
        
        //BreakIDs - write in one line
        std::set<int>::iterator bit = mit->breakIDs.begin();
        for(; bit!=mit->breakIDs.end(); ++bit) {
            fprintf(f,"%d ",*bit);
        }
        fprintf(f,"\n"); //end line
        
        //Section Discontinuities
        std::list<MapDiscont>::iterator lit = mit->discontList.begin();
        for(; lit!=mit->discontList.end(); ++lit) {
            unsigned int t = lit->getTypeNumber(lit->type);
            int p = lit->period;
            nvis::vec2 pos = lit->failurePos;
            fprintf(f,"%d %d %.15f %.15f\n",t,p,pos[0],pos[1]);
        }
    }
    
    ///--------------------------------------------------------------------------------
    
    
    // Close the file.
    fclose(f);
    
    return true; // indicate success
}

} //end topology

#endif
