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


//ManifoldConnection class for autonomous connection generation
//Author:  Wayne Schlei (Purdue University)
//Date:  7/2/2016

#ifndef MANIFOLD_CONNECTION_HPP
#define MANIFOLD_CONNECTION_HPP

#include <vector>
#include <fstream>
#include <iostream>


namespace pmate {

///Manifold Connections templated with deltaV vector type
template <class VEC>
class ManifoldConnection {
public:
    typedef ManifoldConnection<VEC>  SelfType;
    typedef std::pair<int,int>       IntPair;
    ///Constructor
    ManifoldConnection() :
        manifoldDataIdx(IntPair(0,0)),
        uMan(IntPair(0,0)), uTau(0.0),
        sMan(IntPair(0,0)), sTau(0.0),
        isHomoclinic(false), pos(0.0), tof(-1.0),
        depDV(0.0), arrDV(0.0), deltaV(0.0)
    {}
    
    /// ManifoldData indicator (assuming ManifoldData objects are indexed)
    IntPair manifoldDataIdx;
    /// Unstable Manifold [manID,segID] identifier
    IntPair uMan;
    /// Parameter on Unstable Manifold that locates connection point
    double uTau;
    /// Stable Manifold [manID,segID] identifier
    IntPair sMan;
    /// Parameter on Stable Manifold that locates connection point
    double sTau;
    /// Is this a homoclinic connection
    bool isHomoclinic;
    /// Connection point at detected intersection
    VEC pos; //Note: There are likely more map locations for this particular connection.
    
    /// Time of flight from unstable progenitor state to stable progenitor state
    double tof;
    /// Departure Delta-V in position space
    VEC depDV;
    /// Arrival Delta-V in position space
    VEC arrDV;
    /// Augmented Delta-V (at connection point) in position space [HARD CODED TO PLANAR (X,Y)]
    VEC deltaV;
    
    
    ///Write connection data to a given file (as one line)
    void write(FILE* f) const
    {
        //HARD-CODED to 2D delta-V vectors
        fprintf(f,"%d %d %d %d %.15f %d %d %.15f %d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
                manifoldDataIdx.first, manifoldDataIdx.second,
                uMan.first, uMan.second, uTau,
                sMan.first, sMan.second, sTau,
                (isHomoclinic)?1:0, tof, pos[0], pos[1],
                depDV[0], depDV[1], arrDV[0], arrDV[1], deltaV[0], deltaV[1]
               );
    }
    ///Read connection data from a single line in a given file
    void read(FILE* f)
    {
        int ih = 0;
        //HARD-CODED to 2D delta-V vectors
        fscanf(f,"%d %d %d %d %lf %d %d %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &(manifoldDataIdx.first), &(manifoldDataIdx.second),
               &(uMan.first), &(uMan.second), &uTau,
               &(sMan.first), &(sMan.second), &sTau,
               &ih, &tof, &(pos[0]), &(pos[1]),
               &(depDV[0]), &(depDV[1]), &(arrDV[0]), &(arrDV[1]), &(deltaV[0]), &(deltaV[1])
              );
        isHomoclinic = (ih == 1);
    }
    
};

} //end pmate

#endif