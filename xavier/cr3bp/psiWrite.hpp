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


//PSI Write
#ifndef __PSI_WRITE_HPP
#define __PSI_WRITE_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct PointPSI {
    double x, y, z;
    PointPSI(double x=0, double y=0, double z=0) :
        x(x), y(y), z(z) {}
};

//Write PSI Data to file for Avizo
bool PSI_WriteToFile(vector<PointPSI>& points, vector<int>& ids,
                     vector< vector<double> >& data,
                     vector<string>& labels, vector<string>& symbols,
                     const char* filename)
{

    //Check sizes
    int numPts = points.size();
    int numIds = ids.size();
    int numDataPts = data[0].size();
    if (numPts != numIds || numPts != numDataPts) {
        return false;
    }
    int numData = data.size();
    int numSyms = symbols.size();
    int numLabs = labels.size();
    if (numData != numSyms || numData != numLabs) {
        return false;
    }
    
    
    //Open file
    std::ofstream theFile;
    theFile.open( filename );
    
    //Write Header:
    theFile << "# PSI Format 1.0\n#\n";
    //Labels
    theFile << "# column[0] = \"x\"\n";
    theFile << "# column[1] = \"y\"\n";
    theFile << "# column[2] = \"z\"\n";
    theFile << "# column[3] = \"Id\"\n";
    char buffer[200];
    int n=0;
    for (int k=0; k<numLabs; k++) {
        n = sprintf(buffer, "# column[%d] = \"%s\"\n",4+k,labels[k].c_str());
        theFile << buffer;
    }
    
    //Symbols
    theFile << "#\n";
    for (int k=0; k<numSyms; k++) {
        n = sprintf(buffer, "# symbol[%d] = \"%s\"\n",4+k,symbols[k].c_str());
        theFile << buffer;
    }
    
    //Type - assume all floating-point data for now
    theFile << "#\n";
    for (int k=0; k<numSyms; k++) {
        n = sprintf(buffer, "# type[%d] = float\n", 4+k);
        theFile << buffer;
    }
    
    //NumPts p q
    theFile << "\n";
    theFile << numPts << " " << numData << " " << numData << "\n";
    theFile << "1.00 0.00 0.00\n";
    theFile << "0.00 1.00 0.00\n";
    theFile << "0.00 0.00 1.00\n";
    theFile << "\n";
    
    //x y z id data per point
    for (int i=0; i<numPts; i++) {
        theFile << points[i].x << " " << points[i].y << " " << points[i].z << " " << ids[i];
        for (int k=0; k<numData; k++) {
            theFile << " " << data[k][i];
        }
        theFile << "\n";
    }
    //Close file
    theFile.close();
    
    return true;
};


#endif
