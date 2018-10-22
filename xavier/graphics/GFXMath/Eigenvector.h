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


#ifndef EIGENVECTOR_H
#define EIGENVECTOR_H

#include "VectorN.h"

class Eigenvector {
public:
	Eigenvector() { isReal = false; }
	template <class T> Eigenvector(const T eigenvalue, const T *eigenvector, const unsigned int vectorDim) {
		set(eigenvalue, eigenvector, vectorDim);
	}

	template <class T> inline void set(const T eigenvalue, const T *eigenvector, const unsigned int vectorDim) {
		isReal = true;
		value = eigenvalue;

		dim = vectorDim;
		for (unsigned int i=0; i<vectorDim; i++) vec[i] = eigenvector[i];
	}

	template <class T> inline void setEigenvector(const T *eigenvector)
		{ vec[0] = eigenvector[0];  vec[1] = eigenvector[1];  vec[2] = eigenvector[2]; }
	template <class T> inline void setEigenvector(const T val0, const T val1, const T val2)
		{ vec[0] = val0;  vec[1] = val1;  vec[2] = val2; }

	bool isReal;
	double value;

	unsigned int dim;
	double vec[3];
};

#endif