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


#ifndef GRID_3D_REGULAR_SPARSE_I2_H
#define GRID_3D_REGULAR_SPARSE_I2_H

#include "Grid3D_Regular_Base.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Sparse_I2 : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Sparse_I2() { }

	Grid3D_Regular_Sparse_I2(Vector3ui numberOfCells, T outsideValue) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Sparse_I2(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid3D_Regular_Sparse_I2() { clear(false); }

	void clear(bool freePointers) {
		if (data.size() == 0) return;

		/*if (freePointers) {
			for (unsigned int i=0; i<data.size(); i++)
			for (unsigned int j=0; j<data[i].size(); j++)
			for (unsigned int k=0; k<data[i][j].size(); k++)
				if (data[i][j][k] != NULL) { delete data[i][j][k];  data[i][j][k] = NULL; }
		}*/

		for (unsigned int i=0; i<data.size(); i++) {
			for (unsigned int j=0; j<data[i].size(); j++) {
				data[i][j].clear();
				depthIndices[i][j].clear();
			}
			
			data[i].clear();
			colIndices[i].clear();
			depthIndices[i].clear();
		}
		data.clear();
		colIndices.clear();
		depthIndices.clear();
	}

	inline int getType() { return GRID_REGULAR_SPARSE_I2; }

	void rebuildGrid(Vector3ui numberOfCells, T outsideValue) {
		/*Grid3D_Regular_Base::rebuildGrid(numberOfCells, outsideValue);

		if (data.size() > 0) { clear(false); }*/

		Grid3D_Regular_Sparse_Base::rebuildGrid(numberOfCells, outsideValue);

		colIndices.resize(numberOfCells.x);
		depthIndices.resize(numberOfCells.x);
		data.resize(numberOfCells.x);
	}

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int col;
		if (colIndices[x].find(y, &col)) {
			unsigned int depth;
			if (depthIndices[x][col].find(z, &depth)) data[x][col][depth] = val; 
			else {
				depthIndices[x][col].insert(depth, z);
				data[x][col].insert(depth, val);
			}
		}
		else {
			colIndices[x].insert(col, y);

			depthIndices[x].insert(col, SortedStlVector<unsigned int>());
			depthIndices[x][col].add(z);

			data[x].insert(col, StlVector<T>());
			data[x][col].add(val);
		}
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int col;
		if (colIndices[x].find(y, &col)) {
			unsigned int depth;
			if (depthIndices[x][col].find(z, &depth)) return data[x][col][depth]; 
		}
		return outsideVal;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		unsigned int col;
		if (colIndices[x].find(y, &col)) {
			unsigned int depth;
			if (depthIndices[x][col].find(z, &depth)) return &data[x][col][depth]; 
		}
		return NULL;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int col;
		if (colIndices[x].find(y, &col)) {
			unsigned int depth;
			if (depthIndices[x][col].find(z, &depth)) return true; 
		}
		return false;
	}

protected:
	StlVector< SortedStlVector<unsigned int> >colIndices;
	StlVector< StlVector< SortedStlVector<unsigned int> > >depthIndices;
	StlVector< StlVector< StlVector<T> > > data;
};

#endif