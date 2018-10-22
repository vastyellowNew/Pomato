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


#ifndef GRID_3D_REGULAR_SPARSE_VEC_H
#define GRID_3D_REGULAR_SPARSE_VEC_H

#include "Grid3D_Regular_Base.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Sparse_Vec : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Sparse_Vec() { }

	Grid3D_Regular_Sparse_Vec(Vector3ui numberOfCells, T outsideValue) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Sparse_Vec(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid3D_Regular_Sparse_Vec() { clear(false); }

	void clear(bool freePointers) {
		if (data.size() == 0) return;

		/*if (freePointers) {
			for (unsigned int i=0; i<data.size(); i++)
				if (data[i] != NULL) { delete data[i];  data[i] = NULL; }
		}*/

		data.clear();
		indices.clear();
	}

	inline int getType() { return GRID_REGULAR_SPARSE_VEC; }

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int index = getCellIndex(x,y,z), dataIndex;
		if (indices.find(index, &dataIndex)) data[dataIndex] = val;
		else {
			indices.insert(dataIndex, index);
			data.insert(dataIndex, val);
		}
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int dataIndex;
		if (indices.find(getCellIndex(x,y,z), &dataIndex)) return data[dataIndex];
		return outsideVal;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		unsigned int dataIndex;
		if (indices.find(getCellIndex(x,y,z), &dataIndex)) return &data[dataIndex];
		return NULL;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		return indices.find(getCellIndex(x,y,z));
	}

protected:
	SortedStlVector<unsigned int> indices;
	StlVector<T> data;
};

#endif