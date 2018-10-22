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


#ifndef GRID_3D_REGULAR_SPARSE_I1_H
#define GRID_3D_REGULAR_SPARSE_I1_H

#include "Grid3D_Regular_Base.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Sparse_I1 : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Sparse_I1() { }

	Grid3D_Regular_Sparse_I1(Vector3ui numberOfCells, T outsideValue) {
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Sparse_I1(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid3D_Regular_Sparse_I1() { clear(); }

	void clear() {
		if (data.size() == 0) return;

		for (unsigned int i=0; i<data.size(); i++) {
			for (unsigned int j=0; j<data[i].size(); j++) {
				data[i][j].clear();
				depthIndices[i][j].clear();
			}
			
			data[i].clear();
			depthIndices[i].clear();
		}
		data.clear();
		depthIndices.clear();
	}

	inline void deletePointers() {
		for (unsigned int i=0; i<data.size(); i++)
		for (unsigned int j=0; j<data[i].size(); j++) data[i][j].deleteElements();
	}

	inline int getType() const { return GRID_REGULAR_SPARSE_I1; }

	void rebuildGrid(Vector3ui numberOfCells) {
		Grid3D_Regular_Base::rebuildGrid(numberOfCells);

		depthIndices.resize(numberOfCells.x);
		data.resize(numberOfCells.x);

		for (unsigned int i=0; i<depthIndices.size(); i++) {
			depthIndices[i].resize(numberOfCells.y);
			data[i].resize(numberOfCells.y);
		}

		numberOfEntries = 0;
		numberOfDiagnolsSet = 0;
	}

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int depth;
		if (depthIndices[x][y].find(z, &depth)) data[x][y].vec[depth] = val;
		else {
			depthIndices[x][y].insert(depth, z);
			data[x][y].insert(depth, val);
			numberOfEntries++;

			if (x == y && y == z) numberOfDiagnolsSet++;
		}
	}

	inline T set_returnOld(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int depth;
		if (depthIndices[x][y].find(z, &depth)) { T old = data[x][y][depth];  data[x][y].vec[depth] = val;  return old; }
		else {
			depthIndices[x][y].insert(depth, z);
			data[x][y].insert(depth, val);
			numberOfEntries++;

			if (x == y && y == z) numberOfDiagnolsSet++;

			return outsideVal;
		}
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int depth;
		if (depthIndices[x][y].find(z, &depth)) return data[x][y][depth];
		else return outsideVal;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		unsigned int depth;
		if (depthIndices[x][y].find(z, &depth)) return &data[x][y].vec[depth];
		else return NULL;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int depth;
		if (depthIndices[x][y].find(z, &depth)) return true;
		else return false;
	}

	// to do: it seems that if the outside vector type is StlVector, then things like resize stop working
	vector< vector< SortedStlVector<unsigned int> > >depthIndices;
	vector< vector< StlVector<T> > > data;

	unsigned int numberOfEntries;
	unsigned int numberOfDiagnolsSet;

	inline double getAverageNumberOfBytesPerElement() const;

protected:

};

template <class T> inline double Grid3D_Regular_Sparse_I1<T>::getAverageNumberOfBytesPerElement() const { return sizeof(T); }

// I don't know how to do get the template pattern matching to work nicely with this stuff (without making a specialized template class)
//   so instead, when a new type is added, just add a specialization here
template <> inline double Grid3D_Regular_Sparse_I1< vector<int> >::getAverageNumberOfBytesPerElement() const {
	double size=0;
	for (unsigned int i=0; i<data.size(); i++)
	for (unsigned int j=0; j<data[i].size(); j++)
	for (unsigned int k=0; k<data[i][j].size(); k++) { size += (double)data[i][j][k].size()*sizeof(int); }
	return size/(double)data.size();
}
template <> inline double Grid3D_Regular_Sparse_I1< vector<unsigned int> >::getAverageNumberOfBytesPerElement() const {
	double size=0;
	for (unsigned int i=0; i<data.size(); i++)
	for (unsigned int j=0; j<data[i].size(); j++)
	for (unsigned int k=0; k<data[i][j].size(); k++) { size += (double)data[i][j][k].size()*sizeof(unsigned int); }
	return size/(double)data.size();
}

#endif