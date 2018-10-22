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


#ifndef GRID_2D_REGULAR_SPARSE_I1_H
#define GRID_2D_REGULAR_SPARSE_I1_H

#include "Grid2D_Regular_Base.h"

#include "StlVectorUtils.h"

template <class T> class Grid2D_Regular_Sparse_I1 : public Grid2D_Regular_Base<T> {
public:
	Grid2D_Regular_Sparse_I1() { }

	Grid2D_Regular_Sparse_I1(Vector2ui numberOfCells, T outsideValue) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(Vector2d(0,0), Vector2d(1,1));
	}

	Grid2D_Regular_Sparse_I1(Vector2ui numberOfCells, T outsideValue, Vector2d startPosition, Vector2d cellSizes) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid2D_Regular_Sparse_I1() { clear(); }

	void clear() {
		if (data.size() == 0) return;

		for (unsigned int i=0; i<data.size(); i++) {
			data[i].clear();
			columnIndices[i].clear();
		}
		data.clear();
		columnIndices.clear();
	}

	inline void deletePointers() {
		for (unsigned int i=0; i<data.size(); i++) data[i].deleteElements();
	}

	inline int getType() { return GRID2D_REGULAR_SPARSE_I1; }

	void rebuildGrid(Vector2ui numberOfCells, T outsideValue) {
		Grid2D_Regular_Base::rebuildGrid(numberOfCells, outsideValue);

		columnIndices.resize(numberOfCells.x);
		data.resize(numberOfCells.x);

		numberOfEntries = 0;
		numberOfDiagnolsSet = 0;
	}

	inline void set(const unsigned int x, const unsigned int y, const T &val) {
		unsigned int column;
		if (columnIndices[x].find(y, &column)) data[x].vec[column] = val;
		else {
			columnIndices[x].insert(column, y);
			data[x].insert(column, val);
			numberOfEntries++;

			if (x == y) numberOfDiagnolsSet++;
		}
	}

	inline T set_returnOld(const unsigned int x, const unsigned int y, const T &val) {
		unsigned int column;
		if (columnIndices[x].find(y, &column)) { T old = data[x][column];  data[x].vec[column] = val;  return old; }
		else {
			columnIndices[x].insert(column, y);
			data[x].insert(column, val);
			numberOfEntries++;

			if (x == y) numberOfDiagnolsSet++;

			return outsideVal;
		}
	}

	inline T get(const unsigned int x, const unsigned int y) const {
		unsigned int column;
		if (columnIndices[x].find(y, &column)) return data[x][column];
		else return outsideVal;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y) {
		unsigned int column;
		if (columnIndices[x].find(y, &column)) return &data[x].vec[column];
		else return NULL;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y) const {
		unsigned int column;
		if (columnIndices[x].find(y, &column)) return true;
		else return false;
	}

	StlVector< SortedStlVector<unsigned int> >columnIndices;
	StlVector< StlVector<T> > data;

	unsigned int numberOfEntries;
	unsigned int numberOfDiagnolsSet;

	inline double getAverageNumberOfBytesPerElement() const;

protected:

};

template <class T> inline double Grid2D_Regular_Sparse_I1<T>::getAverageNumberOfBytesPerElement() const { return sizeof(T); }

// I don't know how to do get the template pattern matching to work nicely with this stuff (without making a specialized template class)
//   so instead, when a new type is added, just add a specialization here
template <> inline double Grid2D_Regular_Sparse_I1< vector<int> >::getAverageNumberOfBytesPerElement() const {
	double size=0;
	for (unsigned int i=0; i<data.size(); i++)
	for (unsigned int j=0; j<data[i].size(); j++) { size += (double)data[i][j].size()*sizeof(int); }
	return size/(double)data.size();
}
template <> inline double Grid2D_Regular_Sparse_I1< vector<unsigned int> >::getAverageNumberOfBytesPerElement() const {
	double size=0;
	for (unsigned int i=0; i<data.size(); i++)
	for (unsigned int j=0; j<data[i].size(); j++) { size += (double)data[i][j].size()*sizeof(unsigned int); }
	return size/(double)data.size();
}

#endif