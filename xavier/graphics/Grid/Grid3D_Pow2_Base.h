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


#ifndef GRID_3D_POW_2_BASE_H
#define GRID_3D_POW_2_BASE_H

#include "Grid3D_Regular_Base.h"
#include "Math/MathUtils.h"

template <class T> class Grid3D_Pow2_Base : public Grid3D_Regular_Base<T> {
public:
	virtual void rebuildGrid(Vector3ui numberOfCells, T outsideValue) {
		initialNumCells = numberOfCells;

		Vector3ui pow2NumCells(MathUtils::roundUpPowerOf2(numberOfCells.x), MathUtils::roundUpPowerOf2(numberOfCells.y), MathUtils::roundUpPowerOf2(numberOfCells.z));

		Grid3D_Regular_Base::rebuildGrid(pow2NumCells, outsideValue);
	}

	virtual inline int getType() = 0;

	Vector3ui initialNumCells;

	inline unsigned int getCellIndex(const unsigned int x, const unsigned int y, const unsigned int z) const { return MathUtils::interleaveBits(x,y,z); }

	inline unsigned int getCellIndexComponentX(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentX(index); }
	inline unsigned int getCellIndexComponentY(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentY(index); }
	inline unsigned int getCellIndexComponentZ(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentZ(index); }
	inline Vector3ui getCellIndices(const unsigned int index) const { unsigned int x,y,z;  MathUtils::unInterleaveBits(index,x,y,z);  return Vector3ui(x,y,z); }

protected:

};

#endif