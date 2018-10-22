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


#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "graphics/Math/VectorN.h"

template <class T> class BoundingBox {
public:
	BoundingBox() { }

	virtual inline Vector3<T> getCenter() const = 0;
	virtual inline Vector3<T> getMin() const = 0;
	virtual inline Vector3<T> getMax() const = 0;
	virtual inline Vector3<T> getDimensions() const  = 0;
	virtual inline Vector3<T> getHalfDimensions() const  = 0;

	inline float getBoxIntersectionVolume(const Vector3f &boxMin, const Vector3f &boxMax) const {
		Vector3f center = getCenter();
		Vector3f bbHalfDim = getHalfDimensions();

		Vector3f bbMin = center-bbHalfDim;
		Vector3f bbMax = center+bbHalfDim;

		bbMin = bbMin.maxVector(boxMin);
		bbMax = bbMax.minVector(boxMax);

		return (bbMax.x-bbMin.x)*(bbMax.y-bbMin.y)*(bbMax.z-bbMin.z);
	}

	inline bool isPointInside(const Vector3f &pt) const { return pt.isInsideBox(getMin(), getMax()); }

};

template <class T> class AxisAlignedBoundingBox : public BoundingBox<T> {
public:
	AxisAlignedBoundingBox() { }

	template <class S> inline void set(const Vector3<S> &bbCenter, const Vector3<S> &bbDim) {
		center = bbCenter;
		halfDim = bbDim * 0.5;
	}
	
	inline Vector3<T> getCenter() const { return center; }
	inline Vector3<T> getMin() const { return center - halfDim; }
	inline Vector3<T> getMax() const { return center + halfDim; }
	inline Vector3<T> getDimensions() const { return 2.0 * halfDim; }
	inline Vector3<T> getHalfDimensions() const { return halfDim; }

	Vector3<T> center, halfDim;
};

template <class T> class AxisAlignedBoundingBox_MinSize : public BoundingBox<T> {
public:
	AxisAlignedBoundingBox_MinSize() { }

	template <class S> inline void set(const Vector3<S> &bbCenter, const Vector3<S> &bbDim) {
		center = bbCenter;
		halfDim = bbDim.maxValue()*0.5;
	}
	
	inline Vector3<T> getCenter() const { return center; }
	inline Vector3<T> getMin() const { return center - halfDim; }
	inline Vector3<T> getMax() const { return center + halfDim; }
	inline Vector3<T> getDimensions() const { float dim = 2.0 * halfDim;  return Vector3<T>(dim, dim, dim); }
	inline Vector3<T> getHalfDimensions() const { return Vector3<T>(halfDim, halfDim, halfDim); }

	Vector3f center;
	float halfDim;
};

//inline Vector3<T> getCenter() const { }
//inline Vector3<T> getMin() const { }
//inline Vector3<T> getMax() const { }
//inline Vector3<T> getDimensions() const { }
//inline Vector3<T> getHalfDimensions() const { }

#endif