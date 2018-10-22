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


#ifndef __INTERVAL_HPP__
#define __INTERVAL_HPP__

#include <iostream>

namespace xavier {
template<typename T>
struct interval {
    interval() : __min(T(0)), __max(T(0)) {}
    interval(T min, T max) : __min(min), __max(max) {}
    interval(const interval<T>& _int) : __min(_int.__min), __max(_int.__max) {}
    
    bool inside(T val) const
    {
        return (val >= __min && val <= __max);
    }
    
    bool empty() const
    {
        return (__max <= __min);
    }
    
    T length() const
    {
        return std::max(0., __max - __min);
    }
    
    T __min, __max;
};
}


template<typename T>
inline xavier::interval<T> intersect(const xavier::interval<T>& i0, const xavier::interval<T>& i1)
{
    return xavier::interval<T>(std::max(i0.__min, i1.__min),
                               std::min(i0.__max, i1.__max));
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const xavier::interval<T>& i)
{
    os << "[ " << i.__min << ", " << i.__max << " ]";
    return os;
}


#endif


