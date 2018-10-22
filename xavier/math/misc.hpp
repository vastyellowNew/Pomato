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


#ifndef __XAVIER_MATH_MISC_HPP__
#define __XAVIER_MATH_MISC_HPP__

template<typename T>
inline T sign(T a)
{
    return (a < 0) ? -1 : 1;
}

template<typename T>
inline T mindist(T a, T b)
{
    if (a > 0.5*b) {
        return a - b;
    } else if (a < -0.5*b) {
        return a + b;
    } else {
        return a;
    }
}

template<typename I>
inline I ipow(I n, I k)
{
    I result = 1;
    while (k) {
        if (k & 1) {
            result *= n;
        }
        k >>= 1;
        n *= n;
    }
    
    return result;
}

#endif