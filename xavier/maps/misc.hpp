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


#ifndef __MISC_HPP__
#define __MISC_HPP__

#include <vector>
#include <list>
#include <boost/rational.hpp>
#include <math/fixed_vector.hpp>

//Prototypes

namespace xavier {
template<typename T>
inline void push_front(const T& val, std::vector<T>& _in)
{
    if (!_in.size()) {
        _in.push_back(val);
    } else if (nvis::all(val != _in[0])) {
        std::vector<T> _out;
        _out.reserve(_in.size() + 1);
        _out.push_back(val);
        std::copy(_in.begin(), _in.end(), std::back_inserter(_out));
        std::swap(_in, _out);
    }
}

template<typename T>
inline double value(const boost::rational<T>& r)
{
    return (double)r.numerator() / (double)r.denominator();
}

template<typename Iterator>
inline double min_dist(double v, const Iterator& begin, const Iterator& end)
{
    std::list<double> dist;
    for (Iterator it = begin ; it != end ; ++it) {
        dist.push_back(fabs(*it - v));
    }
    return *std::min_element(dist.begin(), dist.end());
}

//Templated secant function
template<typename T>
double secant(const T& v0, const T& v1)
{
    const double threshold = 0.75*M_PI;
    const double tiny = 1.0e-3;
    
    T a = v0 / nvis::norm(v0);
    T b = v1 / nvis::norm(v1);
    if (acos(nvis::inner(a, b)) < threshold) {
        return 0.5;
    }
    if (nvis::norm(v0) < tiny || nvis::norm(v1) < tiny) {
        return 0.5;
    }
    
    T dir = v0 - v1;
    dir /= nvis::norm(dir);
    double f0 = nvis::inner(v0, dir);
    double f1 = nvis::inner(v1, dir);
    
    return -f0/(f1-f0);
}



} // xavier

#endif




