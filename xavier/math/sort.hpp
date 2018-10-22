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


#ifndef __XAVIER_SORT_HPP__
#define __XAVIER_SORT_HPP__

#include <iostream>
#include <complex>
#include <vector>
#include <map>
#include <algorithm>
#include <math/fixed_vector.hpp>

namespace xavier {
template< typename T >
bool LtIndexed(const T& v1, const T& v2)
{
    return (v1.first < v2.first);
}

template< typename T >
bool GtIndexed(const T& v1, const T& v2)
{
    return (v1.first > v2.first);
}

template< typename T >
inline void sort_ids(std::vector< unsigned int >& sorted,
                     const std::vector< T >& values, bool increasing = true)
{
    typedef std::pair< T, unsigned int > indexed_type;
    
    std::vector< indexed_type > tmp(values.size());
    for (unsigned int i = 0 ; i < values.size() ; i++) {
        tmp[i] = indexed_type(values[i], i);
    }
    
    if (increasing) {
        std::sort(tmp.begin(), tmp.end(), LtIndexed< indexed_type >);
    } else {
        std::sort(tmp.begin(), tmp.end(), GtIndexed< indexed_type >);
    }
    sorted.resize(tmp.size());
    for (unsigned int i = 0 ; i < sorted.size() ; i++) {
        sorted[i] = tmp[i].second;
    }
}

template< typename T >
inline unsigned int min_id(const std::vector< T >& values)
{
    return std::distance(values.begin(), std::min_element(values.begin(), values.end()));
}

template< typename T >
inline unsigned int max_id(const std::vector< T >& values)
{
    return std::distance(values.begin(), std::max_element(values.begin(), values.end()));
}
template<typename T1, typename T2>
struct Lt_pair_lexicographic {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1)
    {
        return (p0.first < p1.first) ||
               (p0.first == p1.first && p0.second < p1.second);
    }
};

template<typename T1, typename T2>
struct Lt_pair_first {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1)
    {
        return (p0.first < p1.first);
    }
};

template<typename T1, typename T2>
struct Lt_pair_second {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1)
    {
        return (p0.second < p1.second);
    }
};

template<typename V, int N>
struct Lt_vector_N {
    bool operator()(const V& v0, const V& v1)
    {
        return (v0[N] < v1[N]);
    }
};

struct Lt_fabs {
    bool operator()(double a, double b)
    {
        return fabs(a) < fabs(b);
    }
};

// A simple fixed length container that keeps its elements in the order
// defined by std::less<T>
template<typename T, int N, class Comp = std::less<T> >
class fixed_sorted_vector {
public:
    typedef nvis::fixed_vector<T, N>    vec_type;
    typedef fixed_sorted_vector<T, N>   self_type;
    typedef Comp                        less_than;
    
    fixed_sorted_vector()
        : _data() {}
        
    fixed_sorted_vector(const self_type& sa)
        : _data(sa._data) {}
        
    fixed_sorted_vector(const vec_type& v)
        : _data(v)
    {
        std::sort(_data.begin(), _data.end(), less_than());
    }
    
    fixed_sorted_vector(const T& v)
        : _data(v) {}
        
    fixed_sorted_vector(const T& v0, const T& v1)
        : _data(v0, v1)
    {
        std::sort(_data.begin(), _data.end(), less_than());
    }
    
    fixed_sorted_vector(const T& v0, const T& v1, const T& v2)
        : _data(v0, v1, v2)
    {
        std::sort(_data.begin(), _data.end(), less_than());
    }
    
    fixed_sorted_vector(const T& v0, const T& v1, const T& v2, const T& v3)
        : _data(v0, v1, v2, v3)
    {
        std::sort(_data.begin(), _data.end(), less_than());
    }
    
    fixed_sorted_vector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4)
        : _data(v0, v1, v2, v3, v4)
    {
        std::sort(_data.begin(), _data.end(), less_than());
    }
    
    template<typename V>
    fixed_sorted_vector(const V& v)
    {
        for (int i=0 ; i<N ; ++i) {
            _data[i] = v[i];
        }
        std::sort(_data.begin(), _data.end(), less_than());
    }
    
    bool operator==(const self_type& sa) const
    {
        return nvis::all(_data == sa._data);
    }
    
    bool operator<(const self_type& sa) const
    {
        nvis::lexicographical_order lex;
        return lex(_data, sa._data);
    }
    
    const T& operator[](size_t i) const
    {
        return _data[i];
    }
    
    const T& min() const
    {
        return _data[0];
    }
    
    const T& max() const
    {
        return _data[N-1];
    }
    
private:
    vec_type _data;
};

template< typename T >
void sort(const std::vector< T >& values,
          std::vector< unsigned int >& sorted);
          
template< typename T, typename Order >
void sort(const std::vector< T >& values,
          std::vector< unsigned int >& sorted,
          const Order& order);
          
template< typename T >
bool test(const std::vector< T >& values,
          const std::vector< unsigned int >& sorted);
          
template< typename T >
inline void sort(const std::vector< T >& values,
                 std::vector< unsigned int >& sorted)
{
    sorted.clear();
    typedef std::list< unsigned int >                   list_type;
    list_type::const_iterator                           lit;
    std::map< T, list_type >                            _map;
    typename std::map< T, list_type >::iterator         mit;
    typename std::map< T, list_type >::const_iterator   cmit;
    for (unsigned int i = 0 ; i < values.size() ; i++) {
        mit = _map.find(values[i]);
        if (mit == _map.end()) {
            _map[values[i]].push_back(i);
        } else {
            mit->second.push_back(i);
        }
    }
    
    for (cmit = _map.begin() ; cmit != _map.end() ; ++cmit) {
        for (lit = cmit->second.begin() ; lit != cmit->second.end(); ++lit) {
            sorted.push_back(*lit);
        }
    }
}

template< typename T, typename Order >
inline void sort(const std::vector< T >& values,
                 std::vector< unsigned int >& sorted,
                 const Order& order)
{
    sorted.clear();
    typedef std::list< unsigned int >                           list_type;
    list_type::const_iterator                                   lit;
    std::map< T, list_type, Order >                             _map;
    typename std::map< T, list_type, Order >::iterator          mit;
    typename std::map< T, list_type, Order >::const_iterator    cmit;
    for (unsigned int i = 0 ; i < values.size() ; i++) {
        mit = _map.find(values[i]);
        if (mit == _map.end()) {
            _map[values[i]].push_back(i);
        } else {
            mit->second.push_back(i);
        }
    }
    
    for (cmit = _map.begin() ; cmit != _map.end() ; ++cmit) {
        for (lit = cmit->second.begin() ; lit != cmit->second.end(); ++lit) {
            sorted.push_back(*lit);
        }
    }
}

template< typename T >
bool test(const std::vector< T >& values,
          const std::vector< unsigned int >& sorted)
{
    if (!sorted.size() && !values.size()) {
        return true;
    }
    
    T old, cur;
    old = values[sorted[0]];
    for (unsigned int i = 1 ; i < sorted.size() ; i++) {
        cur = values[sorted[i]];
        if (old > cur) {
            return false;
        }
        
        old = cur;
    }
    
    return true;
}

}


#endif








