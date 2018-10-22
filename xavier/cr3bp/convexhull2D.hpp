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


#ifndef __CR3Bpoint_type_CONVEX_HULL_2D_Hpoint_typepoint_type__
#define __CR3Bpoint_type_CONVEX_HULL_2D_Hpoint_typepoint_type__

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>

namespace orbital {

template<typename P>
class convexhull {
public:
    typedef P                       point_type;
    typedef std::list<int>          ring_type;
    
    convexhull() : __points() {}
    
    void clear()
    {
        __points.clear();
    }
    
    template<typename Iterator>
    void add(const Iterator& begin, const Iterator& end)
    {
        std::copy(begin, end, std::back_inserter(__points));
    }
    
    void compute(ring_type& hull) const
    {
        int N = __points.size();
        std::vector<int> ids(N);
        for (int i=0 ; i<N ; ++i) {
            ids[i]=i;
        }
        
        std::sort(ids.begin(), ids.end(), lexico_more(__points));
        ring_type Lupper;
        grow(Lupper, ids);
        std::cerr << "Lupper: ";
        std::copy(Lupper.begin(), Lupper.end(), std::ostream_iterator<int>(std::cerr, ", "));
        std::cerr << '\n';
        
        std::sort(ids.begin(), ids.end(), lexico_less(__points));
        ring_type Llower;
        grow(Llower, ids);
        std::cerr << "Llower: ";
        std::copy(Llower.begin(), Llower.end(), std::ostream_iterator<int>(std::cerr, ", "));
        std::cerr << '\n';
        
        Llower.pop_back();
        Llower.pop_front();
        Lupper.insert(Lupper.end(), Llower.begin(), Llower.end());
        
        ring_type::iterator lowest = std::min_element(Lupper.begin(), Lupper.end());
        ring_type _back(Lupper.begin(), lowest);
        Lupper.erase(Lupper.begin(), lowest);
        Lupper.insert(Lupper.end(), _back.begin(), _back.end());
        
        hull.clear();
        std::copy(Lupper.begin(), Lupper.end(), std::back_inserter(hull));
    }
    
private:
    std::vector<point_type> __points;
    
    inline double whichSide(const point_type& p0, const point_type& p1, const point_type& p2) const
    {
        point_type v0 = p1-p0;
        point_type v1 = p2-p1;
        return v0[0]*v1[1] - v0[1]*v1[0];
    }
    
    struct lexico_less {
        lexico_less(const std::vector<point_type>& points) : __points(points) {}
        const std::vector<point_type>& __points;
        
        bool operator()(int i0, int i1)
        {
            const point_type& p0 = __points[i0];
            const point_type& p1 = __points[i1];
            return (p0[0] < p1[0] || (p0[0] == p1[0] && p0[1] < p1[1]));
        }
    };
    
    struct lexico_more {
        lexico_more(const std::vector<point_type>& points) : __points(points) {}
        const std::vector<point_type>& __points;
        
        bool operator()(int i0, int i1)
        {
            const point_type& p0 = __points[i0];
            const point_type& p1 = __points[i1];
            return (p0[0] > p1[0] || (p0[0] == p1[0] && p0[1] > p1[1]));
        }
    };
    
    void grow(ring_type& L, const std::vector<int>& ids) const
    {
        L.push_back(ids[0]);
        L.push_back(ids[1]);
        for (int i=2 ; i<ids.size() ; ++i) {
            L.push_back(ids[i]);
            while(L.size() >= 3) {
                std::list<int>::iterator i0, i1, i2;
                i0 = L.end();
                i0--;
                i2 = i0;
                i0--;
                i1 = i0;
                i0--;
                if (whichSide(__points[*i0], __points[*i1], __points[*i2]) <= 0) {
                    L.erase(i1);
                } else {
                    break;
                }
            }
        }
    }
};

}


#endif