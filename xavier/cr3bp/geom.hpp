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


#ifndef __CR3BP_GEOM_HPP__
#define __CR3BP_GEOM_HPP__

#include "convexhull2D.hpp"
#include <vector>
#include <list>

namespace orbital {

template<typename P>
static P innermost(const std::vector<P>& points)
{
    typedef P                      point_type;
    typedef std::list<int>         ring_type;
    typedef std::list<point_type>  polygon_type;
    
    convexhull<point_type> chull;
    chull.add(points.begin(), points.end());
    
    ring_type hull;
    chull.compute(hull);
    
    std::cout << "Convex hull: ";
    std::copy(hull.begin(), hull.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << '\n';
    std::vector<bool> removed(points.size(), false);
    
    int nbremoved=0;
    for (ring_type::const_iterator it=hull.begin() ; it!=hull.end() ; ++it, ++nbremoved) {
        removed[*it] = true;
    }
    std::cout << "nb removed = " << nbremoved << '\n';
    
    std::vector<point_type> inner;
    while (nbremoved < points.size()-2) {
        for (int i=0 ; i<points.size() ; ++i) {
            if (!removed[i]) {
                inner.push_back(points[i]);
            }
        }
        
        chull.clear();
        chull.add(inner.begin(), inner.end());
        hull.clear();
        chull.compute(hull);
        std::cout << "Convex hull: ";
        std::copy(hull.begin(), hull.end(), std::ostream_iterator<int>(std::cout, ", "));
        std::cout << '\n';
        
        for (ring_type::const_iterator it=hull.begin() ; it!=hull.end() ; ++it, ++nbremoved) {
            removed[*it] = true;
        }
        std::cout << "nb removed = " << nbremoved << '\n';
    }
    
    // hull contains the vertices of the innermost convex layer
    point_type innermost(0);
    for (ring_type::const_iterator it=hull.begin() ; it!=hull.end() ; ++it) {
        innermost += inner[*it];
    }
    
    return innermost / hull.size();
}

}

#endif