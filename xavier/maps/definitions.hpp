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


#ifndef __XAVIER_MAPS_DEFINITIONS_HPP__
#define __XAVIER_MAPS_DEFINITIONS_HPP__

#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <boost/rational.hpp>
//#include <boost/regex.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <data/edge.hpp>
#include <data/grid.hpp>
#include <data/raster_data.hpp>
#include <math/angle.hpp>
#include <math/sort.hpp>
#include "metric.hpp"
#include "map_analysis_param.hpp"

namespace xavier {

/// Storing initial sampling data.  "steps" are all the map iterates from the propagation
struct orbit_data {
    orbit_data() : steps(), wn(), periods(), isICValid(true) {}
    typedef std::vector<double>    value_type;
    
    std::vector<nvis::vec2> steps;
    std::vector<double> wn; //winding number set (at end of run)
    std::vector<int> periods;
    /// Check if there are existing iterates
    bool hasSteps() const
    {
        return (bool) steps.size();
    }
    /// Valid if within viable region of motion (inside feasible region or no imaginary velocity)
    bool isICValid;
};

/// External dataset_type
typedef raster_data<orbit_data, double, 2>        dataset_type;

// data type to associate a period to an object type that supports sorting
template<typename T, typename Compare = std::less<T> >
struct tagged : public std::pair<T, int> {
    typedef std::pair<T, int>   pair_type;
    typedef tagged<T, Compare>    self_type;
    
    tagged() : pair_type(T(), 0) {}
    tagged(T _t, int _p) : pair_type(_t, _p) {}
    int& tag()
    {
        return this->second;
    }
    int tag() const
    {
        return this->second;
    }
    T& object()
    {
        return this->first;
    }
    const T& object() const
    {
        return this->first;
    }
    bool operator<(const self_type& other) const
    {
        if (this->second < other.second) {
            return true;
        } else if (this->second > other.second) {
            return false;
        } else {
            Compare c;
            return c(this->first, other.first);
        }
    }
};
typedef edge<nvis::ivec2, nvis::lexicographical_order>            edge_type;
///Associate a sortable period with each edge
typedef tagged<edge_type>                     p_edge_type;
///Associate a sortable period with each cell
typedef tagged<nvis::ivec2, nvis::lexicographical_order>    p_cell_type;


struct rational_surface_suspected : public std::runtime_error {
    explicit rational_surface_suspected(const std::string& __arg)
        : std::runtime_error(__arg) {}
};

struct higher_precision_required : public std::runtime_error {
    explicit higher_precision_required(const std::string& __arg)
        : std::runtime_error(__arg) {}
};

struct rational_surface_found : public std::runtime_error {
    typedef std::pair<ivec_type, ivec_type> edge_type;
    typedef std::pair<edge_type, ivec_type>    edge_point;
    
    explicit rational_surface_found(const std::string& __arg)
        : std::runtime_error(__arg) {}
        
    edge_point xings[2];
    rational_type wn;
};

struct rational_segment {
    typedef std::pair<edge_type, ivec_type>    edge_point;
    edge_point pt[2];
    xavier::rational_type sf;
};
// Less Than operator for 2D bounding boxes
struct Lt_bbox {
    nvis::lexicographical_order Lt;
    
    bool operator()(const nvis::bbox2& b0, const nvis::bbox2& b1) const
    {
        if (Lt(b0.min(), b1.min())) {
            return true;
        } else if (Lt(b1.min(), b0.min())) {
            return false;
        }
        return Lt(b0.max(), b1.max());
    }
};
typedef tagged<nvis::bbox2, xavier::Lt_bbox>            p_box_type;

} // namespace xavier

#endif
