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


#ifndef __ORBITS_HPP__
#define __ORBITS_HPP__

#include <vector>
#include <list>
#include <iostream>
#include <math/fixed_vector.hpp>
#include "period.hpp"

namespace map_analysis {
class orbit {
public:
    orbit() : __points(), __q(-1.) {}
    
    orbit(const std::vector<nvis::vec2>& points, double period = -1)
        : __points(points), __q(period) {}
        
    orbit(const std::vector<nvis::vec2>& points, const std::vector<nvis::vec2>& errors,
          double period = -1)
        : __points(points), __errors(errors), __q(period) {}
        
    orbit(const std::list<nvis::vec2>& points, double period = -1)
        : __points(points.begin(), points.end()), __q(period) {}
        
    orbit(const std::list<nvis::vec2>& points, const std::list<nvis::vec2>& errors, double period = -1)
        : __points(points.begin(), points.end()), __errors(errors.begin(), errors.end()), __q(period) {}
        
    size_t size() const
    {
        return __points.size();
    }
    
    const nvis::vec2& operator[](size_t i) const
    {
        assert(i < __points.size());
        return __points[i];
    }
    
    nvis::vec2& operator[](size_t i)
    {
        assert(i < __points.size());
        return __points[i];
    }
    
    const std::vector<nvis::vec2>& points() const
    {
        return __points;
    }
    
    const std::vector<nvis::vec2>& errors() const
    {
        return __errors;
    }
    
    double period() const
    {
        return __q;
    }
    
    double& period()
    {
        return __q;
    }
    
private:
    std::vector<nvis::vec2>     __points;
    std::vector<nvis::vec2>     __errors;
    double __q;
};

namespace static_data {
extern std::vector<orbit>   central_map_orbits;
}

class point_data {
public:
    point_data() : __orbit(-1), __index(-1) {}
    point_data(size_t which_orbit, size_t idx = -1)
        : __orbit(which_orbit), __index(idx) {}
        
    size_t orbit_id() const
    {
        assert(__orbit < map_analysis::static_data::central_map_orbits.size());
        return __orbit;
    }
    
    void set_orbit(size_t which_orbit)
    {
        __orbit = which_orbit;
    }
    
    double period() const
    {
        assert(__orbit < static_data::central_map_orbits.size());
        return static_data::central_map_orbits[__orbit].period();
    }
    
    const nvis::vec2& pos() const
    {
        assert(__orbit < static_data::central_map_orbits.size());
        return static_data::central_map_orbits[__orbit][__index];
    }
    
    size_t index() const
    {
        return __index;
    }
    
    double scalar_value() const
    {
        return period();
    }
    
    const nvis::vec2& error() const
    {
        const std::vector<nvis::vec2>& errors = static_data::central_map_orbits[__orbit].errors();
        return errors[__index];
    }
    
private:
    size_t  __orbit, __index;
};

inline nvis::vec2 vector_value(const point_data& pt, int period, const metric_type& metric)
{
    const orbit& chain = static_data::central_map_orbits[pt.orbit_id()];
    size_t N = chain.size();
    if (period >= N) {
        std::cerr << "WARNING: attempting to determine " << period
                  << "-map on a chain of length " << N
                  << " seeded at " << static_data::central_map_orbits[pt.orbit_id()][0]
                  << std::endl;
        double large = std::numeric_limits<double>::max();
        return nvis::vec2(large, large);
    }
    int i = pt.index();
    if (i >= period && i + period < N) {
        const nvis::vec2& x0 = chain[i-period];
        const nvis::vec2& x2 = chain[i+period];
        return 0.5*(metric.displacement(x0, x2));
    } else if (i + period < N) {
        const nvis::vec2& x1 = chain[i];
        const nvis::vec2& x2 = chain[i+period];
        return metric.displacement(x1, x2);
    } else if (i >= period) {
        const nvis::vec2& x0 = chain[i-period];
        const nvis::vec2& x1 = chain[i];
        return metric.displacement(x0, x1);
    } else {
        double large = std::numeric_limits<double>::max();
        return nvis::vec2(large, large);
    }
}

inline std::pair<nvis::vec2, nvis::vec2>
vector_and_error_value(const point_data& pt, int period,
                       const metric_type& metric)
{
    const orbit& chain = static_data::central_map_orbits[pt.orbit_id()];
    size_t N = chain.size();
    if (period >= N) {
        std::cerr << "WARNING: attempting to determine " << period
                  << "-map on a chain of length " << N
                  << " seeded at " << static_data::central_map_orbits[pt.orbit_id()][0]
                  << std::endl;
        double large = std::numeric_limits<double>::max();
        return std::make_pair(nvis::vec2(large, large), large);
    }
    int i = pt.index();
    if (i >= period && i + period < N) {
        const nvis::vec2& x0 = chain[i-period];
        const nvis::vec2& x2 = chain[i+period];
        return std::make_pair(0.5*(metric.displacement(x0, x2)), chain.errors()[i-period] + chain.errors()[i+period]);
    } else if (i + period < N) {
        const nvis::vec2& x1 = chain[i];
        const nvis::vec2& x2 = chain[i+period];
        return std::make_pair(metric.displacement(x1, x2), chain.errors()[i] + chain.errors()[i+period]);
    } else if (i >= period) {
        const nvis::vec2& x0 = chain[i-period];
        const nvis::vec2& x1 = chain[i];
        return std::make_pair(metric.displacement(x0, x1), chain.errors()[i-period] + chain.errors()[i]);
    } else {
        double large = std::numeric_limits<double>::max();
        return std::make_pair(nvis::vec2(large, large), large);
    }
}

inline nvis::vec2
vector_value(size_t orbit_id, size_t idx, size_t period, const metric_type& metric)
{
    point_data pd(orbit_id, idx);
    return vector_value(pd, period, metric);
}

inline std::pair<nvis::vec2, nvis::vec2>
vector_and_error_value(size_t orbit_id, size_t idx, size_t period,
                       const metric_type& metric)
{
    point_data pd(orbit_id, idx);
    return vector_and_error_value(pd, period, metric);
}

}


#endif







