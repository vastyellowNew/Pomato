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


#ifndef MAPS_WINDING_HPP
#define MAPS_WINDING_HPP

#include <math/angle.hpp>
#include <map>
#include <set>

namespace xavier {

template<typename Metric, int Dim>
struct periodic_dim_winding_number {

    static const int periodic_dim = Dim;
    
    periodic_dim_winding_number() : _metric() {}
    periodic_dim_winding_number(const periodic_dim_winding_number& wn) : _metric(wn._metric) {}
    
    template<typename Iterator>
    double operator()(const Iterator& begin, const Iterator& end) const
    {
    
        typedef typename Iterator::value_type value_type;
        const value_type& x0 = *begin;
        
        double w = 0;
        double number = 0;
        Iterator it=begin;
        size_t counter;
        for (++it, counter=1 ; it!=end ; ++it, ++counter) {
            const value_type& x = *it;
            double d = nvis::norm(_metric.displacement(x0, x));
            if (d == 0) {
                return (x[periodic_dim] - x0[periodic_dim])/(float)counter;
            } else {
                w += 1./d; // weight is inverse propertional to periodic distance to x0
                number += (x[periodic_dim] - x0[periodic_dim])/(float)counter/d;
            }
        }
        return number/w;
    }
    
    Metric _metric;
};

struct planar_rotation_winding_number {

    planar_rotation_winding_number(bool fwd=true) : _forward(fwd) {}
    planar_rotation_winding_number(const planar_rotation_winding_number& wn) : _forward(wn._forward) {}
    
    template<typename Iterator>
    double operator()(const Iterator& begin, const Iterator& end) const
    {
        typedef typename Iterator::value_type value_type;
        const value_type& x0 = *begin;
        
        double w = 0;
        double r = 0;
        
        // compute centroid
        value_type center(0);
        size_t counter=0;
        for (Iterator it=begin ; it!=end; ++it, ++counter) {
            center += *it;
        }
        center /= (double)counter;
        
        // compute weighted angular speed about this center
        value_type last = x0;
        Iterator it = begin;
        double theta;
        for (++it, counter=1 ; it!=end ; ++it, ++counter) {
            value_type next = *it;
            double dtheta = signed_angle(last-center, next-center);
            if (_forward && dtheta<0) {
                dtheta += 2*M_PI;
            } else if (!_forward && dtheta>0) {
                dtheta -= 2.*M_PI;
            }
            theta += dtheta;
        }
        
        return 2.*M_PI*counter/fabs(theta);
    }
    
    bool _forward;
};

double max_angle(const std::vector<nvis::vec2>& points)
{
    nvis::vec2 center(0);
    for (int i=0 ; i<points.size() ; ++i) {
        center += points[i];
    }
    center /= (double)points.size();
    
    std::set<double> angles;
    for (int i=0 ; i<points.size() ; ++i) {
        double theta = xavier::positive_angle(points[i]-center);
        angles.insert(theta);
    }
    
    double max_theta = -1;
    for (std::set<double>::const_iterator it=angles.begin() ; it!=angles.end() ; ++it) {
        std::set<double>::const_iterator jt = it;
        ++jt;
        if (jt==angles.end()) {
            jt=angles.begin();    // make it circular
        }
        std::set<double>::const_iterator kt = jt;
        ++kt;
        if (kt==angles.end()) {
            kt=angles.begin();    // make it circular
        }
        double theta = *jt - *it;
        if (theta < 0) {
            theta += 2.*M_PI;
        }
        if (theta > max_theta) {
            max_theta = theta;
        }
        theta = *kt - *jt;
        if (theta < 0) {
            theta += 2.*M_PI;
        }
        if (theta > max_theta) {
            max_theta = theta;
        }
        if (theta < 0) {
            theta += 2.*M_PI;
        }
        if (theta > max_theta) {
            max_theta = theta;
        }
    }
    return max_theta;
}

}

#endif