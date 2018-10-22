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


#ifndef __COLORS_HPP__
#define __COLORS_HPP__

#include <math/fixed_vector.hpp>
#include <vector>
#include <map>

namespace xavier {

inline nvis::vec3 hsv2rgb(double hue, double sat, double val)
{
    double chroma = val * sat;
    double h_ = hue / 60.;
    double x = chroma * (1. - fabs(fmod(h_, 2.) - 1));
    
    nvis::vec3 rgb_;
    if(val == 0) {
        rgb_ = nvis::vec3(0, 0, 0);
    } else if(0 <= h_ && h_ < 1) {
        rgb_ = nvis::vec3(chroma, x, 0);
    } else if(1 <= h_ && h_ < 2) {
        rgb_ = nvis::vec3(x, chroma, 0);
    } else if(2 <= h_ && h_ < 3) {
        rgb_ = nvis::vec3(0, chroma, x);
    } else if(3 <= h_ && h_ < 4) {
        rgb_ = nvis::vec3(0, x, chroma);
    } else if(4 <= h_ && h_ < 5) {
        rgb_ = nvis::vec3(x, 0, chroma);
    } else if(5 <= h_ && h_ < 6) {
        rgb_ = nvis::vec3(chroma, 0, x);
    }
    
    double m = val - chroma;
    return rgb_ + nvis::vec3(m, m, m);
}

const nvis::fvec3 red(1, 0, 0);
const nvis::fvec3 green(0, 1, 0);
const nvis::fvec3 blue(0, 0, 1);
const nvis::fvec3 yellow(1, 1, 0);
const nvis::fvec3 cyan(0, 1, 1);
const nvis::fvec3 white(1, 1, 1);
const nvis::fvec3 black(0, 0, 0);
const nvis::fvec3 orange(1, 0.5, 0);
const nvis::fvec3 magenta(1, 0, 1);
// some funky color names after Apple's color editor
const nvis::fvec3 cayenne(0.5, 0, 0);
const nvis::fvec3 midnight(0, 0, 0.5);
const nvis::fvec3 aqua(0, 0.5, 1);
const nvis::fvec3 sea_green(0, 1, 0.5);
const nvis::fvec3 lime(0.5, 1, 0);
const nvis::fvec3 framboise(1, 0, 0.5);
const nvis::fvec3 carnation(1, 0.5, 1);

const nvis::fvec3 rainbow[] = {
    black, midnight, blue, aqua,
    cyan, sea_green, green, lime,
    yellow, orange, red, framboise,
    magenta, carnation, white
};

static void spiral_scale(std::vector<nvis::fvec3>& colors, int n, double minval, int r=1)
{
    double dhue = (double)r*360./(double)n;
    double dval = (1-minval)/(double)n;
    colors.resize(n);
    double hue = 0;
    double val = minval;
    for (int i=0 ; i<n ; ++i) {
        colors[i] = hsv2rgb(hue, 1, val);
        hue += dhue;
        if (hue > 360) {
            hue -= 360;
        }
        val += dval;
    }
}

template<typename T>
struct adaptive_color_map {
    adaptive_color_map(const std::vector<T>& vals, const std::vector<nvis::fvec3>& scale)
        : t(scale.size()), colors(scale)
    {
        std::vector<T> tmp(vals);
        std::sort(tmp.begin(), tmp.end());
        unsigned int step = tmp.size() / (scale.size() - 1);
        for (int i = 0 ; i < t.size() ; ++i) {
            t[i] = tmp[i*step];
        }
        t.back() = tmp.back();
    }
    
    nvis::fvec3 operator()(const T& val) const
    {
        unsigned int bin = std::distance(t.begin(),
                                         std::lower_bound(t.begin(), t.end(), val));
        if (bin > t.size() - 2) {
            bin = t.size() - 2;
        }
        T u = (val - t[bin]) / (t[bin+1] - t[bin]);
        return (1. - u)*colors[bin] + u*colors[bin+1];
    }
    
    std::vector<T>              t;
    std::vector<nvis::fvec3>    colors;
};

template<typename T>
struct discrete_color_map {

    discrete_color_map(const std::vector<T>& cps, const std::vector<nvis::fvec3>& scale)
        : colors(scale)
    {
        assert(cps.size() == colors.size());
        for (int i=0 ; i<cps.size() ; ++i) {
            lookup[cps[i]] = i;
        }
    }
    
    nvis::fvec3 operator()(const T& val) const
    {
        typename std::map<T, int>::const_iterator it = lookup.find(val);
        if (it == lookup.end()) {
            return nvis::fvec3(0,0,0);
        } else {
            return colors[it->second];
        }
    }
    
    std::map<T, int>            lookup;
    std::vector<nvis::fvec3>    colors;
};

template<typename T>
struct fixed_color_map {

    fixed_color_map(const std::vector<T>& cps, const std::vector<nvis::fvec3>& scale)
        : t(cps), colors(scale)
    {
        std::sort(t.begin(), t.end());
        assert(t.size() == colors.size());
    }
    
    nvis::fvec3 operator()(const T& val) const
    {
        int bin = std::distance(t.begin(), std::lower_bound(t.begin(), t.end(), val))-1;
        if (bin<0) {
            bin=0;
        } else if (bin > t.size() - 2) {
            bin = t.size() - 2;
        }
        T u = (val - t[bin]) / (t[bin+1] - t[bin]);
        return (1. - u)*colors[bin] + u*colors[bin+1];
    }
    
    std::vector<T>              t;
    std::vector<nvis::fvec3>    colors;
};

template<typename T>
struct band_color_map {
    band_color_map(const std::vector<T>& cps, const std::vector<nvis::fvec3>& scale)
        : t(cps), colors(scale)
    {
        std::sort(t.begin(), t.end());
        assert(t.size()+2 == colors.size());
    }
    
    nvis::fvec3 operator()(const T& val) const
    {
        unsigned int bin = std::distance(t.begin(),
                                         std::lower_bound(t.begin(), t.end(), val));
        return colors[bin];
    }
    
    std::vector<T>              t;
    std::vector<nvis::fvec3>    colors;
};

} // xavier

#endif
