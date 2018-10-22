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


#ifndef __ANGLE_HPP__
#define __ANGLE_HPP__

#include <math/fixed_vector.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <math.h>

//isnan macro (windows v linux)
#ifdef _WIN32
#include <boost/math/special_functions/fpclassify.hpp>
#define ISNANCALL(X) { boost::math::isnan( X ) }
#define ISINFCALL(X) { boost::math::isinf( X ) }
#else
#define ISNANCALL(X) ( (bool) std::isnan( (double) X ) )
#define ISINFCALL(X) ( (bool) std::isinf( (double) X ) )
#endif

namespace {
inline double sign(const double& x)
{
    return (x >= 0 ? 1 : -1);
}
}

namespace xavier {
const double TWO_PI = 2*M_PI;

inline double signed_angle(const nvis::vec2& x)
{
    double l = nvis::norm(x);
    if (!l) {
        return 0;
    }
    
    double theta = acos(x[0] / l);
    
    if (theta > M_PI) {
        throw std::runtime_error("invalid angle value\n");
    }
    
    if (x[1] < 0) {
        theta *= -1;
    }
    return theta;
}

inline double positive_angle(const nvis::vec2& x)
{
    double theta = signed_angle(x);
    return theta < 0 ? theta + 2.*M_PI : theta;
}

inline double positive_modulo(double angle, double period = 2*M_PI)
{
    return angle >= 0 ? fmod(angle, period) : period + fmod(angle, period);
}

inline double modulo(double angle, double period = 2*M_PI)
{
    return fmod(angle, period);
}

inline double signed_angle(const nvis::vec2& v0, const nvis::vec2& v1)
{
    if (!nvis::norm(v0) || !nvis::norm(v1)) {
        return 0;
    }
    nvis::vec2 w0 = v0 / nvis::norm(v0);
    nvis::vec2 w1 = v1 / nvis::norm(v1);
    double dot = nvis::inner(w0, w1);
    if (dot <= -1) {
        return M_PI;
    } else if (dot >= 1) {
        return 0;
    }
    double theta = acos(dot);
    bool isnanTest = ISNANCALL( theta );
    if ( isnanTest) {
        std::ostringstream os;
        os << "signed angle return nan. v0 = " << v0 << ", v1 = " << v1 << ", w0 = " << w0 << ", w1 = " << w1 << std::endl;
        std::cerr << os.str();
    }
    double det = w0[0] * w1[1] - w0[1] * w1[0];
    return sign(det)*theta;
}

template<int N>
inline double unsigned_angle(const nvis::fixed_vector<double, N>& v0, const nvis::fixed_vector<double, N>& v1)
{
    typedef typename nvis::fixed_vector<double, N>  vec_type;
    
    if (!nvis::norm(v0) || !nvis::norm(v1)) {
        return 0;
    }
    vec_type w0 = v0 / nvis::norm(v0);
    vec_type w1 = v1 / nvis::norm(v1);
    return acos(nvis::inner(w0, w1));
}

inline double to_degrees(double rad)
{
    return rad/M_PI*180.;
}


}

#endif




