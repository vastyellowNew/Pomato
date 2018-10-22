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


#ifndef __XAVIER_MATH_HPP__
#define __XAVIER_MATH_HPP__

#include <complex>
#include <vector>
#include <map>
#include <algorithm>
#include <math/fixed_vector.hpp>

namespace xavier {

inline bool bary_coords(nvis::vec3& b, const nvis::vec3& p,
                        const nvis::vec3& p0, const nvis::vec3& p1, const nvis::vec3& p2)
{
    // bring problem back to 2D
    nvis::vec3 q, e1(p1), e2(p2);
    double x, y, x1, x2, y2, det, delta;
    
    e1 -= p0;
    e2 -= p0;
    x1 = nvis::norm(e1);
    e1 *= 1 / x1; // e1 normalized
    x2 = nvis::inner(e2, e1);
    e2 -= x2 * e1; // e2 orthogonal to e1
    y2 = nvis::norm(e2);
    e2 *= 1 / y2;
    
    q = p - p0;
    x = nvis::inner(q, e1);
    y = nvis::inner(q, e2);
    
    delta = x1 * y2; // NB: y1==0
    b[1] = (x * y2 - y * x2) / delta;
    b[2] = x1 * y / delta;
    b[0] = 1.0 - b[1] - b[2];
    
    return (b[0] >= 0 && b[1] >= 0 && b[2] >= 0);
}


inline int quadratic_equation(double a, double b, double c,
                              std::complex<double> x[2])
// PAR: a*x^2+b*x+c=0.0, x[0..1] gets the complex solutions
// POST: x[0], x[1] contain the solution
// RETURN: number of different solutions
//         (-1,0,1,2) -1 stands for infinite many!
{
    // mario
    if (fabs(a) < 1e-9) {
        // linear equation
        if (fabs(b) < 1e-9) {
            return 0;
        }
        // FIXME: there may be an infinite number, how to deal with that?
        x[0] = std::complex<double>(-c / b, 0.);
        return 1;
    }
    // end mario
    
    double d = b * b - 4 * a * c;
    if (d > 0.0) {
        double q = b >= 0.0 ? -0.5 * (b + sqrt(d)) : -0.5 * (b - sqrt(d));
        x[0] = (q / a); // x[0]=q/a+i*0
        x[1] = (c / q); // x[1]=c/q+i*0
        return 2;
    }
    if (d == 0.0) {
        if (a != 0.0) {
            x[0] = std::complex<double>(-0.5 * b / a);
            // although there is only one solution
            // set x[1] for convenience to the same value as x[0]
            x[1] = x[0];
            return 1;
        } else if (c == 0.0) {
            return -1;
        } else {
            return 0;
        }
    }
    // ASSERT: d<0.0
    std::complex<double> dd = d;
    dd = sqrt(dd);
    dd = (real(b * dd) >= 0) ? dd : -dd;
    std::complex<double> q = -0.5 * (b + dd);
    x[0] = (q / a);
    x[1] = (c / q);
    return 2;
}

inline int cubic_equation(double a3, double a2, double a1, double a0,
                          std::complex<double> x[3])
// PAR: a3*x^3+a2*x^2+a1*x+a0=0.0, x[0..2] gets the solution
// POST: x[0], x[1], x[2] contain the solutions
// RETURN: number of different solutions
//         (-1,0,1,2,3) -1 stands for infinite many!
// REMARK:  ideas taken from "Numerical Recipes in C", p.184/185
{
    if (a3 == 0.0) {
        int n = quadratic_equation(a2, a1, a0, x);
        x[2] = x[1];
        return n;
    }
    a2 /= a3;
    a1 /= a3;
    a0 /= a3;
    double Q = (a2 * a2 - 3.0 * a1) / 9.0;
    double R = (2.0 * a2 * a2 * a2 - 9.0 * a2 * a1 + 27.0 * a0) / 54.0;
    double comp = R * R - Q * Q * Q;
    if (comp >= 0.0) {
        double sgn_R = (R >= 0.0) ? 1.0 : -1.0;
        double A = fabs(R) + sqrt(comp);
        A = pow(A, 1.0 / 3.0);
        A *= (-sgn_R);
        double B = (A != 0.0) ? Q / A : 0.0;
        x[0] = std::complex<double>((A + B) - a2 / 3.0);
        x[1] = std::complex<double>(-0.5 * (A + B) - a2 / 3.0, 0.5 * sqrt(3.0) * (A - B));
        x[2] = std::complex<double>(-0.5 * (A + B) - a2 / 3.0, -0.5 * sqrt(3.0) * (A - B));
        return 3;
    }
    double theta = acos(R / sqrt(Q * Q * Q));
    x[0] = std::complex<double>(-2.0 * sqrt(Q) * cos(theta / 3.0) - a2 / 3);
    x[1] = std::complex<double>(-2.0 * sqrt(Q) * cos((theta + 2 * M_PI) / 3.0) - a2 / 3);
    x[2] = std::complex<double>(-2.0 * sqrt(Q) * cos((theta - 2 * M_PI) / 3.0) - a2 / 3);
    return 3;
}

}

#endif








