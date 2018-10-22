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


#ifndef __RATIONAL_HPP__
#define __RATIONAL_HPP__

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <boost/rational.hpp>

namespace xavier {

template<typename I, typename F>
inline F value(const boost::rational<I>& r)
{
    return (F)r.numerator() / (F)r.denominator();
}

template<typename I, typename F>
inline boost::rational<I> rational_approx(F v, I maxnum, I maxden=0)
{
    // naive method
    typedef I                       _int;
    typedef F                       _float;
    typedef boost::rational<_int>   _rational;

    if (!maxden) {
        maxden = maxnum;
    }
    std::map<_float, _rational> approx;
    for (_int num=0 ; num<=maxnum ; ++num) {
        for (_int den=1 ; den<=maxden ; ++den) {
            _rational r(num, den);
            _float err = fabs(fabs(v)-value<_int, _float>(r));
            approx.insert(std::pair<_float, _rational>(err, r));
        }
    }
    _rational r = approx.begin()->second;
    if (v < 0) {
        r *= _rational(-1, 1);
    }
    return r;
}

/** Continued Fraction algorithm : Improved version for performance (WRS)
 *  Reference: http://en.wikipedia.org/wiki/Continued_fraction Under "infinite continued fractions"
 *  Basically, we compute the best rational approximation based on the denominator.
 *  Default max denominator (md) is 64, but it's best to use p_max in the topology
 *  extraction algorithm.
 */
template<typename I, typename F>
inline boost::rational<I> rational_approx_CF(F f, const I md=64)
{
    typedef boost::rational<I>  _rational;
    /*  a: continued fraction coefficients. */
    long a, h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
    long x, d, n = 1;
    long i, neg = 0;

    //A bad max denominator (md)
    if (md <= 1) {
        return _rational((I) f,1);
    }

    if (f < 0) {
        neg = 1;
        f = -f;
    }

    while (f != floor(f)) {
        n <<= 1;
        f *= 2;
    }
    d = f;

    /* continued fraction and check denominator each step */
    for (i = 0; i < 64; i++) {
        a = n ? d / n : 0;
        if (i && !a) {
            break;
        }

        x = d;
        d = n;
        n = x % n;

        x = a;
        if (k[1] * a + k[0] >= md) {
            x = (md - k[0]) / k[1];
            if (x * 2 >= a || k[1] >= md) {
                i = 65;
            } else {
                break;
            }
        }

        h[2] = x * h[1] + h[0];
        h[0] = h[1];
        h[1] = h[2];
        k[2] = x * k[1] + k[0];
        k[0] = k[1];
        k[1] = k[2];
    }

    try {
        _rational r(neg ? -h[1] : h[1], k[1]);
        //*denom = k[1];
        //*num = neg ? -h[1] : h[1];
        return r;
    }
    catch (boost::bad_rational& e) {
        return _rational((I) (f*md));
    }
}


}



#endif
