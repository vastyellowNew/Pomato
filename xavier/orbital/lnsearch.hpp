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


/** Line Search Method
 *  - Useful for minimizing a scalar function f(x,lambda) wrt
 *    parameter lambda using a 1D line search.
 *
 *  T - Template class that represents the function to optimize
 *  with it's 1st order derivatives in Funcd.
 *  The functor Funcd must be of this form:
 *   template<class T>
 *   struct Funcd { //Derivatives
        double EPS;
        T &func;
        double f;
        Funcd(T &funcc) : EPS(), func(funcc) {}
        ///Function call
        double operator()(vector &x)
            { return f = func(x); }
        ///First-order derivatives
        void df(vector &x, array &df)
     };

     With the actual scalar function as the functor (i.e., T)
     struct Func {
        double operator()(vector &x)
        {//Function evaluation}
     };
 */
#ifndef LNSEARCH_HPP
#define LNSEARCH_HPP

#include <vector>
#include <limits>

using namespace std;

//namespace Numerical Recipes {

namespace orbital {

template<class T>
void lnsearch(const vector<double>& xold, const double fold, vector<double>& g,
              vector<double>& p, vector<double>& x, double& f, const double stpmax,
              bool& check, const double& TOLX, T& func)
{
    const double ALF = 1.0e-4;
    double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
    int i,n=(int)xold.size();
    check = false;
    for (i=0; i<n; i++) {
        sum += p[i]*p[i];
    }
    sum = sqrt(sum);
    if (sum > stpmax)
        for(i=0; i<n; i++) {
            p[i] *= stpmax/sum;    //scale if attempted step is too big
        }
    for (i=0; i<n; i++) {
        slope += g[i]*p[i];
    }
    if (slope >= 0.0) {
        throw("Roundoff problem in lnsearch.");
    }
    test = 0.0;
    for(i=0; i<n; i++) {
        temp = fabs(p[i])/max(fabs(xold[i]),1.0);
        if (temp > test) {
            test = temp;
        }
    }
    alamin = TOLX/test;
    alam=1.0;
    for (;;) {
        for(i=0; i<n; i++) {
            x[i]=xold[i]+alam*p[i];
        }
        try {
            f = func(x);
        } catch(...) {
            throw("Error in Function Call");
        }
        if (alam < alamin) {
            for(i=0; i<n; i++) {
                x[i]=xold[i];
            }
            check = true;
            return;
        } else if (f <= fold+ALF*alam*slope) {
            return;    //Sufficient function decrease
        } else {
            if (alam == 1.0) {
                tmplam = -slope/(2.0*(f-fold-slope));
            } else {
                rhs1=f-fold-alam*slope;
                rhs2=f2-fold-alam2*slope;
                a=(rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a==0.0) {
                    tmplam = -slope/(2.0*b);
                } else {
                    disc = b*b-3.0*a*slope;
                    if (disc < 0.0) {
                        tmplam = 0.5*alam;
                    } else if (b <= 0.0) {
                        tmplam = (-b+sqrt(disc))/(3.0*a);
                    } else {
                        tmplam = -slope/(b+sqrt(disc));
                    }
                }
                if (tmplam>0.5*alam) {
                    tmplam = 0.5*alam;
                }
            }
        }
        alam2 = alam;
        f2 = f;
        alam = max(tmplam,0.1*alam);
    }
}

} //End namespace
#endif