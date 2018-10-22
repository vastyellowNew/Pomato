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


/** Quasi-Newton Method for Multidim Root-finding
 *  Uses a Line Search Method with backtracking for updates
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
#ifndef QUASI_NEWTON_HPP
#define QUASI_NEWTON_HPP

#include <vector>
#include <limits>
#include <iostream>
//Eigen Lib
#include <Eigen/Core>
#include <Eigen/Dense>
//Linesearch from Numerical Recipes
#include <orbital/lnsearch.hpp>

using namespace std;
using namespace Eigen;

//namespace NumericalRecipes {

namespace orbital {
/// Numerical Jacobian via forward difference
template<class T>
struct NRfdjac {
    const double EPS; //Set to approximate sqrt of machine precision
    T& func;
    
    //Init with user-supplied function or functor that returns
    //the vector of functions to be zeroed
    NRfdjac(T& funcc) : EPS(1.e-8), func(funcc) {}
    
    ///Compute the Jacobian and return in df
    void computeDF(vector<double>& x, vector<double>& fvec, MatrixXd& df)
    {
        //computes Jacobian matrix - Assumes already sized to (m,n)
        df.setZero();
        int n = (int)x.size();
        int m = (int)fvec.size();
        vector<double> xh = x;
        for(int j=0; j<n; j++) {
            double temp = xh[j];
            double h = EPS*fabs(temp);
            if (h == 0.0) {
                h = EPS;
            }
            xh[j] = temp + h;
            h = xh[j] - temp;
            vector<double> f(fvec.size(),0.0);
            try {
                func.eval(xh,f);
            } catch(...) {
                throw("Function Call Error");
            }
            xh[j] = temp;
            for (int i=0; i<m; i++) {
                df(i,j) = (f[i] - fvec[i])/h;
            }
        }
        //DF is computed as a MatrixXd object
    }
};

/// Returns f = 0.5*F*F. Also stores value of F in fvec
template <class T>
struct NRfmin {
public :
    vector<double> fvec;
    T& func;
    //Init with user supplied function/functor
    NRfmin(T& funcc) : func(funcc) {}
    double operator() (vector<double>& x)
    {
        //int n = x.size();
        double sum = 0.0;
        try {
            func.eval(x,fvec);
        } catch(...) {
            throw("Function Call Error");
        }
        int m = (int)fvec.size();
        for (int i=0; i<m; i++) {
            sum += fvec[i]*fvec[i];
        }
        return 0.5*sum;
    }
};


/// A QuasiNewton Root-finding method.  Based on Numerical Recipes version.
template<class T>
class QuasiNewton {
public :
    QuasiNewton(T& vecfunc, bool output=false) :
        func(vecfunc),
        analyticDF(false),
        verbose(output)
    {
        MAXITS = 200;
        TOLF = 1.0e-8;
        TOLMIN = 1.0e-12;
        STPMX = 100.0;
        TOLX = numeric_limits<double>::epsilon();
        SLOWTOL = 0.01; //If haven't reduced error by more than 1% 5x
        jobID = 0;
        numIters = 0;
        analyticDF = func.derivsSupplied();
    }
    
    void setMaxIts(const int& it)
    {
        MAXITS = it;
    }
    void setTol(const double& tol)
    {
        TOLF = tol;
    }
    void setJobID(const int& id)
    {
        jobID = id;
    }
    int getNumIters()
    {
        return numIters;
    }
    
    /// Solver function
    void solve(vector<double>& x, double& f, double& f0, bool& check)
    {
        int i, j, its, n=(int)x.size();
        int slowSteps = 0;
        double den,fold,stpmax,sum,temp,test;
        MatrixXd DF;
        vector<double> g(n),p(n),xold(n,0.0);
        VectorXd dX,F;
        //Create objects
        NRfmin<T> fmin(func);
        vector<double>& fvec = fmin.fvec; //Make an alias
        NRfdjac<T> fdjac(func);
        //First call - init
        f = fmin(x);
        int m = (int)fvec.size();
        F.resize(m);
        dX.resize(n); //search direction (dX = p)
        DF.resize(m,n);
        test = 0.0; //Test for initial guess being a root.
        for (j=0; j<m; j++) {
            if (fabs(fvec[j]) > test) {
                test = fabs(fvec[j]);
            }
        }
        if (test < 0.01*TOLF) {
            check = false;
            return;
        }
        //Calculate stpmax for line searches
        sum=0.0;
        for(i=0; i<n; i++) {
            sum += x[i]*x[i];
        }
        stpmax = STPMX*max(sqrt(sum),(double) n);
        
        //Iteration loop
        for (its=0; its<MAXITS; its++) {
            numIters = its;
            //Build derivatives - analytic or numerical
            if (analyticDF) {
                func.df(x,DF);
            } else {
                fdjac.computeDF(x,fvec,DF);
            }
            //Gradient f for line search
            for (i=0; i<n; i++) {
                sum = 0.0;
                for (j=0; j<m; j++) {
                    sum += DF(j,i)*fvec[j];
                }
                g[i] = sum;
            }
            for(i=0; i<n; i++) {
                xold[i]=x[i];
            }
            fold = f;
            for(j=0; j<m; j++) {
                F(j) = -fvec[j];    //RHS for linear eqn
            }
            
            //Solve linear equations with Eigen Lib
            dX = DF.colPivHouseholderQr().solve(F); //Min-Norm if(m<n) or Newton (m=n)
            for(i=0; i<n; i++) {
                p[i] = dX(i);
            }
            
            //Line search on step - returns p,g,f,x
            try {
                lnsearch<NRfmin<T> >(xold,fold,g,p,x,f,stpmax,check,TOLX,fmin);
            } catch(...) {
                if (verbose) {
                    std::cerr << "(" << jobID << ") QN: Iter = " << its+1 << " encountered function call error\n";
                }
                throw("Function call error");
            }
            
            if(its == 0) {
                f0 = fold;
                if (verbose) {
                    std::cerr << "(" << jobID << ") QN: Iter = 0  ||F(X)|| = " << sqrt(2.0*fold) << " \n";
                }
            }
            if (verbose) {
                std::cerr << "(" << jobID << ") QN: Iter = " << its+1 << " ||F(X)|| = " << sqrt(2.0*f) << " \n";
                //Print out x and F
                /*  std::vector<double> tempF;
                  try {
                    func.eval(x,tempF);
                    std::cerr << "(" << jobID << ") QN: Iter = " << its << " X = [ ";
                    for(int i=0;i<n;i++) {
                  std::cerr << x[i] ;
                  if (i<n-1) std::cerr << ", ";
                    }
                    std::cerr << "]\n";
                    std::cerr << "(" << jobID << ") QN: Iter = " << its << " F = [ ";
                    for(int i=0;i<m;i++) {
                  std::cerr << tempF[i] ;
                  if (i<m-1) std::cerr << ", ";
                    }
                    std::cerr << "]\n";
                  } catch(...) {
                    std::cerr << "PoincareMapError\n";
                  }*/
                
            }
            
            //Test for convergence
            test = 0.0;
            for (i=0; i<n; i++) if (fabs(fvec[i]) > test) {
                    test = fabs(fvec[i]);
                }
            if (test < TOLF) {
                check = false;
                return;
            }
            //Check for spurious convergence (grad f = 0)
            if (check) {
                //Relative gradient is (deltaf/f)/(deltax/x) [and deltaf~=(grad f) * deltax]
                test = 0.0;
                den = max(f,0.5*n);
                for (i=0; i<n; i++) {
                    temp = fabs(g[i])*max(fabs(x[i]),1.0)/den;
                    if (temp > test) {
                        test = temp;
                    }
                }
                check = (test<TOLMIN) ? true : false;
                return;
            }
            //Test for convergence on dX
            test = 0.0;
            for (i=0; i<n; i++) {
                temp = (fabs(x[i] - xold[i]))/max(fabs(x[i]),1.0);
                if (temp > test) {
                    test = temp;
                }
            }
            if (test < TOLX) {
                return;
            }
            //Test for extremely slow convergence, and exit if occurring
            if(its>0) {
                if( fabs((f-fold)/fold) < SLOWTOL ) {
                    slowSteps++;
                } else {
                    slowSteps = 0;
                }
                if (slowSteps >= 5) {
                    throw("Change in F is too slow in quasiNewton.solve()");
                }
            }
        } //End main loop
        throw("MAXITS exceeded in quasiNewton.solve()");
    }
    
    bool isJacobianSupplied()
    {
        return analyticDF;
    }
    
private :
    /// F(X) functor class
    T& func;
    int MAXITS;
    double TOLF,TOLMIN,STPMX,TOLX;
    double SLOWTOL;
    int jobID, numIters;
    int rows,cols;
    bool analyticDF;
    bool verbose;
};



}//End namespace
#endif