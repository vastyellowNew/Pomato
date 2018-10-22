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


#ifndef __MAPS_LIB_NEWTON_HPP
#define __MAPS_LIB_NEWTON_HPP

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <vector>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <limits>
#include "metric.hpp"

namespace xavier {
// 1D central difference derivative
template<typename Func, int N>
inline nvis::fixed_vector<double, N> cd1d(const Func& f, const nvis::fixed_vector<double, N>& x, double h,  int dim)
{
    nvis::fixed_vector<double, 2> dx(0);
    dx[dim] = h;
    return 0.5*(f(x+dx)-f(x-dx))/h;
}

template<typename Func, int N>
inline nvis::fixed_matrix<double, N> central_differences(const Func& f, const nvis::fixed_vector<double, N>& x, double h)
{
    typedef nvis::fixed_matrix<double, N> mat_type;
    
    mat_type Jt(0);
    for (int i=0 ; i<N ; ++i) {
        Jt[i] = cd1d<Func, N>(f, x, h, i);
    }
    return nvis::transpose(Jt);
}

// Richardson's extrapolation method
template<typename Func, int N>
inline nvis::fixed_matrix<double, N> richardson(const Func& f, const nvis::fixed_vector<double, N>& x, double h, double t=2.)
{
    typedef nvis::fixed_matrix<double, N> mat_type;
    
    double ht = h/t;
    double t2 = t*t;
    
    mat_type Jh = central_differences<Func, N>(f, x, h);
    mat_type Jht = central_differences<Func, N>(f, x, ht);
    
    mat_type J = 1/(t2-1)*(t2*Jht - Jh);
    return J;
}

extern bool record_newton_steps;
extern std::vector<nvis::vec2> newton_steps;

/// RHS_ONLY_WRAPPER - indicates that there is no analytical Jacobian function (Use Richardson Extrapolation)
template<typename MAP, int N>
struct rhs_only_wrapper {
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::fixed_matrix<double, N>   mat_type;
    typedef rhs_only_wrapper<MAP, N>        self_type;
    
    rhs_only_wrapper(const MAP& pmap, const metric<double, N>& metric, int period)
        :  _p(period), _h(0.1), _map(pmap), _metric(metric) {}
        
    void set_h(double h)
    {
        _h = h;
    }
    
    vec_type operator()(const vec_type& x) const
    {
        vec_type y = _map.map(x, _p);
        return _metric.displacement(x, y);
    }
    
    mat_type jacobian(const vec_type& x, double h=0) const
    {
        if (!h) {
            h = _h;
        }
        //If we have no Jacobian, use richardson (stencil+centralDiff)
        return richardson<self_type, N>(*this, x, h, 10.);
    }
    
    int _p;
    double _h;
    const MAP& _map;
    const metric<double, N> _metric;
};

/// Wrapper function for the Right-Hand Side of Newton Solution (f(x) = P^p(x)-x => Displacement)
template<typename MAP, int N>
struct rhs_wrapper {
    typedef MAP                             map_type;
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::fixed_matrix<double, N>   mat_type;
    typedef rhs_wrapper<MAP, N>             self_type;
    typedef typename MAP::return_type       return_type;
    
    rhs_wrapper(const map_type& pmap, const metric<double, N>& metric, int period, double Jeps=0)
        : _map(pmap), _metric(metric), _p(period), _jeps(Jeps),
          _lastx(std::numeric_limits<double>::max(), 0)
    {}
    
    vec_type operator()(const vec_type& x) const
    {
        if (nvis::all(x == _lastx)) {
            return _lastv;
        } else {
            try {
                return_type rmi = _map.map_complete(x, _p);
                _lastx = x;
                _lastv = _metric.displacement(x, rmi.x);
                _lastJ = rmi.J;
                for (int i=0 ; i<N ; ++i) {
                    _lastJ[i][i] -= 1;    //Subtract Identity matrix
                }
            } catch(std::runtime_error& e) {
                _lastx[0] = std::numeric_limits<double>::max();
                throw;
            }
            return _lastv;
        }
    }
    
    /// Jacobian of Newton Function (df(x)/dx)
    mat_type jacobian(const vec_type& x) const
    {
        if (nvis::all(x == _lastx)) {
            return _lastJ;
        } else {
            try {
                return_type rmi = _map.map_complete(x, _p);
                _lastx = x;
                _lastv = _metric.displacement(x, rmi.x);
                _lastJ = rmi.J;
                for (int i=0 ; i<N ; ++i) {
                    _lastJ[i][i] -= 1;    //Subtract Identity matrix
                }
            } catch(std::runtime_error& e) {
                _lastx[0] = std::numeric_limits<double>::max();
                throw;
            }
            return _lastJ;
        }
    }
    
    mat_type cd_jacobian(const vec_type& x, const double h=0.05) const
    {
        return central_differences<self_type, N>((*this), x, h);
    }
    
    int _p;
    double _jeps;
    const map_type& _map;
    const metric<double, N> _metric;
    mutable vec_type _lastx, _lastv;
    mutable mat_type _lastJ;
};

template<int N>
struct box_constraint {
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::bounding_box<vec_type>    box_type;
    
    box_constraint(const box_type& box) : _box(box) {}
    
    vec_type operator()(const vec_type& from, const vec_type& to) const
    {
        return to;
        assert(_box.inside(from));
        if (_box.inside(to)) {
            return to;
        }
        vec_type mid = 0.5*(from + to);
        return (*this)(from, mid);
    }
    
    double size() const
    {
        return nvis::norm(_box.size());
    }
    
    box_type _box;
};

// line search method used to increase robustness of Newton's method
template<typename RHS, int N>
bool lnsearch(const RHS& rhs, nvis::fixed_vector<double, N>& x, nvis::fixed_vector<double, N>& f,
              const nvis::fixed_vector<double, N>& dd, double maxlength)
{
    typedef nvis::fixed_vector<double, N>   vec_type;
    
    double lambda = 1.0;
    const double alpha = 1e-4;
    
    vec_type xsave = x, fsave = f, d = norm(dd) > maxlength ? dd * maxlength / norm(dd) : dd;
    
    for (unsigned int i = 0; i < 7; ++i) {
        x = xsave + lambda * d;
        f = rhs(x);
        if (norm(f) < (1 - alpha*lambda)*norm(fsave)) {
            return true;
        }
        
        lambda *= 0.5;
    }
    
    return false;
}

// A Line Search with Backtracking given a max step -> Uses a more advanced backtracking
template<typename RHS, int N>
bool lnsearchMinStrategy(const RHS& rhs,
                         nvis::fixed_vector<double, N>& x, nvis::fixed_vector<double, N>& f,
                         const nvis::fixed_vector<double, N>& dd, double maxlength)
{
    // x, f -> original state and constrant vector
    // dd -> Newton search direction (input as p = -J^-1*F(x))
    // This will find a new "x" that is along the Newton search
    //   direction that decreases "f" 'sufficiently'
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::fixed_matrix<double, N>   mat_type;
    
    double lambda = 1.0; //Start with Newton step
    const double alpha = 1e-4;
    const double tolx = std::numeric_limits<double>::epsilon();
    
    vec_type xsave = x;
    vec_type fsave = f;
    double fold = 0.5*nvis::inner(f,f);
    double fscalar;
    vec_type p = norm(dd) > maxlength ? dd*maxlength/norm(dd) : dd;
    
    //Compute gradient information
    mat_type JT = nvis::transpose(rhs.jacobian(x));
    vec_type gradient(0); //Gradient
    // f * J : vector * matrix
    for (int i=0; i<N; i++) {
        gradient[i] = nvis::inner(f,JT[i]);
    }
    double slope = nvis::inner(gradient,p);
    if (slope >= 0.0) {
        std::cerr << "Roundoff problem in lnsearchMinStrategy.\n";
    }
    
    //Compute lambda_min
    double test = 0.0, temp = 0.0;
    for (int i=0; i<N; i++) {
        temp = abs(p[i]) / ( abs(x[i]) > 1.0 ? abs(x[i]) : 1.0);
        if (temp > test) {
            test = temp;
        }
    }
    double lambda_min = tolx / test;
    
    //Iterate to find successful step
    int maxIters = 500, iter = 0;
    double lambda2 = lambda, fscalar2 = fold, tmplam = 0.0; //Initialize loop variables
    for (;;) {
        //Compute new x and f
        x = xsave + lambda*p;
        f = rhs(x);
        fscalar = 0.5 * inner(f,f);
        if (lambda < lambda_min) {
            //Convergence on Delta_x -> calling program should verify converngence
            std::cerr << " lnsearchMS: lambda = " << lambda << " < lambda_min = " << lambda_min << '\n';
            x = xsave;
            return true;
        } else if (iter == maxIters) {
            std::cerr << "Too many iterations in lnsearchMinStrategy.\n";
            x = xsave;
            return false;
        } else if (fscalar <= fold + alpha*lambda*slope) {
            //Sufficient function decrease
            std::cerr << " lnsearchMS: Sufficient Decrease with lambda = " << lambda << '\n';
            return true;
        } else {
            //Backtrack
            if (lambda == 1.0) {
                //First Time
                tmplam = -slope/(2.0*(fscalar-fold-slope));
            } else {
                //Subsequent backtracks
                double rhs1 = fscalar - fold - lambda*slope;
                double rhs2 = fscalar2 - fold - lambda2*slope;
                double a = (rhs1/(lambda*lambda)-rhs2/(lambda2*lambda2))/(lambda-lambda2);
                double b = (-lambda2*rhs1/(lambda*lambda)+lambda*rhs2/(lambda2*lambda2))/(lambda-lambda2);
                if (a==0.0) {
                    tmplam = -slope/(2.0*b);
                } else {
                    double disc = b*b-3.0*a*slope;
                    if(disc<0.0) {
                        tmplam = 0.5*lambda;
                    } else if (b <= 0.0) {
                        tmplam = (-b+sqrt(disc))/(3.0*a);
                    } else {
                        tmplam = -slope/(b+sqrt(disc));
                    }
                }
                //Implied maximum change
                if (tmplam > 0.5*lambda) {
                    tmplam = 0.5*lambda;
                }
            }
        }
        //Update loop variables
        lambda2 = lambda;
        fscalar2 = fscalar;
        iter++;
        //Implied minimum change
        lambda = (tmplam > 0.1*lambda ? tmplam : 0.1*lambda);
    }
}

template<typename RHS, int N>
bool newton(const RHS& rhs, nvis::fixed_vector<double, N>& x,
            double eps, size_t maxiter, double maxlength, bool verbose = false)
{
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::fixed_matrix<double, N>   mat_type;
    
    vec_type d, f; // d is destination, f is rhs
    mat_type J;
    
    std::ostringstream os;
    
    if (verbose) {
        os << "entering Newton at " << x << std::endl;
    }
    
    vec_type best;
    double minnorm = std::numeric_limits<double>::max();
    bool check_improve = false;
    int nb_failed = 0;
    
    unsigned int k;
    try {
        f = rhs(x);
        double dinit = nvis::norm(f);
        for (k = 0; k < maxiter; k++) {
            if (verbose) {
                os <<"newton: k = " << k << ", norm(f) = " << norm(f) << " (eps=" << eps << ") at " << x << '\n';
            }
            if (record_newton_steps) {
                newton_steps.push_back(x);
            }
            
            double _norm = norm(f);
            if (_norm < eps) {
                check_improve = true;
            }
            if (_norm < minnorm) {
                minnorm = _norm;
                best = x;
                if (check_improve) {
                    nb_failed = 0;
                }
            } else if (check_improve) {
                ++nb_failed;
                if (nb_failed > 5) {
                    break;
                }
            }
            
            // determine local search direction
            J = rhs.jacobian(x);
            d = nvis::solve(J, J * x - f);
            
            // do a relaxation linesearch
            // (updates x and f)
            lnsearch<RHS, N>(rhs, x, f, d - x, maxlength);
        }
        
        if (k == maxiter) {
            if (verbose) {
                double dfin = nvis::norm(rhs(x));
                os << "\t\t initial distance = " << dinit
                   << ", final distance = " << dfin << ". failed.\n";
            }
        }
    } catch (std::runtime_error& e) {
        if (verbose) {
            os << e.what() << std::endl;
            os << "\texception caught in Newton. current position is " << x << std::endl;
            std::cerr << os.str() << std::flush;
        }
        return false;
    }
    
    if (k == maxiter) {
        x = best;
    }
    std::cerr << os.str() << std::flush;
    return (minnorm < eps);
}

// variant of Newton's method in which search is confined to a prescribed bounding box
template<typename RHS, int N>
bool newton_in_box(const RHS& rhs, nvis::fixed_vector<double, N>& x,
                   const nvis::bounding_box<nvis::fixed_vector<double, N> >& box,
                   double eps, size_t maxiter, bool verbose = false)
{
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::fixed_matrix<double, N>   mat_type;
    
    vec_type d, f; // d is destination, f is rhs
    mat_type J;
    box_constraint<N> constraint(box);
    double maxlength = 0.5*constraint.size();
    
    std::ostringstream os;
    
    if (verbose) {
        os << "entering Newton in a box at " << x << std::endl;
    }
    vec_type best;
    double minnorm = std::numeric_limits<double>::max();
    bool must_improve = false;
    int nb_failed = 0;
    
    unsigned int k;
    try {
        f = rhs(x);
        double dinit = nvis::norm(f);
        minnorm = dinit;
        for (k = 0; k < maxiter; k++) {
            if (verbose) {
                os << "newton: k = " << k << ", norm(f) = "  << norm(f) << " (eps=" << eps << ") at " << x << std::endl;
            }
            if (record_newton_steps) {
                newton_steps.push_back(x);
            }
            double _norm = norm(f);
            if (_norm < eps) {
                must_improve = true;
            }
            if (_norm < minnorm) {
                minnorm = _norm;
                best = x;
                if (must_improve) {
                    nb_failed = 0;
                }
            } else if (must_improve) {
                ++nb_failed;
                if (nb_failed > 5) {
                    break;
                }
            }
            
            // determine local search direction
            J = rhs.jacobian(x);
            if (verbose) {
                os << "   rhs.jacobian(x) = " << J << std::endl;
            }
            d = nvis::solve(J, J * x - f);
            if (verbose) {
                os << "   Destination (d) = (" << d << ")\n";
            }
            vec_type p = nvis::solve(J,-1.0*f);
            if (verbose) {
                os << "   Step Direction (p) = (" << p << ") or (" << p/norm(p) << ")\n";
            }
            
            // do a relaxation linesearch
            // (updates x and f)
            vec_type save(x);
            //lnsearch<RHS, N>(rhs, x, f, d - x, maxlength);
            bool cvg = lnsearchMinStrategy<RHS, N>(rhs, x, f, p, maxlength);
            if (verbose) {
                os << "   lnsearchMinStrategy Result: converge = " << cvg << ", x = (" << x << "), f = (" << f << "), maxlength = " << maxlength << std::endl;
            }
            if (verbose) {
                os << "     Last vs Previous: xold = (" << save << "), xnew = (" << x << ")\n";
            }
            x = constraint(save, x);
            if (verbose) {
                os << "     Kept: x = (" << x << ")\n";
            }
        }
        
        if (k == maxiter) {
            if (verbose) {
                os << "\t\t initial distance = " << dinit
                   << ", final distance = " << minnorm << ". failed.\n";
            }
        }
    } catch (std::runtime_error& e) {
        if (verbose) {
            os << e.what() << std::endl;
            os << "\texception caught in Newton. current position is " << x << std::endl;
            std::cerr << os.str() << std::flush;
        }
        return false;
    }
    
    if (k == maxiter) {
        x = best;
    }
    std::cerr << os.str() << std::flush;
    return (minnorm < eps);
}

}

#endif

