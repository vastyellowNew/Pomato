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


#ifndef __pcr3bp_hpp
#define __pcr3bp_hpp

#include <cmath>
#include <vector>
#include <complex>
#include <stdexcept>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/bezier_spline.hpp>
#include <maps/rhs_base.hpp>
#include <math/angle.hpp>

namespace orbital {

/*
    Earth-Moon:
        C = 2.96, mu = 1.21505714306e-2
        suggested bounds: -0.3 -0.5 1.2 0.5

    Jupiter-Europa:
        C = 3.000, mu = 2.528017705e-5
        suggested bounds: -1.5 -0.2 -1 0.2
*/

class pcr3bp: public rhs_base<double, 4> {
public:
    typedef rhs_base<double, 4> base_type;
    
    // state variables
    // 0: x, 1: y, 2: dxdt, 3: dydt
    // derivative variables
    // coefficients of 4x4 spatial derivative, in row-wise order
    
    typedef nvis::fixed_vector<double, 5>   vec5;
    
    double yd(double x, double xd) const
    {
    
        std::complex<double> _x(x, 0), _xd(xd, 0);
        std::complex<double> _y = sqrt(2.*(1. - _mu) / sqrt(pow(_x + _mu, 2)) +
                                       2.*_mu / sqrt(pow(_x - 1. + _mu, 2)) + pow(_x, 2) - pow(_xd, 2) - _C);
                                       
        if (_y.imag()) {
            throw std::runtime_error("pcr3bp: imaginary velocity");
        }
        return _y.real();
    }
    
    pcr3bp(double C, double mu) : _bbox(state_type(m_inf), state_type(p_inf)), _C(C), _mu(mu)
    {}
    
    pcr3bp(double C, double mu, const box_type& bbox) : _bbox(bbox), _C(C), _mu(mu)  {}
    
    xstate_type operator()(const double& t, const xstate_type& y) const
    {
        xstate_type r;
        // dx/dt = v
        nvis::subv<0, 4, double, 20>(r) = interpolate(nvis::subv<0, 4, double, 20>(y));
        // dJ/dt = dv/dx * J
        mat_type m = Jacobi(nvis::subv<0, 4, double, 20>(y));
        mat_type J;
        J[0] = nvis::subv< 4, 4, double, 20>(y);
        J[1] = nvis::subv< 8, 4, double, 20>(y);
        J[2] = nvis::subv<12, 4, double, 20>(y);
        J[3] = nvis::subv<16, 4, double, 20>(y);
        
        m = m*J;
        nvis::subv< 4, 4, double, 20>(r) = m[0];
        nvis::subv< 8, 4, double, 20>(r) = m[1];
        nvis::subv<12, 4, double, 20>(r) = m[2];
        nvis::subv<16, 4, double, 20>(r) = m[3];
        
        return r;
    }
    
    state_type value(const double& t, const state_type& y) const
    {
        return interpolate(y);
    }
    
    mat_type derivative(const double& t, const state_type& y) const
    {
        return Jacobi(y);
    }
    
    value_type hamiltonian(const double& t, const xstate_type& y) const
    {
        //The Jacobi Constant
        return JacobiConstant(y);
    }
    
    value_type desired_hamiltonian() const
    {
        return _C;
    }
    
    state_type hamiltonian_first_partials(const double& t, const xstate_type& y) const
    {
        return JacobiPartials(y);
    }
    
    const box_type& bounds() const
    {
        return _bbox;
    }
    
    base_type* clone() const
    {
        return new pcr3bp(_C, _mu, _bbox);
    }
    
private:

    state_type interpolate(const state_type& p) const
    {
    
        const double& x =   p[0];
        const double& y =   p[1];
        const double& xd =  p[2];
        const double& yd =  p[3];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1 + _mu) * (x - 1 + _mu) + y * y);
        double d3   = d*d*d;
        double r3   = r*r*r;
        
        state_type v;
        v[0]    =   xd;
        v[1]    =   yd;
        v[2]    =   2 * yd  + x - (1 - _mu) * (x + _mu) / d3    - _mu * (x - 1 + _mu) / r3;
        v[3]    =  -2 * xd + y - (1 - _mu) * y / d3             - _mu * y / r3;
        
        return v;
    }
    
    /// Linearized equations of motion "A" matrix - Also a spatial derivative of a state (or Jacobian)
    mat_type Jacobi(const state_type& p) const
    {
        const double& x =   p[0];
        const double& y =   p[1];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1.0 + _mu) * (x - 1.0 + _mu) + y * y);
        double d3   = d*d*d;
        double d5   = d3*d*d;
        double r3   = r*r*r;
        double r5   = r3*r*r;
        
        //Partials
        mat_type J(0); //Initialize
        
        double Uxx = 1.0 -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*(x+_mu)*(x+_mu)/d5+3.0*_mu*(x-1.0+_mu)*(x-1.0+_mu)/r5;
        double Uyy = 1.0 -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*y*y/d5+3.0*_mu*y*y/r5;
        double Uxy = 3.0*(1.0-_mu)*(x+_mu)*y/d5+3.0*_mu*(x-1.0+_mu)*y/r5;
        
        J(0,2) = J(1,3) = 1.0;
        J(2,3) = 2.0;
        J(3,2) = -2.0;
        J(2,0) = Uxx;
        J(2,1) = Uxy;
        J(3,0) = Uxy;
        J(3,1) = Uyy;
        
        return J;
    }
    /// Jacobi Constant function
    value_type JacobiConstant(const xstate_type& p) const
    {
    
        const double& x =   p[0];
        const double& y =   p[1];
        const double& xd =  p[2];
        const double& yd =  p[3];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1 + _mu) * (x - 1 + _mu) + y * y);
        
        value_type Ustar, JC;
        Ustar = (1.0-_mu)/d+_mu/r + (x*x+y*y)/2.0;
        JC = 2.0*Ustar - (xd*xd+yd*yd);
        
        return JC;
    }
    /// Jacobi Constant first-order partials (gradient)
    state_type JacobiPartials(const xstate_type& p) const
    {
        const double& x =   p[0];
        const double& y =   p[1];
        const double& xd =  p[2];
        const double& yd =  p[3];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1 + _mu) * (x - 1 + _mu) + y * y);
        double d3   = d*d*d;
        double r3   = r*r*r;
        
        state_type partials;
        
        partials[0]    =  2.0*(x - (1 - _mu) * (x + _mu) / d3   - _mu * (x - 1 + _mu) / r3);
        partials[1]    =  2.0*(y - (1 - _mu) * y / d3 - _mu * y / r3);
        partials[2]    = -2.0*xd;
        partials[3]    = -2.0*yd;
        return partials;
    }
    
    nvis::bounding_box<nvis::fixed_vector<double, 4> > _bbox;
    double _C, _mu;
};

class pcr3bp_reduced {
public:
    typedef double                          value_type;
    typedef nvis::fixed_vector<double, 4>   state_type;
    typedef nvis::fixed_vector<double, 4>   xstate_type; // fake extended state
    typedef nvis::fixed_matrix<double, 4>   mat_type;
    typedef nvis::bounding_box<state_type>  box_type;
    static const int dimension = 4;
    
public:
    // state variables
    // 0: x, 1: y, 2: dxdt, 3: dydt
    // derivative variables
    // coefficients of 4x4 spatial derivative, in row-wise order
    
    typedef nvis::fixed_vector<double, 5>   vec5;
    
    double yd(double x, double xd) const
    {
    
        std::complex<double> _x(x, 0), _xd(xd, 0);
        std::complex<double> _y = sqrt(2.*(1. - _mu) / sqrt(pow(_x + _mu, 2)) +
                                       2.*_mu / sqrt(pow(_x - 1. + _mu, 2)) + pow(_x, 2) - pow(_xd, 2) - _C);
                                       
        if (_y.imag()) {
            throw std::runtime_error("pcr3bp: imaginary velocity");
        }
        return _y.real();
    }
    
    pcr3bp_reduced(double C, double mu) : _bbox(state_type(m_inf), state_type(p_inf)), _C(C), _mu(mu)
    {}
    
    pcr3bp_reduced(double C, double mu, const box_type& bbox) : _bbox(bbox), _C(C), _mu(mu)  {}
    
    xstate_type operator()(const double& t, const xstate_type& y) const
    {
        return interpolate(y);
    }
    
    state_type value(const double& t, const state_type& y) const
    {
        return interpolate(y);
    }
    
    mat_type derivative(const double& t, const state_type& y) const
    {
        return Jacobi(y);
    }
    
    value_type hamiltonian(const double& t, const xstate_type& y) const
    {
        //The Jacobi Constant
        return JacobiConstant(y);
    }
    
    value_type desired_hamiltonian() const
    {
        return _C;
    }
    
    state_type hamiltonian_first_partials(const double& t, const xstate_type& y) const
    {
        return JacobiPartials(y);
    }
    
    const box_type& bounds() const
    {
        return _bbox;
    }
    
    pcr3bp_reduced* clone() const
    {
        return new pcr3bp_reduced(_C, _mu, _bbox);
    }
    
private:

    state_type interpolate(const state_type& p) const
    {
    
        const double& x =   p[0];
        const double& y =   p[1];
        const double& xd =  p[2];
        const double& yd =  p[3];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1 + _mu) * (x - 1 + _mu) + y * y);
        double d3   = d*d*d;
        double r3   = r*r*r;
        
        state_type v;
        v[0]    =   xd;
        v[1]    =   yd;
        v[2]    =   2 * yd  + x - (1 - _mu) * (x + _mu) / d3    - _mu * (x - 1 + _mu) / r3;
        v[3]    =  -2 * xd + y - (1 - _mu) * y / d3             - _mu * y / r3;
        
        return v;
    }
    
    /// Linearized equations of motion "A" matrix - Also a spatial derivative of a state (or Jacobian)
    mat_type Jacobi(const state_type& p) const
    {
        const double& x =   p[0];
        const double& y =   p[1];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1.0 + _mu) * (x - 1.0 + _mu) + y * y);
        double d3   = d*d*d;
        double d5   = d3*d*d;
        double r3   = r*r*r;
        double r5   = r3*r*r;
        
        //Partials
        mat_type J(0); //Initialize
        
        double Uxx = 1.0 -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*(x+_mu)*(x+_mu)/d5+3.0*_mu*(x-1.0+_mu)*(x-1.0+_mu)/r5;
        double Uyy = 1.0 -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*y*y/d5+3.0*_mu*y*y/r5;
        double Uxy = 3.0*(1.0-_mu)*(x+_mu)*y/d5+3.0*_mu*(x-1.0+_mu)*y/r5;
        
        J(0,2) = J(1,3) = 1.0;
        J(2,3) = 2.0;
        J(3,2) = -2.0;
        J(2,0) = Uxx;
        J(2,1) = Uxy;
        J(3,0) = Uxy;
        J(3,1) = Uyy;
        
        return J;
    }
    /// Jacobi Constant function
    value_type JacobiConstant(const xstate_type& p) const
    {
    
        const double& x =   p[0];
        const double& y =   p[1];
        const double& xd =  p[2];
        const double& yd =  p[3];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1 + _mu) * (x - 1 + _mu) + y * y);
        
        value_type Ustar, JC;
        Ustar = (1.0-_mu)/d+_mu/r + (x*x+y*y)/2.0;
        JC = 2.0*Ustar - (xd*xd+yd*yd);
        
        return JC;
    }
    
    /// Jacobi Constant first-order partials (gradient)
    state_type JacobiPartials(const xstate_type& p) const
    {
        const double& x =   p[0];
        const double& y =   p[1];
        const double& xd =  p[2];
        const double& yd =  p[3];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y);
        double r    = sqrt((x - 1 + _mu) * (x - 1 + _mu) + y * y);
        double d3   = d*d*d;
        double r3   = r*r*r;
        
        state_type partials;
        
        partials[0]    =  2.0*(x - (1 - _mu) * (x + _mu) / d3   - _mu * (x - 1 + _mu) / r3);
        partials[1]    =  2.0*(y - (1 - _mu) * y / d3 - _mu * y / r3);
        partials[2]    = -2.0*xd;
        partials[3]    = -2.0*yd;
        return partials;
    }
    
    nvis::bounding_box<nvis::fixed_vector<double, 4> > _bbox;
    double _C, _mu;
};

} // orbital

#endif // __pcr3bp_hpp

























