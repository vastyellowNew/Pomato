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


#ifndef __cr3bp_hpp
#define __cr3bp_hpp

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
#include <maps/mapExceptions.hpp>

namespace {
double p_inf        = std::numeric_limits<double>::max();
double m_inf        = std::numeric_limits<double>::min();
};

namespace orbital {

/*
        Earth-Moon:
                C = 2.96, mu = 1.21505714306e-2
                suggested bounds: -0.4 -2.5 1.1 2.5

        Jupiter-Europa:
                C = 3.000, mu = 2.528017705e-5
                suggested bounds: -1.5 -0.2 -1 0.2
*/


/** Circular Restricted Three-body Problem - (RHS class)
 *    - Throws a singularity error if trying to access a state inside one of the
 *      bodies (which will in turn throw a divide by zero.)
 */
class cr3bp: public rhs_base<double, 6, 2> {
public:

    typedef rhs_base<double, 6, 2> base_type;
    
    // state variables
    // 0: x, 1: y, 2: z, 3: dxdt, 4: dydt, 5: dzdt
    // derivative variables
    // coefficients of 6x6 spatial derivative, in row-wise order
    
    typedef nvis::fixed_vector<double, 7>         vec7;
    
    
    ///Constructor with Jacobi Constant C and mass ratio mu
    cr3bp(const double& C,const double& mu) :
        _bbox(state_type(m_inf),
              state_type(p_inf)), _C(C), _mu(mu),
#if defined(HX_HAS_STD) || defined(C_0X)
        singularityTol(1e-20)
#else
        singularityTol(1e-20), primaries {0.0}, safeDist {0.0}
#endif
    {
        primaries[0] = xstate_type(0.0);
        primaries[1] = xstate_type(0.0);
        primaries[0][0] = -_mu;
        primaries[1][0] = 1.0-_mu;
        safeDist[0] = 0.01;  //Roughly Earth radius in Earth-Moon problem
        safeDist[1] = 0.0045; //Roughly Moon radius in Earth-Moon problem
        singularityTol = 1e10*std::numeric_limits<double>::min();
        
    }
    
    ///Constructor with a bounding box
    cr3bp(const double& C, const double& mu, const box_type& bbox) :
        _bbox(bbox), _C(C), _mu(mu),
#if defined(HX_HAS_STD) || defined(C_0X)
        singularityTol(1e-20)
#else
        singularityTol(1e-20), primaries {0.0}, safeDist {0.0}
#endif
    {
        primaries[0] = xstate_type(0.0);
        primaries[1] = xstate_type(0.0);
        primaries[0][0] = -_mu;
        primaries[1][0] = 1.0-_mu;
        safeDist[0] = 0.01;  //Roughly Earth radius in Earth-Moon problem
        safeDist[1] = 0.0045; //Roughly Moon radius in Earth-Moon problem
        singularityTol = 1e10*std::numeric_limits<double>::min();
    }
    
    //Probably should make a copy constructor!
    
    ///Full-state derivative operator called by the integrator (e.g., nvis::dorpri5)
    xstate_type operator()(const double& t, const xstate_type& y) const
    {
        xstate_type rState(0);
        // dx/dt = v
        nvis::subv<0, 6, double, 42>(rState) = interpolate(nvis::subv<0, 6, double, 42>(y));
        // dJ/dt = dv/dx * J
        mat_type m = Jacobi(nvis::subv<0, 6, double, 42>(y));
        mat_type J(0);
        J[0] = nvis::subv< 6, 6, double, 42>(y);
        J[1] = nvis::subv<12, 6, double, 42>(y);
        J[2] = nvis::subv<18, 6, double, 42>(y);
        J[3] = nvis::subv<24, 6, double, 42>(y);
        J[4] = nvis::subv<30, 6, double, 42>(y);
        J[5] = nvis::subv<36, 6, double, 42>(y);
        
        m = m*J;
        nvis::subv< 6, 6, double, 42>(rState) = m[0];
        nvis::subv<12, 6, double, 42>(rState) = m[1];
        nvis::subv<18, 6, double, 42>(rState) = m[2];
        nvis::subv<24, 6, double, 42>(rState) = m[3];
        nvis::subv<30, 6, double, 42>(rState) = m[4];
        nvis::subv<36, 6, double, 42>(rState) = m[5];
        
        return rState;
    }
    
    ///Derivative operator for just the 6D position/velocity states
    state_type value(const double& t, const state_type& y) const
    {
        return interpolate(y);
    }
    
    ///Return the 6D state derivative (or linear A matrix in xdot=Ax)
    mat_type derivative(const double& t, const state_type& y) const
    {
        return Jacobi(y);
    }
    
    ///Computing the Hamiltonian of the system (CR3BP: JacobiConstant = Hamiltonian)
    value_type hamiltonian(const double& t, const xstate_type& y) const
    {
        //The Jacobi Constant
        return JacobiConstant(y);
    }
    
    ///The desired value of the Hamiltonian to which we are testing
    value_type desired_hamiltonian() const
    {
        return _C;
    }
    
    ///First Derivatives of Hamiltonian with respect to state elements
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
        return new cr3bp(_C, _mu, _bbox);
    }
    
    /// Get the mass parameter value (private access)
    double getMu() const
    {
        return _mu;
    }
    /// Get the Jacobi parameter value (private access)
    double getJacobi() const
    {
        return _C;
    }
    /// Set the Jacobi parameter value (private access)
    void setJacobi(const double value)
    {
        _C = value;
    }
    /// Set the singularity tolerance
    void setSingularityTolerance(const double& tol)
    {
        singularityTol = tol;
    }
    /// Get the singularity tolerance
    double getSingularityTolerance() const
    {
        return singularityTol;
    }
    
    /// Referencing singularities
    const xstate_type* singularities() const
    {
        return primaries;
    }
    /// Number of singularities
    int getNumSingularities() const
    {
        return numSingularities;
    }
    /// Minimum numerically safe distance from singularity (in nd map units sqrt(x^2+xd^2))
    value_type getSingularitySafeDistance(const int& idx) const
    {
        if(idx>=numSingularities) {
            return safeDist[1];
        }
        return safeDist[idx];
    }
    /// Set the numerically safe distance from singularity (in nd map units)
    void setSingularitySafeDistance(const int& idx, const double& r)
    {
        if(idx<numSingularities) {
            safeDist[idx] = r;
        } else {
            throw std::runtime_error("Call to cr3pb.hpp::setSingularitySafeDistance() index invalid!");
        }
    }
    
private:
    ///Evaluate the equations of motion
    state_type interpolate(const state_type& p) const
    {
    
        const double& x  =         p[0];
        const double& y  =         p[1];
        const double& z  =         p[2];
        const double& xd =         p[3];
        const double& yd =         p[4];
        const double& zd =         p[5];
        
        double d        = sqrt((x + _mu) * (x + _mu) + y * y + z * z);
        double r        = sqrt((x - 1.0 + _mu) * (x - 1.0 + _mu) + y * y + z * z);
        double d3       = d*d*d;
        double r3       = r*r*r;
        
        //Throw singularity error if d or r is zero to specified singularity tolerance
        if ( d3 <= singularityTol || r3 <= singularityTol) {
            throw xavier::MapSingularityError<state_type>(p);
        }
        
        state_type v;
        v[0]    =   xd;
        v[1]    =   yd;
        v[2]    =   zd;
        v[3]    =   2.0 * yd  + x - (1.0 - _mu) * (x + _mu) / d3   - _mu * (x - 1.0 + _mu) / r3;
        v[4]    =  -2.0 * xd  + y - (1.0 - _mu) * y / d3           - _mu * y / r3;
        v[5]    =                 - (1.0 - _mu) * z / d3           - _mu * z / r3;
        
        return v;
    }
    
    /// Linearized equations of motion "A" matrix - Also a spatial derivative of a state (or Jacobian)
    mat_type Jacobi(const state_type& p) const
    {
        const double& x =         p[0];
        const double& y =         p[1];
        const double& z =         p[2];
        //const double& xd =         p[3];
        //const double& yd =         p[4];
        //const double& zd =         p[5];
        
        double d    = sqrt((x + _mu) * (x + _mu) + y * y + z * z);
        double r    = sqrt((x - 1.0 + _mu) * (x - 1.0 + _mu) + y * y + z * z);
        double d3   = d*d*d;
        double d5   = d3*d*d;
        double r3   = r*r*r;
        double r5   = r3*r*r;
        
        //Partials
        mat_type J(0); //Initialize
        
        double Uxx = 1.0 -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*(x+_mu)*(x+_mu)/d5+3.0*_mu*(x-1.0+_mu)*(x-1.0+_mu)/r5;
        double Uyy = 1.0 -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*y*y/d5+3.0*_mu*y*y/r5;
        double Uzz = -(1.0-_mu)/d3-_mu/r3+3.0*(1.0-_mu)*z*z/d5+3.0*_mu*z*z/r5;
        double Uxy = 3.0*(1.0-_mu)*(x+_mu)*y/d5+3.0*_mu*(x-1.0+_mu)*y/r5;
        double Uxz = 3.0*(1.0-_mu)*(x+_mu)*z/d5+3.0*_mu*(x-1.0+_mu)*z/r5;
        double Uyz = 3.0*(1.0-_mu)*y*z/d5+3.0*_mu*y*z/r5;
        
        J(0,3) = J(1,4) = J(2,5) = 1.0;
        J(3,4) = 2.0;
        J(4,3) = -2.0;
        J(3,0) = Uxx;
        J(3,1) = Uxy;
        J(3,2) = Uxz;
        J(4,0) = Uxy;
        J(4,1) = Uyy;
        J(4,2) = Uyz;
        J(5,0) = Uxz;
        J(5,1) = Uyz;
        J(5,2) = Uzz;
        
        return J;
    }
    /// Jacobi Constant function
    value_type JacobiConstant(const xstate_type& p) const
    {
    
        const double& x  =         p[0];
        const double& y  =         p[1];
        const double& z  =         p[2];
        const double& xd =         p[3];
        const double& yd =         p[4];
        const double& zd =         p[5];
        
        double d        = sqrt((x + _mu) * (x + _mu) + y * y + z * z);
        double r        = sqrt((x - 1.0 + _mu) * (x - 1.0 + _mu) + y * y + z * z);
        
        
        value_type Ustar, JC;
        Ustar = (1.0-_mu)/d+_mu/r + (x*x+y*y)/2.0;
        JC = 2.0*Ustar - (xd*xd+yd*yd+zd*zd);
        
        return JC;
    }
    /// Jacobi Constant first-order partials (gradient)
    state_type JacobiPartials(const xstate_type& p) const
    {
        const double& x  =         p[0];
        const double& y  =         p[1];
        const double& z  =         p[2];
        const double& xd =         p[3];
        const double& yd =         p[4];
        const double& zd =         p[5];
        
        double d        = sqrt((x + _mu) * (x + _mu) + y * y + z * z);
        double r        = sqrt((x - 1.0 + _mu) * (x - 1.0 + _mu) + y * y + z * z);
        double d3       = d*d*d;
        double r3       = r*r*r;
        
        state_type partials;
        
        partials[0]    =  2.0*(x - (1.0 - _mu) * (x + _mu) / d3   - _mu * (x - 1.0 + _mu) / r3);
        partials[1]    =  2.0*(y - (1.0 - _mu) * y / d3 - _mu * y / r3);
        partials[2]    =  2.0*(-(1.0 - _mu) * z / d3 - _mu * z / r3);
        partials[3]    = -2.0*xd;
        partials[4]    = -2.0*yd;
        partials[5]    = -2.0*zd;
        return partials;
    }
    
    ///Bounding box - holding bounds of analysis
    nvis::bounding_box<nvis::fixed_vector<double, 6> > _bbox;
    ///CR3BP parameters
    double _C, _mu;
    ///Tolerance on singularity
    double singularityTol;
    /// Singularities
    xstate_type primaries[2];
    /// Singularity safe distance
    value_type safeDist[2];
};

} // orbital

#endif // __cr3bp_hpp