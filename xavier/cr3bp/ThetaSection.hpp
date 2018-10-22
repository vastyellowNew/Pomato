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


#ifndef __CR3BP_THETA_SECTION_HPP__
#define __CR3BP_THETA_SECTION_HPP__

#include <maps/section.hpp>
#include <cstdio>

namespace orbital {

/// Planar Hyperplane defined by an angle of rotation from the x-axis in rotating frame
template<typename RHS, int N1, int N2>
class ThetaSection : public section<double, N2, 2> {
public:
    typedef RHS                              rhs_type;
    typedef typename rhs_type::xstate_type   xstate_type;
    typedef typename rhs_type::state_type    state_type;
    
    static const int state_dim =  N1;
    static const int space_dim =  N1/2;
    
    typedef section<double, N2, 2>           base_type;
    typedef typename base_type::lvec_type    lvec_type;
    typedef typename base_type::lmat_type    lmat_type;
    typedef typename base_type::lgMatrix     lgMatrix;
    typedef typename base_type::gvec_type    gvec_type;
    typedef typename base_type::lbox_type    lbox_type;
    typedef typename base_type::value_type   value_type;
    typedef ThetaSection<rhs_type, N1, N2> self_type;
    
    bool verbose;
    /// The angle defining the hyperplane wrt x-axis
    double theta0;
    ///Option to indicate r*thetaDot direction
    bool isPositive;
    
    ThetaSection(const rhs_type& rhs) :
        verbose(false), theta0(0.0), isPositive(true), _rhs(rhs), _ct(1.0), _st(0.0), _counter(0), _dist(0)
    {}
    ThetaSection(const self_type& other) :
        verbose(false), theta0(other.theta0), isPositive(other.isPositive), _rhs(other._rhs),
        _ct(cos(theta0)), _st(sin(theta0)), _counter(0), _dist(0)
    {}
    
    const lbox_type& bounds() const
    {
        return _bounds;
    }
    
    lbox_type& bounds()
    {
        return _bounds;
    }
    
    void set_rhs(const rhs_type& rhs)
    {
        _rhs = rhs;
    }
    
    /// Set the angle defining the hyperplane (rad)
    void setTheta(const double& angle)
    {
        theta0 = angle;
        _ct = cos(theta0);
        _st = sin(theta0);
    }
    /// Set the angle defining the hyperplane (deg)
    void setThetaDeg(const double& angleDeg)
    {
        theta0 = angleDeg*M_PI/180.0;
        _ct = cos(theta0);
        _st = sin(theta0);
    }
    /// Get theta in radians
    double getTheta() const
    {
        return theta0;
    }
    /// Get theta in degrees
    double getThetaDeg() const
    {
        return theta0*180.0/M_PI;
    }
    
    
    std::pair<lvec_type, lmat_type> project(const gvec_type& x) const
    {
        std::pair<lvec_type, lmat_type> rout;
        //Section coords are r,rdot
        value_type r = sqrt(x[0]*x[0]+x[1]*x[1]);
        value_type rdot = x[space_dim]*_ct + x[space_dim+1]*_st;
        rout.first[0] = r;
        rout.first[1] = rdot;
        //These are NOT correct, but I don't want to gather the derivatives from state variables:
        //May be best to compute these numerically
        rout.second(0,0) = x[state_dim];                            //dr/dr0
        rout.second(0,1) = x[state_dim+space_dim];                  //dr/drdot0
        rout.second(1,0) = x[(space_dim+1)*state_dim];              //drdot/dr0
        rout.second(1,1) = x[(space_dim+1)*state_dim+space_dim];    //drdot/drdot0
        return rout;
    }
    
    /// Get the local projection of global STM matrix
    lgMatrix localProjection(const gvec_type& x) const
    {
        lgMatrix localSTM(0);
        //These are NOT corrext:
        for(int i=0; i<state_dim; i++) {
            int stmID = state_dim + i;          //First row = dr/dX(r,theta)
            localSTM(0,i) = x[stmID];
            stmID = state_dim + 3*state_dim + i; //Third row = drdot/dX(r,theta)
            localSTM(1,i) = x[stmID];
        }
        return localSTM;
    }
    
    gvec_type unproject(const lvec_type& x) const
    {
        gvec_type gv(0);
        gv[0] = x[0]*_ct;
        gv[1] = x[0]*_st;
        value_type vt = get_Vt(x[0],x[1]);
        gv[space_dim  ] = x[1]*_ct - vt*_st;
        gv[space_dim+1] = x[1]*_st + vt*_ct;
        // initialize flow map Jacobian to Identity
        for (int i=0 ; i<state_dim ; ++i) {
            gv[state_dim*(i+1)+i] = 1;
        }
        return gv;
    }
    
    value_type distance(const gvec_type& x) const
    {
        value_type dist = x[1]*_ct - x[0]*_st;
        return (isPositive)? dist : -dist;
    }
    
    /// Transverse velocity component (perpendicular to section)
    value_type transverseVelocity( const lvec_type& x) const
    {
        //lvec_type sx = project(x).first;
        //return get_yd(sx[0],sx[1]);
        return get_Vt(x[0],x[1]);
    }
    
    gvec_type distance_first_partials(const gvec_type& x) const
    {
        ///Gradient (cartesian) of hyperplane equation (theta0-theta_desired)
        gvec_type gv(0);
        gv[0] = -_st;
        gv[1] = _ct;
        return gv;
    }
    
    state_type mapToState(const lvec_type& x) const
    {
        state_type v(0);
        v[0] = x[0]*_ct;
        v[1] = x[0]*_st;
        value_type vt = get_Vt(x[0],x[1]);
        v[space_dim  ] = x[1]*_ct - vt*_st;
        v[space_dim+1] = x[1]*_st + vt*_ct;
        return v;
    }
    
    ///Mirror theorem
    lvec_type mirror(const lvec_type& x) const
    {
        //Fix!!!
        return lvec_type(x[0],-x[1]);
    }
    
    /// Symmetry test
    bool isSymmetric() const
    {
        return false;
    }
    
    /** Compute the r*thetaDot component on the section given (r,rdot)
         *  Note:  Throws a std::runtime_error if x,xd values lead to an
         *  imaginary value.
         */
    double get_Vt(double r, double rdot) const
    {
        //z=0 on section
        double mu = _rhs.getMu();
        double x = r*_ct;
        double y = r*_st;
        double d  = sqrt((x + mu) * (x + mu) + y*y);
        double r1  = sqrt((x - 1 + mu) * (x - 1 + mu) + y*y);
        
        value_type Ustar, vtsqr;
        Ustar = (1.0-mu)/d + mu/r1 + (x*x+y*y)/2.0;
        vtsqr = 2.0*Ustar - _rhs.getJacobi() - rdot*rdot;
        
        if (vtsqr < 0) {
            throw std::runtime_error("ThetaSection (CR3BP): imaginary velocity for Vt");
        }
        return ((isPositive)? sqrt(vtsqr) : -sqrt(vtsqr));
    }
    
private:
    lbox_type           _bounds;
    rhs_type            _rhs;
    double              _ct,_st;
    
    mutable size_t      _counter;
    mutable double      _dist;
    mutable gvec_type   _last;
};



}

#endif
