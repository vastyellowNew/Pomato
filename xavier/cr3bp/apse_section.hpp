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


/** Apse Section (Periapsis or Apoapsis) in CR3BP
 *  Author:  Wayne Schlei (Purdue University)
 *  Notes: Ode to Amanda Haapala (Chalk)
 *         This is for a PLANAR section (in x-y plane of rotating frame)!
 *         Also, it's hard to reproduce propagation in interactive setting (Avizo)
 *         without a direction flag given an (x,y) projection. (Periapsis v. Apoapsis)
 */
#ifndef __CR3BP_APSE_SECTION_HPP__
#define __CR3BP_APSE_SECTION_HPP__

#include <maps/section.hpp>
#include <cmath>

namespace orbital {

/// Section for detecting apsidies in CR3BP
template<typename RHS, int N1, int N2>
class apse_section : public section<double, N2, 2> {
public :
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
    typedef apse_section<rhs_type, N1, N2>   self_type;
    
    bool verbose;
    ///Option on type: Periapsis (true) or Apoapsis (false)
    bool isPeriapsis;
    ///Option on seed direction
    bool seedPrograde;
    
    ///Option for which body to evalute apsidies about
    enum ApseRef {
        BARYCENTER = 0,//Barycenter (0,0)
        P1,            //Primary 1 (-mu,0)
        P2             //Primary 2 (1-mu,0)
    } referenceBody;
    
    apse_section(const rhs_type& rhs) : verbose(false),
        isPeriapsis(true), seedPrograde(true), referenceBody(BARYCENTER),
        _rhs(rhs), _counter(0), _dist(0) {}
    apse_section(const self_type& other) : verbose(false),
        isPeriapsis(other.isPeriapsis), seedPrograde(other.seedPrograde),
        referenceBody(other.referenceBody), _bounds(other._bounds),
        _rhs(other._rhs), _counter(0), _dist(0) {}
        
    const lbox_type& bounds() const
    {
        return _bounds;
    }
    
    lbox_type& bounds()
    {
        return _bounds;
    }
    
    /// Projection of total state leads to (x,y) plot
    std::pair<lvec_type, lmat_type> project(const gvec_type& x) const
    {
        std::pair<lvec_type, lmat_type> r;
        r.first[0] = x[0];
        r.first[1] = x[1];
        //r.first[2] = x[2]; //Assume as zero
        //Map Jacobian
        r.second(0,0) = x[state_dim];
        r.second(0,1) = x[state_dim+1];
        r.second(1,0) = x[2*state_dim];
        r.second(1,1) = x[2*state_dim+1];
        return r;
    }
    
    gvec_type unproject(const lvec_type& x) const
    {
        gvec_type gv(0);
        state_type s = mapToState(x);
        for (int i=0 ; i<state_dim ; ++i) {
            gv[i] = s[i];
        }
        // initialize flow map Jacobian to Identity
        for (int i=0 ; i<state_dim ; ++i) {
            gv[state_dim*(i+1)+i] = 1;
        }
        return gv;
    }
    
    /// Get the local projection of global STM matrix
    lgMatrix localProjection(const gvec_type& x) const
    {
        lgMatrix localSTM(0);
        //These aren't right... ???
        for(int i=0; i<state_dim; i++) {
            int stmID = state_dim + i;           //First row = dx/dX
            localSTM(0,i) = x[stmID];
            stmID = state_dim + 1*state_dim + i; //Third row = dy/dX
            localSTM(1,i) = x[stmID];
        }
        return localSTM;
    }
    
    /// Distance from section ( dot(r,v) = 0 OR rdot = 0 )
    value_type distance(const gvec_type& xx) const
    {
        value_type rdot = 0.0;
        double l = get_l();
        //Planar r,theta
        const double& x = xx[0];
        const double& y = xx[1]; //double &z = xx[2];
        const double& xd = xx[3];
        const double& yd = xx[4]; //double &zd = xx[5];
        rdot = (x-l)*xd + y*yd;
        return (isPeriapsis) ? rdot : -rdot;
    }
    
    /// Transverse velocity component (perpendicular to section)
    value_type transverseVelocity( const lvec_type& xIn) const
    {
        //In apse map, rdot = 0 (with respect to referenceBody)
        //So velocity through section is rdoubledot (acceleration)
        value_type rdd = 0.0;
        value_type l = get_l();
        //State References
        state_type xx = mapToState(xIn);
        double& x = xx[0];
        double& y = xx[1]; //double &z = xx[2];
        double& xd = xx[3];
        double& yd = xx[4]; //double &zd = xx[5];
        //Acceleration
        state_type xxDer = _rhs.value(0,xx); //state derivatives
        double& xdd = xxDer[3];
        double& ydd = xxDer[4]; //double &zdd = xxDer[5];
        
        //Staight from Amanda
        double r = sqrt( (x-l)*(x-l) + y*y );
        double rdot = (x-l)*xd + y*yd;
        rdd = (xd*xd+yd*yd + (x-l)*xdd + y*ydd)/r - (rdot*rdot)/(r*r*r);
        return rdd;
    }
    
    /// First partial derivatives of section equation : d(rdot)/dx_i
    gvec_type distance_first_partials(const gvec_type& xx) const
    {
        gvec_type gv(0);
        value_type l = get_l();
        //State References
        const double& x = xx[0];
        const double& y = xx[1]; //const double &z = xx[2];
        const double& xd = xx[3];
        const double& yd = xx[4]; //const double &zd = xx[5];
        //drdot / dx
        gv[0] = xd;
        //drdot / dy
        gv[1] = yd;
        //drdot / dz
        //gv[2] = ;
        //drdot / dxdot
        gv[3] = x-l;
        //drdot / dydot
        gv[4] = y;
        //drdot / dzdot
        //gv[5] = ;
        return gv;
    }
    
    /// Get the 'x' modifier based on refernceBody
    double get_l() const
    {
        double l = 0.0;
        switch(referenceBody) {
            case P1 :
                l = - _rhs.getMu();
                break;
            case P2 :
                l =  1.0 - _rhs.getMu();
                break;
            default :
                break;
        }
        return l;
    }
    
    /// Convert a state ON the map (2D) to a full 6D state vector (rho)
    state_type mapToState(const lvec_type& xx) const
    {
        state_type sv(0);
        sv[0] = xx[0]; //x
        sv[1] = xx[1]; //y
        //Velocity determined by Jacobi constant value in RHS and referenceBody
        double vMag = get_v(xx[0],xx[1]);
        nvis::vec3 r(xx[0],xx[1],0.0), v(0.0);
        double l = get_l();
        r[0] = xx[0] - l;
        r /= nvis::norm(r);
        
        //Assume CCW for seed, but flip if necessary (seedPrograde)
        //Note: HARD-CODED FOR 2D CASE!
        v[0] = -r[1]*vMag;
        v[1] = r[0]*vMag;
        v[2] = r[2]*vMag;
        
        if (!seedPrograde) {
            //Flip the velocity direction to CW
            v[0] *= -1.0;
            v[1] *= -1.0;
            //v[2] *= -1.0;??
        }
        //Note on Direction:  Should only be applied on seeding initial conditions.
        //Techinally, lvec_type should have a direction flag.  We lose some info
        // converting from an (x,y) map state to (x,y,xdot,ydot)!
        
        //Set based on velocity
        sv[3] = v[0];
        sv[4] = v[1];
        //sv[5] = v[2];
        return sv;
    }
    
    ///Mirror theorem
    lvec_type mirror(const lvec_type& x) const
    {
        return lvec_type(x[0],-x[1]);
    }
    
    /// Symmetry test
    bool isSymmetric() const
    {
        return true;
    }
    
    /** Compute the velocity magnitude on the section dot(r,v)=0 given (x,y)
      *  Note:  Throws a std::runtime_error if x,y values lead to an
      *  imaginary velocity.
      */
    double get_v(double x, double y) const
    {
        //apse section
        double mu = _rhs.getMu();
        double d  = sqrt((x + mu) * (x + mu) + y*y);
        double r  = sqrt((x - 1 + mu) * (x - 1 + mu) + y*y );
        
        value_type Ustar, vsqr;
        Ustar = (1.0-mu)/d+mu/r + (x*x+y*y)/2.0;
        vsqr = 2.0*Ustar - _rhs.getJacobi();
        
        if (vsqr < 0) {
            throw std::runtime_error("cr3bp: imaginary velocity");
        }
        return sqrt(vsqr);
    }
    
    /// Evaluate the dot product given that (r,v) are stored in state_type object
    inline double dot(const state_type& state) const
    {
        double val = 0.0;
        for (int i=0; i<space_dim; i++) {
            val+= state[i]*state[space_dim+i];
        }
        return val;
    }
    
    /// Look up state derivatives given a local section coordinate
    state_type getStateDerivatives(const lvec_type& x) const
    {
        state_type y = mapToState(x);
        return _rhs.value(0.0,y);
    }
private:
    lbox_type           _bounds;
    const rhs_type      _rhs;
    
    mutable size_t      _counter;
    mutable double      _dist;
    
    
};

} //end orbital


#endif //__CR3BP_APSE_SECTION_HPP__
