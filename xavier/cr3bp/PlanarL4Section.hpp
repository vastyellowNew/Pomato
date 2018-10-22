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


#ifndef __CR3BP_PLANARL4_SECTION_HPP__
#define __CR3BP_PLANARL4_SECTION_HPP__

#include <maps/section.hpp>
#define YL4 0.86602540378444
namespace orbital {

template<typename RHS, int N1, int N2>
class planar_section_L4 : public section<double, N2, 2> {
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
    typedef planar_section_L4<rhs_type, N1, N2> self_type;
    
    bool verbose;
    lvec_type seed;
    ///Option to indicate ydot direction
    bool isPositive;
    
    planar_section_L4(const rhs_type& rhs) :
        verbose(false), isPositive(true), _rhs(rhs), _counter(0), _dist(0) {}
    planar_section_L4(const self_type& other) :
        verbose(false), isPositive(other.isPositive), _rhs(other._rhs), _counter(0), _dist(0) {}
        
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
    
    
    std::pair<lvec_type, lmat_type> project(const gvec_type& x) const
    {
        std::pair<lvec_type, lmat_type> r;
        r.first[0] = x[0];
        r.first[1] = x[space_dim];
        r.second(0,0) = x[state_dim];
        r.second(0,1) = x[state_dim+space_dim];
        r.second(1,0) = x[(space_dim+1)*state_dim];
        r.second(1,1) = x[(space_dim+1)*state_dim+space_dim];
        return r;
    }
    
    gvec_type unproject(const lvec_type& x) const
    {
        gvec_type gv(0);
        gv[0] = x[0];
        gv[1] = YL4; //Non-zero y
        gv[space_dim  ] = x[1];
        gv[space_dim+1] = get_yd(x[0], x[1]);
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
        for(int i=0; i<state_dim; i++) {
            int stmID = state_dim + i; //First row = dx/dX
            localSTM(0,i) = x[stmID];
            stmID = state_dim + 3*state_dim + i; //Third row = dxdot/dX
            localSTM(1,i) = x[stmID];
        }
        return localSTM;
    }
    
    /// Section is y=y_L4
    value_type distance(const gvec_type& x) const
    {
        return (isPositive)? (x[1]-YL4) : (YL4-x[1]);
    }
    
    /// Transverse velocity component (perpendicular to section)
    value_type transverseVelocity( const lvec_type& x) const
    {
        //lvec_type sx = project(x).first;
        //return get_yd(sx[0],sx[1]);
        return get_yd(x[0],x[1]);
    }
    
    gvec_type distance_first_partials(const gvec_type& x) const
    {
        gvec_type gv(0);
        gv[1] = 1;
        return gv;
    }
    
    state_type mapToState(const lvec_type& x) const
    {
        state_type v(0);
        v[0] = x[0];
        v[1] = YL4; //Non-zero y
        v[space_dim] = x[1];
        v[space_dim+1] = get_yd(x[0], x[1]);
        return v;
    }
    
    ///Mirror theorem - doesn't apply!
    lvec_type mirror(const lvec_type& x) const
    {
        return lvec_type(x[0],-x[1]);
    }
    
    /// Symmetry test
    bool isSymmetric() const
    {
        return false;
    }
    
    /** Compute the ydot component on the section y=0 given (x,xdot)
         *  Note:  Throws a std::runtime_error if x,xd values lead to an
         *  imaginary yd value.
         */
    double get_yd(double x, double xd) const
    {
        //z=0 on section
        double y = YL4;
        double mu = _rhs.getMu();
        double d  = sqrt((x + mu) * (x + mu) + y*y );
        double r  = sqrt((x - 1 + mu) * (x - 1 + mu) +y*y );
        
        value_type Ustar, ydsqr;
        Ustar = (1.0-mu)/d+mu/r + (x*x+y*y)/2.0;
        ydsqr = 2.0*Ustar - _rhs.getJacobi() - xd*xd;
        
        if (ydsqr < 0) {
            throw std::runtime_error("cr3bp: imaginary velocity");
        }
        return ((isPositive)? sqrt(ydsqr) : -sqrt(ydsqr));
    }
    
private:
    lbox_type           _bounds;
    rhs_type            _rhs;
    
    mutable size_t      _counter;
    mutable double      _dist;
    mutable gvec_type   _last;
};


}

#endif
