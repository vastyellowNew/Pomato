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


#ifndef __SECTION_HPP__
#define __SECTION_HPP__

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <map>

namespace orbital {

template<typename T, size_t N, size_t M>
class section {
public:
    typedef T                                  value_type;
    typedef nvis::fixed_vector<T, N>           gvec_type;   // global coordinates
    typedef nvis::fixed_matrix<T, N, N>        gmat_type;   // global coordinates
    typedef nvis::fixed_vector<T, M>           lvec_type;   // local coordinates
    typedef nvis::fixed_matrix<T, M, M>        lmat_type;   // local coordinates
    typedef nvis::fixed_matrix<T, M, N>        lgMatrix;    // local projection of global differentials
    typedef nvis::bounding_box<lvec_type>      lbox_type;   // local coordinates
    
    static const int global_dimension = N;
    static const int local_dimension = M;
    
    /// Project global dimension state to a section (local) state and matrix (if STM is available)
    virtual std::pair<lvec_type, lmat_type> project(const gvec_type& x) const = 0;
    
    /// Transform local state (on section) to global state
    virtual gvec_type  unproject(const lvec_type& x) const = 0;
    
    /// Get the local projection of global STM matrix
    virtual lgMatrix localProjection(const gvec_type& x) const = 0;
    
    /// Distance from section
    virtual value_type distance(const gvec_type& x) const = 0;
    
    /// Transverse velocity component (perpendicular to section)
    virtual value_type transverseVelocity( const lvec_type& x) const = 0;
    
    /// First-order spatial derivative of distance from section
    virtual gvec_type  distance_first_partials(const gvec_type& x) const = 0;
    
    /// Application of Mirror Theorem on the section (i.e., if -t, what happens to state)
    virtual lvec_type  mirror(const lvec_type& x) const = 0;
    /// Symmetry test:  "return true" if mirror() is to be used
    virtual bool isSymmetric() const = 0;
    
    /// Section bounds
    const lbox_type& local_bounds() const
    {
        return _lbox;
    }
    /// Section bounds
    lbox_type& local_bounds()
    {
        return _lbox;
    }
    
    
protected:
    lbox_type    _lbox;
};


// Other 'section' type objects that could be useful (need work to conform to standards)

/// Frame class : for analysis within a rectangular section
template<typename T, int N, int M>
class frame {
public:
    typedef nvis::fixed_vector<T, N>     gvec_type;   // global coordinates
    typedef nvis::fixed_matrix<T, N, N>  gmat_type;   // global coordinates
    typedef nvis::fixed_vector<T, M>     lvec_type;   // local coordinates
    typedef nvis::fixed_matrix<T, M, M>  lmat_type;   // local coordinates
    typedef frame<T, N, M>               self_type;
    
    
    frame(const std::vector<gvec_type>& basis, const gvec_type& origin)
        : _basis(basis.begin(), basis.end()), _origin(origin)
    {
        assert(_basis.size() == M);
    }
    
    
    frame(const gvec_type basis[], const gvec_type& origin) : _basis(M), _origin(origin)
    {
        for (int i=0 ; i<M ; ++i) {
            _basis[i] = basis[i];
        }
    }
    
    
    frame(const self_type& f)
        : _basis(f.basis.begin(), f.basis.end()), _origin(f._origin) {}
        
        
    lvec_type project(const gvec_type& x) const
    {
        gvec_type y = x-_origin;
        lvec_type r;
        for (int i=0 ; i<M ; ++i) {
            r[i] = nvis::inner(y, _basis[i]);
        }
        return r;
    }
    
    
    lmat_type project(const gmat_type& m) const
    {
        lmat_type r(0);
        for (int i=0 ; i<M ; ++i) {
            gvec_type v = m*_basis[i];
            for (int j=0 ; j<M ; ++j) {
                r(i,j) = nvis::inner(v, _basis[j]);
            }
        }
        return r;
    }
    
    
    gvec_type unproject(const lvec_type& x) const
    {
        gvec_type y = _origin;
        for (int i=0 ; i<M ; ++i) {
            y += x[i] * _basis[i];
        }
        return y;
    }
    
    
    gvec_type origin() const
    {
        return _origin;
    }
    
    
private:
    std::vector<gvec_type>    _basis;
    gvec_type                 _origin;
};


/// Flat Section class : utilizes a 'frame' object for ease of use
template<typename T, int N, int M>
class flat_section : public section<T, N, M> {
public:
    typedef T value_type;
    typedef frame<T, N, M>           frame_type;
    typedef flat_section<T, N, M>     self_type;
    typedef typename frame_type::lvec_type lvec_type;
    typedef typename frame_type::gvec_type gvec_type;
    typedef typename frame_type::lmat_type lmat_type;
    typedef typename frame_type::gmat_type gmat_type;
    
    
    flat_section(const frame_type& frame)
        : _frame(frame) {}
        
        
    lvec_type project(const gvec_type& x) const
    {
        return _frame.project(x);
    }
    
    
    lmat_type project(const gmat_type& m) const
    {
        return _frame.project(m);
    }
    
    
    gvec_type unproject(const lvec_type& x) const
    {
        return _frame.unproject(x);
    }
    
    
    // this distance is not (cannot be) signed!
    value_type distance(const gvec_type& x) const
    {
        return nvis::norm(x-unproject(project(x)));
    }
    
    
    // this 1D distance is signed
    value_type distance(const gvec_type& x, const gvec_type& normal) const
    {
        return nvis::inner(x-unproject(project(x)), normal);
    }
    
    
    self_type* clone() const
    {
        return new self_type(_frame);
    }
    
    
private:
    frame_type _frame;
};


/// Hyperplane Section:  N-1 dimensional hyperplane
template<typename T, int N>
class hyperplanar_section : public section<T, N, N-1> {
public:
    typedef T value_type;
    typedef frame<T, N, N-1>           plane_type;
    typedef hyperplanar_section<T, N>   self_type;
    typedef typename plane_type::lvec_type lvec_type;
    typedef typename plane_type::gvec_type gvec_type;
    typedef typename plane_type::lmat_type lmat_type;
    typedef typename plane_type::gmat_type gmat_type;
    
    
    hyperplanar_section(const plane_type& plane, const gvec_type& normal)
        : _plane(plane), _normal(normal) {}
        
        
    // this distance is signed
    value_type distance(const gvec_type& x) const
    {
        return nvis::inner(_normal, x-_plane.origin());
    }
    
    
    lvec_type project(const gvec_type& x) const
    {
        return _plane.project(x);
    }
    
    
    lmat_type project(const gmat_type& m) const
    {
        return _plane.project(m);
    }
    
    
    gvec_type unproject(const lvec_type& x) const
    {
        return _plane.unproject(x);
    }
    
    
    self_type* clone() const
    {
        return new self_type(_plane, _normal);
    }
    
    
private:
    plane_type _plane;
    gvec_type  _normal;
};





}

#endif
