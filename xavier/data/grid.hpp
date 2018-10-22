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


#ifndef __GRID_HPP__
#define __GRID_HPP__

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <boost/static_assert.hpp>
#include <stdexcept>
#include <maps/metric.hpp>

namespace xavier {
template<typename T, int N>
struct id_set {
    typedef T                               value_type;
    typedef nvis::fixed_vector<T, N>        vec_type;
    
    id_set() : __ids(-1) {}
    id_set(const value_type ids[N])
    {
        for (int i = 0 ; i < N ; ++i) {
            __ids[i] = ids[i];
        }
        std::sort(__ids.begin(), __ids.end());
    }
    id_set(value_type i0, value_type i1)
    {
        BOOST_STATIC_ASSERT(N == 2);
        if (i0 < i1) {
            __ids[0] = i0;
            __ids[1] = i1;
        } else {
            __ids[0] = i1;
            __ids[1] = i0;
        }
    }
    id_set(value_type i0, value_type i1, value_type i2)
    {
        BOOST_STATIC_ASSERT(N == 3);
        __ids[0] = i0;
        __ids[1] = i1;
        __ids[2] = i2;
        std::sort(__ids.begin(), __ids.end());
    }
    id_set(value_type i0, value_type i1, value_type i2, value_type i3)
    {
        BOOST_STATIC_ASSERT(N == 4);
        __ids[0] = i0;
        __ids[1] = i1;
        __ids[2] = i2;
        __ids[3] = i3;
        std::sort(__ids.begin(), __ids.end());
    }
    id_set(const vec_type& ids) : __ids(ids)
    {
        std::sort(__ids.begin(), __ids.end());
    }
    
    vec_type __ids;
};

struct Lt_id_set {
    template<typename T, int N>
    bool operator()(const id_set<T, N>& f0, const id_set<T, N>& f1)
    {
        nvis::lexicographical_order Lt;
        return Lt(f0.__ids, f1.__ids);
    }
};

typedef id_set<int, 2>        id2;
typedef id_set<int, 3>        id3;
typedef id_set<int, 4>        id4;

/*const nvis::ivec4 voxel_faces[6] = {
    nvis::ivec4(0, 1, 2, 3), nvis::ivec4(4, 5, 6, 7),
    nvis::ivec4(0, 1, 5, 4), nvis::ivec4(1, 2, 6, 5),
    nvis::ivec4(2, 3, 7, 6), nvis::ivec4(3, 0, 4, 7)
};

const nvis::ivec3 voxel_vertices[8] = {
    nvis::ivec3(0, 0, 0), nvis::ivec3(1, 0, 0),
    nvis::ivec3(1, 1, 0), nvis::ivec3(0, 1, 0),
    nvis::ivec3(0, 0, 1), nvis::ivec3(1, 0, 1),
    nvis::ivec3(1, 1, 1), nvis::ivec3(0, 1, 1)
};*/

template<typename Type, int N>
class grid {
    static Type __modulo(Type a, Type b)
    {
        Type r = fmod(a, b - 1);
        if (!r) {
            return 0;
        } else {
            return a >= 0 ? r : b + r - 1;
        }
    }
    
public:
    typedef Type                                     data_type;
    typedef int                                      index_type;
    typedef nvis::fixed_vector<index_type, N>        ivec_type;
    typedef nvis::fixed_vector<data_type, N>         vec_type;
    typedef nvis::fixed_vector<data_type, N>         vertex_type;
    typedef nvis::fixed_vector<bool, N>              bvec_type;
    typedef nvis::bounding_box<vertex_type>          bounds_type;
    typedef grid<data_type, N>                       self_type;
    typedef metric<data_type, N>                     metric_type;
    
    grid(const ivec_type& size, const vec_type& spacing)
        : __size(size), __spacing(spacing), __origin(0), __periodic(false),
          __verbose(false)
    {
        __metric.bounds() = bounds();
        __metric.periodic() = __periodic;
    }
    
    grid(const ivec_type& size, const bounds_type& _bounds)
        : __size(size), __periodic(false), __verbose(false)
    {
        __spacing = _bounds.size() / vec_type(size - ivec_type(1));
        __origin = _bounds.min();
        __metric.bounds() = _bounds;
        __metric.periodic() = __periodic;
    }
    
    grid(const ivec_type& size, const vec_type& spacing, const bvec_type& periodic)
        : __size(size), __spacing(spacing), __origin(0), __periodic(periodic),
          __verbose(false)
    {
        __metric.bounds() = bounds();
        __metric.periodic() = __periodic;
    }
    
    grid(const ivec_type& size, const vec_type& spacing, const vec_type& origin)
        : __size(size), __spacing(spacing), __origin(origin), __periodic(false),
          __verbose(false)
    {
        __metric.bounds() = bounds();
        __metric.periodic() = __periodic;
    }
    
    grid(const ivec_type& size, const vec_type& spacing, const vec_type& origin, const bvec_type& periodic)
        : __size(size), __spacing(spacing), __origin(origin), __periodic(periodic),
          __verbose(false)
    {
        __metric.bounds() = bounds();
        __metric.periodic() = __periodic;
    }
    
    void setNewBounds(const ivec_type& size, const bounds_type& _bounds)
    {
        __size = size;
        __spacing = _bounds.size() / vec_type(size - ivec_type(1));
        __origin = _bounds.min();
        __metric.bounds() = _bounds;
    }
    
    index_type size() const
    {
        index_type s = __size[0];
        for (int i = 1 ; i < N ; ++i) {
            s *= __size[i];
        }
        return s;
    }
    
    const ivec_type& dimensions() const
    {
        return __size;
    }
    
    static int dimension()
    {
        return N;
    }
    
    void verbose(bool v) const
    {
        __verbose = v;
    }
    
    vertex_type operator()(index_type i, index_type j) const
    {
        BOOST_STATIC_ASSERT(N == 2);
        return __origin + vec_type(imodulo(ivec_type(i, j)))*__spacing;
    }
    
    vertex_type operator()(index_type i, index_type j, index_type k) const
    {
        BOOST_STATIC_ASSERT(N == 3);
        return __origin + vec_type(imodulo(ivec_type(i, j, k)))*__spacing;
    }
    
    vertex_type operator()(const ivec_type& ids) const
    {
        return __origin + vec_type(imodulo(ids))*__spacing;
    }
    
    index_type index(index_type i, index_type j) const
    {
        BOOST_STATIC_ASSERT(N == 2);
        ivec_type c = modulo(ivec_type(i, j));
        return c[0] + c[1]*__size[0];
    }
    
    index_type index(index_type i, index_type j, index_type k) const
    {
        BOOST_STATIC_ASSERT(N == 3);
        ivec_type c = imodulo(ivec_type(i, j, k));
        return c[0] + __size[0]*(c[1] + __size[1]*c[2]);
    }
    
    index_type index(const ivec_type& c) const
    {
        ivec_type __c = imodulo(c);
        index_type idx = __c[0];
        index_type k = __size[0];
        for (int i = 1 ; i < N ; ++i) {
            idx += k * __c[i];
            k *= __size[i];
        }
        return idx;
    }
    
    ivec_type coordinates(index_type idx) const
    {
        index_type tmp = idx;
        ivec_type r;
        for (int i = 0 ; i < N - 1 ; ++i) {
            r[i] = tmp % __size[i];
            tmp /= __size[i];
        }
        r[N-1] = tmp;
        return r;
    }
    
    std::pair<ivec_type, vec_type> local_coordinates(const vec_type& x) const
    {
        vec_type y(x);
        y -= __origin;
        y /= __spacing;
        ivec_type c;
        for (int i = 0 ; i < N ; ++i) {
            c[i] = floor(y[i]);
        }
        y -= vec_type(c);
        c = imodulo(c);
        if (__verbose) {
            std::cout << "grid::local_coordinates: x = " << x << ", "
                      << " y = " << y << ", c = " << c << '\n';
        }
        return std::make_pair(c, y);
    }
    
    ivec_type imodulo(const ivec_type& c) const
    {
        if (!nvis::any(__periodic)) {
            return c;
        }
        ivec_type pc(c);
        for (int i = 0 ; i < N ; ++i) {
            if (__periodic[i]) {
                pc[i] = __modulo(c[i], __size[i]);
                // std::cerr << pc[i] << " = __modulo(" << c[i] << ", " << __size[i] << ")\n";
            } else if (c[i] < 0 || c[i] >= __size[i]) {
                // std::cerr << "c[" << i << "] = " << c[i] << " is out of range\n";
                throw std::runtime_error("grid::imodulo(): invalid coordinates");
            }
        }
        return pc;
    }
    
    vec_type dmodulo(const vec_type& x) const
    {
        vec_type y = __metric.modulo(x);
        if (!bounds().inside(y)) {
            throw std::runtime_error("grid::dmodulo(): invalid coordinates");
        }
        return y;
    }
    
    const metric_type& get_metric() const
    {
        return __metric;
    }
    
    bool on_boundary(const ivec_type& c) const
    {
        for (int i = 0 ; i < N ; ++i)
            if (!__periodic[i] && (c[i] == 0 || c[i] == __size[i] - 1)) {
                return true;
            }
        return false;
    }
    
    bvec_type which_boundary(const ivec_type& c) const
    {
        bvec_type r(false);
        for (int i = 0 ; i < N ; ++i)
            if (!__periodic[i] && (c[i] == 0 || c[i] == __size[i] - 1)) {
                r[i] = true;
            }
        return r;
    }
    
    data_type step(index_type i) const
    {
        return __spacing[i];
    }
    
    const vec_type& spacing() const
    {
        return __spacing;
    }
    
    bounds_type bounds() const
    {
        bounds_type b;
        b.min() = __origin;
        b.max() = b.min() + vec_type(__size - ivec_type(1)) * __spacing;
        return b;
    }
    
    std::list<index_type> neighbors(index_type idx) const
    {
        ivec_type c = coordinates(idx);
        std::list<index_type> r;
        for (int i = 0 ; i < N ; ++i) {
            index_type l = c[i];
            ivec_type previous(c), next(c);
            --previous[l];
            ++next[l];
            try {
                previous = modulo(previous);
                r.push_back(index(previous));
            } catch (...) {}
            try {
                next = modulo(next);
                r.push_back(index(next));
            } catch (...) {}
        }
        return r;
    }
    
private:
    ivec_type            __size, __offset;
    vec_type             __spacing, __origin;
    bvec_type            __periodic;
    metric_type          __metric;
    mutable bool         __verbose;
};

}


#endif






























