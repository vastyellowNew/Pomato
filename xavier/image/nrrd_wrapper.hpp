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


#ifndef __NRRD_WRAPPER_HPP__
#define __NRRD_WRAPPER_HPP__

#include <teem/nrrd.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <math/bounding_box.hpp>
#include <math/fixed_vector.hpp>

namespace xavier {

template<typename T>
inline bool invalid(T v)
{
    return (std::isnan(v) || std::isinf(v));
}

template<int N> inline
nvis::bounding_box<nvis::fixed_vector<double, N> > bounds(const Nrrd* nrrd)
{
    typedef nvis::fixed_vector<double, N>   vec_type;
    typedef nvis::bounding_box<vec_type>    bounds_type;
    
    bounds_type b;
    for (int i = 0 ; i < N ; ++i) {
        const NrrdAxisInfo& axis = nrrd->axis[nrrd->dim-N+i];
        double step = axis.spacing;
        if (invalid(step)) {
            step = 1;
        }
        double width = (axis.size - 1) * step;
        if (!invalid(axis.min)) {
            b.min()[i] = axis.min;
            b.max()[i] = axis.min + width;
        } else if (!invalid(axis.max)) {
            b.max()[i] = axis.max;
            b.min()[i] = axis.max - width;
        } else {
            b.min()[i] = 0;
            b.max()[i] = width;
        }
        if (axis.center == nrrdCenterCell) {
            b.min()[i] -= 0.5 * step;
            b.max()[i] += 0.5 * step;
        }
    }
    return b;
}

template<int N> inline
nvis::fixed_vector<double, N> step(const Nrrd* nrrd)
{
    nvis::fixed_vector<double, N> s;
    for (int i = 0 ; i < N ; ++i) {
        s[i] = nrrd->axis[nrrd->dim-N+i].spacing;
        if (invalid(s[i])) {
            s[i] = 1;
        }
    }
    return s;
}

template<typename T>
struct nrrd_data_wrapper {
    nrrd_data_wrapper(const Nrrd* nrrd) : __nrrd(nrrd) {}
    T operator[](size_t i) const
    {
        switch (__nrrd->type) {
            case nrrdTypeChar:
                return ((char*)__nrrd->data)[i];
            case nrrdTypeUChar:
                return ((unsigned char*)__nrrd->data)[i];
            case nrrdTypeShort:
                return ((short*)__nrrd->data)[i];
            case nrrdTypeUShort:
                return ((unsigned short*)__nrrd->data)[i];
            case nrrdTypeInt:
                return ((int*)__nrrd->data)[i];
            case nrrdTypeUInt:
                return ((unsigned int*)__nrrd->data)[i];
            case nrrdTypeLLong:
                return ((long int*)__nrrd->data)[i];
            case nrrdTypeULLong:
                return ((unsigned long int*)__nrrd->data)[i];
            case nrrdTypeFloat:
                return ((float*)__nrrd->data)[i];
            case nrrdTypeDouble:
                return ((double*)__nrrd->data)[i];
            default:
                std::cerr << "unrecognized data type\n";
                throw std::runtime_error("unrecognized data type\n");
        }
    }
    T operator()(size_t i) const
    {
        return (*this)[i];
    }
    
    const Nrrd* __nrrd;
};

struct nrrd_params {

    nrrd_params() : data_type(nrrdTypeUnknown) {}
    
    int data_type;
    std::vector<double> mins, spc;
    std::vector<size_t> dims;
    std::vector<bool> cell_based;
};

template<typename T1, typename T2>
inline void to_vector(std::vector<T1>& vals, const void* data, size_t size)
{
    vals.resize(size);
    for (size_t i = 0 ; i < size ; ++i) {
        vals[i] = ((T2*)data)[i];
    }
}

template<typename T>
inline void to_vector(std::vector<T>& vals, const Nrrd* nin)
{
    size_t size = 1;
    for (size_t i = 0 ; i < nin->dim ; ++i) {
        size *= nin->axis[i].size;
    }
    // std::cerr << "size = " << size << std::endl;
    
    vals.clear();
    switch (nin->type) {
        case nrrdTypeChar:
            return to_vector<T, char>(vals, nin->data, size);
        case nrrdTypeUChar:
            return to_vector<T, unsigned char>(vals, nin->data, size);
        case nrrdTypeShort:
            return to_vector<T, short>(vals, nin->data, size);
        case nrrdTypeUShort:
            return to_vector<T, unsigned short>(vals, nin->data, size);
        case nrrdTypeInt:
            return to_vector<T, int>(vals, nin->data, size);
        case nrrdTypeUInt:
            return to_vector<T, unsigned int>(vals, nin->data, size);
        case nrrdTypeLLong:
            return to_vector<T, long int>(vals, nin->data, size);
        case nrrdTypeULLong:
            return to_vector<T, unsigned long int>(vals, nin->data, size);
        case nrrdTypeFloat:
            return to_vector<T, float>(vals, nin->data, size);
        case nrrdTypeDouble:
            return to_vector<T, double>(vals, nin->data, size);
        default:
            std::cerr << "unrecognized data type\n";
            throw std::runtime_error("unrecognized data type\n");
    }
}

template<typename T1, typename T2>
inline T1* to_array(const void* data, size_t size)
{
    T1* array = (T1*)calloc(size, sizeof(T1));
    for (int i = 0 ; i < size ; ++i) {
        array[i] = ((T2*)data)[i];
    }
    return array;
}

template<typename T>
inline T* to_array(const Nrrd* nin)
{
    size_t size = 1;
    for (int i = 0 ; i < nin->dim ; ++i) {
        size *= nin->axis[i].size;
    }
    
    switch (nin->type) {
        case nrrdTypeChar:
            return to_array<T, char>(nin->data, size);
        case nrrdTypeUChar:
            return to_array<T, unsigned char>(nin->data, size);
        case nrrdTypeShort:
            return to_array<T, short>(nin->data, size);
        case nrrdTypeUShort:
            return to_array<T, unsigned short>(nin->data, size);
        case nrrdTypeInt:
            return to_array<T, int>(nin->data, size);
        case nrrdTypeUInt:
            return to_array<T, unsigned int>(nin->data, size);
        case nrrdTypeLLong:
            return to_array<T, long int>(nin->data, size);
        case nrrdTypeULLong:
            return to_array<T, unsigned long int>(nin->data, size);
        case nrrdTypeFloat:
            return to_array<T, float>(nin->data, size);
        case nrrdTypeDouble:
            return to_array<T, double>(nin->data, size);
        default:
            std::cerr << "unrecognized data type\n";
            throw std::runtime_error("unrecognized data type\n");
    }
}

template<int N>
nvis::bounding_box< nvis::fixed_vector<double, N> >
inline compute_bounds(const Nrrd* nrrd)
{
    nvis::bounding_box<nvis::fixed_vector<double, N> > box;
    for (int i = 0 ; i < N ; ++i) {
        const NrrdAxisInfo& axis = nrrd->axis[nrrd->dim-N+i];
        box.min()[i] = axis.min;
        if (std::isnan(box.min()[i])) {
            box.min()[i] = 0.;
        }
        double step = axis.spacing;
        if (std::isnan(step)) {
            step = 1.;
        }
        box.max()[i] = box.min()[i] + step * (axis.size - 1);
    }
    return box;
}

inline Nrrd* readNrrd(const std::string& filename)
{
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    return nin;
}

inline void writeNrrd(void* data, const std::string& filename, int data_type, const std::vector< size_t >& dims,
                      const std::vector< double >& spacing)
{
    Nrrd* nout = nrrdNew();
    
    if (nrrdWrap_nva(nout, data, data_type, dims.size(), &dims[0])) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    if (spacing.size() == dims.size()) {
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&spacing[0]);
    }
    if (nrrdSave(filename.c_str(), nout, NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
}

inline void writeNrrd(void* data, const std::string& filename, int data_type, const std::vector< size_t >& dims)
{
    Nrrd* nout = nrrdNew();
    
    if (nrrdWrap_nva(nout, data, data_type, dims.size(), &dims[0])) {
        std::cerr << "writeNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    std::vector<double> spacing(dims.size());
    std::fill(spacing.begin(), spacing.end(), 1);
    if (spacing.size() == dims.size()) {
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&spacing[0]);
    }
    if (nrrdSave(filename.c_str(), nout, NULL)) {
        std::cerr << "writeNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
}

inline void writeNrrd(void* data, const std::string& filename, const nrrd_params& params)
{
    if (params.data_type == nrrdTypeUnknown) {
        std::cerr << "writeNrrd(): data type is undefined. Data will not be saved.\n";
        return;
    } else if (!params.dims.size()) {
        std::cerr << "writeNrrd(): dimensions were not provided. Data will not be saved.\n";
        return;
    }
    size_t dim = params.dims.size();
    
    Nrrd* nout = nrrdNew();
    
    if (nrrdWrap_nva(nout, data, params.data_type, dim, &params.dims[0])) {
        std::cerr << "writeNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    if (params.spc.size() == dim) {
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&params.spc[0]);
    }
    if (params.mins.size() == dim) {
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, (const void*)&params.mins[0]);
    }
    if (params.cell_based.size() == dim) {
        int cent[dim];
        for (size_t i = 0 ; i < dim ; ++i) {
            cent[i] = (params.cell_based[i] ? nrrdCenterCell : nrrdCenterNode);
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, (const void*)cent);
    }
    if (nrrdSave(filename.c_str(), nout, NULL)) {
        std::cerr << "writeNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
}

}

#endif











































