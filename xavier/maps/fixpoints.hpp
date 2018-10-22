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


//Detection of Fixed-points on the map and storage
//Authors:  Xavier Tricoche and Wayne Schlei
//Date:  6/15/2013

#ifndef __MAPS_LIB_FIXPOINT_HPP__
#define __MAPS_LIB_FIXPOINT_HPP__

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <cstdio>
#include <exception>
#include <iostream>
#include "newton.hpp"
#include <math/angle.hpp>
#include <map>

namespace xavier {

/// Fixed point structure:  Stores the Eigen vals/vecs, type, iterates, and isolation classifier
struct fixpoint {
    typedef nvis::fixed_vector<double,6>   StateType; //Need to template!
    typedef StateType                      vec6;
    typedef nvis::vec2                     VecType;
    ///Base constructor:  Creates a default fixpoint at origin with negative time
    fixpoint() : pos(0.0), saddle(false), K(-1), orbitID(-1), t(-1.0), timePeriod(-1.0),
#if defined(HX_HAS_STD) || defined(C_0X)
        si(0.0), otherSI(0.0), isolated(true) {}
#else
        si(0.0), otherSI(0.0), isolated(true),
        evec{0.0}, fullEvec{0.0}, eval{0.0} {}
#endif
        
    ///Copy Constructor:  Build a fixpoint from another
    fixpoint(const fixpoint& fp)
        : pos(fp.pos), saddle(fp.saddle), K(fp.K),
          orbitID(fp.orbitID), t(fp.t), timePeriod(fp.timePeriod),
#if defined(HX_HAS_STD) || defined(C_0X)
          si(fp.si), otherSI(fp.otherSI), isolated(fp.isolated)
#else
          si(fp.si), otherSI(fp.otherSI), isolated(fp.isolated),
          evec{0.0}, fullEvec{0.0}, eval{0.0}
#endif
    {
        evec[0] = fp.evec[0];
        evec[1] = fp.evec[1];
        fullEvec[0] = fp.fullEvec[0];
        fullEvec[1] = fp.fullEvec[1];
        eval[0] = fp.eval[0];
        eval[1] = fp.eval[1];
    }
    
    /// Map Coordinate of fixed point
    VecType pos;
    /// Is the fixed point of type: saddle?
    bool saddle;
    /// Eigenvectors (0-stable, 1-unstable)
    VecType evec[2];
    /// Full-state eigenvectors (0-stable, 1-unstable)
    StateType fullEvec[2];
    /// Eigenvalues (0-stable, 1-unstable - if saddle)
    double eval[2];
    /// Period parameter (K or p - number of iterates)
    unsigned int K;
    /// Orbit Index:  represents the parent orbit of this fixed point
    int orbitID;
    /// Orbital Period (time)
    double t, timePeriod;
    /// Stability indexes (si->In-Plane, otherSI->Out-of-Plane)
    double si, otherSI;
    /// Isolated orbit after filtering? [Not really used]
    bool isolated;
};

/// Map quadrant structure
struct map_quad {

    map_quad() {}
    
    nvis::vec2 pos(int i) const
    {
        switch (i) {
            case 0:
                return _bounds.min();
            case 1:
                return nvis::vec2(_bounds.max()[0], _bounds.min()[1]);
            case 2:
                return _bounds.max();
            case 3:
                return nvis::vec2(_bounds.min()[0], _bounds.max()[1]);
			default:
			throw std::runtime_error("invalid value in map_quad::pos(int) const");
        }
    }
    
    nvis::bbox2& bounds()
    {
        return _bounds;
    }
    
    nvis::vec2& val(int i)
    {
        return _v[i];
    }
    
    nvis::bbox2         _bounds;
    nvis::vec2          _v[4];
};


/// Not sure:  Position displacement check
struct Lt_pos_epsilon {
    Lt_pos_epsilon(double epsilon) : _eps(epsilon), _Lt() {}
    
    bool operator()(const nvis::vec2& x0, const nvis::vec2& x1) const
    {
        if (nvis::norm(x0-x1) < _eps) {
            return false;
        } else {
            return _Lt(x0, x1);
        }
    }
    
    double _eps;
    nvis::lexicographical_order _Lt;
};

/// Class to check if a fixpoint is close to another and at a higher period for filtering
struct is_highperiod_duplicate {
    is_highperiod_duplicate() : eps(1e-4) {}
    is_highperiod_duplicate(const double& epsilon) : eps(epsilon) {}
    
    bool operator()(fixpoint& first, fixpoint& second)
    {
        if (nvis::norm(first.pos - second.pos) < eps) {
            //Check that first_period is a multiple (or same) as second_period
            if ( first.K % second.K == 0 ) {
                return true; //unique->removes first
            } else {
                return false;
            }
        } else {
            return false;
        }
    }
    
private :
    double eps;
};

/// Class predicate to sort with high periods in front, low periods in back
template<typename T>
bool fpPeriodCompare(const T& fp1, const T& fp2)
{
    //Check if period is higher (opposite of less than call)
    if (fp1.K > fp2.K) {
        return true;
    } else if (fp1.K < fp2.K) {
        return false;
    } else { //fp1.K = fp2.K, sort by position
        nvis::lexicographical_order _Lt;
        return _Lt(fp1.pos,fp2.pos);
    }
}


//Function Declarations:
//Note:  Implementations of all functions moved to fixpoints.cpp to avoid redefinitions
///Nullspace of 2x2 matrix
inline nvis::vec2 nullspace(const nvis::mat2& A);

///Eigen values/vectors of 2x2 matrix
inline bool eigen(nvis::vec2 evecs[2], double evals[2], const nvis::mat2& J);

///Operator for printing fixpoint to ostream
std::ostream& operator<<(std::ostream& os, const fixpoint& fp);

///Test if two fixed points are the same type
bool similar(const xavier::fixpoint& fp0, const xavier::fixpoint& fp1);

///Write BASIC vector of xavier::fixpoint conditions to a file (p x0 x1)
template<typename FPVEC>
bool writeFixpointVecBasic(const char* filename, const FPVEC& fpv)
{
    FILE* f = fopen(filename,"w");
    if (!f) {
        std::cerr << "Write Basic Fixpoint Data FAILED!\n";
        return false;
    }
    
    int numFPs = (int) fpv.size();
    fprintf(f,"# Basic Fixed Point Information\n");
    fprintf(f,"%d\n",numFPs);
    for(int n=0; n<numFPs; n++) {
        int p = fpv[n].K;
        nvis::vec2 x = fpv[n].pos;
        //Write the BASIC line (p,x0,x1)
        fprintf(f,"%d %.15f %.15f\n",p,x[0],x[1]);
    }
    fclose(f);
    return true;
}

///Read BASIC vector of xavier::fixpoint conditions to a file (p x0 x1)
template<typename FP>
bool readFixpointVecBasic(const char* filename, std::vector<FP>& fpv)
{
    FILE* f = fopen(filename,"r");
    if (!f) {
        std::cerr << "Read Basic Fixpoint Data FAILED!\n";
        return false;
    }
    char buf[80];
    fgets(buf, 80, f);
    
    //Read size
    int numFPs;
    fscanf(f,"%d",numFPs);
    
    //Read all the fixed points
    for (int n=0; n<numFPs; n++) {
        int p;
        nvis::vec2 pos;
        fscanf(f,"%d %lf %lf",&p,&(pos[0]),&(pos[1]));
        FP afp;
        afp.K = p;
        afp.pos = pos;
        fpv.push_back(afp);
    }
    fclose(f);
    return true;
}

///Setting a fixpoint object based on input J, p, and x (performs linear analysis)
void linear_analysis(const nvis::mat2& J, unsigned int period,
                     const nvis::vec2& x, xavier::fixpoint& fp);
                     
/*  TEMPLATE FUNCTIONS : PROTOTYPES */
/// Sets a fixpoint with linear analysis (Useful for comparing coarse and fine datasets)
template<typename MAP>
bool linear_analysis(const MAP& map, unsigned int period, const metric<double, 2>& metric,
                     const nvis::vec2& x, xavier::fixpoint& fp,
                     double hmin, double hmax);
                     
/// Split a bounding box into 4 quadrants (quadtree?)
void split_box(const nvis::bbox2& box, std::vector<nvis::bbox2>& subs);

/// Search to find the quadrant
template<typename RHS>
void search(const RHS& rhs, map_quad& quad, int depth, bool must_proceed,
            std::vector<nvis::vec2>& found);
            
/// Find the seed point for fixed-point computation from a cell
template<typename RHS>
bool find_seed(const RHS& rhs, const metric<double, 2>& metric, const nvis::bbox2& bounds,
               nvis::vec2& first_guess, int depth);
               
///Fixed-Point search using Newton's line search method (specifically with Richardson't extrapolation method
/// for the Jacobian Matrix at the current guess)
template<typename MAP>
bool meta_newton(const MAP& pmap, const metric<double, 2>& metric, const nvis::bbox2& bounds,
                 const nvis::vec2& first_guess, int depth,
                 int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                 double eps, double Jeps=0, bool verbose=false,
                 size_t maxiter=50, bool prefound=false);
                 
                 
/// Newton iteration scheme for Standard Map -> Will assume map repeats itself in a particular dimension
template<typename MAP>
bool meta_newton_stdmap(const MAP& pmap, const metric<double, 2>& metric, const nvis::bbox2& bounds,
                        const nvis::vec2& first_guess, int depth,
                        int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                        double eps, bool verbose=false,
                        size_t maxiter=50, bool prefound=false);
                        
///Fixed-Point search using quasi-Newton method (A simpler approach than meta_newton)
template<typename MAP>
bool simple_newton(const MAP& pmap, const metric<double, 2>& metric, const nvis::bbox2& bounds,
                   const nvis::vec2& first_guess, int depth,
                   int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                   double eps, double Jeps=0, bool verbose=false,
                   size_t maxiter=50, bool prefound=false);
                   
/// Linear analysis for chains?
template<typename MAP>
bool linear_chain_analysis(const MAP& pmap, const metric<double, 2>& metric,
                           const nvis::vec2& first_guess, int period,
                           std::vector<fixpoint>& fps, double maxlength,
                           double eps, double Jeps=0, bool verbose=false,
                           size_t maxiter=50);
}

/* TEMPLATE FUNCTIONS : SOURCE IMPLEMENTATION */
#include <maps/fixpoints.tpp>


#endif
