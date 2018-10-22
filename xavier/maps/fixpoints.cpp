// fixpoints.tpp
//  Author: Wayne Schlei
//  Date:  6/15/2013
// Purpose:  Store functions relevant to fixpoints.hpp.
// Aviod redefinitions when used in multiple classes.

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <iostream>
#include "newton.hpp"
#include <math/angle.hpp>
#include <map>
//API
#include <maps/fixpoints.hpp>
using namespace xavier;

//Need to have these compiled to avoid multiple definitions

namespace xavier {

extern bool record_search_steps;
extern std::vector<std::pair<nvis::vec2, nvis::vec2> > search_steps;

///Nullspace of 2x2 matrix
inline nvis::vec2 nullspace(const nvis::mat2& A)
{
    nvis::vec2 row1(A[0]);
    nvis::vec2 row2(A[1]);
    nvis::vec2 e;
    if (nvis::norm(row1) > nvis::norm(row2)) {
        e = nvis::vec2(-row1[1], row1[0]);
    } else {
        e = nvis::vec2(-row2[1], row2[0]);
    }
    e /= nvis::norm(e);
    return e;
}

///Eigen values/vectors of 2x2
inline bool eigen(nvis::vec2 evecs[2], double evals[2], const nvis::mat2& J)
{
    double tr = nvis::trace(J);
    double det = nvis::det(J);
    double delta = tr * tr - 4 * det;
    
    if (delta >= 0) { // real eigenvalues: saddle type
        evals[0] = 0.5 * (tr - sqrt(delta)); // small eigenvalue first
        evals[1] = 0.5 * (tr + sqrt(delta));
        
        nvis::mat2 A(J);
        A[0][0] -= evals[0];
        A[1][1] -= evals[0];
        evecs[0] = nullspace(A);
        
        A = J;
        A[0][0] -= evals[1];
        A[1][1] -= evals[1];
        evecs[1] = nullspace(A);
        
        return true;
    }
    
    return false;
}

///Operator for printing fixpoint to ostream
std::ostream& operator<<(std::ostream& os, const fixpoint& fp)
{
    os << "fixed point: pos=" << fp.pos
       << ", period=" << fp.K
       << ", " << (fp.saddle ? "saddle" : "center");
    if (fp.saddle) {
        os << ", ev0=" << fp.evec[0] << ", ev1=" << fp.evec[1];
    }
    
    return os;
}

///Test if two fixed points are the same type
bool similar(const xavier::fixpoint& fp0, const xavier::fixpoint& fp1)
{
    if ((fp0.saddle && !fp1.saddle) || (fp1.saddle && !fp0.saddle)) {
        return false;
    } else if (!fp0.saddle) {
        return true;
    } else if (fabs(nvis::inner(fp0.evec[0], fp1.evec[0])) < 0.95) {
        return false;
    } else if (fabs(nvis::inner(fp0.evec[1], fp1.evec[1])) < 0.95) {
        return false;
    } else {
        return true;
    }
}

///Setting a fixpoint object?
void linear_analysis(const nvis::mat2& J, unsigned int period,
                     const nvis::vec2& x, xavier::fixpoint& fp)
{
    fp.pos = x;
    fp.K = period;
    fp.saddle = eigen(fp.evec, fp.eval, J);
}




/// Split the box into 4 quadrants
void split_box(const nvis::bbox2& box, std::vector<nvis::bbox2>& subs)
{
    nvis::vec2 center = box.center();
    const nvis::vec2& min = box.min();
    const nvis::vec2& max = box.max();
    nvis::bbox2 sub(min, center);
    subs.push_back(sub);
    sub.min() = center;
    sub.max() = max;
    subs.push_back(sub);
    sub.min() = nvis::vec2(center[0], min[1]);
    sub.max() = nvis::vec2(max[0], center[1]);
    subs.push_back(sub);
    sub.min() = nvis::vec2(min[0], center[1]);
    sub.max() = nvis::vec2(center[0], max[1]);
    subs.push_back(sub);
}



} //end xavier