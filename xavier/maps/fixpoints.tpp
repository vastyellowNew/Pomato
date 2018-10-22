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


/// fixpoints.tpp - Source implementation for fixpoints.hpp functions with templates
//  Wayne Schlei & Xavier Tricoche
//  7/14/2013
#ifndef __FIXPOINT_TPP__
#define __FIXPOINT_TPP__

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <iostream>
#include "newton.hpp"
#include <math/angle.hpp>
#include <map>
//API
#include <maps/fixpoints.hpp>


namespace xavier {

extern bool record_search_steps;
extern std::vector<std::pair<nvis::vec2, nvis::vec2> > search_steps;

/// Useful for comparing coarse and fine datasets
template<typename MAP>
bool linear_analysis(const MAP& map, unsigned int period, const metric<double, 2>& metric,
                     const nvis::vec2& x, fixpoint& fp,
                     double hmin, double hmax)
{
    rhs_only_wrapper<MAP, 2> rhs(map, metric, period);
    nvis::mat2 J_coarse = rhs.jacobian(x, hmax);
    nvis::mat2 J_fine = rhs.jacobian(x, hmin);
    xavier::fixpoint fp_coarse, fp_fine;
    linear_analysis(J_coarse, period, x, fp_coarse);
    linear_analysis(J_fine, period, x, fp_fine);
    std::ostringstream os;
    os << "fp_coarse = " << fp_coarse << ", fp_fine = " << fp_fine << std::endl;
    std::cerr << os.str();
    if (similar(fp_coarse, fp_fine)) {
        fp = fp_fine;
        return true;
    } else {
        nvis::mat2 J_mid = rhs.jacobian(x, 0.5*(hmin + hmax));
        xavier::fixpoint fp_mid;
        linear_analysis(J_mid, period, x, fp_mid);
        if (similar(fp_coarse, fp_mid)) {
            fp = fp_mid;
        } else if (similar(fp_mid, fp_fine)) {
            fp = fp_fine;
        } else {
            fp = fp_mid;
        }
        return false;
    }
}

/// Search to find the quadrant
template<typename RHS>
void search(const RHS& rhs, map_quad& quad, int depth, bool must_proceed,
            std::vector<nvis::vec2>& found)
{

    double theta = 0;
    bool split = false;
    for (int e=0 ; e<4 && !split ; ++e) {
        double dtheta = signed_angle(quad.val(e), quad.val((e+1)%4));
        if (dtheta > 2./3.*M_PI) {
            split = true;
        } else {
            theta += dtheta;
        }
    }
    
    if (!split) {
        long int idx = lrint(0.5*theta/M_PI);
        if (idx == 0 && !must_proceed) {
            return;
        }
    }
    
    if (depth == 0) {
        found.push_back(0.5*(quad.pos(0) + quad.pos(2)));
        return;
    }
    
    nvis::vec2 midv[5];
    try {
        for (int i=0 ; i<4 ; ++i) {
            nvis::vec2 x = 0.5*(quad.pos(i) + quad.pos((i+1)%4));
            midv[i] = rhs(x);
            if (record_search_steps) {
                search_steps.push_back(std::pair<nvis::vec2, nvis::vec2>(x, midv[i]));
            }
        }
        midv[4] = rhs(0.5*(quad.pos(0) + quad.pos(2)));
    } catch(std::runtime_error& err) {
        return;
    }
    if (record_search_steps) {
        search_steps.push_back(std::pair<nvis::vec2, nvis::vec2>(0.5*(quad.pos(0) + quad.pos(2)), midv[4]));
    }
    
    map_quad subquads[4];
    subquads[0].bounds().min() = quad.pos(0);
    subquads[0].bounds().max() = 0.5*(quad.pos(0) + quad.pos(2));
    subquads[0].val(0) = quad.val(0);
    subquads[0].val(1) = midv[0];
    subquads[0].val(2) = midv[4];
    subquads[0].val(3) = midv[3];
    
    subquads[1].bounds().min() = 0.5*(quad.pos(0) + quad.pos(1));
    subquads[1].bounds().max() = 0.5*(quad.pos(1) + quad.pos(2));
    subquads[1].val(0) = midv[0];
    subquads[1].val(1) = quad.val(1);
    subquads[1].val(2) = midv[1];
    subquads[1].val(3) = midv[4];
    
    subquads[2].bounds().min() = 0.5*(quad.pos(0) + quad.pos(2));
    subquads[2].bounds().max() = quad.pos(2);
    subquads[2].val(0) = midv[4];
    subquads[2].val(1) = midv[1];
    subquads[2].val(2) = quad.val(2);
    subquads[2].val(3) = midv[2];
    
    subquads[3].bounds().min() = 0.5*(quad.pos(0) + quad.pos(3));
    subquads[3].bounds().max() = 0.5*(quad.pos(2) + quad.pos(3));
    subquads[3].val(0) = midv[3];
    subquads[3].val(1) = midv[4];
    subquads[3].val(2) = midv[2];
    subquads[3].val(3) = quad.val(3);
    
    search(rhs, subquads[0], depth-1, false, found);
    search(rhs, subquads[1], depth-1, false, found);
    search(rhs, subquads[2], depth-1, false, found);
    search(rhs, subquads[3], depth-1, false, found);
}

/// Find the seed point for fixed-point computation from a cell
template<typename RHS>
bool find_seed(const RHS& rhs, const metric<double, 2>& metric, const nvis::bbox2& bounds,
               nvis::vec2& first_guess, int depth)
{
    std::map<double, nvis::vec2> norm_to_location;
    
    std::vector<nvis::bbox2> cur_res;
    cur_res.push_back(bounds);
    
    try {
        for (int d=0 ; d<depth ; ++d) {
            std::vector<nvis::bbox2> next_res;
            for (int i=0 ; i<cur_res.size() ; ++i) {
                nvis::vec2 x = cur_res[i].center();
                nvis::vec2 v = rhs(x);
                norm_to_location.insert(std::pair<double, nvis::vec2>(nvis::norm(v), x));
                split_box(cur_res[i], next_res);
            }
            std::swap(next_res, cur_res);
        }
    } catch(...) {
        return false;
    }
    
    first_guess = norm_to_location.begin()->second;
    return true;
}

///Fixed-Point search using Newton's line search method (specifically with Richardson't extrapolation method
/// for the Jacobian Matrix at the current guess)
template<typename MAP>
bool meta_newton(const MAP& pmap, const metric<double, 2>& metric, const nvis::bbox2& bounds,
                 const nvis::vec2& first_guess, int depth,
                 int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                 double eps, double Jeps, bool verbose,
                 size_t maxiter, bool prefound)
{

    typedef rhs_only_wrapper<MAP, 2> rhs_type;
    
    rhs_type rhs(pmap, metric, period);
    iterates.resize(period);
    nvis::vec2 x=first_guess, f;
    try {
        f = rhs(x);
    } catch(...) {
        return false;
    }
    
    nvis::mat2 Jinit = rhs.jacobian(x);
    double dinit = nvis::det(Jinit);
    
    if (nvis::norm(f) > eps) {
        if (verbose) {
            std::cerr << "norm at seed point (" << first_guess << ") = " << nvis::norm(f) << " is too large\n";
        }
        bool ok = find_seed(rhs, metric, bounds, x, depth);
        if (!ok) {
            if (verbose) {
                std::cerr << "unable to find seed\n";
            }
            return false;
        } else if (verbose) {
            std::cerr << "initial Jacobian was: " << Jinit << " with determinant " << dinit << '\n';
            std::cerr << "improved seed found\n"
                      << "norm at new seed (" << x << ") is " << nvis::norm(rhs(x)) << '\n';
            for (int k=0 ; k<8 ; ++k) {
                double h = 0.1 / (double)(1 << k);
                nvis::mat2 Jfin = rhs.jacobian(x, h);
                double dfin = nvis::det(Jfin);
                std::cerr << "final Jacobian (h=" << h << ") is " << Jfin << " with determinant " << dfin << '\n';
            }
            
            nvis::mat2 Jrichardson = richardson<rhs_type, 2>(rhs, x, 0.1, 10.);
            std::cerr << "Richardson gives us " << Jrichardson << " with determinant "
                      << nvis::det(Jrichardson) << '\n';
                      
        }
        try {
            nvis::vec2 g = rhs(x);
            if (nvis::norm(f) < nvis::norm(g)) {
                x = first_guess;
            }
        } catch(...) {
            if (verbose) {
                std::cerr << "exception caught in metanewton\n";
            }
            return false;
        }
    }
    
    bool found = newton_in_box<rhs_type,2>(rhs, x, bounds, eps, maxiter, verbose);
    if (!found && !prefound) {
        return false;
    }
    if (!bounds.inside(x)) {
        // something has to be done here
    } else {
    }
    
    try {
        if (found) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP, 2> jrhs(pmap, metric, period, Jeps);
            iterates[0] = x;
            for (int i=1 ; i<period ; ++i) {
                iterates[i] = metric.modulo(amap->map(iterates[i-1], 1));
            }
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
        } else if (prefound) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP, 2> jrhs(pmap, metric, period, Jeps);
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            iterates.clear(); // we won't trust iterates from this location
        }
    } catch(...) {
        return false;
    }
    
    return true;
}

/// Newton iteration scheme for Standard Map -> Will assume map repeats itself in a particular dimension
template<typename MAP>
bool meta_newton_stdmap(const MAP& pmap, const metric<double, 2>& metric, const nvis::bbox2& bounds,
                        const nvis::vec2& first_guess, int depth,
                        int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                        double eps, bool verbose,
                        size_t maxiter, bool prefound)
{

    std::ostringstream os;
    
    if (verbose) {
        os << "meta-newton: guess = " << first_guess << ", period = " << period << '\n';
    }
    
    typedef rhs_wrapper<MAP, 2> rhs_type;
    
    rhs_type rhs(pmap, metric, period);
    iterates.resize(period);
    nvis::vec2 x=first_guess, f;
    try {
        f = rhs(x);
    } catch(...) {
        if (verbose) {
            os << "unable to compute map at seed point" << std::endl << std::flush;
            std::cerr << os.str();
        }
        return false;
    }
    
    if (nvis::norm(f) > eps) {
        bool ok = find_seed(rhs, metric, bounds, x, depth);
        if (!ok) {
            return false;
        }
        try {
            nvis::vec2 g = rhs(x);
            if (nvis::norm(f) < nvis::norm(g)) {
                x = first_guess;
            }
        } catch(...) {
            return false;
        }
    }
    
    bool found = newton_in_box<rhs_type,2>(rhs, x, bounds, eps, maxiter, verbose);
    if (!found && !prefound) {
        if (verbose) {
            os << "FAILED" << std::endl;
            std::cerr << os.str();
        }
        return false;
    }
    if (!bounds.inside(x)) {
        if (verbose) {
            os << "left boundaries" << std::endl;
        }
    } else {
        if (verbose) {
            os << "SUCCESSFUL" << std::endl;
        }
    }
    
    if (verbose && found) {
        os << "\n\nfound a critical point of period " << period << " at " << x << '\n'
           << "processing iterates\n" << std::flush;
    }
    
    try {
        if (found) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP, 2> jrhs(pmap, metric, period, 0);
            iterates[0] = x;
            for (int i=1 ; i<period ; ++i) {
                iterates[i] = metric.modulo(amap->map(iterates[i-1], 1));
                if (verbose) {
                    os << "iterate #" << i << " from " << x << " is " << iterates[i] << std::endl << std::flush;
                }
            }
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            if (verbose) {
                os << "found " << fp << std::endl << std::flush;
            }
        } else if (prefound) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP, 2> jrhs(pmap, metric, period, 0);
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            iterates.clear(); // we won't trust iterates from this location
        }
    } catch(...) {
        std::cerr << os.str();
        return false;
    }
    
    std::cerr << os.str();
    return true;
}


///Fixed-Point search using Newton's line search method (A simpler approach than meta_newton - also different jacobian format)
template<typename MAP>
bool simple_newton(const MAP& pmap, const metric<double, 2>& metric, const nvis::bbox2& bounds,
                   const nvis::vec2& first_guess, int depth,
                   int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                   double eps, double Jeps, bool verbose,
                   size_t maxiter, bool prefound)
{

    typedef rhs_wrapper<MAP, 2> rhs_type;
    
    rhs_type rhs(pmap, metric, period);
    iterates.resize(period);
    nvis::vec2 x=first_guess, f;
    try {
        f = rhs(x);
    } catch(...) {
        return false;
    }
    
    nvis::mat2 Jinit = rhs.jacobian(x);
    double dinit = nvis::det(Jinit);
    
    //If given first guess is inadequate
    if (nvis::norm(f) > eps) {
        if (verbose) {
            std::cerr << "norm at seed point (" << first_guess << ") = " << nvis::norm(f) << " is too large\n";
        }
        //Try to find a better seed point with a brute-force subsample
        bool ok = find_seed(rhs, metric, bounds, x, depth);
        if (!ok) {
            if (verbose) {
                std::cerr << "unable to find seed\n";
            }
            return false;
        } else if (verbose) {
            std::cerr << "initial Jacobian was: " << Jinit << " with determinant " << dinit << '\n';
            std::cerr << "improved seed found\n"
                      << "  norm at new seed (" << x << ") is " << nvis::norm(rhs(x)) << '\n';
            //Recompute Jacobian
            nvis::mat2 Jnew = rhs.jacobian(x);
            double dnew = nvis::det(Jnew);
            std::cerr << "  new Jacobian: " << Jnew << " with determinant " << dnew << '\n';
        }
        
        //If new seed isn't any better, just take the original
        try {
            nvis::vec2 g = rhs(x);
            if (nvis::norm(f) < nvis::norm(g)) {
                x = first_guess;
            }
        } catch(...) {
            if (verbose) {
                std::cerr << "exception caught in simple_newton\n";
            }
            return false;
        }
    }
    
    //Run the newton search
    bool found = newton_in_box<rhs_type,2>(rhs, x, bounds, eps, maxiter, verbose);
    if (!found && !prefound) {
        return false;
    }
    if (!bounds.inside(x)) {
        // If not inside the bounds
        // Do something .. has yet to be done -> Skip
    } else {
    }
    
    try {
        if (found) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP, 2> jrhs(pmap, metric, period, Jeps);
            iterates[0] = x;
            for (int i=1 ; i<period ; ++i) {
                iterates[i] = amap->map(iterates[i-1], 1);
            }
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
        } else if (prefound) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP, 2> jrhs(pmap, metric, period, Jeps);
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            iterates.clear(); // we won't trust iterates from this location
        }
    } catch(...) {
        return false;
    }
    
    return true;
}


template<typename MAP>
bool linear_chain_analysis(const MAP& pmap, const metric<double, 2>& metric,
                           const nvis::vec2& first_guess, int period,
                           std::vector<fixpoint>& fps, double maxlength,
                           double eps, double Jeps, bool verbose,
                           size_t maxiter)
{
    rhs_wrapper<MAP, 2> rhs(pmap, metric, period, Jeps);
    fps.resize(period);
    
    nvis::vec2 x = first_guess;
    
    std::ostringstream os;
    
    bool saddle;
    for (int i=0 ; i<period ; ++i) {
        // find precise location of this fixed point
        bool found = newton(rhs, x, eps, maxiter, maxlength, verbose);
        if (!found) {
            return false;
        }
        
        fixpoint& fp = fps[i];
        nvis::mat2 J = rhs.jacobian(x);
        linear_analysis(J, period, x, fp);
        // check type consistency
        if (!i) {
            saddle = (fp.saddle ? 1 : 0);
        } else if ((saddle && !fp.saddle) || (!saddle && fp.saddle)) {
            os << "incompatible linear types for period " << period << " - giving up" << std::endl;
            std::cerr << os.str();
            return false;
        }
        if (i<period-1) {
            x = metric.modulo(pmap.map(x, 1));
        }
    }
    
    return true;
}

} //xavier

#endif